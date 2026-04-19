//! End-to-end alignment driver.
//!
//! Pipeline per PLAN.md §3.1:
//! ```text
//! target ──► PosTable ──┐
//!                       ├─► seed_search ──► DiagHash ──► HSP ──► chain ──► anchor ──► gapped_extend ──► Record
//! query  ──► seeds  ────┘
//! ```
//! `seed_search` and `chain` are optional steps: `--nogapped` skips from HSP
//! straight to Record, `--nochain` skips the chain filter. Everything runs
//! under a rayon work-stealing pool fanning out over `(target × strand)`
//! pairs.

use std::sync::Mutex;

use rayon::prelude::*;

use crate::anchor::choose_anchor;
use crate::chain::best_chain;
use crate::edit_script::{EditOp, EditScript};
use crate::gapped_extend::{extend as gapped_extend, GappedParams};
use crate::hsp::{Hsp, HspParams};
use crate::output::{Record, Strand};
use crate::pos_table::PosTable;
use crate::scoring::ScoringMatrix;
use crate::seed_search::{search, SearchParams};
use crate::seeds::SeedPattern;
use crate::sequences::{PackedSeq, Sequence};
use crate::tweener::{interpolate, TweenerConfig};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandSpec {
    Plus,
    Minus,
    Both,
}

impl StrandSpec {
    fn strands(self) -> &'static [Strand] {
        match self {
            StrandSpec::Plus => &[Strand::Plus],
            StrandSpec::Minus => &[Strand::Minus],
            StrandSpec::Both => &[Strand::Plus, Strand::Minus],
        }
    }
}

/// Public configuration for `run()`.
#[derive(Debug, Clone)]
pub struct Config {
    pub pattern: SeedPattern,
    pub matrix: ScoringMatrix,
    pub step: usize,
    pub hsp: HspParams,
    pub gapped: GappedParams,
    pub strand: StrandSpec,
    /// Maximum seed-word multiplicity before a word is dropped as "hot". 0
    /// disables the filter.
    pub max_word_count: u32,
    /// Honor soft-masking (lowercase `acgt`) when picking seed positions.
    /// `true` matches upstream lastz's default; set to `false` for
    /// `--nomasking`.
    pub respect_masking: bool,
    /// Run gapped affine extension after HSP. When false, the pipeline
    /// stops at HSPs (equivalent to upstream `--nogapped`).
    pub gapped_enabled: bool,
    /// Run chaining between HSPs and gapped extension. When false,
    /// every HSP is extended independently (equivalent to `--nochain`).
    pub chain_enabled: bool,
    /// Sliding window width for anchor selection inside each HSP.
    pub anchor_window: u32,
    /// Transition substitutions tolerated per seed (upstream
    /// `--transition`/`--notransition`; 0 = none, 1 = one, 2 = two).
    pub transitions: u8,
    /// Within-target chunking. When the target's length exceeds
    /// `chunk.primary_bp`, it is split into overlapping chunks of that
    /// size plus `chunk.halo_bp` of overlap on each side, and the
    /// pipeline is rayon-parallelised over `(chunk × query × strand)`
    /// triples. `primary_bp = 0` (the default) disables chunking.
    pub chunk: ChunkConfig,
    /// When `Some`, apply KegAlign's Shannon-entropy gate to every HSP:
    /// drop HSPs whose target-slice entropy, scaled against the raw
    /// score, is below this threshold. Opt-in via `--entropy`.
    pub entropy_threshold: Option<i32>,
    /// Inter-alignment interpolation (upstream `tweener.c`). When `Some`,
    /// runs after chaining on each `(target, query, strand)` triple with
    /// the supplied tweener parameters.
    pub tweener: Option<TweenerConfig>,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            pattern: SeedPattern::twelve_of_nineteen(),
            matrix: ScoringMatrix::hoxd70(),
            step: 1,
            hsp: HspParams::default(),
            gapped: GappedParams::default(),
            strand: StrandSpec::Both,
            max_word_count: 0,
            respect_masking: true,
            gapped_enabled: true,
            chain_enabled: false,
            anchor_window: 31,
            transitions: 1,
            entropy_threshold: None,
            chunk: ChunkConfig::default(),
            tweener: None,
        }
    }
}

/// Within-target chunking configuration.
#[derive(Debug, Clone, Copy)]
pub struct ChunkConfig {
    /// Primary region length per chunk. `0` disables chunking entirely —
    /// the whole target is one chunk.
    pub primary_bp: u32,
    /// Overlap on each side of each primary region. Must be large enough
    /// to contain any single alignment that starts near a chunk boundary,
    /// otherwise alignments crossing boundaries will be truncated by
    /// `y_drop`. Typical value: `~50_000` bp for HOXD70 defaults.
    pub halo_bp: u32,
}

impl Default for ChunkConfig {
    fn default() -> Self {
        Self { primary_bp: 0, halo_bp: 0 }
    }
}

/// One target chunk produced by `chunk_target`. Coordinates in `seq` are
/// chunk-local; `global_start` maps chunk-local position 0 to its index
/// in the parent target sequence. HSPs emitted while aligning against
/// this chunk need `global_start` added to their `t_start`.
#[derive(Debug, Clone)]
struct TargetChunk {
    parent_index: usize,
    global_start: u32,
    seq: PackedSeq,
    /// Chunk-local coordinate cut-off beyond which an HSP's anchor is
    /// considered to live in the next chunk's primary region and should
    /// be suppressed. Equal to `primary_bp + halo_bp` for non-final
    /// chunks; `seq.len()` for the final chunk (so nothing is dropped).
    primary_end_local: u32,
    /// Chunk-local offset where this chunk's primary region starts.
    /// HSPs anchored before this offset live in the previous chunk's
    /// primary region and should be suppressed to avoid double-emission.
    /// Equal to `halo_bp` for non-first chunks; `0` for the first.
    primary_start_local: u32,
}

/// Split a target sequence into overlapping chunks. Returns a single
/// whole-sequence chunk when `cfg.primary_bp == 0` or the target is
/// shorter than one primary-plus-halo window.
fn chunk_target(
    parent_index: usize,
    target: &Sequence,
    cfg: &ChunkConfig,
) -> Vec<TargetChunk> {
    let t_len = target.seq.len() as u32;
    if cfg.primary_bp == 0 || t_len <= cfg.primary_bp + cfg.halo_bp {
        return vec![TargetChunk {
            parent_index,
            global_start: 0,
            seq: target.seq.clone(),
            primary_start_local: 0,
            primary_end_local: t_len,
        }];
    }

    let mut chunks = Vec::new();
    let mut primary_start = 0u32;
    while primary_start < t_len {
        let primary_end = (primary_start + cfg.primary_bp).min(t_len);
        let chunk_lo = primary_start.saturating_sub(cfg.halo_bp);
        let chunk_hi = (primary_end + cfg.halo_bp).min(t_len);
        let seq = target.seq.slice_to_new(chunk_lo as usize, chunk_hi as usize);
        let primary_start_local = primary_start - chunk_lo;
        let primary_end_local = primary_end - chunk_lo;
        chunks.push(TargetChunk {
            parent_index,
            global_start: chunk_lo,
            seq,
            primary_start_local,
            primary_end_local,
        });
        if primary_end == t_len {
            break;
        }
        primary_start = primary_end;
    }
    chunks
}

type PerGroup = (usize, usize, Strand, Vec<Record>);

/// Run the full pipeline and return all records.
pub fn run(targets: &[Sequence], queries: &[Sequence], config: &Config) -> Vec<Record> {
    // Build the chunk list up front so we can parallelise over it with
    // rayon. When chunking is disabled this is a single chunk per target.
    let chunks: Vec<TargetChunk> = targets
        .iter()
        .enumerate()
        .flat_map(|(ti, t)| chunk_target(ti, t, &config.chunk))
        .collect();

    let per_group: Mutex<Vec<PerGroup>> = Mutex::new(Vec::new());

    chunks.par_iter().for_each(|chunk| {
        let target = &targets[chunk.parent_index];
        let mut table = PosTable::build(&chunk.seq, &config.pattern, config.step);
        if config.max_word_count > 0 {
            table.prune_hot_words(config.max_word_count);
        }
        let chunk_ascii = chunk.seq.to_ascii();
        let target_len = target.seq.len() as u32;

        for (qi, query) in queries.iter().enumerate() {
            let plus_ascii = query.seq.to_ascii();
            let rc_cache = query.seq.reverse_complement();
            let rc_ascii = rc_cache.to_ascii();

            for &strand in config.strand.strands() {
                let (qseq, qascii) = match strand {
                    Strand::Plus => (&query.seq, &plus_ascii),
                    Strand::Minus => (&rc_cache, &rc_ascii),
                };
                let hsps = search(
                    &table,
                    &chunk.seq,
                    qseq,
                    &config.matrix,
                    &SearchParams {
                        step: config.step,
                        hsp: config.hsp,
                        transitions: config.transitions,
                        entropy_threshold: config.entropy_threshold,
                    },
                );
                if hsps.is_empty() {
                    continue;
                }

                // Drop HSPs whose anchor (chunk-local `t_start`) lives in
                // a neighbouring chunk's primary region. The surviving
                // chunk will still emit this HSP from its own primary-
                // region anchor, so dropping here avoids double-counting
                // without losing alignments.
                let mut kept: Vec<Hsp> = hsps
                    .into_iter()
                    .filter(|h| {
                        h.t_start >= chunk.primary_start_local
                            && h.t_start < chunk.primary_end_local
                    })
                    .collect();

                if config.chain_enabled {
                    kept = best_chain(&kept);
                }

                let mut recs: Vec<Record> = Vec::with_capacity(kept.len());
                for hsp in kept {
                    let record = if config.gapped_enabled {
                        build_gapped_record(
                            hsp,
                            &chunk.seq,
                            qseq,
                            &chunk_ascii,
                            qascii,
                            &config.matrix,
                            &config.gapped,
                            config.anchor_window,
                            target,
                            query,
                            strand,
                            chunk.global_start,
                            target_len,
                        )
                    } else {
                        Some(build_ungapped_record(
                            hsp,
                            &chunk_ascii,
                            qascii,
                            target,
                            query,
                            strand,
                            chunk.global_start,
                            target_len,
                        ))
                    };
                    if let Some(r) = record {
                        recs.push(r);
                    }
                }

                let recs = if let Some(tween_cfg) = &config.tweener {
                    interpolate(
                        recs,
                        target,
                        // tweener re-scans inter-chain gaps against the
                        // full target / query (not the chunk), because
                        // gaps between chain members can cross chunk
                        // boundaries. Pass full sequences + ASCII.
                        &target.seq.to_ascii(),
                        query,
                        qseq,
                        qascii,
                        strand,
                        &config.matrix,
                        config,
                        tween_cfg,
                    )
                } else {
                    recs
                };

                if !recs.is_empty() {
                    per_group
                        .lock()
                        .unwrap()
                        .push((chunk.parent_index, qi, strand, recs));
                }
            }
        }
    });

    let mut grouped = per_group.into_inner().unwrap();
    grouped.sort_by_key(|(ti, qi, strand, _)| (*ti, *qi, *strand as u8));

    let mut out = Vec::new();
    for (_, _, _, mut recs) in grouped {
        recs.sort_by_key(|r| (r.t_start, r.q_start, r.t_span, r.q_span));
        // Different ungapped HSPs on the same diagonal can extend to the
        // same gapped alignment — happens especially with transitions
        // enabled, since one homologous region produces many seed hits.
        // Drop such duplicates based on span + coordinates. Scores
        // necessarily match for duplicates since the DP is deterministic.
        recs.dedup_by(|a, b| {
            a.t_start == b.t_start
                && a.q_start == b.q_start
                && a.t_span == b.t_span
                && a.q_span == b.q_span
                && a.query_strand == b.query_strand
        });
        out.extend(recs);
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn build_gapped_record(
    hsp: Hsp,
    target_seq: &PackedSeq,
    query_seq: &PackedSeq,
    target_ascii: &[u8],
    query_ascii: &[u8],
    matrix: &ScoringMatrix,
    params: &GappedParams,
    anchor_window: u32,
    target: &Sequence,
    query: &Sequence,
    strand: Strand,
    // Chunk-local → global offset on the target axis; 0 when chunking
    // disabled. `target_len` is the parent target length (distinct from
    // `target_seq.len()` when target_seq is a chunk).
    chunk_global_start: u32,
    target_len: u32,
) -> Option<Record> {
    let anchor = choose_anchor(hsp, target_seq, query_seq, matrix, anchor_window);
    let aln = gapped_extend(anchor, target_seq, query_seq, matrix, params)?;

    let t_bases = target_ascii[aln.t_start as usize..aln.t_end() as usize].to_vec();
    let q_bases = query_ascii[aln.q_start as usize..aln.q_end() as usize].to_vec();

    Some(Record {
        target_name: target.name.clone(),
        target_len,
        query_name: query.name.clone(),
        query_len: query.seq.len() as u32,
        query_strand: strand,
        t_start: aln.t_start + chunk_global_start,
        q_start: aln.q_start,
        t_span: aln.t_len,
        q_span: aln.q_len,
        score: aln.score,
        script: aln.script,
        target_bases: t_bases,
        query_bases: q_bases,
    })
}

#[allow(clippy::too_many_arguments)]
fn build_ungapped_record(
    hsp: Hsp,
    target_ascii: &[u8],
    query_ascii: &[u8],
    target: &Sequence,
    query: &Sequence,
    strand: Strand,
    chunk_global_start: u32,
    target_len: u32,
) -> Record {
    let t0 = hsp.t_start as usize;
    let t1 = t0 + hsp.length as usize;
    let q0 = hsp.q_start as usize;
    let q1 = q0 + hsp.length as usize;
    let mut script = EditScript::new();
    script.push(EditOp::Match, hsp.length);
    Record {
        target_name: target.name.clone(),
        target_len,
        query_name: query.name.clone(),
        query_len: query.seq.len() as u32,
        query_strand: strand,
        t_start: hsp.t_start + chunk_global_start,
        q_start: hsp.q_start,
        t_span: hsp.length,
        q_span: hsp.length,
        score: hsp.score,
        script,
        target_bases: target_ascii[t0..t1].to_vec(),
        query_bases: query_ascii[q0..q1].to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequences::PackedSeq;

    fn seq(name: &str, bases: &[u8]) -> Sequence {
        Sequence { name: name.into(), seq: PackedSeq::from_ascii(bases) }
    }

    #[test]
    fn gapped_pipeline_recovers_single_indel() {
        // Target: GATTACACATGGCATGTCGA  (20 bp)
        // Query:  GATTAC_ACATGGCATGTCGA has a deletion relative to target.
        // Easier: query with one extra base.
        let targets = vec![seq(
            "t",
            b"AAAAAAAAAAGATTACACATGGCATGTCGAAAAAAAAAA",
        )];
        let queries = vec![seq("q", b"GATTACAACATGGCATGTCGA")]; // extra A at position 6
        let cfg = Config {
            pattern: SeedPattern::solid(8),
            strand: StrandSpec::Plus,
            hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
            gapped: GappedParams { y_drop: 3_000, gapped_threshold: 500 },
            gapped_enabled: true,
            chain_enabled: false,
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        assert!(!recs.is_empty(), "no records");
        // Find the record that covers the full query
        let r = recs.iter().max_by_key(|r| r.q_span).unwrap();
        assert_eq!(r.q_span, 21);
        assert_eq!(r.t_span, 20);
        let cigar = r.script.to_cigar();
        assert!(cigar.contains("1I"), "cigar was {cigar}");
    }

    #[test]
    fn nogapped_mode_emits_single_match_block() {
        let targets = vec![seq(
            "t",
            b"AAAAAAAAAAAGATTACACATGGCATGTCGAAAAAAAAAAA",
        )];
        let queries = vec![seq("q", b"GATTACACATGGCATGTCGA")];
        let cfg = Config {
            pattern: SeedPattern::solid(8),
            strand: StrandSpec::Plus,
            hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
            gapped_enabled: false,
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        assert_eq!(recs.len(), 1);
        assert_eq!(recs[0].script.to_cigar(), "20M");
    }

    #[test]
    fn chain_enabled_drops_conflicting_hsp() {
        // Two planted homology blocks separated by a stretch of Ns long
        // enough that x-drop (910) at ambig_score (-100 per N) can't
        // bridge them. 12 Ns = -1200 penalty, well past the threshold.
        let t_bytes: Vec<u8> = b"AAAAAAAAAAA"
            .iter()
            .chain(b"GATTACACATGGCATG".iter())
            .chain(b"AAAAAAAAAAAA".iter()) // 12-bp low-complexity separator
            .chain(b"CTGACCGAATGCATCA".iter())
            .chain(b"AAAAAAAAAAA".iter())
            .copied()
            .collect();
        let queries = vec![seq(
            "q",
            b"GATTACACATGGCATGNNNNNNNNNNNNCTGACCGAATGCATCA", // 12 Ns
        )];
        let targets = vec![seq("t", &t_bytes)];
        let cfg = Config {
            pattern: SeedPattern::solid(10),
            strand: StrandSpec::Plus,
            hsp: HspParams { x_drop: 910, hsp_threshold: 500 },
            gapped: GappedParams { y_drop: 3_000, gapped_threshold: 500 },
            gapped_enabled: false,
            chain_enabled: true,
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        // The chain should cover both distinct HSPs.
        assert!(recs.len() >= 2, "got {recs:#?}");
        // Records are sorted by t_start, so the second one must start after
        // the first one ends (strict increase in chain).
        for pair in recs.windows(2) {
            assert!(pair[0].t_end() <= pair[1].t_start);
            assert!(pair[0].q_end() <= pair[1].q_start);
        }
    }
}
