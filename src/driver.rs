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
use crate::sequences::Sequence;

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
    /// Run gapped affine extension after HSP. When false, the pipeline
    /// stops at HSPs (equivalent to upstream `--nogapped`).
    pub gapped_enabled: bool,
    /// Run chaining between HSPs and gapped extension. When false,
    /// every HSP is extended independently (equivalent to `--nochain`).
    pub chain_enabled: bool,
    /// Sliding window width for anchor selection inside each HSP.
    pub anchor_window: u32,
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
            gapped_enabled: true,
            chain_enabled: false,
            anchor_window: 31,
        }
    }
}

type PerGroup = (usize, usize, Strand, Vec<Record>);

/// Run the full pipeline and return all records.
pub fn run(targets: &[Sequence], queries: &[Sequence], config: &Config) -> Vec<Record> {
    let per_group: Mutex<Vec<PerGroup>> = Mutex::new(Vec::new());

    targets.par_iter().enumerate().for_each(|(ti, target)| {
        let mut table = PosTable::build(&target.seq, &config.pattern, config.step);
        if config.max_word_count > 0 {
            table.prune_hot_words(config.max_word_count);
        }
        let target_ascii = target.seq.to_ascii();

        for (qi, query) in queries.iter().enumerate() {
            let rc_cache = query.seq.reverse_complement();
            let rc_ascii = rc_cache.to_ascii();

            for &strand in config.strand.strands() {
                let (qseq, qascii) = match strand {
                    Strand::Plus => (&query.seq, &query.seq.to_ascii()),
                    Strand::Minus => (&rc_cache, &rc_ascii),
                };
                let hsps = search(
                    &table,
                    &target.seq,
                    qseq,
                    &config.matrix,
                    &SearchParams { step: config.step, hsp: config.hsp },
                );
                if hsps.is_empty() {
                    continue;
                }

                let kept = if config.chain_enabled {
                    best_chain(&hsps)
                } else {
                    hsps
                };

                let mut recs: Vec<Record> = Vec::with_capacity(kept.len());
                for hsp in kept {
                    let record = if config.gapped_enabled {
                        build_gapped_record(
                            hsp,
                            &target.seq,
                            qseq,
                            &target_ascii,
                            qascii,
                            &config.matrix,
                            &config.gapped,
                            config.anchor_window,
                            target,
                            query,
                            strand,
                        )
                    } else {
                        Some(build_ungapped_record(
                            hsp,
                            &target_ascii,
                            qascii,
                            target,
                            query,
                            strand,
                        ))
                    };
                    if let Some(r) = record {
                        recs.push(r);
                    }
                }

                if !recs.is_empty() {
                    per_group.lock().unwrap().push((ti, qi, strand, recs));
                }
            }
        }
    });

    let mut grouped = per_group.into_inner().unwrap();
    grouped.sort_by_key(|(ti, qi, strand, _)| (*ti, *qi, *strand as u8));

    let mut out = Vec::new();
    for (_, _, _, mut recs) in grouped {
        recs.sort_by_key(|r| (r.t_start, r.q_start));
        out.extend(recs);
    }
    out
}

#[allow(clippy::too_many_arguments)]
fn build_gapped_record(
    hsp: Hsp,
    target_seq: &crate::sequences::PackedSeq,
    query_seq: &crate::sequences::PackedSeq,
    target_ascii: &[u8],
    query_ascii: &[u8],
    matrix: &ScoringMatrix,
    params: &GappedParams,
    anchor_window: u32,
    target: &Sequence,
    query: &Sequence,
    strand: Strand,
) -> Option<Record> {
    let anchor = choose_anchor(hsp, target_seq, query_seq, matrix, anchor_window);
    let aln = gapped_extend(anchor, target_seq, query_seq, matrix, params)?;

    let t_bases = target_ascii[aln.t_start as usize..aln.t_end() as usize].to_vec();
    let q_bases = query_ascii[aln.q_start as usize..aln.q_end() as usize].to_vec();

    Some(Record {
        target_name: target.name.clone(),
        target_len: target.seq.len() as u32,
        query_name: query.name.clone(),
        query_len: query.seq.len() as u32,
        query_strand: strand,
        t_start: aln.t_start,
        q_start: aln.q_start,
        t_span: aln.t_len,
        q_span: aln.q_len,
        score: aln.score,
        script: aln.script,
        target_bases: t_bases,
        query_bases: q_bases,
    })
}

fn build_ungapped_record(
    hsp: Hsp,
    target_ascii: &[u8],
    query_ascii: &[u8],
    target: &Sequence,
    query: &Sequence,
    strand: Strand,
) -> Record {
    let t0 = hsp.t_start as usize;
    let t1 = t0 + hsp.length as usize;
    let q0 = hsp.q_start as usize;
    let q1 = q0 + hsp.length as usize;
    let mut script = EditScript::new();
    script.push(EditOp::Match, hsp.length);
    Record {
        target_name: target.name.clone(),
        target_len: target.seq.len() as u32,
        query_name: query.name.clone(),
        query_len: query.seq.len() as u32,
        query_strand: strand,
        t_start: hsp.t_start,
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
        // Two non-overlapping HSPs plus a third that conflicts on query.
        // With --chain the third should be dropped.
        let t_bytes: Vec<u8> = b"AAAAA"
            .iter()
            .chain(b"GATTACACATGGCATG".iter())
            .chain(b"AAAAAAAA".iter())
            .chain(b"CTGACCGAATGCATCA".iter())
            .chain(b"AAAAA".iter())
            .copied()
            .collect();
        let targets = vec![seq("t", &t_bytes)];
        let queries = vec![seq(
            "q",
            b"GATTACACATGGCATGNNNNNNNNCTGACCGAATGCATCA",
        )];
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
        assert!(recs.len() >= 2);
        // Records are sorted by t_start, so the second one must start after
        // the first one ends (strict increase in chain).
        for pair in recs.windows(2) {
            assert!(pair[0].t_end() <= pair[1].t_start);
            assert!(pair[0].q_end() <= pair[1].q_start);
        }
    }
}
