//! End-to-end alignment driver.
//!
//! Builds one `PosTable` per target sequence, then fans out query × strand
//! combinations across a `rayon` thread pool. Per the PLAN.md §3.2 concurrency
//! model, the `PosTable` and scoring matrix are shared immutably (`Arc`-free
//! because rayon's `par_iter` only needs a `&`); each worker keeps its own
//! `DiagHash` inside `seed_search::search`.

use std::sync::Mutex;

use rayon::prelude::*;

use crate::hsp::HspParams;
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

/// Public configuration for `run()`. Mirrors the common lastz flag set.
#[derive(Debug, Clone)]
pub struct Config {
    pub pattern: SeedPattern,
    pub matrix: ScoringMatrix,
    pub step: usize,
    pub hsp: HspParams,
    pub strand: StrandSpec,
    /// Maximum seed-word multiplicity before a word is dropped as "hot". 0
    /// disables the filter.
    pub max_word_count: u32,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            pattern: SeedPattern::twelve_of_nineteen(),
            matrix: ScoringMatrix::hoxd70(),
            step: 1,
            hsp: HspParams::default(),
            strand: StrandSpec::Both,
            max_word_count: 0,
        }
    }
}

/// Align every `(target, query, strand)` triple and return the ungapped
/// alignment records. Output order is stable: records are grouped by target,
/// then query, then strand, then by the internal HSP sort order from
/// `seed_search::search`.
pub fn run(
    targets: &[Sequence],
    queries: &[Sequence],
    config: &Config,
) -> Vec<Record> {
    // Per-target indexed build. Targets alignments are independent so we
    // parallelize the outer loop.
    let all: Mutex<Vec<(usize, usize, Strand, Vec<crate::hsp::Hsp>, bool)>> =
        Mutex::new(Vec::new());

    targets.par_iter().enumerate().for_each(|(ti, target)| {
        let mut table = PosTable::build(&target.seq, &config.pattern, config.step);
        if config.max_word_count > 0 {
            table.prune_hot_words(config.max_word_count);
        }

        for (qi, query) in queries.iter().enumerate() {
            let rc_cache = query.seq.reverse_complement();

            for &strand in config.strand.strands() {
                let qseq = match strand {
                    Strand::Plus => &query.seq,
                    Strand::Minus => &rc_cache,
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
                let is_rc = matches!(strand, Strand::Minus);
                all.lock().unwrap().push((ti, qi, strand, hsps, is_rc));
            }
        }
    });

    // Assemble Records in a stable order.
    let mut grouped = all.into_inner().unwrap();
    grouped.sort_by_key(|&(ti, qi, strand, _, _)| (ti, qi, strand as u8));

    let mut records = Vec::new();
    for (ti, qi, strand, hsps, is_rc) in grouped {
        let target = &targets[ti];
        let query = &queries[qi];
        let q_rc = if is_rc { Some(query.seq.reverse_complement()) } else { None };
        let qseq = match strand {
            Strand::Plus => &query.seq,
            Strand::Minus => q_rc.as_ref().unwrap(),
        };
        let target_ascii = target.seq.to_ascii();
        let query_ascii = qseq.to_ascii();

        for hsp in hsps {
            let t0 = hsp.t_start as usize;
            let t1 = t0 + hsp.length as usize;
            let q0 = hsp.q_start as usize;
            let q1 = q0 + hsp.length as usize;
            records.push(Record {
                target_name: target.name.clone(),
                target_len: target.seq.len() as u32,
                query_name: query.name.clone(),
                query_len: query.seq.len() as u32,
                query_strand: strand,
                hsp,
                target_bases: target_ascii[t0..t1].to_vec(),
                query_bases: query_ascii[q0..q1].to_vec(),
            });
        }
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn end_to_end_plus_strand_finds_embedded_match() {
        let targets = vec![Sequence {
            name: "t".into(),
            seq: crate::sequences::PackedSeq::from_ascii(
                b"AAAAAAAAAAAGATTACACATGGCATGTCGAAAAAAAAAAA",
            ),
        }];
        let queries = vec![Sequence {
            name: "q".into(),
            seq: crate::sequences::PackedSeq::from_ascii(b"GATTACACATGGCATGTCGA"),
        }];
        let cfg = Config {
            pattern: SeedPattern::solid(8),
            strand: StrandSpec::Plus,
            hsp: HspParams { x_drop: 500, hsp_threshold: 500 },
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        assert_eq!(recs.len(), 1, "got {recs:#?}");
        assert_eq!(recs[0].hsp.t_start, 11);
        assert_eq!(recs[0].hsp.q_start, 0);
        assert_eq!(recs[0].hsp.length, 20);
        assert_eq!(recs[0].query_strand, Strand::Plus);
    }

    #[test]
    fn finds_reverse_complement_match_on_minus_strand() {
        // Target embeds the reverse complement of the query; the plus-strand
        // scan finds nothing, but the minus-strand scan hits.
        let query_ascii: &[u8] = b"GATTACACATGGCATGTCGA";
        let rc_of_query = crate::sequences::PackedSeq::from_ascii(query_ascii)
            .reverse_complement()
            .to_ascii();
        let mut target_bytes: Vec<u8> = b"AAAAAAAAAAA".to_vec();
        target_bytes.extend_from_slice(&rc_of_query);
        target_bytes.extend_from_slice(b"AAAAAAAAAAA");
        let targets = vec![Sequence {
            name: "t".into(),
            seq: crate::sequences::PackedSeq::from_ascii(&target_bytes),
        }];
        let queries = vec![Sequence {
            name: "q".into(),
            seq: crate::sequences::PackedSeq::from_ascii(query_ascii),
        }];
        let cfg = Config {
            pattern: SeedPattern::solid(8),
            strand: StrandSpec::Minus,
            hsp: HspParams { x_drop: 500, hsp_threshold: 500 },
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        assert!(!recs.is_empty(), "expected an RC hit, got {recs:#?}");
        assert!(recs.iter().all(|r| r.query_strand == Strand::Minus));
    }

    #[test]
    fn both_strands_scans_both_orientations() {
        // Plus-strand match in one region; the minus-strand scan sees the
        // same target bases reversed/complemented and finds nothing
        // (non-palindromic), so we expect exactly plus-strand records.
        let targets = vec![Sequence {
            name: "t".into(),
            seq: crate::sequences::PackedSeq::from_ascii(
                b"AAAAAAAAAAAGATTACACATGGCATGTCGAAAAAAAAAAA",
            ),
        }];
        let queries = vec![Sequence {
            name: "q".into(),
            seq: crate::sequences::PackedSeq::from_ascii(b"GATTACACATGGCATGTCGA"),
        }];
        let cfg = Config {
            pattern: SeedPattern::solid(8),
            strand: StrandSpec::Both,
            hsp: HspParams { x_drop: 500, hsp_threshold: 500 },
            ..Config::default()
        };
        let recs = run(&targets, &queries, &cfg);
        assert!(recs.iter().any(|r| r.query_strand == Strand::Plus));
        assert!(!recs.iter().any(|r| r.query_strand == Strand::Minus));
    }
}
