//! Inter-alignment interpolation ("tweening").
//!
//! After a best-chain filter leaves gaps between adjacent chain members,
//! upstream `tweener.c` re-runs the seed → HSP → extend pipeline on each
//! inter-chain window with a more sensitive seed pattern (typically
//! `match12` or `match10`). Any new high-scoring alignment found in the
//! window is inserted back into the chain.
//!
//! This is a sensitivity lever, not a performance lever: it recovers
//! alignments the default seed missed. The returned records are in
//! ascending `t_start` order, merged with the input.

use crate::anchor::choose_anchor;
use crate::chain::best_chain;
use crate::driver::Config;
use crate::gapped_extend::extend as gapped_extend;
use crate::hsp::{Hsp, HspParams};
use crate::output::{Record, Strand};
use crate::pos_table::PosTable;
use crate::scoring::ScoringMatrix;
use crate::seed_search::{search, SearchParams};
use crate::seeds::SeedPattern;
use crate::sequences::{PackedSeq, Sequence};

/// Parameters for the tweener pass.
#[derive(Debug, Clone)]
pub struct TweenerConfig {
    /// Pattern used for the re-scan. Defaults to a solid match-12 when
    /// omitted — denser than the default `12of19` so more seeds hit.
    pub pattern: SeedPattern,
    /// Minimum gap size (in target coordinates) below which interpolation
    /// is skipped. Very small gaps are rarely worth the extra work.
    pub min_gap: u32,
    /// Upper bound on gap size to interpolate into. Runaway windows
    /// dominated by low-complexity repeats can otherwise blow up
    /// pos_table builds.
    pub max_gap: u32,
    /// HSP threshold relaxation factor. New HSPs must score at least
    /// `baseline_hsp_threshold / relax_factor` — matching upstream's
    /// "tweener is more lenient" behaviour.
    pub relax_factor: i32,
}

impl Default for TweenerConfig {
    fn default() -> Self {
        Self {
            pattern: SeedPattern::solid(12),
            min_gap: 100,
            max_gap: 100_000,
            relax_factor: 2,
        }
    }
}

/// Walk adjacent pairs of chain records; for each gap that satisfies the
/// tweener's bounds, re-scan the corresponding (target, query) window with
/// a denser seed and append any new high-scoring alignments. Input and
/// output records are for one `(target, query, strand)` triple.
pub fn interpolate(
    mut records: Vec<Record>,
    target: &Sequence,
    target_ascii: &[u8],
    query: &Sequence,
    query_packed: &PackedSeq,
    query_ascii: &[u8],
    strand: Strand,
    matrix: &ScoringMatrix,
    driver_cfg: &Config,
    tween_cfg: &TweenerConfig,
) -> Vec<Record> {
    if records.len() < 2 {
        return records;
    }
    records.sort_by_key(|r| (r.t_start, r.q_start));

    let relaxed_hsp_thresh =
        (driver_cfg.hsp.hsp_threshold / tween_cfg.relax_factor.max(1)).max(100);
    let hsp_params = HspParams {
        x_drop: driver_cfg.hsp.x_drop,
        hsp_threshold: relaxed_hsp_thresh,
    };    let relaxed_gapped_thresh =
        (driver_cfg.gapped.gapped_threshold / tween_cfg.relax_factor.max(1)).max(100);
    let gapped_params = crate::gapped_extend::GappedParams {
        y_drop: driver_cfg.gapped.y_drop,
        gapped_threshold: relaxed_gapped_thresh,
    };

    let mut extras: Vec<Record> = Vec::new();
    for pair in records.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        let t_gap_start = a.t_end();
        let t_gap_end = b.t_start;
        let q_gap_start = a.q_end();
        let q_gap_end = b.q_start;

        if t_gap_end <= t_gap_start || q_gap_end <= q_gap_start {
            continue;
        }
        let t_gap = t_gap_end - t_gap_start;
        let q_gap = q_gap_end - q_gap_start;
        if t_gap < tween_cfg.min_gap || q_gap < tween_cfg.min_gap {
            continue;
        }
        if t_gap > tween_cfg.max_gap || q_gap > tween_cfg.max_gap {
            continue;
        }

        let t_sub = subseq(&target.seq, t_gap_start as usize, t_gap_end as usize);
        let q_sub = subseq(query_packed, q_gap_start as usize, q_gap_end as usize);

        let mut table = PosTable::build(&t_sub, &tween_cfg.pattern, 1);
        if driver_cfg.max_word_count > 0 {
            table.prune_hot_words(driver_cfg.max_word_count);
        }
        let hsps = search(
            &table,
            &t_sub,
            &q_sub,
            matrix,
            &SearchParams {
                step: 1,
                hsp: hsp_params,
                transitions: driver_cfg.transitions,
                entropy_threshold: driver_cfg.entropy_threshold,
            },
        );
        if hsps.is_empty() {
            continue;
        }

        let kept: Vec<Hsp> = if driver_cfg.chain_enabled {
            best_chain(&hsps)
        } else {
            hsps
        };

        for local_hsp in kept {
            let anchor = choose_anchor(
                local_hsp,
                &t_sub,
                &q_sub,
                matrix,
                driver_cfg.anchor_window,
            );
            let Some(aln) = gapped_extend(anchor, &t_sub, &q_sub, matrix, &gapped_params)
            else {
                continue;
            };

            // Translate local coordinates back to global.
            let global_t_start = t_gap_start + aln.t_start;
            let global_q_start = q_gap_start + aln.q_start;
            let t_bases = target_ascii
                [global_t_start as usize..(global_t_start + aln.t_len) as usize]
                .to_vec();
            let q_bases = query_ascii
                [global_q_start as usize..(global_q_start + aln.q_len) as usize]
                .to_vec();
            extras.push(Record {
                target_name: target.name.clone(),
                target_len: target.seq.len() as u32,
                query_name: query.name.clone(),
                query_len: query.seq.len() as u32,
                query_strand: strand,
                t_start: global_t_start,
                q_start: global_q_start,
                t_span: aln.t_len,
                q_span: aln.q_len,
                score: aln.score,
                script: aln.script,
                target_bases: t_bases,
                query_bases: q_bases,
            });
        }
    }

    records.extend(extras);
    records.sort_by_key(|r| (r.t_start, r.q_start));
    // Drop exact duplicates introduced by the re-scan.
    records.dedup_by(|a, b| {
        a.t_start == b.t_start
            && a.q_start == b.q_start
            && a.t_span == b.t_span
            && a.q_span == b.q_span
    });
    records
}

fn subseq(seq: &PackedSeq, start: usize, end: usize) -> PackedSeq {
    let mut out = PackedSeq::with_capacity(end - start);
    for i in start..end {
        let byte = if seq.is_valid(i) {
            let c = crate::dna::decode_base(seq.code(i));
            if seq.is_masked(i) { c.to_ascii_lowercase() } else { c }
        } else {
            b'N'
        };
        out.push_ascii(byte);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::driver::StrandSpec;
    use crate::gapped_extend::GappedParams;

    fn matrix() -> ScoringMatrix {
        ScoringMatrix::hoxd70()
    }

    fn minimal_driver_cfg() -> Config {
        Config {
            pattern: SeedPattern::solid(19),
            matrix: matrix(),
            step: 1,
            hsp: HspParams { x_drop: 910, hsp_threshold: 3_000 },
            gapped: GappedParams { y_drop: 9_400, gapped_threshold: 3_000 },
            strand: StrandSpec::Plus,
            max_word_count: 0,
            respect_masking: true,
            gapped_enabled: true,
            chain_enabled: true,
            anchor_window: 31,
            transitions: 1,
            entropy_threshold: None,
            chunk: crate::driver::ChunkConfig::default(),
            tweener: None,
        }
    }

    #[test]
    fn no_gap_means_passthrough() {
        let recs = vec![];
        let t = Sequence { name: "t".into(), seq: PackedSeq::from_ascii(b"A") };
        let q = Sequence { name: "q".into(), seq: PackedSeq::from_ascii(b"A") };
        let out = interpolate(
            recs,
            &t,
            b"A",
            &q,
            &q.seq,
            b"A",
            Strand::Plus,
            &matrix(),
            &minimal_driver_cfg(),
            &TweenerConfig::default(),
        );
        assert!(out.is_empty());
    }

    #[test]
    fn finds_missed_middle_alignment() {
        // Scenario: original scan used match19 and missed a small 15bp
        // homologous block in the middle. Tweener's match12 should find it.
        //
        // We construct fake records at start + end of a target, then ask
        // the tweener to look for anything in the middle.
        let target_bytes: Vec<u8> = b"GATTACACATGGCATGTCGA" // 20bp — first aln
            .iter()
            .chain(b"NNNNNNNNNNNNNNNNNNNN".iter()) // 20bp spacer
            .chain(b"AGCTAGCTAGCTAGC".iter()) // 15bp missed match
            .chain(b"NNNNNNNNNNNNNNNNNNNN".iter()) // 20bp spacer
            .chain(b"GATTACACATGGCATGTCGA".iter()) // 20bp — last aln
            .copied()
            .collect();
        let query_bytes = target_bytes.clone(); // identity
        let target = Sequence {
            name: "t".into(),
            seq: PackedSeq::from_ascii(&target_bytes),
        };
        let query = Sequence {
            name: "q".into(),
            seq: PackedSeq::from_ascii(&query_bytes),
        };
        let target_ascii = target.seq.to_ascii();
        let query_ascii = query.seq.to_ascii();

        // Synthesize two chain records at target[0..20] and target[75..95].
        let make_rec = |t_start: u32, t_end: u32| -> Record {
            let mut script = crate::edit_script::EditScript::new();
            script.push(crate::edit_script::EditOp::Match, t_end - t_start);
            Record {
                target_name: target.name.clone(),
                target_len: target.seq.len() as u32,
                query_name: query.name.clone(),
                query_len: query.seq.len() as u32,
                query_strand: Strand::Plus,
                t_start,
                q_start: t_start,
                t_span: t_end - t_start,
                q_span: t_end - t_start,
                score: 1_900,
                script,
                target_bases: target_ascii[t_start as usize..t_end as usize].to_vec(),
                query_bases: query_ascii[t_start as usize..t_end as usize].to_vec(),
            }
        };
        let recs = vec![make_rec(0, 20), make_rec(75, 95)];

        let mut cfg = minimal_driver_cfg();
        // Pretend the outer scan would never have hit a 15bp block.
        cfg.hsp.hsp_threshold = 5_000;
        cfg.gapped.gapped_threshold = 5_000;

        // With relax_factor = 2, thresholds drop to 2500 each; a 15bp
        // perfect run scores ~1400 under HOXD70 so still doesn't reach the
        // relaxed bar. Force relaxation all the way down.
        let tween = TweenerConfig {
            pattern: SeedPattern::solid(12),
            min_gap: 10,
            max_gap: 10_000,
            relax_factor: 100, // aggressive for test: thresholds become 100
        };

        let out = interpolate(
            recs,
            &target,
            &target_ascii,
            &query,
            &query.seq,
            &query_ascii,
            Strand::Plus,
            &matrix(),
            &cfg,
            &tween,
        );

        // We should have at least one extra record covering the middle
        // block (inside target[20..75]).
        let middle_found = out
            .iter()
            .any(|r| r.t_start >= 20 && r.t_end() <= 75 && r.t_span >= 15);
        assert!(middle_found, "tweener missed the middle block: {out:#?}");
    }

    #[test]
    fn respects_max_gap() {
        // Make the gap huge; tweener should skip it entirely.
        let t = Sequence { name: "t".into(), seq: PackedSeq::from_ascii(b"A") };
        let q = Sequence { name: "q".into(), seq: PackedSeq::from_ascii(b"A") };

        let r1 = Record {
            target_name: "t".into(),
            target_len: 100,
            query_name: "q".into(),
            query_len: 100,
            query_strand: Strand::Plus,
            t_start: 0,
            q_start: 0,
            t_span: 10,
            q_span: 10,
            score: 0,
            script: crate::edit_script::EditScript::new(),
            target_bases: vec![],
            query_bases: vec![],
        };
        let r2 = Record {
            t_start: 99,
            q_start: 99,
            t_span: 1,
            q_span: 1,
            ..r1.clone()
        };
        // max_gap = 5 means the 89bp gap is skipped.
        let tween = TweenerConfig { max_gap: 5, ..TweenerConfig::default() };
        let out = interpolate(
            vec![r1, r2],
            &t,
            b"",
            &q,
            &q.seq,
            b"",
            Strand::Plus,
            &matrix(),
            &minimal_driver_cfg(),
            &tween,
        );
        // Just the two originals; no new interpolation.
        assert_eq!(out.len(), 2);
    }
}
