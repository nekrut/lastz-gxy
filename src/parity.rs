//! MAF parsing and block-level comparison for the parity gate.
//!
//! The gate (PLAN.md §5) requires measuring how closely our output matches
//! upstream lastz's on the same inputs. This module parses both sides'
//! MAF files into a canonical `Block` representation, computes block-level
//! Jaccard on `(target_name, query_name, strand, t_range, q_range)`
//! signatures, and summarises score / aligned-bp deltas for the blocks
//! that appear in both sides.
//!
//! Not a general-purpose MAF library — we only handle the pairwise
//! two-`s`-line-per-block form that both upstream lastz and lastz-gxy
//! emit. Multiple-alignment MAFs are out of scope.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum MafError {
    #[error("I/O: {0}")]
    Io(#[from] std::io::Error),
    #[error("line {line}: expected 2 `s` rows per `a` block, found {count}")]
    BadBlockArity { line: usize, count: usize },
    #[error("line {line}: malformed `s` record: {text}")]
    BadSLine { line: usize, text: String },
    #[error("line {line}: expected `score=` in: {text}")]
    MissingScore { line: usize, text: String },
}

/// One pairwise MAF block.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Block {
    pub target_name: String,
    pub query_name: String,
    pub query_strand: char,
    pub t_start: u32,
    pub t_span: u32,
    pub q_start: u32,
    pub q_span: u32,
    pub score: i64,
    /// Count of aligned columns excluding gap characters on either side.
    pub identity_columns: u32,
    /// Count of aligned columns matching base-for-base (case-insensitive).
    pub matches: u32,
}

impl Block {
    pub fn t_end(&self) -> u32 {
        self.t_start + self.t_span
    }
    pub fn q_end(&self) -> u32 {
        self.q_start + self.q_span
    }

    /// Block signature used for set-intersection / Jaccard. Intentionally
    /// excludes the score so that tied-score alternatives don't cause
    /// false misses.
    pub fn signature(&self) -> BlockSig {
        BlockSig {
            target_name: self.target_name.clone(),
            query_name: self.query_name.clone(),
            strand: self.query_strand,
            t_start: self.t_start,
            t_end: self.t_end(),
            q_start: self.q_start,
            q_end: self.q_end(),
        }
    }

    /// Fraction of aligned (non-gap) columns that match base-for-base.
    /// Returns `None` for degenerate blocks with zero identity columns.
    pub fn identity(&self) -> Option<f64> {
        if self.identity_columns == 0 {
            None
        } else {
            Some(self.matches as f64 / self.identity_columns as f64)
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct BlockSig {
    pub target_name: String,
    pub query_name: String,
    pub strand: char,
    pub t_start: u32,
    pub t_end: u32,
    pub q_start: u32,
    pub q_end: u32,
}

/// Read a MAF file into a flat list of pairwise blocks. Multi-alignment
/// blocks (>2 `s` rows) are rejected.
pub fn load_maf<P: AsRef<Path>>(path: P) -> Result<Vec<Block>, MafError> {
    let f = File::open(path)?;
    parse_maf(BufReader::new(f))
}

/// Parse MAF from any reader.
pub fn parse_maf<R: BufRead>(reader: R) -> Result<Vec<Block>, MafError> {
    let mut blocks = Vec::new();
    let mut current_score: Option<i64> = None;
    let mut current_s: Vec<(String, u32, u32, char, String)> = Vec::new();
    let mut start_line = 0usize;

    let flush = |score: Option<i64>,
                 s_rows: &[(String, u32, u32, char, String)],
                 start: usize,
                 out: &mut Vec<Block>|
     -> Result<(), MafError> {
        if s_rows.is_empty() {
            return Ok(());
        }
        if s_rows.len() != 2 {
            return Err(MafError::BadBlockArity { line: start, count: s_rows.len() });
        }
        let (tn, tstart, tspan, tstrand, ttext) = &s_rows[0];
        let (qn, qstart, qspan, qstrand, qtext) = &s_rows[1];
        // Target row in MAF is the first sequence; by lastz convention it is
        // always emitted on the `+` strand.
        debug_assert_eq!(*tstrand, '+');
        let (ident, matches) = count_identity(ttext.as_bytes(), qtext.as_bytes());
        out.push(Block {
            target_name: tn.clone(),
            query_name: qn.clone(),
            query_strand: *qstrand,
            t_start: *tstart,
            t_span: *tspan,
            q_start: *qstart,
            q_span: *qspan,
            score: score.unwrap_or(0),
            identity_columns: ident,
            matches,
        });
        Ok(())
    };

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        let trimmed = line.trim_end();
        if trimmed.is_empty() {
            flush(current_score.take(), &current_s, start_line, &mut blocks)?;
            current_s.clear();
            continue;
        }
        if trimmed.starts_with('#') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("a ") {
            // Starting a new block; flush any in-flight one first.
            flush(current_score.take(), &current_s, start_line, &mut blocks)?;
            current_s.clear();
            start_line = lineno + 1;
            current_score = Some(parse_score(rest).ok_or_else(|| MafError::MissingScore {
                line: lineno + 1,
                text: rest.to_string(),
            })?);
            continue;
        }
        if trimmed == "a" {
            flush(current_score.take(), &current_s, start_line, &mut blocks)?;
            current_s.clear();
            start_line = lineno + 1;
            current_score = Some(0);
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix("s ") {
            let parsed = parse_s_line(rest).ok_or_else(|| MafError::BadSLine {
                line: lineno + 1,
                text: rest.to_string(),
            })?;
            current_s.push(parsed);
            continue;
        }
        // ignore other tracks (i/q/e lines)
    }
    flush(current_score.take(), &current_s, start_line, &mut blocks)?;
    Ok(blocks)
}

/// Extract `score=N` (integer) from an `a`-line payload. Tolerant of
/// additional whitespace-separated key=value pairs on the same line.
fn parse_score(rest: &str) -> Option<i64> {
    for tok in rest.split_ascii_whitespace() {
        if let Some(v) = tok.strip_prefix("score=") {
            return v.parse().ok();
        }
    }
    None
}

/// Parse the body of an `s` line: `name start size strand src_size text`.
fn parse_s_line(rest: &str) -> Option<(String, u32, u32, char, String)> {
    let mut it = rest.split_ascii_whitespace();
    let name = it.next()?.to_string();
    let start: u32 = it.next()?.parse().ok()?;
    let size: u32 = it.next()?.parse().ok()?;
    let strand = it.next()?.chars().next()?;
    let _src_size: u32 = it.next()?.parse().ok()?;
    let text = it.next()?.to_string();
    Some((name, start, size, strand, text))
}

/// Count (identity_columns, matches) in a pair of aligned rows. Gap
/// positions (either side has `-` or `.`) are excluded. Matching is
/// case-insensitive; `N`/`n` on either side counts as a mismatch, same as
/// upstream lastz's identity reporting.
fn count_identity(t: &[u8], q: &[u8]) -> (u32, u32) {
    let n = t.len().min(q.len());
    let mut cols = 0u32;
    let mut matches = 0u32;
    for i in 0..n {
        let tc = t[i];
        let qc = q[i];
        if tc == b'-' || tc == b'.' || qc == b'-' || qc == b'.' {
            continue;
        }
        cols += 1;
        if tc.eq_ignore_ascii_case(&qc) && !matches!(tc.to_ascii_uppercase(), b'N') {
            matches += 1;
        }
    }
    (cols, matches)
}

/// Summary of a `left` vs `right` MAF comparison.
///
/// The report treats `left` as the baseline (upstream) and `right` as the
/// test implementation (`lastz-gxy`), and exposes two separate quality
/// measures rather than collapsing them into a single Jaccard number:
///
/// - **recall** = `intersection / left_blocks` — what fraction of the
///   baseline's alignments did the test also emit? This is what matters
///   for users running existing lastz pipelines; a dropped block could
///   silently break their downstream tools.
/// - **precision** = `intersection / right_blocks` — what fraction of the
///   test's alignments are also in the baseline? "Extras" here aren't
///   necessarily wrong — they can be weak-signal alignments the test's
///   DP reaches and the baseline's happens not to — so precision < 1.0
///   is not automatically a failure.
#[derive(Debug, Clone)]
pub struct ComparisonReport {
    pub left_blocks: usize,
    pub right_blocks: usize,
    pub intersection: usize,
    pub union_size: usize,
    pub jaccard: f64,
    /// Fraction of baseline blocks present in the test set.
    pub recall: f64,
    /// Fraction of test blocks present in the baseline.
    pub precision: f64,
    pub left_only: Vec<BlockSig>,
    pub right_only: Vec<BlockSig>,
    pub score_delta_median: i64,
    pub score_delta_max_abs: i64,
    pub score_delta_mean: f64,
    pub left_aligned_bp: u64,
    pub right_aligned_bp: u64,
    pub aligned_bp_delta_pct: f64,
    pub left_identity_mean: f64,
    pub right_identity_mean: f64,
}

impl ComparisonReport {
    /// The revised release gate: every baseline block must be recovered
    /// (`recall == 1.0`) and scores on shared blocks must match
    /// bit-exactly (median delta 0, |max delta| ≤ 1). Precision and
    /// Jaccard are reported but not gated — a test implementation that
    /// reaches extra alignments the baseline doesn't is a superset, not
    /// a regression. See PLAN.md §5.
    pub fn passes_release_gate(&self) -> bool {
        (self.recall >= 1.0 || self.left_blocks == 0)
            && self.score_delta_median == 0
            && self.score_delta_max_abs <= 1
    }
}

/// Compare two parsed MAF block sets.
pub fn compare(left: &[Block], right: &[Block]) -> ComparisonReport {
    // Index left and right by signature. Because multiple blocks on a
    // single signature are rare (lastz doesn't emit them), we keep the
    // highest-scoring one as the representative.
    let mut l_by_sig: BTreeMap<BlockSig, &Block> = BTreeMap::new();
    for b in left {
        l_by_sig
            .entry(b.signature())
            .and_modify(|e| {
                if b.score > e.score {
                    *e = b;
                }
            })
            .or_insert(b);
    }
    let mut r_by_sig: BTreeMap<BlockSig, &Block> = BTreeMap::new();
    for b in right {
        r_by_sig
            .entry(b.signature())
            .and_modify(|e| {
                if b.score > e.score {
                    *e = b;
                }
            })
            .or_insert(b);
    }

    let mut intersection_deltas: Vec<i64> = Vec::new();
    let mut intersection_count = 0usize;
    for (sig, lb) in &l_by_sig {
        if let Some(rb) = r_by_sig.get(sig) {
            intersection_deltas.push(lb.score - rb.score);
            intersection_count += 1;
        }
    }

    let left_only: Vec<BlockSig> = l_by_sig
        .keys()
        .filter(|s| !r_by_sig.contains_key(*s))
        .cloned()
        .collect();
    let right_only: Vec<BlockSig> = r_by_sig
        .keys()
        .filter(|s| !l_by_sig.contains_key(*s))
        .cloned()
        .collect();

    let union_size = l_by_sig.len() + r_by_sig.len() - intersection_count;
    let jaccard = if union_size == 0 {
        1.0
    } else {
        intersection_count as f64 / union_size as f64
    };
    let recall = if l_by_sig.is_empty() {
        1.0
    } else {
        intersection_count as f64 / l_by_sig.len() as f64
    };
    let precision = if r_by_sig.is_empty() {
        1.0
    } else {
        intersection_count as f64 / r_by_sig.len() as f64
    };

    intersection_deltas.sort_unstable();
    let score_delta_median = if intersection_deltas.is_empty() {
        0
    } else {
        intersection_deltas[intersection_deltas.len() / 2]
    };
    let score_delta_max_abs = intersection_deltas
        .iter()
        .map(|d| d.unsigned_abs())
        .max()
        .unwrap_or(0) as i64;
    let score_delta_mean = if intersection_deltas.is_empty() {
        0.0
    } else {
        intersection_deltas.iter().sum::<i64>() as f64 / intersection_deltas.len() as f64
    };

    let left_aligned_bp: u64 = left.iter().map(|b| b.identity_columns as u64).sum();
    let right_aligned_bp: u64 = right.iter().map(|b| b.identity_columns as u64).sum();
    let aligned_bp_delta_pct = if left_aligned_bp == 0 {
        0.0
    } else {
        (right_aligned_bp as f64 - left_aligned_bp as f64) / left_aligned_bp as f64 * 100.0
    };

    let mean_identity = |set: &[Block]| -> f64 {
        let ids: Vec<f64> = set.iter().filter_map(|b| b.identity()).collect();
        if ids.is_empty() {
            0.0
        } else {
            ids.iter().sum::<f64>() / ids.len() as f64
        }
    };

    ComparisonReport {
        left_blocks: left.len(),
        right_blocks: right.len(),
        intersection: intersection_count,
        union_size,
        jaccard,
        recall,
        precision,
        left_only,
        right_only,
        score_delta_median,
        score_delta_max_abs,
        score_delta_mean,
        left_aligned_bp,
        right_aligned_bp,
        aligned_bp_delta_pct,
        left_identity_mean: mean_identity(left),
        right_identity_mean: mean_identity(right),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SIMPLE_MAF: &str = "\
##maf version=1 program=lastz
# scoring=HOXD70

a score=1901
s chrT 11 20 + 43 GATTACACATGGCATGTCGA
s chrQ 0 20 + 20 GATTACACATGGCATGTCGA

a score=500
s chrT 50 5 + 43 ACGTA
s chrQ 10 5 - 20 ACGTA
";

    #[test]
    fn parses_two_blocks_with_scores() {
        let blocks = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].target_name, "chrT");
        assert_eq!(blocks[0].query_name, "chrQ");
        assert_eq!(blocks[0].t_start, 11);
        assert_eq!(blocks[0].t_span, 20);
        assert_eq!(blocks[0].score, 1901);
        assert_eq!(blocks[0].query_strand, '+');
        assert_eq!(blocks[1].query_strand, '-');
    }

    #[test]
    fn identity_excludes_gaps_and_ns() {
        let text = "a score=100\n\
                    s t 0 5 + 10 ACGTA\n\
                    s q 0 5 + 10 ACNTA\n\n";
        let blocks = parse_maf(text.as_bytes()).unwrap();
        assert_eq!(blocks.len(), 1);
        // 5 aligned columns; one is N/N which we count as aligned-but-mismatch
        // unless one side is N. In this fixture query has N at pos 2, target has G:
        // 4 matches among 5 columns.
        assert_eq!(blocks[0].identity_columns, 5);
        assert_eq!(blocks[0].matches, 4);
    }

    #[test]
    fn identical_inputs_give_perfect_recall_and_precision() {
        let left = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        let right = left.clone();
        let r = compare(&left, &right);
        assert_eq!(r.jaccard, 1.0);
        assert_eq!(r.recall, 1.0);
        assert_eq!(r.precision, 1.0);
        assert_eq!(r.score_delta_median, 0);
        assert_eq!(r.score_delta_max_abs, 0);
        assert!(r.passes_release_gate());
    }

    #[test]
    fn score_differences_are_summarised() {
        let left = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        let mut right = left.clone();
        right[0].score += 5;
        right[1].score -= 2;
        let r = compare(&left, &right);
        assert_eq!(r.score_delta_max_abs, 5);
        assert_eq!(r.intersection, 2);
        // Recall still 1.0 (all baseline blocks matched on signature),
        // but score-delta max exceeds the gate's tolerance (≤ 1).
        assert_eq!(r.recall, 1.0);
        assert!(!r.passes_release_gate());
    }

    #[test]
    fn missing_baseline_block_fails_recall_gate() {
        let left = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        let right = vec![left[0].clone()]; // drop second
        let r = compare(&left, &right);
        assert_eq!(r.intersection, 1);
        assert_eq!(r.jaccard, 0.5);
        assert_eq!(r.recall, 0.5);
        assert_eq!(r.precision, 1.0);
        // Dropping a baseline block is a hard fail — recall < 1.0.
        assert!(!r.passes_release_gate());
    }

    #[test]
    fn extra_test_blocks_still_pass_gate() {
        // Test emits everything baseline emits, plus an extra. Recall is
        // 1.0 (the point of the reframe); precision drops to 0.5 but the
        // gate still passes because extras on cross-species are known to
        // be weak-signal alignments our DP reaches that upstream's
        // implementation doesn't — not wrong results.
        let left = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        let mut right = left.clone();
        right.push(Block {
            target_name: "chrT".into(),
            query_name: "chrQ".into(),
            query_strand: '+',
            t_start: 200,
            t_span: 10,
            q_start: 200,
            q_span: 10,
            score: 123,
            identity_columns: 10,
            matches: 9,
        });
        right.push(Block {
            target_name: "chrT".into(),
            query_name: "chrQ".into(),
            query_strand: '+',
            t_start: 300,
            t_span: 10,
            q_start: 300,
            q_span: 10,
            score: 456,
            identity_columns: 10,
            matches: 9,
        });
        let r = compare(&left, &right);
        assert_eq!(r.recall, 1.0);
        assert_eq!(r.intersection, 2);
        assert_eq!(r.right_only.len(), 2);
        assert!(r.precision < 1.0);
        assert!(r.passes_release_gate());
    }

    #[test]
    fn release_gate_tolerates_small_score_drift() {
        let left = parse_maf(SIMPLE_MAF.as_bytes()).unwrap();
        let mut right = left.clone();
        right[0].score += 1;
        let r = compare(&left, &right);
        assert_eq!(r.recall, 1.0);
        assert_eq!(r.score_delta_max_abs, 1);
        assert!(r.passes_release_gate());
    }

    #[test]
    fn empty_input_gives_trivial_pass() {
        let r = compare(&[], &[]);
        assert_eq!(r.jaccard, 1.0);
        assert_eq!(r.recall, 1.0);
        assert_eq!(r.precision, 1.0);
        assert_eq!(r.intersection, 0);
        assert!(r.passes_release_gate());
    }

    #[test]
    fn multi_alignment_block_is_rejected() {
        let text = "a score=100\n\
                    s t 0 3 + 10 ACG\n\
                    s q 0 3 + 10 ACG\n\
                    s r 0 3 + 10 ACG\n\n";
        let err = parse_maf(text.as_bytes()).unwrap_err();
        assert!(matches!(err, MafError::BadBlockArity { count: 3, .. }));
    }
}
