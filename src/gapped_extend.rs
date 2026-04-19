//! Scalar 3-state affine-gap x-drop extension.
//!
//! Given an anchor at `(t_anchor, q_anchor)` that is known to be a high-
//! scoring column, extend leftward and rightward with the standard
//! Gotoh-style 3-state affine recurrence:
//!
//! ```text
//! M[i,j] = s(t[i-1], q[j-1]) + max( M[i-1,j-1], X[i-1,j-1], Y[i-1,j-1] )
//! X[i,j] = max( M[i-1,j] - (O+E), X[i-1,j] - E )   // gap in query
//! Y[i,j] = max( M[i,j-1] - (O+E), Y[i,j-1] - E )   // gap in target
//! ```
//!
//! Extension terminates along each direction when `best_seen - current_best`
//! exceeds `y_drop` — equivalent to the upstream `--ydrop` flag.
//!
//! The returned `GappedAlignment` spans the entire aligned region (left
//! extension + anchor column + right extension) with one merged edit script.
//! This is the scalar reference implementation; Phase 3 bolts on striped
//! SIMD lanes while keeping this module as the parity oracle.

use crate::anchor::Anchor;
use crate::edit_script::{EditOp, EditScript};
use crate::scoring::ScoringMatrix;
use crate::sequences::PackedSeq;

/// Parameters for gapped extension. Defaults match upstream HOXD70.
#[derive(Debug, Clone, Copy)]
pub struct GappedParams {
    pub y_drop: i32,
    pub gapped_threshold: i32,
}

impl Default for GappedParams {
    fn default() -> Self {
        Self {
            // `--ydrop` default in upstream lastz for HOXD70.
            y_drop: 9_400,
            // `--gappedthresh` default.
            gapped_threshold: 3_000,
        }
    }
}

const NEG_INF: i32 = i32::MIN / 4;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum State {
    M,
    X,
    Y,
}

/// The product of gapped extension.
#[derive(Debug, Clone)]
pub struct GappedAlignment {
    pub t_start: u32,
    pub q_start: u32,
    pub t_len: u32,
    pub q_len: u32,
    pub score: i32,
    pub script: EditScript,
    pub source: Anchor,
}

impl GappedAlignment {
    #[inline]
    pub fn t_end(&self) -> u32 {
        self.t_start + self.t_len
    }
    #[inline]
    pub fn q_end(&self) -> u32 {
        self.q_start + self.q_len
    }
}

/// Extend an anchor into a gapped alignment. Returns `None` if the total
/// score does not reach `params.gapped_threshold`.
pub fn extend(
    anchor: Anchor,
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    params: &GappedParams,
) -> Option<GappedAlignment> {
    // Anchor column score first (single matched position).
    let anchor_t = anchor.t_pos as usize;
    let anchor_q = anchor.q_pos as usize;
    let anchor_col_score = column_score(target, query, matrix, anchor_t, anchor_q);

    // Right extension: slice the target and query beginning one base past
    // the anchor and DP them head-on.
    let right = extend_one_side(
        slice_codes(target, anchor_t + 1, target.len()),
        slice_codes(query, anchor_q + 1, query.len()),
        matrix,
        params,
    );

    // Left extension: reverse both input strings, run the same DP, then
    // reverse the resulting script back so positions line up.
    let mut t_left = slice_codes(target, 0, anchor_t);
    let mut q_left = slice_codes(query, 0, anchor_q);
    t_left.reverse();
    q_left.reverse();
    let left_raw = extend_one_side(t_left, q_left, matrix, params);

    let total_score =
        anchor_col_score.saturating_add(left_raw.score).saturating_add(right.score);
    if total_score < params.gapped_threshold {
        return None;
    }

    let mut script = EditScript::new();
    // Left script is built from the reversed inputs, so we emit it in
    // reverse to restore the true left-to-right order.
    script.extend(&left_raw.script.reversed());
    script.push(EditOp::Match, 1); // anchor column
    script.extend(&right.script);

    let t_start = anchor.t_pos - left_raw.t_span;
    let q_start = anchor.q_pos - left_raw.q_span;
    let t_len = left_raw.t_span + 1 + right.t_span;
    let q_len = left_raw.q_span + 1 + right.q_span;

    Some(GappedAlignment {
        t_start,
        q_start,
        t_len,
        q_len,
        score: total_score,
        script,
        source: anchor,
    })
}

/// Copy the (code, is_valid) pairs of `seq[start..end]` into a small
/// in-memory slice. We keep the validity bit so invalid positions reuse the
/// matrix's `ambig_score` without per-cell branching on bit twiddles.
fn slice_codes(seq: &PackedSeq, start: usize, end: usize) -> Vec<(u8, bool)> {
    (start..end)
        .map(|i| (seq.code(i), seq.is_valid(i)))
        .collect()
}

#[inline]
fn column_score(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    t: usize,
    q: usize,
) -> i32 {
    if !target.is_valid(t) || !query.is_valid(q) {
        matrix.ambig_score
    } else {
        matrix.score(target.code(t), query.code(q))
    }
}

#[inline]
fn pair_score(t: (u8, bool), q: (u8, bool), matrix: &ScoringMatrix) -> i32 {
    if t.1 && q.1 {
        matrix.score(t.0, q.0)
    } else {
        matrix.ambig_score
    }
}

struct OneSide {
    score: i32,
    t_span: u32,
    q_span: u32,
    script: EditScript,
}

/// Fill the 3-state affine DP for one direction (forward only).
/// Caller is responsible for reversing the input slices and output script
/// when doing a left extension.
fn extend_one_side(
    target: Vec<(u8, bool)>,
    query: Vec<(u8, bool)>,
    matrix: &ScoringMatrix,
    params: &GappedParams,
) -> OneSide {
    let n = target.len();
    let m = query.len();

    if n == 0 || m == 0 {
        return OneSide {
            score: 0,
            t_span: 0,
            q_span: 0,
            script: EditScript::new(),
        };
    }

    let gap_open_ext = matrix.gap_open + matrix.gap_extend; // cost of opening a gap
    let gap_ext = matrix.gap_extend; // cost of extending an existing gap
    let y_drop = params.y_drop;

    // Adaptive band radius. A single Y-state chain (gap in target) can
    // propagate a score of `best - (gap_open_ext + k*gap_extend)` through
    // `k` cells before y-drop pruning fires. The band at any row can
    // therefore grow by up to `(y_drop - gap_open_ext) / gap_extend + 1`
    // columns in either direction relative to the previous row's band.
    // We use that as the per-row expansion radius — tight enough that
    // full-matrix behaviour is preserved (cells outside the band would be
    // pruned as below-threshold anyway), loose enough to not truncate
    // any legitimate alignment path.
    //
    // A fixed +1 was tried in a previous pass and regressed parity 14→7
    // on cat × pig; the computed radius is ~300 for HOXD70 defaults.
    // See PLAN.md §3.4 for the failure-mode writeup.
    let band_radius: usize = if gap_ext > 0 {
        ((y_drop.saturating_sub(gap_open_ext) / gap_ext).max(0) + 2) as usize
    } else {
        m + 1 // unbounded — fall back to full-matrix
    };

    // Band-allocated DP storage. Each row stores `band_width` cells,
    // with `row_lo[i]` mapping local index 0 to global j-coordinate
    // `row_lo[i]`. Allocation is `(n+1) * band_width` instead of the
    // `(n+1) * (m+1)` full-matrix scheme that came before — for a
    // 30 kbp x 30 kbp extension that's ~430 MB instead of ~21 GB,
    // so wall-time on the sars-cov-* fixture drops from ~3 min to
    // a few seconds (PLAN.md §8 open perf bug, now closed).
    //
    // band_width must contain any single row's compute window, which
    // expands by `band_radius` on each side of the previous row's
    // alive band. Worst-case width is `2 * band_radius` (alive band)
    // plus `2 * band_radius` (per-side expansion) plus a small slack
    // for the off-by-one. For tiny m, fall back to full row width.
    let band_width = (4 * band_radius + 4).min(m + 1);
    let mut row_lo: Vec<usize> = vec![0; n + 1];
    let size = (n + 1) * band_width;

    let mut mm = vec![NEG_INF; size];
    let mut xx = vec![NEG_INF; size];
    let mut yy = vec![NEG_INF; size];
    let mut m_from = vec![State::M; size];
    let mut x_from = vec![State::M; size];
    let mut y_from = vec![State::M; size];

    // Per-cell bounds checks against the band are inlined into the
    // DP inner loop. The check is a pair of `usize` comparisons
    // against the hoisted `lo_prev` / `hi_prev_excl` and the current
    // row's `lo_this`; see the loop body below.

    mm[0] = 0;

    let mut best_score: i32 = 0;
    let mut best_cell = (0usize, 0usize, State::M);

    // Band bounds for the previous row. Row 0's band is the single seed
    // cell (0, 0); subsequent rows expand to [lo - band_radius, hi + band_radius].
    let mut lo: usize = 0;
    let mut hi: usize = 0;

    // Row 0 (the seed cell) lives at local index 0 of row 0's storage.
    row_lo[0] = 0;

    for i in 0..=n {
        let j_lo = lo.saturating_sub(band_radius);
        let j_hi = (hi + band_radius).min(m).min(j_lo + band_width - 1);

        // Position row i's storage so j_lo lands at local index 0.
        // The compute window [j_lo, j_hi] fits inside the band_width-
        // wide slice by construction.
        row_lo[i] = j_lo;

        // Hoist per-row invariants out of the inner-cell loop. The
        // index arithmetic inside the loop becomes `j - lo_this`
        // instead of `j - row_lo[i]`, and the prev-row bounds reduce
        // to compares against two locals (lo_prev, hi_prev_excl).
        let lo_this = j_lo;
        let lo_prev = if i >= 1 { row_lo[i - 1] } else { 0 };
        let hi_prev_excl = lo_prev + band_width;
        let base_this = i * band_width;
        let base_prev = if i >= 1 { (i - 1) * band_width } else { 0 };

        for j in j_lo..=j_hi {
            if i == 0 && j == 0 {
                continue;
            }

            let mut best = NEG_INF;
            let here = base_this + (j - lo_this);

            // Saturating arithmetic is overkill — NEG_INF is
            // `i32::MIN / 4`, so `NEG_INF - gap_open_ext` = ~-537M
            // and cumulative subtraction over the whole extension
            // (bounded by O(sequence_len)) never reaches `i32::MIN`.

            // Read the three prev-row neighbours (`(i-1, j)` and
            // `(i-1, j-1)`) with a single bounds check each, inlined
            // from the old closure-based helper.
            let prev_in_band = |j: usize| j >= lo_prev && j < hi_prev_excl;

            // X: gap in query (target advances, query does not). Requires i>=1.
            if i >= 1 {
                let (prev_m, prev_x) = if prev_in_band(j) {
                    let off = base_prev + (j - lo_prev);
                    (mm[off], xx[off])
                } else {
                    (NEG_INF, NEG_INF)
                };
                let from_m = prev_m - gap_open_ext;
                let from_x = prev_x - gap_ext;
                let (v, from) = if from_m >= from_x {
                    (from_m, State::M)
                } else {
                    (from_x, State::X)
                };
                xx[here] = v;
                x_from[here] = from;
                if v > best {
                    best = v;
                }
            }

            // Y: gap in target (query advances, target does not). Requires j>=1.
            if j >= 1 {
                // In-current-row neighbour at (i, j-1). Always in band
                // because j-1 >= lo_this (unless j == lo_this, which
                // reads the cell below lo_this → NEG_INF).
                let (prev_m, prev_y) = if j - 1 >= lo_this {
                    let off = base_this + (j - 1 - lo_this);
                    (mm[off], yy[off])
                } else {
                    (NEG_INF, NEG_INF)
                };
                let from_m = prev_m - gap_open_ext;
                let from_y = prev_y - gap_ext;
                let (v, from) = if from_m >= from_y {
                    (from_m, State::M)
                } else {
                    (from_y, State::Y)
                };
                yy[here] = v;
                y_from[here] = from;
                if v > best {
                    best = v;
                }
            }

            // M: match/mismatch column. Requires i>=1 && j>=1.
            if i >= 1 && j >= 1 {
                let s = pair_score(target[i - 1], query[j - 1], matrix);
                let (from_m, from_x, from_y) = if prev_in_band(j - 1) {
                    let off = base_prev + (j - 1 - lo_prev);
                    (mm[off], xx[off], yy[off])
                } else {
                    (NEG_INF, NEG_INF, NEG_INF)
                };
                let (prev, from) = if from_m >= from_x && from_m >= from_y {
                    (from_m, State::M)
                } else if from_x >= from_y {
                    (from_x, State::X)
                } else {
                    (from_y, State::Y)
                };
                let v = prev + s;
                mm[here] = v;
                m_from[here] = from;
                if v > best {
                    best = v;
                }
            }

            if best > best_score {
                best_score = best;
                let state = if mm[here] == best {
                    State::M
                } else if xx[here] == best {
                    State::X
                } else {
                    State::Y
                };
                best_cell = (i, j, state);
            }
        }

        // X-drop prune: within the cells we actually touched this row,
        // any whose best state is more than `y_drop` below `best_score`
        // is dead. Shrink the band to the surviving window; cells
        // outside `[j_lo, j_hi]` weren't allocated in this row and
        // don't need touching.
        let threshold = best_score.saturating_sub(y_drop);
        let mut new_lo: Option<usize> = None;
        let mut new_hi = 0usize;
        for j in j_lo..=j_hi {
            let here = i * band_width + (j - row_lo[i]);
            let alive = mm[here].max(xx[here]).max(yy[here]);
            if alive >= threshold {
                if new_lo.is_none() {
                    new_lo = Some(j);
                }
                new_hi = j;
            } else {
                // Kill this cell so downstream transitions don't revive
                // it via saturating arithmetic.
                mm[here] = NEG_INF;
                xx[here] = NEG_INF;
                yy[here] = NEG_INF;
            }
        }
        if new_lo.is_none() {
            // Entire row is dead → extension cannot improve.
            break;
        }
        lo = new_lo.unwrap();
        hi = new_hi;
        if lo == hi && i > 0 {
            let here = i * band_width + (lo - row_lo[i]);
            let alive = mm[here].max(xx[here]).max(yy[here]);
            if alive < best_score.saturating_sub(y_drop) {
                break;
            }
        }
    }

    // Traceback from best cell. Cell (i, j) lives at flat index
    // `i * band_width + (j - row_lo[i])`; the band-allocated layout
    // doesn't change the traceback algorithm, only how we read the
    // stored from-state for each cell.
    let (mut i, mut j, mut state) = best_cell;
    let mut script = EditScript::new();
    let here = |i: usize, j: usize| i * band_width + (j - row_lo[i]);
    while i > 0 || j > 0 {
        match state {
            State::M => {
                script.push(EditOp::Match, 1);
                let prev = m_from[here(i, j)];
                i -= 1;
                j -= 1;
                state = prev;
            }
            State::X => {
                // Gap in query: target advanced, query did not.
                script.push(EditOp::DeleteQuery, 1);
                let prev = x_from[here(i, j)];
                i -= 1;
                state = prev;
            }
            State::Y => {
                // Gap in target: query advanced, target did not.
                script.push(EditOp::InsertQuery, 1);
                let prev = y_from[here(i, j)];
                j -= 1;
                state = prev;
            }
        }
    }
    // Script was emitted tip-to-root, i.e. right-to-left. Flip.
    let script = script.reversed();

    OneSide {
        score: best_score,
        t_span: best_cell.0 as u32,
        q_span: best_cell.1 as u32,
        script,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hsp::Hsp;

    fn matrix() -> ScoringMatrix {
        ScoringMatrix::hoxd70()
    }

    fn anchor_from_hsp(hsp: Hsp, offset: u32) -> Anchor {
        Anchor {
            t_pos: hsp.t_start + offset,
            q_pos: hsp.q_start + offset,
            source: hsp,
        }
    }

    #[test]
    fn perfect_match_extension_stays_ungapped() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGTACGT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 16, score: 0 };
        let anchor = anchor_from_hsp(hsp, 8);
        let params = GappedParams { y_drop: 500, gapped_threshold: 100 };
        let aln = extend(anchor, &seq, &seq, &matrix(), &params).unwrap();
        assert_eq!(aln.t_start, 0);
        assert_eq!(aln.q_start, 0);
        assert_eq!(aln.t_len, 16);
        assert_eq!(aln.q_len, 16);
        // All columns should be M.
        assert_eq!(aln.script.to_cigar(), "16M");
    }

    #[test]
    fn single_insertion_is_recovered() {
        // Target missing one base relative to query: best alignment has a
        // 1bp gap in the target (InsertQuery).
        let target = PackedSeq::from_ascii(b"ACGTACGTACGT");
        let query = PackedSeq::from_ascii(b"ACGTACGGTACGT"); // extra G at pos 7
        let hsp = Hsp { t_start: 0, q_start: 0, length: 6, score: 0 };
        let anchor = anchor_from_hsp(hsp, 3);
        let params = GappedParams { y_drop: 1_000, gapped_threshold: 100 };
        let aln = extend(anchor, &target, &query, &matrix(), &params).unwrap();
        assert_eq!(aln.t_len, 12);
        assert_eq!(aln.q_len, 13);
        let cigar = aln.script.to_cigar();
        // Must contain exactly one 1I somewhere; the rest matches.
        assert!(cigar.contains("1I"), "cigar was {cigar}");
        assert_eq!(aln.script.target_span(), 12);
        assert_eq!(aln.script.query_span(), 13);
    }

    #[test]
    fn single_deletion_is_recovered() {
        // Query missing one base relative to target.
        let target = PackedSeq::from_ascii(b"ACGTACGGTACGT"); // extra G at pos 7
        let query = PackedSeq::from_ascii(b"ACGTACGTACGT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 6, score: 0 };
        let anchor = anchor_from_hsp(hsp, 3);
        let params = GappedParams { y_drop: 1_000, gapped_threshold: 100 };
        let aln = extend(anchor, &target, &query, &matrix(), &params).unwrap();
        assert_eq!(aln.t_len, 13);
        assert_eq!(aln.q_len, 12);
        let cigar = aln.script.to_cigar();
        assert!(cigar.contains("1D"), "cigar was {cigar}");
    }

    #[test]
    fn below_threshold_returns_none() {
        let seq = PackedSeq::from_ascii(b"ACGT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 4, score: 0 };
        let anchor = anchor_from_hsp(hsp, 1);
        let params = GappedParams { y_drop: 500, gapped_threshold: 100_000 };
        assert!(extend(anchor, &seq, &seq, &matrix(), &params).is_none());
    }

    #[test]
    fn extension_stops_at_diverging_tail() {
        // Matching 24-bp core on diagonal 0, flanked on both sides by
        // strong mismatches (T vs A). y-drop should stop the extension
        // close to the core/flank boundary.
        let target = PackedSeq::from_ascii(b"TTTTTTTTACGTACGTACGTACGTACGTACGTTTTTTTTT");
        let query = PackedSeq::from_ascii(b"AAAAAAAAACGTACGTACGTACGTACGTACGTAAAAAAAA");
        // ^ Note: both strings are 40 bp. Target[8..32] and query[8..32]
        // are both `ACGTACGTACGTACGTACGTACGT` aligned on diagonal 0.
        let hsp = Hsp { t_start: 8, q_start: 8, length: 24, score: 0 };
        let anchor = anchor_from_hsp(hsp, 12);
        let params = GappedParams { y_drop: 300, gapped_threshold: 500 };
        let aln = extend(anchor, &target, &query, &matrix(), &params).unwrap();
        assert!(aln.t_len >= 24 && aln.t_len <= 32, "t_len={}", aln.t_len);
    }

    #[test]
    fn banded_dp_accommodates_large_gap_within_y_drop() {
        // Construct an alignment with a gap near the start so the band
        // has to grow fast on the right side to keep up. The gap size is
        // chosen to fit comfortably within the default band radius
        // (y_drop/gap_extend ≈ 313); a naïve "+1 per row" band would
        // truncate it.
        let core = b"GATTACACATGGCATGTCGAATGCCTAGCATGTCA"; // 35 bp
        let mut target = core.to_vec();
        target.extend_from_slice(core); // 70 bp matching run
        let mut query = Vec::new();
        query.extend_from_slice(core);
        // Insert a 50-bp insertion in the query between the two core runs;
        // fits within band_radius but way beyond any fixed +1/row scheme.
        query.extend_from_slice(b"AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAG");
        query.extend_from_slice(core);
        let t = PackedSeq::from_ascii(&target);
        let q = PackedSeq::from_ascii(&query);

        let hsp = Hsp { t_start: 0, q_start: 0, length: 10, score: 0 };
        let anchor = anchor_from_hsp(hsp, 5);
        let params = GappedParams { y_drop: 9_400, gapped_threshold: 500 };
        let aln = extend(anchor, &t, &q, &matrix(), &params).unwrap();

        // Alignment should reach the far side of the insertion — its
        // target span covers both cores (~70 bp) and the query span is
        // ~120 bp (70 core + 50 insertion).
        assert!(aln.t_len >= 60, "t_len = {}; banded DP truncated", aln.t_len);
        assert!(aln.q_len >= 110, "q_len = {}; banded DP truncated", aln.q_len);
        // Script must encode the long Y-state run (InsertQuery).
        assert!(aln.script.to_cigar().contains("I"));
    }

    #[test]
    fn script_spans_equal_reported_lengths() {
        // Property: whatever the shape of the alignment, the script's
        // per-side spans must match the reported t_len / q_len.
        let target = PackedSeq::from_ascii(b"AAAGATTACACATGGCATGTCGATTTT");
        let query = PackedSeq::from_ascii(b"GATACACATGGCAGTCGA"); // small indels
        let hsp = Hsp { t_start: 3, q_start: 0, length: 6, score: 0 };
        let anchor = anchor_from_hsp(hsp, 3);
        let params = GappedParams { y_drop: 5_000, gapped_threshold: 100 };
        let aln = extend(anchor, &target, &query, &matrix(), &params).unwrap();
        assert_eq!(aln.script.target_span(), aln.t_len);
        assert_eq!(aln.script.query_span(), aln.q_len);
    }
}
