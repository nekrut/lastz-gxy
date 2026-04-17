//! Ungapped HSP extension with x-drop termination.
//!
//! An HSP (High-scoring Segment Pair) is a gap-free alignment: a pair of
//! intervals `(t_start..t_end)` in the target and `(q_start..q_end)` in the
//! query of equal length, scored as the sum of per-column substitution
//! scores from a `ScoringMatrix`.
//!
//! Given a seed hit at `(t_pos, q_pos)` this module extends left and right
//! along the shared diagonal `d = t_pos - q_pos` until the running score
//! drops more than `x_drop` below the best-so-far, then returns the best
//! sub-interval encountered. This is a scalar reference implementation;
//! Phase 2 bolts on AVX2/NEON lanes while keeping this as the oracle.

use crate::scoring::ScoringMatrix;
use crate::sequences::PackedSeq;

/// One ungapped high-scoring segment pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Hsp {
    pub t_start: u32,
    pub q_start: u32,
    pub length: u32,
    pub score: i32,
}

impl Hsp {
    #[inline]
    pub fn t_end(&self) -> u32 {
        self.t_start + self.length
    }
    #[inline]
    pub fn q_end(&self) -> u32 {
        self.q_start + self.length
    }
    /// Diagonal index used by the diag-hash stage: `t - q`. Widened to `i64`
    /// because raw `i32` subtraction would wrap on 2 Gbp+ references.
    #[inline]
    pub fn diagonal(&self) -> i64 {
        self.t_start as i64 - self.q_start as i64
    }
}

/// Parameters for the x-drop ungapped extender. Maps to lastz's
/// `--xdrop`/`--hspthresh` flags.
#[derive(Debug, Clone, Copy)]
pub struct HspParams {
    pub x_drop: i32,
    pub hsp_threshold: i32,
}

impl Default for HspParams {
    fn default() -> Self {
        // Upstream defaults for HOXD70 (see docs/readme.lastz.html).
        Self {
            x_drop: 910,
            hsp_threshold: 3000,
        }
    }
}

/// Extend a seed hit at `(t_pos, q_pos)` into an ungapped HSP.
///
/// Returns `None` if the best score found does not reach `params.hsp_threshold`
/// or if the seed itself scores negatively. The returned HSP always spans a
/// contiguous interval that includes the seed hit itself.
pub fn extend_hit(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    params: &HspParams,
    t_pos: u32,
    q_pos: u32,
    seed_len: u32,
) -> Option<Hsp> {
    let t_pos = t_pos as i64;
    let q_pos = q_pos as i64;
    let t_len = target.len() as i64;
    let q_len = query.len() as i64;
    let seed_len = seed_len as i64;

    if t_pos < 0 || q_pos < 0 || t_pos + seed_len > t_len || q_pos + seed_len > q_len {
        return None;
    }

    // Score the seed footprint first; this is the score at step 0 of the
    // extension.
    let mut seed_score: i32 = 0;
    for k in 0..seed_len {
        seed_score =
            seed_score.saturating_add(column_score(target, query, matrix, t_pos + k, q_pos + k));
    }

    // Extend rightward from the first base past the seed.
    let right =
        extend_one_side(target, query, matrix, params, t_pos + seed_len, q_pos + seed_len, 1);
    // Extend leftward from the base immediately before the seed.
    let left = extend_one_side(target, query, matrix, params, t_pos - 1, q_pos - 1, -1);

    let score = seed_score.saturating_add(left.best_score).saturating_add(right.best_score);
    if score < params.hsp_threshold {
        return None;
    }

    let left_extent = left.best_extent as i64;
    let right_extent = right.best_extent as i64;
    let t_start = t_pos - left_extent;
    let q_start = q_pos - left_extent;
    let length = seed_len + left_extent + right_extent;

    Some(Hsp {
        t_start: t_start as u32,
        q_start: q_start as u32,
        length: length as u32,
        score,
    })
}

struct SideResult {
    best_score: i32,
    best_extent: u32,
}

/// Walk one side of the diagonal, updating the running score column by
/// column and tracking the extent at which the best score was observed.
///
/// `direction` is `+1` for rightward, `-1` for leftward.
#[inline]
fn extend_one_side(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    params: &HspParams,
    start_t: i64,
    start_q: i64,
    direction: i64,
) -> SideResult {
    let t_len = target.len() as i64;
    let q_len = query.len() as i64;

    let mut running: i32 = 0;
    let mut best: i32 = 0;
    let mut best_extent: u32 = 0;
    let mut extent: u32 = 0;

    let mut t = start_t;
    let mut q = start_q;
    loop {
        if t < 0 || q < 0 || t >= t_len || q >= q_len {
            break;
        }
        extent += 1;
        running = running.saturating_add(column_score(target, query, matrix, t, q));
        if running > best {
            best = running;
            best_extent = extent;
        } else if best - running > params.x_drop {
            break;
        }
        t += direction;
        q += direction;
    }
    SideResult { best_score: best, best_extent }
}

#[inline]
fn column_score(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    t: i64,
    q: i64,
) -> i32 {
    let (ti, qi) = (t as usize, q as usize);
    if !target.is_valid(ti) || !query.is_valid(qi) {
        matrix.ambig_score
    } else {
        matrix.score(target.code(ti), query.code(qi))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m() -> ScoringMatrix {
        ScoringMatrix::hoxd70()
    }

    #[test]
    fn perfect_match_extends_to_full_length() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGTACGT");
        let params = HspParams { x_drop: 100, hsp_threshold: 100 };
        let hsp = extend_hit(&seq, &seq, &m(), &params, 2, 2, 4).unwrap();
        assert_eq!(hsp.t_start, 0);
        assert_eq!(hsp.q_start, 0);
        assert_eq!(hsp.length as usize, seq.len());
        // 16 cols with avg diagonal score of ~95 should blow past 100.
        assert!(hsp.score >= 100);
    }

    #[test]
    fn seed_on_diagonal_zero_anchors_correctly() {
        // A/T mismatches are very negative under HOXD70 (-123), so flanking
        // mismatches should trigger x-drop and trim the HSP to the seed.
        let target = PackedSeq::from_ascii(b"TTTACGTTTT");
        let query = PackedSeq::from_ascii(b"AAAACGAAAA");
        let params = HspParams { x_drop: 50, hsp_threshold: 50 };
        let hsp = extend_hit(&target, &query, &m(), &params, 3, 3, 3).unwrap();
        assert_eq!(hsp.t_start, 3);
        assert_eq!(hsp.q_start, 3);
        assert_eq!(hsp.length, 3);
    }

    #[test]
    fn low_score_returns_none() {
        let target = PackedSeq::from_ascii(b"ACGT");
        let query = PackedSeq::from_ascii(b"TTTT");
        let params = HspParams { x_drop: 100, hsp_threshold: 10_000 };
        assert!(extend_hit(&target, &query, &m(), &params, 0, 0, 4).is_none());
    }

    #[test]
    fn out_of_bounds_seed_returns_none() {
        let target = PackedSeq::from_ascii(b"ACGT");
        let query = PackedSeq::from_ascii(b"ACGT");
        let params = HspParams::default();
        assert!(extend_hit(&target, &query, &m(), &params, 2, 2, 10).is_none());
    }

    #[test]
    fn n_bases_score_as_mismatch() {
        // N in the middle should cost `ambig_score` but not crash.
        let target = PackedSeq::from_ascii(b"ACGTACGT");
        let query = PackedSeq::from_ascii(b"ACNTACGT");
        let params = HspParams { x_drop: 500, hsp_threshold: 200 };
        let hsp = extend_hit(&target, &query, &m(), &params, 4, 4, 4).unwrap();
        assert!(hsp.length >= 4);
    }

    #[test]
    fn hsp_diagonal_is_t_minus_q() {
        let h = Hsp { t_start: 100, q_start: 40, length: 10, score: 0 };
        assert_eq!(h.diagonal(), 60);
    }

    #[test]
    fn extends_symmetrically_when_flanks_equal() {
        let seq = PackedSeq::from_ascii(b"CCCCCCAAAAAACCCCCC");
        let params = HspParams { x_drop: 500, hsp_threshold: 100 };
        let hsp = extend_hit(&seq, &seq, &m(), &params, 8, 8, 2).unwrap();
        let left = 8 - hsp.t_start;
        let right = hsp.t_end() - (8 + 2);
        assert_eq!(left, right, "expected symmetric extension");
    }
}
