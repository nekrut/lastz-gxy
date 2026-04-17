//! Anchor selection for gapped extension.
//!
//! Upstream lastz (`segment.c`) slides a fixed-width window over each HSP and
//! picks the midpoint of the window with the highest local score. The anchor
//! becomes the seed for the 2D affine-gap extension that follows. For parity
//! we use the same sliding-window midpoint rule; with a window of 1 or on
//! short HSPs the anchor is simply the HSP midpoint.

use crate::hsp::Hsp;
use crate::scoring::ScoringMatrix;
use crate::sequences::PackedSeq;

/// The single (target, query) position that `gapped_extend` will extend
/// outward from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchor {
    pub t_pos: u32,
    pub q_pos: u32,
    pub source: Hsp,
}

/// Pick an anchor from `hsp` using a sliding window of `window` columns. The
/// anchor is the (target,query) midpoint of the highest-scoring window.
/// `window` is clamped to the HSP length.
pub fn choose_anchor(
    hsp: Hsp,
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    window: u32,
) -> Anchor {
    let length = hsp.length;
    assert!(length >= 1, "empty HSP cannot be anchored");

    let w = window.clamp(1, length) as usize;
    let len = length as usize;

    // Column scores along the HSP's diagonal.
    let mut cols: Vec<i32> = Vec::with_capacity(len);
    for k in 0..len {
        let ti = hsp.t_start as usize + k;
        let qi = hsp.q_start as usize + k;
        let s = if !target.is_valid(ti) || !query.is_valid(qi) {
            matrix.ambig_score
        } else {
            matrix.score(target.code(ti), query.code(qi))
        };
        cols.push(s);
    }

    // Rolling window sum.
    let mut best_start = 0usize;
    let mut best_score: i32 = cols[..w].iter().sum();
    let mut window_score = best_score;
    for start in 1..=len - w {
        window_score = window_score - cols[start - 1] + cols[start + w - 1];
        if window_score > best_score {
            best_score = window_score;
            best_start = start;
        }
    }

    let mid_offset = (best_start + w / 2) as u32;
    Anchor {
        t_pos: hsp.t_start + mid_offset,
        q_pos: hsp.q_start + mid_offset,
        source: hsp,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix() -> ScoringMatrix {
        ScoringMatrix::hoxd70()
    }

    #[test]
    fn anchor_midpoint_on_uniform_hsp() {
        // All matches → any window scores identically, the first (and only
        // tied) window wins so the anchor is w/2 into the HSP.
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGTACGT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 16, score: 0 };
        let a = choose_anchor(hsp, &seq, &seq, &matrix(), 4);
        assert_eq!(a.t_pos, 2); // best_start=0, w=4, mid_offset=2
        assert_eq!(a.q_pos, 2);
    }

    #[test]
    fn anchor_finds_high_scoring_center() {
        // Flanks are mismatches (A/T), center is matches (CCCC/CCCC). A
        // 4-wide window centered over the C-run maximises score.
        let target = PackedSeq::from_ascii(b"AACCCCAA");
        let query = PackedSeq::from_ascii(b"TTCCCCTT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 8, score: 0 };
        let a = choose_anchor(hsp, &target, &query, &matrix(), 4);
        // The 4-wide window with start=2 contains CCCC/CCCC and wins;
        // midpoint offset = 2 + 2 = 4.
        assert_eq!(a.t_pos, 4);
        assert_eq!(a.q_pos, 4);
    }

    #[test]
    fn window_longer_than_hsp_is_clamped() {
        let seq = PackedSeq::from_ascii(b"ACGT");
        let hsp = Hsp { t_start: 0, q_start: 0, length: 4, score: 0 };
        let a = choose_anchor(hsp, &seq, &seq, &matrix(), 1000);
        // Window clamped to 4; best_start forced to 0, midpoint = 2.
        assert_eq!(a.t_pos, 2);
    }

    #[test]
    fn anchor_preserves_source_hsp() {
        let seq = PackedSeq::from_ascii(b"AAAAAAAAAAAAACGTAAAAAAAAAAA");
        let hsp = Hsp { t_start: 13, q_start: 13, length: 4, score: 777 };
        let a = choose_anchor(hsp, &seq, &seq, &matrix(), 2);
        assert_eq!(a.source, hsp);
    }
}
