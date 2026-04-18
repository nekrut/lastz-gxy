//! Shannon-entropy gate on HSP base composition.
//!
//! KegAlign (and other repeat-sensitive aligners) filter HSPs by combining
//! alignment score with base-composition entropy: homopolymeric or
//! dinucleotide-repeat HSPs get dropped regardless of raw score.
//!
//! For an HSP over target positions `[t_start, t_end)`, compute
//! `H = -Σ p_i log2(p_i)` where `p_i` are the frequencies of `A`, `C`, `G`,
//! `T` in the target slice (`N` / invalid bases are excluded from the
//! count). `H` ranges from 0 (pure homopolymer) to 2 (balanced 4-letter).
//!
//! An HSP passes the gate when
//!
//! ```text
//! score * (H / 2) >= entropy_threshold
//! ```
//!
//! With `entropy_threshold = hsp_threshold` this exactly subsumes the
//! default score gate at balanced-base regions and only tightens for
//! low-complexity stretches. The gate is opt-in (upstream lastz does not
//! apply it by default), behind the `--entropy` CLI flag.

use crate::sequences::PackedSeq;

/// Shannon entropy (in bits) of the A/C/G/T frequency over
/// `seq[start..end]`. Returns `0.0` when the window has no valid ACGT bases.
pub fn shannon_entropy(seq: &PackedSeq, start: usize, end: usize) -> f64 {
    let mut counts = [0u32; 4];
    let mut total = 0u32;
    for i in start..end.min(seq.len()) {
        if seq.is_valid(i) {
            counts[seq.code(i) as usize] += 1;
            total += 1;
        }
    }
    if total == 0 {
        return 0.0;
    }
    let t = total as f64;
    let mut h = 0.0f64;
    for &c in &counts {
        if c == 0 {
            continue;
        }
        let p = c as f64 / t;
        h -= p * p.log2();
    }
    h
}

/// Does the HSP (target interval + score) pass the entropy gate? The gate
/// scales the HSP's score by `H / 2` (`H_max = 2` bits for ACGT) and
/// compares against `threshold`. An HSP with perfectly balanced bases
/// passes iff `score >= threshold`; a homopolymer is always rejected.
#[inline]
pub fn passes_entropy_gate(
    target: &PackedSeq,
    t_start: u32,
    t_end: u32,
    score: i32,
    threshold: i32,
) -> bool {
    let h = shannon_entropy(target, t_start as usize, t_end as usize);
    let scaled = (score as f64) * (h / 2.0);
    scaled >= threshold as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn balanced_composition_reaches_max_entropy() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGTACGT");
        let h = shannon_entropy(&seq, 0, 16);
        assert!((h - 2.0).abs() < 1e-9, "H={h}");
    }

    #[test]
    fn homopolymer_has_zero_entropy() {
        let seq = PackedSeq::from_ascii(b"AAAAAAAAAA");
        assert_eq!(shannon_entropy(&seq, 0, 10), 0.0);
    }

    #[test]
    fn dinucleotide_repeat_has_entropy_of_one() {
        let seq = PackedSeq::from_ascii(b"ATATATATAT");
        let h = shannon_entropy(&seq, 0, 10);
        // Two bases at 50/50 → 1 bit.
        assert!((h - 1.0).abs() < 1e-9, "H={h}");
    }

    #[test]
    fn invalid_bases_excluded_from_count() {
        let seq = PackedSeq::from_ascii(b"NNNNACGTACGT");
        let h_all = shannon_entropy(&seq, 0, seq.len());
        // Ns skipped; remaining ACGT-repeating window gives H = 2.
        assert!((h_all - 2.0).abs() < 1e-9, "H={h_all}");
    }

    #[test]
    fn gate_rejects_homopolymer_even_at_high_score() {
        let seq = PackedSeq::from_ascii(b"AAAAAAAAAA");
        assert!(!passes_entropy_gate(&seq, 0, 10, 10_000, 1_000));
    }

    #[test]
    fn gate_accepts_balanced_block_at_threshold() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGT");
        // Balanced 4-letter: score-scaling factor = H/2 = 1.
        assert!(passes_entropy_gate(&seq, 0, 12, 1_000, 1_000));
        // One below threshold.
        assert!(!passes_entropy_gate(&seq, 0, 12, 999, 1_000));
    }

    #[test]
    fn gate_partially_relaxes_dinucleotide_at_high_score() {
        // ATAT... dinuc, H=1 → scaling 0.5. Need score ≥ 2*threshold.
        let seq = PackedSeq::from_ascii(b"ATATATATAT");
        assert!(!passes_entropy_gate(&seq, 0, 10, 1_000, 1_000));
        assert!(passes_entropy_gate(&seq, 0, 10, 2_000, 1_000));
    }
}
