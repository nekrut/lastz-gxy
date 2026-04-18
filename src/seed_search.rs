//! Query-side seed iteration that drives HSP production.
//!
//! Given a reference `PosTable` and a query `PackedSeq`, this module walks
//! every query window at the configured stride, looks up the seed word
//! (plus optional transition variants), and for each reference hit
//! decides — via `DiagHash` — whether to extend into an ungapped HSP.
//! The resulting HSPs are returned in the order they were emitted (which,
//! because hits flow in query-order, is also diagonal-grouped).

use crate::diag_hash::DiagHash;
use crate::hsp::{extend_hit, Hsp, HspParams};
use crate::pos_table::PosTable;
use crate::scoring::ScoringMatrix;
use crate::seeds::SeedExtractor;
use crate::sequences::PackedSeq;

/// Parameters controlling one `search` invocation.
#[derive(Debug, Clone, Copy)]
pub struct SearchParams {
    pub step: usize,
    pub hsp: HspParams,
    /// Number of transition substitutions tolerated per seed. `0` is
    /// upstream's `--notransition`; `1` is the default (`--transition`);
    /// `2` is `--transition=2`.
    pub transitions: u8,
    /// When `Some(threshold)`, reject HSPs whose entropy-scaled score
    /// (`score * H / 2`, `H` = Shannon entropy of the target slice) is
    /// below `threshold`. KegAlign's low-complexity filter, opt-in.
    pub entropy_threshold: Option<i32>,
}

impl Default for SearchParams {
    fn default() -> Self {
        Self {
            step: 1,
            hsp: HspParams::default(),
            transitions: 1,
            entropy_threshold: None,
        }
    }
}

/// Scan every valid query seed window against the reference index and emit
/// ungapped HSPs that reach `params.hsp.hsp_threshold`.
pub fn search(
    table: &PosTable,
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    params: &SearchParams,
) -> Vec<Hsp> {
    let pattern = table.pattern();
    let seed_len = pattern.len() as u32;
    let weight = pattern.weight();
    let extractor = SeedExtractor::new(pattern);

    let mut out = Vec::new();
    let mut diag = DiagHash::new();

    // Precompute the bit offsets of each care position within the packed
    // seed word. Position 0 is the leftmost care bit (most significant in
    // the word). `shift[i]` is the bit offset of care-position `i`.
    let shifts: Vec<u32> = (0..weight).map(|i| 2 * (weight as u32 - 1 - i as u32)).collect();

    // Collect lookups for the base word and all 1-transition variants.
    // Transition: A↔G (codes 0↔2), C↔T (codes 1↔3). XOR with 0b10 at a
    // care-position's bit pair toggles the transition. Two transitions
    // require nested expansion; we cap at `transitions=2` per upstream.
    let variants_cap = match params.transitions {
        0 => 1,
        1 => 1 + weight,
        _ => 1 + weight + weight * weight.saturating_sub(1) / 2,
    };
    let mut variants: Vec<u64> = Vec::with_capacity(variants_cap);

    for (q_pos, word) in extractor.iter_with_stride(query, params.step) {
        variants.clear();
        variants.push(word);
        if params.transitions >= 1 {
            for &s in &shifts {
                variants.push(word ^ (0b10u64 << s));
            }
        }
        if params.transitions >= 2 {
            for i in 0..shifts.len() {
                for j in (i + 1)..shifts.len() {
                    variants.push(word ^ (0b10u64 << shifts[i]) ^ (0b10u64 << shifts[j]));
                }
            }
        }

        for &v in &variants {
            let hits = table.lookup(v);
            for &t_pos in hits {
                if !diag.admit(t_pos, q_pos) {
                    continue;
                }
                let Some(hsp) = extend_hit(target, query, matrix, &params.hsp, t_pos, q_pos, seed_len)
                else {
                    continue;
                };
                if let Some(thresh) = params.entropy_threshold {
                    if !crate::entropy::passes_entropy_gate(
                        target,
                        hsp.t_start,
                        hsp.t_end(),
                        hsp.score,
                        thresh,
                    ) {
                        continue;
                    }
                }
                diag.mark_covered(hsp.t_end(), hsp.q_end());
                out.push(hsp);
            }
        }
    }

    // Upstream emits HSPs in (diagonal, t_start) order for downstream
    // chaining. We sort here so the output is stable regardless of the
    // interior hit order.
    out.sort_unstable_by_key(|h| (h.diagonal(), h.t_start));
    out.dedup();
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seeds::SeedPattern;

    fn m() -> ScoringMatrix {
        ScoringMatrix::hoxd70()
    }

    #[test]
    fn finds_exact_match_at_expected_coords() {
        // Non-repetitive payload so the correct alignment is the only
        // diagonal with a high-scoring HSP.
        let target = PackedSeq::from_ascii(b"AAAAAGATTACACATGGCATGAAAAA");
        let query = PackedSeq::from_ascii(b"GATTACACATGGCATG");
        let pat = SeedPattern::solid(8);
        let table = PosTable::build(&target, &pat, 1);
        let hsps = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 500, hsp_threshold: 1000 },
                transitions: 0,
                entropy_threshold: None,
            },
        );
        assert_eq!(hsps.len(), 1, "expected one HSP, got {hsps:#?}");
        let h = hsps[0];
        assert_eq!(h.t_start, 5);
        assert_eq!(h.q_start, 0);
        assert_eq!(h.length, 16);
    }

    #[test]
    fn transition_variants_find_off_by_one_hits() {
        // Target has ACGTA. Query has ACATA (G→A is a transition). With
        // --notransition the 5-mer seed should miss; with --transition
        // (default) it should hit.
        let target = PackedSeq::from_ascii(b"AAAAACGTAAAAA");
        let query = PackedSeq::from_ascii(b"ACATA");
        let pat = SeedPattern::solid(5);
        let table = PosTable::build(&target, &pat, 1);

        let no_trans = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 500, hsp_threshold: 50 },
                transitions: 0,
                entropy_threshold: None,
            },
        );
        assert!(no_trans.is_empty(), "unexpected hits without transitions");

        let with_trans = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 500, hsp_threshold: 50 },
                transitions: 1,
                entropy_threshold: None,
            },
        );
        assert!(!with_trans.is_empty(), "expected at least one transition hit");
    }

    #[test]
    fn no_match_yields_no_hsps() {
        let target = PackedSeq::from_ascii(b"AAAAAAAAAA");
        let query = PackedSeq::from_ascii(b"CCCCCCCCCC");
        let pat = SeedPattern::solid(5);
        let table = PosTable::build(&target, &pat, 1);
        let hsps = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 100, hsp_threshold: 100 },
                transitions: 0,
                entropy_threshold: None,
            },
        );
        assert!(hsps.is_empty());
    }

    #[test]
    fn diag_hash_collapses_seeds_inside_one_hsp() {
        // A long homologous non-repetitive region has many overlapping
        // seeds on the same diagonal; DiagHash should collapse them so only
        // one HSP per diagonal survives.
        let target = PackedSeq::from_ascii(b"AAAAAGATTACACATGGCATGTCGAAAAAA");
        let query = PackedSeq::from_ascii(b"GATTACACATGGCATGTCGA");
        let pat = SeedPattern::solid(7);
        let table = PosTable::build(&target, &pat, 1);
        let hsps = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 500, hsp_threshold: 500 },
                transitions: 0,
                entropy_threshold: None,
            },
        );
        let unique_diagonals: std::collections::BTreeSet<_> =
            hsps.iter().map(|h| h.diagonal()).collect();
        assert_eq!(
            hsps.len(),
            unique_diagonals.len(),
            "expected at most one HSP per diagonal, got {hsps:#?}"
        );
        // The true alignment must be among the reported HSPs.
        assert!(hsps.iter().any(|h| h.t_start == 5 && h.q_start == 0 && h.length == 20));
    }

    #[test]
    fn two_matches_on_different_diagonals() {
        // Query appears at two distinct non-repetitive target locations, so
        // the two HSPs live on different diagonals.
        let target = PackedSeq::from_ascii(
            b"AAAAGATTACACATGGAAAAAAAAAAAAGATTACACATGGTTTTT",
        );
        let query = PackedSeq::from_ascii(b"GATTACACATGG");
        let pat = SeedPattern::solid(8);
        let table = PosTable::build(&target, &pat, 1);
        let hsps = search(
            &table,
            &target,
            &query,
            &m(),
            &SearchParams {
                step: 1,
                hsp: HspParams { x_drop: 500, hsp_threshold: 500 },
                transitions: 0,
                entropy_threshold: None,
            },
        );
        let diagonals: std::collections::BTreeSet<_> =
            hsps.iter().map(|h| h.diagonal()).collect();
        assert!(diagonals.len() >= 2, "expected ≥2 distinct diagonals, got {hsps:#?}");
    }
}
