//! Query-side seed iteration that drives HSP production.
//!
//! Given a reference `PosTable` and a query `PackedSeq`, this module walks
//! every query window at the configured stride, looks up the seed word, and
//! for each reference hit decides — via `DiagHash` — whether to extend into
//! an ungapped HSP. The resulting HSPs are returned in the order they were
//! emitted (which, because hits flow in query-order, is also
//! diagonal-grouped).

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
}

impl Default for SearchParams {
    fn default() -> Self {
        Self {
            step: 1,
            hsp: HspParams::default(),
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
    let extractor = SeedExtractor::new(pattern);

    let mut out = Vec::new();
    let mut diag = DiagHash::new();

    for (q_pos, word) in extractor.iter_with_stride(query, params.step) {
        let hits = table.lookup(word);
        for &t_pos in hits {
            if !diag.admit(t_pos, q_pos) {
                continue;
            }
            let Some(hsp) = extend_hit(target, query, matrix, &params.hsp, t_pos, q_pos, seed_len)
            else {
                continue;
            };
            diag.mark_covered(hsp.t_end(), hsp.q_end());
            out.push(hsp);
        }
    }

    // Upstream emits HSPs in (diagonal, t_start) order for downstream
    // chaining. We sort here so the output is stable regardless of the
    // interior hit order.
    out.sort_unstable_by_key(|h| (h.diagonal(), h.t_start));
    // A single extended HSP may have consumed multiple seed hits that the
    // diag-hash admitted before it was written out (since admission happens
    // before extension). Drop exact duplicates.
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
            },
        );
        assert_eq!(hsps.len(), 1, "expected one HSP, got {hsps:#?}");
        let h = hsps[0];
        assert_eq!(h.t_start, 5);
        assert_eq!(h.q_start, 0);
        assert_eq!(h.length, 16);
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
            },
        );
        let diagonals: std::collections::BTreeSet<_> =
            hsps.iter().map(|h| h.diagonal()).collect();
        assert!(diagonals.len() >= 2, "expected ≥2 distinct diagonals, got {hsps:#?}");
    }
}
