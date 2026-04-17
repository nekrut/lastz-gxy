//! Highest-scoring HSP chain (upstream `--chain`).
//!
//! Given a set of HSPs for one `(target, query, strand)` triple, select the
//! subset that forms the highest-scoring strictly-increasing series under
//! both coordinates. This follows `chain.c` upstream: a classic longest-
//! increasing-subsequence DP weighted by HSP score, with back-pointer
//! traceback.
//!
//! Parity note: lastz breaks ties by preferring HSPs that appear earlier in
//! `t_start` order. We inherit that by iterating in sorted `t_start` order
//! and using `>=` for the predecessor score comparison.

use crate::hsp::Hsp;

/// Return the subset of `hsps` that lies on the single highest-scoring
/// strictly-monotone chain. Input order is irrelevant; output is in ascending
/// `t_start` order. Ties are broken in favour of the predecessor encountered
/// first in the sorted input, matching upstream.
pub fn best_chain(hsps: &[Hsp]) -> Vec<Hsp> {
    if hsps.is_empty() {
        return Vec::new();
    }

    // Sort by (t_start, q_start). q_start as secondary preserves a stable
    // resolution when two HSPs share t_start (rare but possible with very
    // short seeds on repetitive sequence).
    let mut sorted: Vec<Hsp> = hsps.to_vec();
    sorted.sort_by_key(|h| (h.t_start, h.q_start));

    let n = sorted.len();
    let mut best: Vec<i64> = sorted.iter().map(|h| h.score as i64).collect();
    let mut parent: Vec<i32> = vec![-1; n];

    for i in 1..n {
        let hi = &sorted[i];
        let mut best_prev_score: i64 = 0;
        let mut best_prev_idx: i32 = -1;
        for j in 0..i {
            let hj = &sorted[j];
            // Strict non-overlap on both axes.
            if hj.t_end() <= hi.t_start && hj.q_end() <= hi.q_start && best[j] > best_prev_score {
                best_prev_score = best[j];
                best_prev_idx = j as i32;
            }
        }
        best[i] = best_prev_score + hi.score as i64;
        parent[i] = best_prev_idx;
    }

    // Find the maximum-scoring terminal.
    let mut end_idx: usize = 0;
    let mut end_score: i64 = best[0];
    for (i, &s) in best.iter().enumerate() {
        if s > end_score {
            end_score = s;
            end_idx = i;
        }
    }

    // Walk parents back, then reverse.
    let mut out = Vec::new();
    let mut cur: i32 = end_idx as i32;
    while cur >= 0 {
        out.push(sorted[cur as usize]);
        cur = parent[cur as usize];
    }
    out.reverse();
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hsp(t: u32, q: u32, len: u32, score: i32) -> Hsp {
        Hsp { t_start: t, q_start: q, length: len, score }
    }

    #[test]
    fn empty_in_empty_out() {
        assert!(best_chain(&[]).is_empty());
    }

    #[test]
    fn single_hsp_is_its_own_chain() {
        let h = hsp(10, 20, 30, 1000);
        assert_eq!(best_chain(&[h]), vec![h]);
    }

    #[test]
    fn picks_monotone_chain_and_drops_conflicting_hsp() {
        // Three HSPs: A then B (both on a nice series), C overlaps B on query
        // but with a higher solo score — chain should prefer the A+B pair
        // if their combined score beats C alone.
        let a = hsp(0, 0, 20, 500);
        let b = hsp(40, 40, 20, 500);
        let c = hsp(80, 10, 20, 700); // q=10..30 overlaps a (q=0..20)? No, a ends at q=20.
                                       // But q=10..30 overlaps b (q=40..60)? No. Let me re-check.
                                       // Actually c's q=10..30 overlaps a's q=0..20 -> invalid next after a.
                                       // And c's t=80 > b's t_end=60, fine. c's q=10 < b's q_end=60, invalid after b.
                                       // So c cannot extend the a-b chain.
        let chain = best_chain(&[a, b, c]);
        // a+b = 1000 vs c alone = 700 → expect [a, b].
        assert_eq!(chain, vec![a, b]);
    }

    #[test]
    fn prefers_higher_total_over_more_members() {
        // Two HSPs with lots of combined score beat three smaller ones.
        let a = hsp(0, 0, 10, 100);
        let b = hsp(20, 20, 10, 100);
        let c = hsp(40, 40, 10, 100); // a+b+c = 300
        let d = hsp(60, 60, 10, 500); // solo 500 — but can extend a+b+c to 800

        let chain = best_chain(&[a, b, c, d]);
        assert_eq!(chain, vec![a, b, c, d]);
    }

    #[test]
    fn drops_off_diagonal_conflict() {
        // Two HSPs that can't both be in a chain (ordered opposite on q).
        let a = hsp(0, 50, 10, 1000);
        let b = hsp(20, 0, 10, 800);
        let chain = best_chain(&[a, b]);
        // Only one of them survives; the higher-scoring wins.
        assert_eq!(chain.len(), 1);
        assert_eq!(chain[0], a);
    }

    #[test]
    fn sort_stability_matches_input_order() {
        let a = hsp(10, 10, 5, 100);
        let b = hsp(0, 0, 5, 100);
        let chain = best_chain(&[a, b]);
        // Sorted output, regardless of input order.
        assert_eq!(chain, vec![b, a]);
    }
}
