//! Per-worker diagonal deduplication.
//!
//! When many seed hits fall on the same diagonal inside an already-extended
//! HSP, upstream lastz avoids re-extending every one of them. Each diagonal
//! remembers the rightmost target position already covered by an HSP; any
//! seed hit on that diagonal with `t_pos <= covered_end` is suppressed.
//!
//! This deliberately mirrors the logic in upstream `diag_hash.c`. Correctness
//! (and parity) requires that we feed seed hits in non-decreasing `t_pos`
//! order, which is what `pos_table::lookup` already provides when queried in
//! query-position order.

use rustc_hash::FxHashMap;

/// Diagonal tracker: remembers, per diagonal `d = t - q`, the smallest
/// `t_pos` beyond which new seed hits on that diagonal are "fresh" and need
/// to be extended.
#[derive(Debug, Default, Clone)]
pub struct DiagHash {
    /// `d -> next_t_pos_to_consider`.
    next: FxHashMap<i64, u32>,
}

impl DiagHash {
    pub fn new() -> Self {
        Self::default()
    }

    /// Should the seed hit at `(t_pos, q_pos)` be extended?
    ///
    /// Returns `true` iff there is no prior HSP on this diagonal whose right
    /// edge already covers `t_pos`.
    #[inline]
    pub fn admit(&self, t_pos: u32, q_pos: u32) -> bool {
        let d = t_pos as i64 - q_pos as i64;
        self.next.get(&d).map_or(true, |&n| t_pos >= n)
    }

    /// Mark diagonal `d = t_end - q_end` as covered out to `t_end` (inclusive
    /// of the final matched column, so the next non-dominated hit must have
    /// `t_pos >= t_end`).
    #[inline]
    pub fn mark_covered(&mut self, t_end: u32, q_end: u32) {
        let d = t_end as i64 - q_end as i64;
        self.next
            .entry(d)
            .and_modify(|e| {
                if t_end > *e {
                    *e = t_end;
                }
            })
            .or_insert(t_end);
    }

    /// Clear the table; useful between query sequences on the same worker.
    pub fn clear(&mut self) {
        self.next.clear();
    }

    pub fn len(&self) -> usize {
        self.next.len()
    }
    pub fn is_empty(&self) -> bool {
        self.next.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn admits_first_hit_on_empty_diagonal() {
        let h = DiagHash::new();
        assert!(h.admit(10, 5));
    }

    #[test]
    fn suppresses_dominated_hit_on_same_diagonal() {
        let mut h = DiagHash::new();
        h.mark_covered(100, 50); // diagonal 50, covered out to t=100
        assert!(!h.admit(80, 30)); // diagonal 50, t=80 dominated
        assert!(h.admit(200, 150)); // still diagonal 50 but past the cover
    }

    #[test]
    fn neighbouring_diagonals_are_independent() {
        let mut h = DiagHash::new();
        h.mark_covered(100, 50);
        // Diagonal 51 is untouched.
        assert!(h.admit(100, 49));
    }

    #[test]
    fn mark_covered_monotone_on_same_diagonal() {
        let mut h = DiagHash::new();
        h.mark_covered(200, 100);
        h.mark_covered(150, 50); // same diagonal (100), earlier t_end
        // Should not regress: 200 still the dominating edge.
        assert!(!h.admit(199, 99));
        assert!(h.admit(200, 100));
    }
}
