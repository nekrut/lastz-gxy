//! Alignment edit scripts.
//!
//! An `EditScript` is a run-length-encoded sequence of edit operations
//! describing how target bases pair with query bases. It's CIGAR-compatible:
//!
//! - `Match` — both advance (covers both match and substitution columns;
//!   following upstream lastz, we do not distinguish `=`/`X` here)
//! - `InsertQuery` — query advances, target does not (a gap in the target
//!   relative to the query, CIGAR `I`)
//! - `DeleteQuery` — target advances, query does not (a gap in the query
//!   relative to the target, CIGAR `D`)
//!
//! Run-length encoding is a direct port of `edit_script.c` upstream.

use std::fmt;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    Match,
    InsertQuery,
    DeleteQuery,
}

impl EditOp {
    #[inline]
    pub fn cigar_char(self) -> char {
        match self {
            EditOp::Match => 'M',
            EditOp::InsertQuery => 'I',
            EditOp::DeleteQuery => 'D',
        }
    }

    #[inline]
    pub fn advances_target(self) -> bool {
        matches!(self, EditOp::Match | EditOp::DeleteQuery)
    }

    #[inline]
    pub fn advances_query(self) -> bool {
        matches!(self, EditOp::Match | EditOp::InsertQuery)
    }
}

/// A run-length-encoded edit script.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct EditScript {
    runs: Vec<(EditOp, u32)>,
}

impl EditScript {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_capacity(n: usize) -> Self {
        Self { runs: Vec::with_capacity(n) }
    }

    pub fn is_empty(&self) -> bool {
        self.runs.is_empty()
    }

    pub fn runs(&self) -> &[(EditOp, u32)] {
        &self.runs
    }

    /// Append `n` copies of `op`, merging with the trailing run when they
    /// share an op.
    pub fn push(&mut self, op: EditOp, n: u32) {
        if n == 0 {
            return;
        }
        if let Some(last) = self.runs.last_mut() {
            if last.0 == op {
                last.1 += n;
                return;
            }
        }
        self.runs.push((op, n));
    }

    /// Append another script in-place; preserves run merging at the join.
    pub fn extend(&mut self, other: &EditScript) {
        for &(op, n) in &other.runs {
            self.push(op, n);
        }
    }

    /// Reverse the script's run order — used when concatenating the
    /// left-extension (which is built back-to-front) with the right.
    pub fn reversed(&self) -> Self {
        Self {
            runs: self.runs.iter().rev().copied().collect(),
        }
    }

    /// Total column count (`Match + InsertQuery + DeleteQuery` lengths).
    pub fn columns(&self) -> u32 {
        self.runs.iter().map(|(_, n)| *n).sum()
    }

    /// Bases consumed on the target side.
    pub fn target_span(&self) -> u32 {
        self.runs
            .iter()
            .filter(|(op, _)| op.advances_target())
            .map(|(_, n)| *n)
            .sum()
    }

    /// Bases consumed on the query side.
    pub fn query_span(&self) -> u32 {
        self.runs
            .iter()
            .filter(|(op, _)| op.advances_query())
            .map(|(_, n)| *n)
            .sum()
    }

    /// Render as a CIGAR string (e.g. `20M3I7M`).
    pub fn to_cigar(&self) -> String {
        let mut out = String::with_capacity(self.runs.len() * 3);
        for (op, n) in &self.runs {
            use std::fmt::Write as _;
            let _ = write!(out, "{n}{c}", n = n, c = op.cigar_char());
        }
        out
    }
}

impl fmt::Display for EditScript {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_cigar())
    }
}

/// Render a pair of gapped alignment strings from the script plus the
/// underlying target/query base slices. Gaps are written as `-`.
pub fn render_aligned_bases(
    script: &EditScript,
    target: &[u8],
    query: &[u8],
) -> (Vec<u8>, Vec<u8>) {
    let mut t_out = Vec::with_capacity(script.columns() as usize);
    let mut q_out = Vec::with_capacity(script.columns() as usize);
    let mut ti = 0usize;
    let mut qi = 0usize;
    for &(op, n) in script.runs() {
        for _ in 0..n {
            match op {
                EditOp::Match => {
                    t_out.push(target[ti]);
                    q_out.push(query[qi]);
                    ti += 1;
                    qi += 1;
                }
                EditOp::InsertQuery => {
                    t_out.push(b'-');
                    q_out.push(query[qi]);
                    qi += 1;
                }
                EditOp::DeleteQuery => {
                    t_out.push(target[ti]);
                    q_out.push(b'-');
                    ti += 1;
                }
            }
        }
    }
    (t_out, q_out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn push_merges_adjacent_runs() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 3);
        s.push(EditOp::Match, 2);
        s.push(EditOp::InsertQuery, 1);
        s.push(EditOp::InsertQuery, 2);
        s.push(EditOp::Match, 1);
        assert_eq!(s.runs(), &[(EditOp::Match, 5), (EditOp::InsertQuery, 3), (EditOp::Match, 1)]);
    }

    #[test]
    fn spans_count_correctly() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 10);
        s.push(EditOp::InsertQuery, 3); // query has 3 extra bases, target static
        s.push(EditOp::Match, 5);
        s.push(EditOp::DeleteQuery, 2); // target has 2 extra bases
        s.push(EditOp::Match, 4);

        assert_eq!(s.target_span(), 10 + 5 + 2 + 4);
        assert_eq!(s.query_span(), 10 + 3 + 5 + 4);
        assert_eq!(s.columns(), 10 + 3 + 5 + 2 + 4);
    }

    #[test]
    fn cigar_roundtrip() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 20);
        s.push(EditOp::InsertQuery, 3);
        s.push(EditOp::Match, 7);
        assert_eq!(s.to_cigar(), "20M3I7M");
    }

    #[test]
    fn reversed_runs_flip_order() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 3);
        s.push(EditOp::DeleteQuery, 2);
        s.push(EditOp::Match, 5);
        let r = s.reversed();
        assert_eq!(
            r.runs(),
            &[(EditOp::Match, 5), (EditOp::DeleteQuery, 2), (EditOp::Match, 3)]
        );
    }

    #[test]
    fn render_aligned_bases_inserts_gap_chars() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 3);
        s.push(EditOp::InsertQuery, 2);
        s.push(EditOp::Match, 2);
        let (t, q) = render_aligned_bases(&s, b"AAACC", b"AAAGGCC");
        assert_eq!(t, b"AAA--CC");
        assert_eq!(q, b"AAAGGCC");
    }

    #[test]
    fn render_aligned_bases_handles_delete() {
        let mut s = EditScript::new();
        s.push(EditOp::Match, 3);
        s.push(EditOp::DeleteQuery, 1);
        s.push(EditOp::Match, 1);
        let (t, q) = render_aligned_bases(&s, b"AAAGC", b"AAAC");
        assert_eq!(t, b"AAAGC");
        assert_eq!(q, b"AAA-C");
    }
}
