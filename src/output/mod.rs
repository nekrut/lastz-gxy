//! Alignment output writers.
//!
//! Phase 2 unifies gapped and ungapped alignments into a single `Record`
//! that carries an `EditScript`. An ungapped HSP is encoded as a script of
//! one `Match` run; a gapped alignment is the full DP traceback.

pub mod maf;
pub mod paf;

use crate::edit_script::EditScript;

/// One emitted alignment block. The raw `target_bases` / `query_bases` are
/// the un-gapped base slices; writers use `script` to render them with `-`
/// in the correct positions.
#[derive(Debug, Clone)]
pub struct Record {
    pub target_name: String,
    pub target_len: u32,
    pub query_name: String,
    pub query_len: u32,
    pub query_strand: Strand,
    pub t_start: u32,
    pub q_start: u32,
    pub t_span: u32,
    pub q_span: u32,
    pub score: i32,
    pub script: EditScript,
    pub target_bases: Vec<u8>,
    pub query_bases: Vec<u8>,
}

impl Record {
    pub fn t_end(&self) -> u32 {
        self.t_start + self.t_span
    }
    pub fn q_end(&self) -> u32 {
        self.q_start + self.q_span
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn as_char(self) -> char {
        match self {
            Strand::Plus => '+',
            Strand::Minus => '-',
        }
    }
}
