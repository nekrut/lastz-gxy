//! Alignment output writers.
//!
//! Phase 1 ships MAF (ungapped HSPs only; gapped blocks arrive in Phase 2).
//! AXT, SAM, and PAF land alongside gapped extension in Phase 2.

pub mod maf;

use crate::hsp::Hsp;

/// A single emitted alignment block. Phase 1 populates only the HSP variant;
/// Phase 2 introduces a `Gapped` variant carrying an edit script.
#[derive(Debug, Clone)]
pub struct Record {
    pub target_name: String,
    pub target_len: u32,
    pub query_name: String,
    pub query_len: u32,
    pub query_strand: Strand,
    pub hsp: Hsp,
    pub target_bases: Vec<u8>,
    pub query_bases: Vec<u8>,
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
