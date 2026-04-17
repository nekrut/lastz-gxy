//! Reference and query sequence loading.
//!
//! v1 supports ASCII FASTA via `noodles-fasta`. 2bit and HSX are tracked for
//! Phase 2 (PLAN.md §3.3). Internally sequences are held as `PackedSeq`: a
//! 2-bit-encoded nucleotide buffer plus a parallel `valid` bitset so `N`
//! regions can be skipped without polluting the nucleotide stream.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

use crate::dna::{decode_base, encode_base};

/// One sequence from a FASTA: name + packed nucleotide buffer.
#[derive(Debug, Clone)]
pub struct Sequence {
    pub name: String,
    pub seq: PackedSeq,
}

impl Sequence {
    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// 2-bit-packed DNA with a parallel validity bitset.
///
/// Layout is little-endian within each byte: position `i` of the nucleotide
/// code lives in bits `(i%4)*2 .. (i%4)*2+2` of `codes[i/4]`. Validity lives
/// in `valid[i/8]`, bit `i%8`.
#[derive(Debug, Clone)]
pub struct PackedSeq {
    codes: Vec<u8>,
    valid: Vec<u8>,
    len: usize,
}

impl PackedSeq {
    pub fn new() -> Self {
        Self {
            codes: Vec::new(),
            valid: Vec::new(),
            len: 0,
        }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            codes: Vec::with_capacity((n + 3) / 4),
            valid: Vec::with_capacity((n + 7) / 8),
            len: 0,
        }
    }

    /// Build a `PackedSeq` from an ASCII nucleotide slice.
    pub fn from_ascii(bases: &[u8]) -> Self {
        let mut out = Self::with_capacity(bases.len());
        for &b in bases {
            out.push_ascii(b);
        }
        out
    }

    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    pub fn push_ascii(&mut self, b: u8) {
        let (code, valid) = encode_base(b);
        let bit_off = (self.len & 3) * 2;
        if bit_off == 0 {
            self.codes.push(0);
        }
        let last = self.codes.last_mut().unwrap();
        *last |= (code & 0b11) << bit_off;

        let vbit = self.len & 7;
        if vbit == 0 {
            self.valid.push(0);
        }
        if valid {
            let vlast = self.valid.last_mut().unwrap();
            *vlast |= 1 << vbit;
        }
        self.len += 1;
    }

    /// 2-bit code at position `i`. The value is meaningful only when
    /// `is_valid(i)` is true.
    #[inline]
    pub fn code(&self, i: usize) -> u8 {
        let bit_off = (i & 3) * 2;
        (self.codes[i >> 2] >> bit_off) & 0b11
    }

    #[inline]
    pub fn is_valid(&self, i: usize) -> bool {
        (self.valid[i >> 3] >> (i & 7)) & 1 != 0
    }

    /// Decode back to uppercase ASCII. Invalid positions become `N`.
    pub fn to_ascii(&self) -> Vec<u8> {
        (0..self.len)
            .map(|i| if self.is_valid(i) { decode_base(self.code(i)) } else { b'N' })
            .collect()
    }

    /// Reverse-complement into a new `PackedSeq`.
    pub fn reverse_complement(&self) -> Self {
        let mut out = Self::with_capacity(self.len);
        for i in (0..self.len).rev() {
            if self.is_valid(i) {
                let c = (!self.code(i)) & 0b11;
                out.push_ascii(decode_base(c));
            } else {
                out.push_ascii(b'N');
            }
        }
        out
    }
}

impl Default for PackedSeq {
    fn default() -> Self {
        Self::new()
    }
}

/// Load all sequences from a FASTA file on disk.
pub fn load_fasta<P: AsRef<Path>>(path: P) -> Result<Vec<Sequence>, SeqError> {
    let file = File::open(path.as_ref()).map_err(|e| SeqError::Io {
        path: path.as_ref().display().to_string(),
        source: e,
    })?;
    parse_fasta(BufReader::new(file))
}

/// Parse a FASTA from any buffered reader. Whitespace is skipped; only the
/// first whitespace-delimited token of each header line is kept as the
/// sequence name, matching samtools/noodles convention.
pub fn parse_fasta<R: BufRead>(reader: R) -> Result<Vec<Sequence>, SeqError> {
    let mut seqs = Vec::new();
    let mut cur_name: Option<String> = None;
    let mut cur_seq = PackedSeq::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| SeqError::Io {
            path: format!("line {}", lineno + 1),
            source: e,
        })?;
        if line.is_empty() {
            continue;
        }
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(name) = cur_name.take() {
                seqs.push(Sequence {
                    name,
                    seq: std::mem::take(&mut cur_seq),
                });
            }
            let name = rest
                .split_whitespace()
                .next()
                .ok_or(SeqError::EmptyHeader { line: lineno + 1 })?
                .to_string();
            cur_name = Some(name);
        } else {
            for &b in line.as_bytes() {
                if b.is_ascii_whitespace() {
                    continue;
                }
                cur_seq.push_ascii(b);
            }
        }
    }
    if let Some(name) = cur_name {
        seqs.push(Sequence { name, seq: cur_seq });
    }
    Ok(seqs)
}

#[derive(Debug, Error)]
pub enum SeqError {
    #[error("I/O error ({path}): {source}")]
    Io {
        path: String,
        #[source]
        source: std::io::Error,
    },
    #[error("line {line}: FASTA header '>' without a sequence name")]
    EmptyHeader { line: usize },
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packed_roundtrip_ascii() {
        let input = b"ACGTACGTNACG";
        let packed = PackedSeq::from_ascii(input);
        assert_eq!(packed.len(), input.len());
        let out = packed.to_ascii();
        assert_eq!(out, b"ACGTACGTNACG");
    }

    #[test]
    fn packed_reports_invalid_for_n() {
        let packed = PackedSeq::from_ascii(b"ANCGT");
        assert!(packed.is_valid(0));
        assert!(!packed.is_valid(1));
        assert!(packed.is_valid(2));
    }

    #[test]
    fn packed_reverse_complement_is_involution() {
        let s = PackedSeq::from_ascii(b"ACGTACGT");
        let rc = s.reverse_complement();
        assert_eq!(rc.to_ascii(), b"ACGTACGT");
        let rcrc = rc.reverse_complement();
        assert_eq!(rcrc.to_ascii(), s.to_ascii());
    }

    #[test]
    fn packed_reverse_complement_random() {
        let s = PackedSeq::from_ascii(b"ACCGTTNAG");
        let rc = s.reverse_complement();
        assert_eq!(rc.to_ascii(), b"CTNAACGGT");
    }

    #[test]
    fn parse_simple_fasta() {
        let text = ">chr1 something\nACGT\nACGT\n>chr2\nNNNN\n";
        let seqs = parse_fasta(text.as_bytes()).unwrap();
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0].name, "chr1");
        assert_eq!(seqs[0].seq.to_ascii(), b"ACGTACGT");
        assert_eq!(seqs[1].name, "chr2");
        assert_eq!(seqs[1].seq.to_ascii(), b"NNNN");
        assert!(!seqs[1].seq.is_valid(0));
    }

    #[test]
    fn parse_fasta_handles_blank_lines_and_case() {
        let text = ">a\n\nacgt\nACGT\n";
        let seqs = parse_fasta(text.as_bytes()).unwrap();
        assert_eq!(seqs.len(), 1);
        assert_eq!(seqs[0].seq.to_ascii(), b"ACGTACGT");
    }
}
