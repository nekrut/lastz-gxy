//! Reference and query sequence loading.
//!
//! v1 supports ASCII FASTA via `noodles-fasta`. 2bit and HSX are tracked for
//! Phase 2 (PLAN.md §3.3). Internally sequences are held as `PackedSeq`: a
//! 2-bit-encoded nucleotide buffer, a parallel `valid` bitset so `N` regions
//! can be skipped without polluting the nucleotide stream, and a parallel
//! `mask` bitset recording soft-masked positions (lowercase ACGT in the
//! input FASTA). Soft-masked positions are still valid for extension but
//! are excluded from seeding — the same semantics upstream lastz uses
//! without `--nomasking`.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use thiserror::Error;

use crate::dna::{decode_base, encode_base};

/// 256-entry unpack table: byte `b` maps to the four 2-bit codes packed
/// into it, one per output byte of the resulting `[u8; 4]`. Used by the
/// batch extractors on `PackedSeq` to avoid a scalar shift-and-mask
/// for every position during SIMD HSP extension.
const UNPACK_LUT: [[u8; 4]; 256] = {
    let mut lut = [[0u8; 4]; 256];
    let mut i = 0usize;
    while i < 256 {
        lut[i][0] = (i & 0b11) as u8;
        lut[i][1] = ((i >> 2) & 0b11) as u8;
        lut[i][2] = ((i >> 4) & 0b11) as u8;
        lut[i][3] = ((i >> 6) & 0b11) as u8;
        i += 1;
    }
    lut
};

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

/// 2-bit-packed DNA with parallel `valid` and `mask` bitsets.
///
/// - `codes`: nucleotide packed 4 per byte (little-endian within each byte).
/// - `valid[i/8] bit i%8`: is set when position `i` is unambiguous ACGT.
/// - `mask[i/8] bit i%8`:  is set when position `i` was soft-masked in the
///   input (lowercase `acgt`). `mask` is purely a seeding signal; it does
///   not affect scoring or extension.
#[derive(Debug, Clone)]
pub struct PackedSeq {
    codes: Vec<u8>,
    valid: Vec<u8>,
    mask: Vec<u8>,
    len: usize,
}

impl PackedSeq {
    pub fn new() -> Self {
        Self {
            codes: Vec::new(),
            valid: Vec::new(),
            mask: Vec::new(),
            len: 0,
        }
    }

    pub fn with_capacity(n: usize) -> Self {
        Self {
            codes: Vec::with_capacity((n + 3) / 4),
            valid: Vec::with_capacity((n + 7) / 8),
            mask: Vec::with_capacity((n + 7) / 8),
            len: 0,
        }
    }

    /// Build a `PackedSeq` from an ASCII nucleotide slice, preserving case
    /// as a soft-mask signal (lowercase `acgt` → masked).
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
            self.mask.push(0);
        }
        if valid {
            let vlast = self.valid.last_mut().unwrap();
            *vlast |= 1 << vbit;
        }
        // Soft-mask: ACGT that arrived as lowercase.
        if valid && b.is_ascii_lowercase() {
            let mlast = self.mask.last_mut().unwrap();
            *mlast |= 1 << vbit;
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

    /// Was position `i` soft-masked (lowercase) in the input?
    #[inline]
    pub fn is_masked(&self, i: usize) -> bool {
        (self.mask[i >> 3] >> (i & 7)) & 1 != 0
    }

    /// Clear all soft-mask bits. Used when the caller wants to bypass
    /// masking for a run (the `--nomasking` CLI flag).
    pub fn clear_masks(&mut self) {
        for byte in &mut self.mask {
            *byte = 0;
        }
    }

    /// Batch-extract up to 16 consecutive codes into `out[..n]` starting
    /// at position `start`. Returns the number of codes written, which
    /// is `min(16, self.len() - start)` (so the last chunk near EOF
    /// writes fewer than 16). Positions beyond the sequence are left
    /// unchanged in `out`.
    ///
    /// Uses a precomputed 256-entry LUT keyed by the packed byte so each
    /// source byte produces 4 output codes in a single table lookup plus
    /// a `u32::to_le_bytes` memcpy — materially faster than 16 separate
    /// `code()` calls on realistic HSP extensions. Handles the
    /// unaligned-start case by shifting the first load left by the
    /// intra-byte offset. Correctness is proptested vs per-position
    /// `code()` in `sequences::tests`.
    #[inline]
    pub fn codes_batch(&self, start: usize, out: &mut [u8; 16]) -> usize {
        let n = 16.min(self.len.saturating_sub(start));
        if n == 0 {
            return 0;
        }
        let off = start & 3;
        let mut base = start >> 2;
        let mut filled = 0usize;

        if off != 0 {
            // Emit the tail of the first byte until we hit a 4-aligned
            // boundary or run out of codes.
            let byte = self.codes[base];
            for k in off..4.min(off + n) {
                let code = (byte >> (k * 2)) & 0b11;
                out[filled] = code;
                filled += 1;
            }
            base += 1;
        }

        // Emit 4 codes at a time via the LUT while at least 4 remain.
        while filled + 4 <= n {
            let bytes = UNPACK_LUT[self.codes[base] as usize];
            out[filled..filled + 4].copy_from_slice(&bytes);
            base += 1;
            filled += 4;
        }

        // Tail: fewer than 4 codes left — use code() to finish.
        while filled < n {
            let pos = start + filled;
            let byte = self.codes[pos >> 2];
            out[filled] = (byte >> ((pos & 3) * 2)) & 0b11;
            filled += 1;
        }

        n
    }

    /// Batch-extract the validity bits for up to 16 consecutive positions
    /// starting at `start`. Returns a `u16` with bit `i` set when
    /// `start + i` is valid and `i < n`, where `n = min(16, len - start)`.
    /// Bits `>= n` are zero (so invalid positions beyond EOF read as
    /// invalid, which matches how callers want to treat padding).
    #[inline]
    pub fn valid_mask_batch(&self, start: usize) -> u16 {
        let n = 16.min(self.len.saturating_sub(start));
        if n == 0 {
            return 0;
        }
        let mut mask = 0u16;
        let mut i = 0usize;
        while i < n {
            let pos = start + i;
            let byte = self.valid[pos >> 3];
            // How many contiguous bits from `byte` lie within `start..start+n`?
            let bit_off = pos & 7;
            let take = (n - i).min(8 - bit_off);
            // `take` is in 1..=8; `1u16 << take` avoids the u8 overflow
            // when `take == 8` (byte-aligned full-byte read).
            let chunk = ((byte >> bit_off) as u16) & ((1u16 << take) - 1);
            mask |= chunk << i;
            i += take;
        }
        mask
    }

    /// Batch-extract the mask bits for up to 16 positions starting at
    /// `start`. Returns a `u16` in the same shape as `valid_mask_batch`:
    /// bit `i` set when `start + i` is soft-masked.
    #[inline]
    pub fn mask_bits_batch(&self, start: usize) -> u16 {
        let n = 16.min(self.len.saturating_sub(start));
        if n == 0 {
            return 0;
        }
        let mut out = 0u16;
        let mut i = 0usize;
        while i < n {
            let pos = start + i;
            let byte = self.mask[pos >> 3];
            let bit_off = pos & 7;
            let take = (n - i).min(8 - bit_off);
            // `take` is in 1..=8; `1u16 << take` avoids the u8 overflow
            // when `take == 8` (byte-aligned full-byte read).
            let chunk = ((byte >> bit_off) as u16) & ((1u16 << take) - 1);
            out |= chunk << i;
            i += take;
        }
        out
    }

    /// Decode back to uppercase ASCII. Invalid positions become `N`.
    pub fn to_ascii(&self) -> Vec<u8> {
        (0..self.len)
            .map(|i| if self.is_valid(i) { decode_base(self.code(i)) } else { b'N' })
            .collect()
    }

    /// Copy `self[start..end]` into a fresh `PackedSeq`, preserving code,
    /// validity, and mask bits. Used by within-target chunking so each
    /// chunk owns a standalone buffer with chunk-local coordinates.
    pub fn slice_to_new(&self, start: usize, end: usize) -> Self {
        assert!(start <= end && end <= self.len, "slice_to_new out of range");
        let mut out = Self::with_capacity(end - start);
        for i in start..end {
            // Respect validity explicitly so the destination pushes an `N`
            // where the source had one (which clears both the valid and
            // mask bits in a single push).
            if !self.is_valid(i) {
                out.push_ascii(b'N');
            } else {
                let code = self.code(i);
                let byte = decode_base(code);
                let byte = if self.is_masked(i) {
                    byte.to_ascii_lowercase()
                } else {
                    byte
                };
                out.push_ascii(byte);
            }
        }
        out
    }

    /// Reverse-complement into a new `PackedSeq`. The `mask` bitset is
    /// carried over so position `len-1-i` in the output is masked iff
    /// position `i` in the input was masked.
    pub fn reverse_complement(&self) -> Self {
        let mut out = Self::with_capacity(self.len);
        for i in (0..self.len).rev() {
            let was_masked = self.is_masked(i);
            let byte = if self.is_valid(i) {
                let c = (!self.code(i)) & 0b11;
                let ch = decode_base(c);
                if was_masked { ch.to_ascii_lowercase() } else { ch }
            } else {
                b'N'
            };
            out.push_ascii(byte);
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
        // First four (lowercase) are soft-masked, latter four are not.
        let s = &seqs[0].seq;
        for i in 0..4 {
            assert!(s.is_masked(i), "pos {i} should be masked");
        }
        for i in 4..8 {
            assert!(!s.is_masked(i), "pos {i} should not be masked");
        }
    }

    #[test]
    fn packed_records_soft_mask_for_lowercase() {
        let packed = PackedSeq::from_ascii(b"ACgtAC");
        assert!(!packed.is_masked(0));
        assert!(!packed.is_masked(1));
        assert!(packed.is_masked(2));
        assert!(packed.is_masked(3));
        assert!(!packed.is_masked(4));
        assert!(!packed.is_masked(5));
    }

    #[test]
    fn clear_masks_wipes_the_bitset() {
        let mut packed = PackedSeq::from_ascii(b"ACgtAC");
        assert!(packed.is_masked(2));
        packed.clear_masks();
        for i in 0..packed.len() {
            assert!(!packed.is_masked(i));
        }
    }

    #[test]
    fn slice_to_new_preserves_codes_validity_and_mask() {
        let s = PackedSeq::from_ascii(b"ACgtNACGT");
        let sl = s.slice_to_new(2, 7); // "gtNAC"
        assert_eq!(sl.len(), 5);
        assert_eq!(sl.to_ascii(), b"GTNAC");
        assert!(sl.is_masked(0)); // 'g' → mask carried
        assert!(sl.is_masked(1)); // 't' → mask carried
        assert!(!sl.is_valid(2)); // 'N' → invalid carried
        assert!(!sl.is_masked(3)); // 'A' → unmasked
        assert!(!sl.is_masked(4)); // 'C' → unmasked
    }

    fn make_random_packed(len: usize, seed: u64) -> PackedSeq {
        let mut s = seed;
        let bytes: Vec<u8> = (0..len)
            .map(|_| {
                s = s.wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                let r = (s >> 33) & 0x1F;
                match r {
                    0 => b'N',   // inject some Ns
                    1 => b'a',   // and some lowercase
                    2 => b'g',
                    _ => b"ACGT"[(r as usize) & 3],
                }
            })
            .collect();
        PackedSeq::from_ascii(&bytes)
    }

    #[test]
    fn codes_batch_matches_per_position_api() {
        // Build a mixed-case sequence long enough to exercise aligned
        // and unaligned starts, multiple full LUT loads, and a short
        // tail.
        let seq = make_random_packed(200, 42);
        let mut out = [0u8; 16];
        for start in 0..(seq.len()) {
            let n = seq.codes_batch(start, &mut out);
            let expected_n = 16.min(seq.len() - start);
            assert_eq!(n, expected_n, "wrong n at start={start}");
            for i in 0..n {
                assert_eq!(
                    out[i],
                    seq.code(start + i),
                    "code mismatch at start={start}, i={i}"
                );
            }
        }
    }

    #[test]
    fn valid_mask_batch_matches_per_position_api() {
        let seq = make_random_packed(200, 7);
        for start in 0..(seq.len()) {
            let mask = seq.valid_mask_batch(start);
            let n = 16.min(seq.len() - start);
            for i in 0..16 {
                let bit = (mask >> i) & 1 != 0;
                let expected = i < n && seq.is_valid(start + i);
                assert_eq!(bit, expected, "validity mismatch at start={start}, i={i}");
            }
        }
    }

    #[test]
    fn mask_bits_batch_matches_per_position_api() {
        let seq = make_random_packed(200, 13);
        for start in 0..(seq.len()) {
            let mask = seq.mask_bits_batch(start);
            let n = 16.min(seq.len() - start);
            for i in 0..16 {
                let bit = (mask >> i) & 1 != 0;
                let expected = i < n && seq.is_masked(start + i);
                assert_eq!(bit, expected, "mask mismatch at start={start}, i={i}");
            }
        }
    }

    #[test]
    fn batch_extract_past_end_returns_zero() {
        let seq = PackedSeq::from_ascii(b"ACGT");
        let mut out = [0u8; 16];
        assert_eq!(seq.codes_batch(4, &mut out), 0);
        assert_eq!(seq.valid_mask_batch(4), 0);
        assert_eq!(seq.mask_bits_batch(4), 0);
    }

    #[test]
    fn reverse_complement_carries_mask_over() {
        let packed = PackedSeq::from_ascii(b"ACgtAC"); // mask at 2,3
        let rc = packed.reverse_complement();
        // RC of "ACgtAC" is "GTacGT" (complement + reverse; mask moves too).
        // Original mask positions 2,3 → RC positions 2,3 (len 6).
        assert!(rc.is_masked(2));
        assert!(rc.is_masked(3));
        for i in [0, 1, 4, 5] {
            assert!(!rc.is_masked(i), "pos {i} of RC unexpectedly masked");
        }
    }
}
