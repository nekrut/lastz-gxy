//! Spaced seed patterns and seed-word packing.
//!
//! A seed pattern is a sequence of `1` (care) and `0` (don't-care) positions
//! scanned along a sequence. Each window of length `pattern.len()` that is
//! fully valid (ACGT) contributes one **seed word**: the 2-bit concatenation
//! of nucleotide codes at the "care" positions, read left-to-right, with the
//! leftmost care position occupying the most significant bits.
//!
//! `12of19` and `19of20` are the two patterns used by default in upstream
//! lastz; they are carried over byte-for-byte from `docs/seed_patterns.html`
//! in upstream.

use std::fmt;

use crate::sequences::PackedSeq;

/// A seed pattern: positions that are "care" (`true`) vs "don't-care"
/// (`false`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedPattern {
    mask: Vec<bool>,
    name: String,
}

impl SeedPattern {
    pub fn new(name: impl Into<String>, mask: Vec<bool>) -> Self {
        assert!(!mask.is_empty(), "seed pattern must be non-empty");
        assert!(
            mask.len() <= 32,
            "seed pattern length {} exceeds 32-position limit",
            mask.len()
        );
        Self { mask, name: name.into() }
    }

    /// Parse a textual pattern like `"1110100110010101111"` (lastz style)
    /// into a `SeedPattern`.
    pub fn parse(name: impl Into<String>, spec: &str) -> Result<Self, SeedPatternError> {
        let mut mask = Vec::with_capacity(spec.len());
        for c in spec.chars() {
            match c {
                '1' => mask.push(true),
                '0' => mask.push(false),
                _ => return Err(SeedPatternError::BadChar(c)),
            }
        }
        if mask.is_empty() {
            return Err(SeedPatternError::Empty);
        }
        Ok(Self::new(name, mask))
    }

    /// The 12-of-19 pattern. Default for lastz mammalian alignment.
    /// `1110100110010101111`
    pub fn twelve_of_nineteen() -> Self {
        Self::parse("12of19", "1110100110010101111").unwrap()
    }

    /// The 19-of-20 "match" pattern: every position is care except one.
    /// Used by lastz when high sensitivity with large seed length is needed.
    /// `11111111101111111111`
    pub fn nineteen_of_twenty() -> Self {
        Self::parse("19of20", "11111111101111111111").unwrap()
    }

    /// A solid k-mer pattern (all care positions).
    pub fn solid(k: usize) -> Self {
        Self::new(format!("match{k}"), vec![true; k])
    }

    pub fn len(&self) -> usize {
        self.mask.len()
    }
    pub fn is_empty(&self) -> bool {
        self.mask.is_empty()
    }
    pub fn name(&self) -> &str {
        &self.name
    }
    pub fn mask(&self) -> &[bool] {
        &self.mask
    }

    /// Number of "care" positions; also the seed word's bit-width / 2.
    pub fn weight(&self) -> usize {
        self.mask.iter().filter(|&&c| c).count()
    }

    /// Bit width of a seed word produced by this pattern (2 × weight).
    pub fn bit_width(&self) -> u32 {
        2 * self.weight() as u32
    }

    /// Number of distinct seed words this pattern can produce.
    pub fn num_words(&self) -> u64 {
        1u64 << self.bit_width()
    }
}

/// Iterator-style seed-word extractor over a `PackedSeq`.
///
/// Each call to `next_word_at(pos)` computes the seed word starting at `pos`,
/// or returns `None` if any care position in the window is outside the
/// sequence or hits an invalid base. A future SIMD implementation replaces
/// this body without changing the interface.
pub struct SeedExtractor<'a> {
    pattern: &'a SeedPattern,
    care_positions: Vec<usize>,
}

impl<'a> SeedExtractor<'a> {
    pub fn new(pattern: &'a SeedPattern) -> Self {
        let care_positions = pattern
            .mask()
            .iter()
            .enumerate()
            .filter_map(|(i, &c)| if c { Some(i) } else { None })
            .collect();
        Self { pattern, care_positions }
    }

    pub fn pattern(&self) -> &SeedPattern {
        self.pattern
    }

    /// Compute the seed word at `pos` in `seq`, or `None` if the window is
    /// out of range, contains any invalid base at a care position, or
    /// hits a soft-masked position at any of the pattern's care positions.
    /// This is the conservative reading of upstream's default — empirically
    /// it gives the highest parity on the `pseudocat × pseudopig` fixture
    /// (relaxing to seed-start-only drops Jaccard ~10 %).
    #[inline]
    pub fn word_at(&self, seq: &PackedSeq, pos: usize) -> Option<u64> {
        let end = pos.checked_add(self.pattern.len())?;
        if end > seq.len() {
            return None;
        }
        let mut word: u64 = 0;
        for &offset in &self.care_positions {
            let p = pos + offset;
            if !seq.is_valid(p) || seq.is_masked(p) {
                return None;
            }
            word = (word << 2) | seq.code(p) as u64;
        }
        Some(word)
    }

    /// Yield `(pos, seed_word)` for every window in `seq` at the given
    /// stride. Windows that straddle invalid bases are skipped silently.
    pub fn iter_with_stride<'b>(
        &'b self,
        seq: &'b PackedSeq,
        stride: usize,
    ) -> impl Iterator<Item = (u32, u64)> + 'b {
        let pat_len = self.pattern.len();
        let max_start = seq.len().saturating_sub(pat_len);
        (0..=max_start)
            .step_by(stride)
            .filter_map(move |pos| self.word_at(seq, pos).map(|w| (pos as u32, w)))
    }
}

#[derive(Debug, thiserror::Error)]
pub enum SeedPatternError {
    #[error("seed pattern character '{0}' is not '0' or '1'")]
    BadChar(char),
    #[error("seed pattern is empty")]
    Empty,
}

impl fmt::Display for SeedPattern {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for &c in &self.mask {
            f.write_str(if c { "1" } else { "0" })?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn twelve_of_nineteen_weights_are_correct() {
        let p = SeedPattern::twelve_of_nineteen();
        assert_eq!(p.len(), 19);
        assert_eq!(p.weight(), 12);
        assert_eq!(p.bit_width(), 24);
        assert_eq!(p.num_words(), 1 << 24);
    }

    #[test]
    fn nineteen_of_twenty_weights_are_correct() {
        let p = SeedPattern::nineteen_of_twenty();
        assert_eq!(p.len(), 20);
        assert_eq!(p.weight(), 19);
    }

    #[test]
    fn parse_round_trip() {
        let p = SeedPattern::parse("custom", "1101").unwrap();
        assert_eq!(format!("{p}"), "1101");
        assert_eq!(p.weight(), 3);
    }

    #[test]
    fn parse_rejects_bad_chars() {
        assert!(matches!(
            SeedPattern::parse("x", "11X0"),
            Err(SeedPatternError::BadChar('X'))
        ));
    }

    #[test]
    fn solid_seed_word_packs_care_positions() {
        // Solid 4-mer on "ACGT" should give (A,C,G,T) = (0,1,2,3)
        // packed as 0b00_01_10_11 = 0x1B = 27.
        let seq = PackedSeq::from_ascii(b"ACGT");
        let pat = SeedPattern::solid(4);
        let ext = SeedExtractor::new(&pat);
        assert_eq!(ext.word_at(&seq, 0), Some(0x1B));
    }

    #[test]
    fn spaced_seed_skips_dont_care() {
        // Pattern "1001" picks positions 0 and 3. On "ACGT" that's A,T = 0,3
        // packed as 0b00_11 = 0x3.
        let pat = SeedPattern::parse("test", "1001").unwrap();
        let ext = SeedExtractor::new(&pat);
        let seq = PackedSeq::from_ascii(b"ACGT");
        assert_eq!(ext.word_at(&seq, 0), Some(0b0011));
    }

    #[test]
    fn masked_bases_block_seed_windows() {
        let pat = SeedPattern::solid(3);
        let ext = SeedExtractor::new(&pat);
        // Middle base is lowercase → masked.
        let seq = PackedSeq::from_ascii(b"AcG");
        assert_eq!(ext.word_at(&seq, 0), None);
    }

    #[test]
    fn mask_only_matters_at_care_positions() {
        // Pattern "101" places a don't-care at offset 1. A masked base at
        // offset 1 should NOT block the seed.
        let pat = SeedPattern::parse("test", "101").unwrap();
        let ext = SeedExtractor::new(&pat);
        let seq = PackedSeq::from_ascii(b"AcG");
        assert!(ext.word_at(&seq, 0).is_some());
    }

    #[test]
    fn word_at_returns_none_past_end() {
        let pat = SeedPattern::solid(4);
        let ext = SeedExtractor::new(&pat);
        let seq = PackedSeq::from_ascii(b"ACG");
        assert_eq!(ext.word_at(&seq, 0), None);
    }

    #[test]
    fn iter_with_stride_skips_invalid_windows() {
        let pat = SeedPattern::solid(3);
        let ext = SeedExtractor::new(&pat);
        let seq = PackedSeq::from_ascii(b"ACGNNACGT");
        let positions: Vec<u32> = ext.iter_with_stride(&seq, 1).map(|(p, _)| p).collect();
        // Windows at 0 (ACG) are fine; 1 (CGN), 2 (GNN), 3 (NNA), 4 (NAC) all
        // contain N. Then 5 (ACG), 6 (CGT) are fine. 7 is off the end.
        assert_eq!(positions, vec![0, 5, 6]);
    }

    #[test]
    fn iter_with_stride_honours_stride() {
        let pat = SeedPattern::solid(2);
        let ext = SeedExtractor::new(&pat);
        let seq = PackedSeq::from_ascii(b"ACGTACGT");
        let positions: Vec<u32> = ext.iter_with_stride(&seq, 3).map(|(p, _)| p).collect();
        assert_eq!(positions, vec![0, 3, 6]);
    }
}
