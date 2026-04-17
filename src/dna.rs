//! 2-bit DNA packing and related primitives.
//!
//! Encoding: A=0, C=1, G=2, T=3. N / ambiguous bases are flagged in a parallel
//! "valid" bitset; their nucleotide code is arbitrary (we use 0) and callers
//! must consult `is_valid()` before trusting the code at a position.

/// Encode one ASCII base into a 2-bit code + validity flag.
///
/// Returns `(code, is_valid)`. `code` is in `0..=3`. `is_valid` is `false` for
/// anything that isn't an unambiguous ACGT (upper or lower case).
#[inline]
pub fn encode_base(b: u8) -> (u8, bool) {
    match b {
        b'A' | b'a' => (0, true),
        b'C' | b'c' => (1, true),
        b'G' | b'g' => (2, true),
        b'T' | b't' | b'U' | b'u' => (3, true),
        _ => (0, false),
    }
}

/// Decode a 2-bit code back to its ASCII uppercase letter.
#[inline]
pub fn decode_base(code: u8) -> u8 {
    match code & 0b11 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

/// Complement a 2-bit code: A<->T, C<->G.
#[inline]
pub fn complement2(code: u8) -> u8 {
    (!code) & 0b11
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn acgt_roundtrip() {
        for (b, want) in [(b'A', 0), (b'C', 1), (b'G', 2), (b'T', 3)] {
            let (code, valid) = encode_base(b);
            assert!(valid);
            assert_eq!(code, want);
            assert_eq!(decode_base(code), b);
        }
    }

    #[test]
    fn case_insensitive() {
        assert_eq!(encode_base(b'a'), (0, true));
        assert_eq!(encode_base(b't'), (3, true));
    }

    #[test]
    fn ambiguous_bases_are_invalid() {
        for &b in b"NnRrYyMmKkSsWwBbDdHhVv?-" {
            let (_, valid) = encode_base(b);
            assert!(!valid, "{} should be invalid", b as char);
        }
    }

    #[test]
    fn complement_is_involution() {
        for c in 0..4u8 {
            assert_eq!(complement2(complement2(c)), c);
        }
        assert_eq!(complement2(0), 3); // A -> T
        assert_eq!(complement2(1), 2); // C -> G
    }
}
