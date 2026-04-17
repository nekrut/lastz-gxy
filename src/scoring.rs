//! Scoring matrices and affine-gap parameters.
//!
//! v1 ships the default HOXD70 matrix compiled in. `ScoringMatrix::from_lastz_text`
//! parses the same `--scores` file format upstream uses (see
//! [`docs/scoring.html`](https://github.com/lastz/lastz/blob/master/docs/scoring.html)
//! in upstream lastz), so existing matrix files drop in unchanged.

use thiserror::Error;

use crate::dna::encode_base;

/// HOXD70 scoring matrix (Chiaromonte, Yap, Miller 2002). This is lastz's
/// default for mammalian-scale alignment and is the baseline for every parity
/// test in this crate.
///
/// Indices are [ref][query] using the 2-bit encoding (A=0, C=1, G=2, T=3).
const HOXD70: [[i32; 4]; 4] = [
    //    A     C     G     T
    [91, -114, -31, -123],  // A
    [-114, 100, -125, -31], // C
    [-31, -125, 100, -114], // G
    [-123, -31, -114, 91],  // T
];

/// Default HOXD70 affine-gap parameters as upstream lastz ships them.
pub const HOXD70_GAP_OPEN: i32 = 400;
pub const HOXD70_GAP_EXTEND: i32 = 30;

/// A 4x4 (A,C,G,T) substitution matrix plus affine-gap parameters. Scores for
/// unknown / ambiguous bases fall back to `ambig_score`.
#[derive(Debug, Clone)]
pub struct ScoringMatrix {
    sub: [[i32; 4]; 4],
    pub gap_open: i32,
    pub gap_extend: i32,
    /// Score returned when either base is not ACGT. Upstream's default is
    /// "treat as mismatch at the most negative on-diagonal score" which we
    /// encode as the minimum off-diagonal score in the matrix.
    pub ambig_score: i32,
}

impl ScoringMatrix {
    /// The upstream HOXD70 default.
    pub fn hoxd70() -> Self {
        let mut min_off_diag = i32::MAX;
        for i in 0..4 {
            for j in 0..4 {
                if i != j && HOXD70[i][j] < min_off_diag {
                    min_off_diag = HOXD70[i][j];
                }
            }
        }
        Self {
            sub: HOXD70,
            gap_open: HOXD70_GAP_OPEN,
            gap_extend: HOXD70_GAP_EXTEND,
            ambig_score: min_off_diag,
        }
    }

    /// Score a pair of 2-bit-encoded bases.
    #[inline]
    pub fn score(&self, ref_code: u8, query_code: u8) -> i32 {
        self.sub[(ref_code & 0b11) as usize][(query_code & 0b11) as usize]
    }

    /// Score a pair of ASCII bases; returns `ambig_score` if either side is
    /// not unambiguously ACGT.
    #[inline]
    pub fn score_ascii(&self, ref_b: u8, query_b: u8) -> i32 {
        let (rc, rv) = encode_base(ref_b);
        let (qc, qv) = encode_base(query_b);
        if rv && qv {
            self.score(rc, qc)
        } else {
            self.ambig_score
        }
    }

    /// Maximum substitution score in the matrix. Used by the HSP extender to
    /// bound x-drop termination.
    pub fn max_score(&self) -> i32 {
        *self.sub.iter().flatten().max().unwrap()
    }

    /// Minimum substitution score in the matrix.
    pub fn min_score(&self) -> i32 {
        *self.sub.iter().flatten().min().unwrap()
    }

    /// True iff every substitution score and the ambiguous-score fall
    /// inside `i8` range (`-128..=127`). The SIMD HSP extender needs this
    /// to safely pack the 4×4 matrix into a single `__m128i` lookup table.
    pub fn is_i8_safe(&self) -> bool {
        let in_range = |v: i32| v >= i8::MIN as i32 && v <= i8::MAX as i32;
        self.sub.iter().flatten().all(|&v| in_range(v)) && in_range(self.ambig_score)
    }

    /// Pack the 4×4 substitution matrix as a `[u8; 16]` lookup table indexed
    /// by `(ref_code << 2) | query_code`. Bytes are `i8`-encoded in a `u8`
    /// cast (lossless because the values are already in `i8` range — see
    /// [`is_i8_safe`]).
    pub fn pack_i8_lut(&self) -> [u8; 16] {
        debug_assert!(self.is_i8_safe(), "matrix is not i8-safe");
        let mut lut = [0u8; 16];
        for t in 0..4u8 {
            for q in 0..4u8 {
                let idx = ((t << 2) | q) as usize;
                lut[idx] = (self.sub[t as usize][q as usize] as i8) as u8;
            }
        }
        lut
    }

    /// Parse a lastz-style `--scores` file. The format is a leading header
    /// line of column bases followed by one row per base:
    ///
    /// ```text
    ///     A    C    G    T
    /// A  91 -114  -31 -123
    /// C -114  100 -125  -31
    /// G  -31 -125  100 -114
    /// T -123  -31 -114   91
    /// O = 400
    /// E = 30
    /// ```
    ///
    /// Lines beginning with `#` are comments; blank lines are skipped.
    pub fn from_lastz_text(text: &str) -> Result<Self, ScoringParseError> {
        let mut sub = [[0i32; 4]; 4];
        let mut sub_seen = [[false; 4]; 4];
        let mut col_order: Option<[usize; 4]> = None;
        let mut gap_open: Option<i32> = None;
        let mut gap_extend: Option<i32> = None;

        for (lineno, raw) in text.lines().enumerate() {
            let line = raw.split('#').next().unwrap().trim();
            if line.is_empty() {
                continue;
            }

            // Key=value form: "O = 400", "E = 30", "open = 400", etc.
            if let Some((key, value)) = line.split_once('=') {
                let k = key.trim().to_ascii_lowercase();
                let v: i32 = value
                    .trim()
                    .parse()
                    .map_err(|_| ScoringParseError::BadNumber {
                        line: lineno + 1,
                        text: value.trim().to_string(),
                    })?;
                match k.as_str() {
                    "o" | "open" | "gap_open" => gap_open = Some(v),
                    "e" | "extend" | "gap_extend" => gap_extend = Some(v),
                    _ => {
                        return Err(ScoringParseError::UnknownKey {
                            line: lineno + 1,
                            key: k,
                        });
                    }
                }
                continue;
            }

            let tokens: Vec<&str> = line.split_whitespace().collect();

            // Header: four ACGT tokens, no leading row label.
            if col_order.is_none() && tokens.len() == 4 && tokens.iter().all(|t| t.len() == 1) {
                let mut order = [0usize; 4];
                for (i, t) in tokens.iter().enumerate() {
                    let (code, valid) = encode_base(t.as_bytes()[0]);
                    if !valid {
                        return Err(ScoringParseError::BadHeaderBase {
                            line: lineno + 1,
                            base: t.as_bytes()[0] as char,
                        });
                    }
                    order[i] = code as usize;
                }
                col_order = Some(order);
                continue;
            }

            // Matrix row: <base> s00 s01 s02 s03
            if tokens.len() == 5 {
                let order = col_order.ok_or(ScoringParseError::MissingHeader)?;
                let row_base = tokens[0].as_bytes()[0];
                let (row_code, valid) = encode_base(row_base);
                if !valid {
                    return Err(ScoringParseError::BadHeaderBase {
                        line: lineno + 1,
                        base: row_base as char,
                    });
                }
                for (i, tok) in tokens[1..].iter().enumerate() {
                    let v: i32 = tok.parse().map_err(|_| ScoringParseError::BadNumber {
                        line: lineno + 1,
                        text: (*tok).to_string(),
                    })?;
                    let col = order[i];
                    sub[row_code as usize][col] = v;
                    sub_seen[row_code as usize][col] = true;
                }
                continue;
            }

            return Err(ScoringParseError::UnexpectedLine {
                line: lineno + 1,
                text: line.to_string(),
            });
        }

        if sub_seen.iter().flatten().any(|&s| !s) {
            return Err(ScoringParseError::IncompleteMatrix);
        }

        let mut min_off_diag = i32::MAX;
        for i in 0..4 {
            for j in 0..4 {
                if i != j && sub[i][j] < min_off_diag {
                    min_off_diag = sub[i][j];
                }
            }
        }

        Ok(Self {
            sub,
            gap_open: gap_open.unwrap_or(HOXD70_GAP_OPEN),
            gap_extend: gap_extend.unwrap_or(HOXD70_GAP_EXTEND),
            ambig_score: min_off_diag,
        })
    }
}

impl Default for ScoringMatrix {
    fn default() -> Self {
        Self::hoxd70()
    }
}

#[derive(Debug, Error)]
pub enum ScoringParseError {
    #[error("line {line}: unrecognised base '{base}' in matrix header")]
    BadHeaderBase { line: usize, base: char },
    #[error("line {line}: could not parse '{text}' as integer")]
    BadNumber { line: usize, text: String },
    #[error("line {line}: unknown key '{key}'")]
    UnknownKey { line: usize, key: String },
    #[error("line {line}: unexpected line '{text}'")]
    UnexpectedLine { line: usize, text: String },
    #[error("matrix header (ACGT column labels) missing")]
    MissingHeader,
    #[error("matrix is missing one or more cells")]
    IncompleteMatrix,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hoxd70_is_i8_safe_and_packs_correctly() {
        let m = ScoringMatrix::hoxd70();
        assert!(m.is_i8_safe());
        let lut = m.pack_i8_lut();
        assert_eq!(lut[(0 << 2) | 0] as i8, 91);
        assert_eq!(lut[(1 << 2) | 1] as i8, 100);
        assert_eq!(lut[(0 << 2) | 3] as i8, -123);
        assert_eq!(lut[(1 << 2) | 2] as i8, -125);
    }

    #[test]
    fn oversized_matrix_is_not_i8_safe() {
        let text = "\
                A    C    G    T\n\
            A   200 -114  -31 -123\n\
            C  -114  100 -125  -31\n\
            G   -31 -125  100 -114\n\
            T  -123  -31 -114   91\n";
        let m = ScoringMatrix::from_lastz_text(text).unwrap();
        assert!(!m.is_i8_safe());
    }

    #[test]
    fn hoxd70_diagonal_matches_upstream() {
        let m = ScoringMatrix::hoxd70();
        assert_eq!(m.score(0, 0), 91); // A/A
        assert_eq!(m.score(1, 1), 100); // C/C
        assert_eq!(m.score(2, 2), 100); // G/G
        assert_eq!(m.score(3, 3), 91); // T/T
    }

    #[test]
    fn hoxd70_transitions_softer_than_transversions() {
        let m = ScoringMatrix::hoxd70();
        // Transitions (A<->G, C<->T) score better than transversions.
        assert!(m.score(0, 2) > m.score(0, 1)); // A/G > A/C
        assert!(m.score(1, 3) > m.score(1, 2)); // C/T > C/G
    }

    #[test]
    fn hoxd70_symmetric() {
        let m = ScoringMatrix::hoxd70();
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.score(i, j), m.score(j, i), "asymmetric at ({i},{j})");
            }
        }
    }

    #[test]
    fn ambig_score_is_negative() {
        let m = ScoringMatrix::hoxd70();
        assert_eq!(m.score_ascii(b'N', b'A'), m.ambig_score);
        assert!(m.ambig_score < 0);
    }

    #[test]
    fn parses_hoxd70_text() {
        let text = "\
            # comment
                A    C    G    T\n\
            A   91 -114  -31 -123\n\
            C -114  100 -125  -31\n\
            G  -31 -125  100 -114\n\
            T -123  -31 -114   91\n\
            O = 400\n\
            E = 30\n";
        let m = ScoringMatrix::from_lastz_text(text).unwrap();
        let ref_m = ScoringMatrix::hoxd70();
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.score(i, j), ref_m.score(i, j));
            }
        }
        assert_eq!(m.gap_open, HOXD70_GAP_OPEN);
        assert_eq!(m.gap_extend, HOXD70_GAP_EXTEND);
    }

    #[test]
    fn parses_reordered_columns() {
        // Same matrix but columns in T,G,C,A order.
        let text = "\
                T    G    C    A\n\
            T   91 -114  -31 -123\n\
            G -114  100 -125  -31\n\
            C  -31 -125  100 -114\n\
            A -123  -31 -114   91\n";
        let m = ScoringMatrix::from_lastz_text(text).unwrap();
        let ref_m = ScoringMatrix::hoxd70();
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m.score(i, j), ref_m.score(i, j), "mismatch at ({i},{j})");
            }
        }
    }

    #[test]
    fn rejects_incomplete_matrix() {
        let text = "\
                A    C    G    T\n\
            A   91 -114  -31 -123\n";
        assert!(matches!(
            ScoringMatrix::from_lastz_text(text),
            Err(ScoringParseError::IncompleteMatrix)
        ));
    }
}
