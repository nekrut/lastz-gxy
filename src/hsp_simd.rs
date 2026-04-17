//! AVX2 implementation of ungapped x-drop HSP extension.
//!
//! Same algorithm as the scalar path in `hsp::extend_hit`, but column scores
//! are looked up 16-at-a-time via `PSHUFB` against a 16-byte `[i8; 16]`
//! packing of the 4×4 substitution matrix. The per-column prefix-sum and
//! x-drop check remain scalar (they are inherently sequential); the win
//! comes from batching the matrix lookup and the validity blend.
//!
//! Gated on:
//! - `target_arch = "x86_64"` at compile time
//! - `is_x86_feature_detected!("avx2")` at runtime
//! - `ScoringMatrix::is_i8_safe()` (values must fit in `i8`)
//!
//! Differential parity vs the scalar path is enforced by a proptest
//! (`hsp::tests::simd_matches_scalar`).

#![cfg(target_arch = "x86_64")]

use std::arch::x86_64::*;

use crate::hsp::{Hsp, HspParams, SideResult};
use crate::scoring::ScoringMatrix;
use crate::sequences::PackedSeq;

/// SIMD version of `hsp::extend_hit`. Caller guarantees AVX2 is available
/// and `matrix.is_i8_safe()`.
pub(crate) unsafe fn extend_hit_avx2(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    params: &HspParams,
    t_pos: u32,
    q_pos: u32,
    seed_len: u32,
) -> Option<Hsp> {
    let t_pos = t_pos as i64;
    let q_pos = q_pos as i64;
    let t_len = target.len() as i64;
    let q_len = query.len() as i64;
    let seed_len = seed_len as i64;

    if t_pos < 0 || q_pos < 0 || t_pos + seed_len > t_len || q_pos + seed_len > q_len {
        return None;
    }

    let lut_bytes = matrix.pack_i8_lut();
    let matrix_lut = _mm_loadu_si128(lut_bytes.as_ptr() as *const __m128i);
    let ambig_i8 = matrix.ambig_score as i8;

    // Score the seed footprint (small, scalar is fine).
    let mut seed_score: i32 = 0;
    for k in 0..seed_len {
        seed_score = seed_score.saturating_add(column_score_scalar(
            target,
            query,
            matrix,
            t_pos + k,
            q_pos + k,
        ));
    }

    let right = extend_side(
        target,
        query,
        matrix_lut,
        ambig_i8,
        matrix.ambig_score,
        params.x_drop,
        t_pos + seed_len,
        q_pos + seed_len,
        1,
    );
    let left = extend_side(
        target,
        query,
        matrix_lut,
        ambig_i8,
        matrix.ambig_score,
        params.x_drop,
        t_pos - 1,
        q_pos - 1,
        -1,
    );

    let score = seed_score
        .saturating_add(left.best_score)
        .saturating_add(right.best_score);
    if score < params.hsp_threshold {
        return None;
    }

    let t_start = (t_pos - left.best_extent as i64) as u32;
    let q_start = (q_pos - left.best_extent as i64) as u32;
    let length = (seed_len + left.best_extent as i64 + right.best_extent as i64) as u32;

    Some(Hsp { t_start, q_start, length, score })
}

#[inline]
fn column_score_scalar(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix: &ScoringMatrix,
    t: i64,
    q: i64,
) -> i32 {
    let (ti, qi) = (t as usize, q as usize);
    if !target.is_valid(ti) || !query.is_valid(qi) {
        matrix.ambig_score
    } else {
        matrix.score(target.code(ti), query.code(qi))
    }
}

#[allow(clippy::too_many_arguments)]
#[target_feature(enable = "avx2,sse4.1,ssse3")]
unsafe fn extend_side(
    target: &PackedSeq,
    query: &PackedSeq,
    matrix_lut: __m128i,
    ambig_i8: i8,
    ambig_i32: i32,
    x_drop: i32,
    start_t: i64,
    start_q: i64,
    direction: i64,
) -> SideResult {
    let t_len = target.len() as i64;
    let q_len = query.len() as i64;

    let max_extent: usize = if direction > 0 {
        if start_t >= t_len || start_q >= q_len {
            0
        } else {
            ((t_len - start_t).min(q_len - start_q)) as usize
        }
    } else if direction < 0 {
        if start_t < 0 || start_q < 0 {
            0
        } else {
            ((start_t + 1).min(start_q + 1)) as usize
        }
    } else {
        0
    };

    if max_extent == 0 {
        return SideResult { best_score: 0, best_extent: 0 };
    }

    let ambig_vec = _mm_set1_epi8(ambig_i8);
    let low4_mask = _mm_set1_epi8(0x0F);

    let mut running: i32 = 0;
    let mut best: i32 = 0;
    let mut best_extent: u32 = 0;
    let mut extent: u32 = 0;

    let mut offset = 0usize;
    while offset < max_extent {
        let chunk = (max_extent - offset).min(16);

        // Scalar extract up to 16 (code, validity) pairs; positions beyond
        // `chunk` are padding and marked invalid (their score contributes
        // `ambig_score` but they'll never be reached because the inner
        // scalar loop below only walks `chunk` steps).
        let mut t_codes = [0u8; 16];
        let mut q_codes = [0u8; 16];
        let mut invalid = [0u8; 16];
        for i in 0..chunk {
            let t = start_t + direction * (offset + i) as i64;
            let q = start_q + direction * (offset + i) as i64;
            let ti = t as usize;
            let qi = q as usize;
            t_codes[i] = target.code(ti);
            q_codes[i] = query.code(qi);
            if !target.is_valid(ti) || !query.is_valid(qi) {
                invalid[i] = 0xFF;
            }
        }
        for i in chunk..16 {
            invalid[i] = 0xFF;
        }

        let t_vec = _mm_loadu_si128(t_codes.as_ptr() as *const __m128i);
        let q_vec = _mm_loadu_si128(q_codes.as_ptr() as *const __m128i);
        // idx = (t << 2) | q.  `_mm_slli_epi16` shifts 16-bit lanes — the
        // low 4 bits we care about are untouched by the cross-byte
        // spillover, and the mask afterwards drops the garbage.
        let t_shifted = _mm_slli_epi16(t_vec, 2);
        let idx = _mm_or_si128(t_shifted, q_vec);
        let idx = _mm_and_si128(idx, low4_mask);
        let scores = _mm_shuffle_epi8(matrix_lut, idx);

        let invalid_vec = _mm_loadu_si128(invalid.as_ptr() as *const __m128i);
        let blended = _mm_blendv_epi8(scores, ambig_vec, invalid_vec);

        let mut deltas = [0i8; 16];
        _mm_storeu_si128(deltas.as_mut_ptr() as *mut __m128i, blended);

        for i in 0..chunk {
            extent += 1;
            let s = deltas[i] as i32;
            // Note: `deltas[i]` is already ambig_i32 at invalid positions.
            let _ = ambig_i32;
            running = running.saturating_add(s);
            if running > best {
                best = running;
                best_extent = extent;
            } else if best.saturating_sub(running) > x_drop {
                return SideResult { best_score: best, best_extent };
            }
        }

        offset += chunk;
    }

    SideResult { best_score: best, best_extent }
}
