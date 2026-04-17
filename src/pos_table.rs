//! Reference-side seed position index.
//!
//! Laid out as struct-of-arrays per PLAN.md §3.3:
//!
//! - `positions`: `Vec<u32>` — absolute reference positions, packed by seed
//!   word in ascending word order, then ascending position within each word.
//! - `starts`: `Vec<u32>` — `starts[w]` is the offset of word `w`'s run in
//!   `positions`; `starts[w+1] - starts[w]` is its multiplicity.
//!
//! This replaces upstream's linked-list-per-bucket (`pos_table.c`) layout and
//! gives the HSP stage contiguous, prefetch-friendly reads once the hot
//! bucket is located. `starts` is `num_words + 1` entries so the terminator
//! lets us compute any bucket's slice without a branch.
//!
//! "Radix bucketing" from the plan is implicit: indexing by seed word *is*
//! the radix, and because words are densely indexed (every word up to
//! `num_words = 1 << bit_width`), lookup is a single offset load.

use crate::seeds::{SeedExtractor, SeedPattern};
use crate::sequences::PackedSeq;

/// A built seed-position index for one reference sequence.
#[derive(Debug, Clone)]
pub struct PosTable {
    pattern: SeedPattern,
    positions: Vec<u32>,
    starts: Vec<u32>,
}

impl PosTable {
    /// Build an index over `seq` using `pattern` at stride `step`.
    ///
    /// `step` mirrors upstream's `--step=N` flag: only positions
    /// `0, step, 2*step, ...` are indexed. `step=1` indexes every window.
    pub fn build(seq: &PackedSeq, pattern: &SeedPattern, step: usize) -> Self {
        assert!(step >= 1, "step must be >= 1");
        let extractor = SeedExtractor::new(pattern);
        let num_words = pattern.num_words() as usize;

        // Two-pass radix bucket build. First pass: count occurrences per
        // word. Second pass: scatter positions into their buckets.
        let mut counts = vec![0u32; num_words];
        let mut hits: Vec<(u32, u32)> = extractor
            .iter_with_stride(seq, step)
            .map(|(pos, word)| {
                counts[word as usize] += 1;
                (pos, word as u32)
            })
            .collect();

        // Prefix-sum into `starts`. `starts[num_words]` is the total count.
        let mut starts = vec![0u32; num_words + 1];
        let mut acc: u32 = 0;
        for (i, c) in counts.iter().enumerate() {
            starts[i] = acc;
            acc = acc.checked_add(*c).expect("pos_table overflow (> 4 Gbp)");
        }
        starts[num_words] = acc;

        // Scatter. Reuse `counts` as per-bucket write cursors by rewinding
        // to the prefix-sum base (avoids a second allocation).
        let mut write_cursor = counts; // now holds remaining count per bucket
        let mut out = vec![0u32; acc as usize];

        // We need cursors that start at `starts[w]` and advance. Rebuild:
        for (i, cur) in write_cursor.iter_mut().enumerate() {
            *cur = starts[i];
        }
        // `hits` is in position-ascending order because it was produced by
        // iterating the sequence left-to-right; scattering preserves that.
        for &(pos, word) in &hits {
            let slot = &mut write_cursor[word as usize];
            out[*slot as usize] = pos;
            *slot += 1;
        }
        hits.clear();

        Self {
            pattern: pattern.clone(),
            positions: out,
            starts,
        }
    }

    pub fn pattern(&self) -> &SeedPattern {
        &self.pattern
    }

    /// Number of distinct seed words the underlying pattern supports.
    pub fn num_words(&self) -> usize {
        self.pattern.num_words() as usize
    }

    /// Total number of indexed positions across all words.
    pub fn total_positions(&self) -> usize {
        self.positions.len()
    }

    /// Reference positions that produced seed word `word`, in ascending
    /// order. Returns an empty slice if `word` has no hits.
    #[inline]
    pub fn lookup(&self, word: u64) -> &[u32] {
        let w = word as usize;
        debug_assert!(w < self.num_words(), "word {w} out of range");
        let start = self.starts[w] as usize;
        let end = self.starts[w + 1] as usize;
        &self.positions[start..end]
    }

    /// Drop words whose multiplicity exceeds `limit`. This is the "hot seed"
    /// filter lastz controls with `--maxwordcount`; it prevents the HSP
    /// stage from being dominated by repetitive k-mers.
    pub fn prune_hot_words(&mut self, limit: u32) {
        // Compact: for each word, if count > limit, drop its positions.
        // Rebuild starts/positions in place.
        let num_words = self.num_words();
        let mut new_positions: Vec<u32> =
            Vec::with_capacity(self.positions.len().min(1024));
        let mut new_starts = vec![0u32; num_words + 1];
        let mut acc: u32 = 0;
        for w in 0..num_words {
            new_starts[w] = acc;
            let slice = {
                let s = self.starts[w] as usize;
                let e = self.starts[w + 1] as usize;
                &self.positions[s..e]
            };
            if (slice.len() as u32) <= limit {
                new_positions.extend_from_slice(slice);
                acc += slice.len() as u32;
            }
        }
        new_starts[num_words] = acc;
        self.positions = new_positions;
        self.starts = new_starts;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn indexes_solid_kmer() {
        let seq = PackedSeq::from_ascii(b"ACGTACGT");
        let pat = SeedPattern::solid(3);
        let table = PosTable::build(&seq, &pat, 1);

        // Windows: pos 0 ACG, 1 CGT, 2 GTA, 3 TAC, 4 ACG, 5 CGT.
        // So "ACG" -> [0, 4], "CGT" -> [1, 5].
        let ext = SeedExtractor::new(&pat);
        let acg = ext.word_at(&seq, 0).unwrap();
        let cgt = ext.word_at(&seq, 1).unwrap();
        let gta = ext.word_at(&seq, 2).unwrap();
        let tac = ext.word_at(&seq, 3).unwrap();

        assert_eq!(table.lookup(acg), &[0, 4]);
        assert_eq!(table.lookup(cgt), &[1, 5]);
        assert_eq!(table.lookup(gta), &[2]);
        assert_eq!(table.lookup(tac), &[3]);
    }

    #[test]
    fn positions_are_ascending_per_word() {
        let seq = PackedSeq::from_ascii(b"ACGACGACGACG");
        let pat = SeedPattern::solid(3);
        let table = PosTable::build(&seq, &pat, 1);
        for w in 0..table.num_words() {
            let slice = table.lookup(w as u64);
            for pair in slice.windows(2) {
                assert!(pair[0] < pair[1]);
            }
        }
    }

    #[test]
    fn stride_skips_positions() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGTACGT");
        let pat = SeedPattern::solid(3);
        let table_step1 = PosTable::build(&seq, &pat, 1);
        let table_step3 = PosTable::build(&seq, &pat, 3);

        assert!(table_step3.total_positions() < table_step1.total_positions());

        let acg = SeedExtractor::new(&pat).word_at(&seq, 0).unwrap();
        assert_eq!(table_step3.lookup(acg), &[0, 12]);
    }

    #[test]
    fn skips_windows_with_n() {
        let seq = PackedSeq::from_ascii(b"ACGNACG");
        let pat = SeedPattern::solid(3);
        let table = PosTable::build(&seq, &pat, 1);
        let acg = SeedExtractor::new(&pat).word_at(&seq, 0).unwrap();
        // Only position 0 and 4 produce "ACG"; windows containing N skipped.
        assert_eq!(table.lookup(acg), &[0, 4]);
    }

    #[test]
    fn prune_hot_words_drops_repeats() {
        let seq = PackedSeq::from_ascii(b"AAAAAAAAAAACG");
        let pat = SeedPattern::solid(3);
        let mut table = PosTable::build(&seq, &pat, 1);
        let aaa = SeedExtractor::new(&pat).word_at(&seq, 0).unwrap();
        assert!(table.lookup(aaa).len() >= 9);
        table.prune_hot_words(4);
        assert!(table.lookup(aaa).is_empty());
        // But ACG survives (multiplicity 1).
        let acg = SeedExtractor::new(&pat).word_at(&seq, 10).unwrap();
        assert_eq!(table.lookup(acg), &[10]);
    }

    #[test]
    fn starts_array_is_monotone() {
        let seq = PackedSeq::from_ascii(b"ACGTACGTACGT");
        let pat = SeedPattern::solid(4);
        let table = PosTable::build(&seq, &pat, 1);
        for pair in table.starts.windows(2) {
            assert!(pair[0] <= pair[1]);
        }
        assert_eq!(
            table.starts[table.num_words()] as usize,
            table.total_positions()
        );
    }
}
