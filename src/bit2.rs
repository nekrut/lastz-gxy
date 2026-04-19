//! UCSC [2bit](https://genome.ucsc.edu/FAQ/FAQformat.html#format7) input
//! reader.
//!
//! 2bit is the reference-sequence format most UCSC / Galaxy-scale lastz
//! pipelines run against: a tiny binary header plus 2-bit-packed DNA
//! with overlay arrays encoding `N`-runs and soft-masked intervals. We
//! parse it into the same `Sequence` type `load_fasta` produces so the
//! rest of the pipeline is format-agnostic.
//!
//! Binary layout (little-endian assumed; big-endian files are flagged as
//! an error rather than transparently byteswapped — none of the common
//! 2bit producers emit them):
//!
//! ```text
//! header:           u32 magic=0x1A412743, u32 version=0, u32 seq_count, u32 reserved
//! per-sequence index record:
//!     u8 name_size, [name_size] bytes name, u32 offset
//! per-sequence data at offset:
//!     u32 dna_size
//!     u32 n_block_count
//!     [n_block_count] u32 n_block_starts, [n_block_count] u32 n_block_sizes
//!     u32 mask_block_count
//!     [mask_block_count] u32 mask_block_starts, [mask_block_count] u32 mask_block_sizes
//!     u32 reserved
//!     [ceil(dna_size/4)] bytes packed DNA, 4 bases per byte,
//!         high-order bits first, code: T=00, C=01, A=10, G=11
//! ```
//!
//! This is *different* from our internal `PackedSeq` encoding (A=00,
//! C=01, G=10, T=11). We decode to ASCII with `N` overlays applied from
//! the N-block arrays and lowercase from the mask-block arrays, then
//! reuse `PackedSeq::from_ascii` so case → soft-mask mapping is shared
//! with the FASTA path.

use std::fs::File;
use std::io::Read;
use std::path::Path;

use thiserror::Error;

use crate::sequences::{PackedSeq, Sequence};

const TWOBIT_MAGIC_LE: u32 = 0x1A41_2743;
const TWOBIT_MAGIC_BE: u32 = 0x4327_411A;

#[derive(Debug, Error)]
pub enum Bit2Error {
    #[error("I/O: {0}")]
    Io(#[from] std::io::Error),
    #[error("not a 2bit file (magic = 0x{found:08x})")]
    BadMagic { found: u32 },
    #[error("unsupported 2bit version {0}")]
    BadVersion(u32),
    #[error("truncated 2bit file (needed {need} bytes at offset {at})")]
    Truncated { at: usize, need: usize },
    #[error("sequence name is not UTF-8")]
    NonUtf8Name,
}

/// Read every sequence from a 2bit file into a `Vec<Sequence>`. The file
/// is read in full rather than mmapped because each sequence is decoded
/// into an owned `PackedSeq` anyway — the mmap'd bytes would outlive
/// nothing.
pub fn load_2bit<P: AsRef<Path>>(path: P) -> Result<Vec<Sequence>, Bit2Error> {
    let mut file = File::open(path)?;
    let mut buf = Vec::new();
    file.read_to_end(&mut buf)?;
    parse_2bit(&buf)
}

/// Parse a 2bit file from its byte representation. Exposed for
/// round-trip testing; most callers want `load_2bit`.
pub fn parse_2bit(bytes: &[u8]) -> Result<Vec<Sequence>, Bit2Error> {
    // Peek at the raw magic bytes to detect endianness. UCSC-native LE
    // files store it as 0x1A 0x41 0x27 0x43 (LE u32 = 0x1A412743);
    // big-endian files store the reversed bytes (read as LE u32 =
    // 0x4327411A). Older pipelines still emit BE.
    if bytes.len() < 4 {
        return Err(Bit2Error::Truncated { at: 0, need: 4 });
    }
    let first_u32_le = u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
    let big_endian = match first_u32_le {
        TWOBIT_MAGIC_LE => false,
        TWOBIT_MAGIC_BE => true,
        other => return Err(Bit2Error::BadMagic { found: other }),
    };

    let mut r = Cursor::new(bytes, big_endian);
    let _magic = r.u32()?;
    let version = r.u32()?;
    if version != 0 {
        return Err(Bit2Error::BadVersion(version));
    }
    let seq_count = r.u32()? as usize;
    let _reserved = r.u32()?;

    // Index: (name, offset) per sequence.
    let mut index: Vec<(String, u32)> = Vec::with_capacity(seq_count);
    for _ in 0..seq_count {
        let name_size = r.u8()? as usize;
        let name_bytes = r.bytes(name_size)?.to_vec();
        let name = String::from_utf8(name_bytes).map_err(|_| Bit2Error::NonUtf8Name)?;
        let offset = r.u32()?;
        index.push((name, offset));
    }

    let mut out = Vec::with_capacity(seq_count);
    for (name, offset) in index {
        let seq = parse_sequence_record(bytes, offset as usize, big_endian)?;
        out.push(Sequence { name, seq });
    }
    Ok(out)
}

fn parse_sequence_record(bytes: &[u8], offset: usize, big_endian: bool) -> Result<PackedSeq, Bit2Error> {
    let mut r = Cursor::new(bytes, big_endian);
    r.seek(offset)?;
    let dna_size = r.u32()? as usize;
    let n_count = r.u32()? as usize;
    let n_starts: Vec<u32> = (0..n_count).map(|_| r.u32()).collect::<Result<_, _>>()?;
    let n_sizes: Vec<u32> = (0..n_count).map(|_| r.u32()).collect::<Result<_, _>>()?;
    let m_count = r.u32()? as usize;
    let m_starts: Vec<u32> = (0..m_count).map(|_| r.u32()).collect::<Result<_, _>>()?;
    let m_sizes: Vec<u32> = (0..m_count).map(|_| r.u32()).collect::<Result<_, _>>()?;
    let _reserved = r.u32()?;

    let packed_len = (dna_size + 3) / 4;
    let packed = r.bytes(packed_len)?.to_vec();

    // Decode packed DNA into an ASCII buffer. 2bit packs 4 bases per
    // byte high-order-first with codes T=00, C=01, A=10, G=11.
    const ACGT: [u8; 4] = [b'T', b'C', b'A', b'G'];
    let mut ascii = vec![b'A'; dna_size];
    for i in 0..dna_size {
        let byte = packed[i >> 2];
        let shift = 6 - (i & 3) * 2;
        let code = ((byte >> shift) & 0b11) as usize;
        ascii[i] = ACGT[code];
    }

    // Apply N-block overlays (these override any ACGT the packed-DNA
    // byte claimed; 2bit writers leave a placeholder code in that byte).
    for (s, sz) in n_starts.iter().zip(n_sizes.iter()) {
        let lo = *s as usize;
        let hi = lo + *sz as usize;
        for b in &mut ascii[lo..hi.min(dna_size)] {
            *b = b'N';
        }
    }

    // Apply soft-mask overlays (lowercase so PackedSeq::from_ascii sets
    // the mask bit, matching the FASTA path).
    for (s, sz) in m_starts.iter().zip(m_sizes.iter()) {
        let lo = *s as usize;
        let hi = lo + *sz as usize;
        for b in &mut ascii[lo..hi.min(dna_size)] {
            if *b != b'N' {
                *b = b.to_ascii_lowercase();
            }
        }
    }

    Ok(PackedSeq::from_ascii(&ascii))
}

/// Minimal cursor over a byte slice. We hand-roll this instead of
/// pulling in `byteorder` to keep the dependency tree small.
struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
    big_endian: bool,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8], big_endian: bool) -> Self {
        Self { buf, pos: 0, big_endian }
    }

    fn seek(&mut self, to: usize) -> Result<(), Bit2Error> {
        if to > self.buf.len() {
            return Err(Bit2Error::Truncated { at: self.pos, need: to - self.pos });
        }
        self.pos = to;
        Ok(())
    }

    fn bytes(&mut self, n: usize) -> Result<&'a [u8], Bit2Error> {
        if self.pos + n > self.buf.len() {
            return Err(Bit2Error::Truncated { at: self.pos, need: n });
        }
        let out = &self.buf[self.pos..self.pos + n];
        self.pos += n;
        Ok(out)
    }

    fn u8(&mut self) -> Result<u8, Bit2Error> {
        Ok(self.bytes(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, Bit2Error> {
        let bs = self.bytes(4)?;
        let arr = [bs[0], bs[1], bs[2], bs[3]];
        Ok(if self.big_endian {
            u32::from_be_bytes(arr)
        } else {
            u32::from_le_bytes(arr)
        })
    }
}

/// Build a minimal 2bit byte stream from one or more `(name, ASCII)`
/// pairs. Used in tests to exercise the parser round-trip without
/// requiring `faToTwoBit` on the host.
#[cfg(test)]
fn build_2bit(seqs: &[(&str, &[u8])]) -> Vec<u8> {
    let mut out = Vec::new();
    // Header.
    out.extend_from_slice(&TWOBIT_MAGIC_LE.to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes()); // version
    out.extend_from_slice(&(seqs.len() as u32).to_le_bytes());
    out.extend_from_slice(&0u32.to_le_bytes()); // reserved

    // Compute index sizes so we can lay out sequence-record offsets.
    let mut index_size = 0usize;
    for (name, _) in seqs {
        index_size += 1 + name.len() + 4;
    }
    let mut seq_start = 16 + index_size; // header is 16 bytes

    // Emit index with placeholder offsets we'll back-patch below.
    let mut offset_patch_positions = Vec::new();
    for (name, _) in seqs {
        out.push(name.len() as u8);
        out.extend_from_slice(name.as_bytes());
        offset_patch_positions.push(out.len());
        out.extend_from_slice(&0u32.to_le_bytes()); // placeholder
    }

    // Emit sequence records, patching offsets.
    for ((_, bases), patch_pos) in seqs.iter().zip(offset_patch_positions.iter()) {
        let offset = out.len();
        out[*patch_pos..*patch_pos + 4].copy_from_slice(&(offset as u32).to_le_bytes());

        // Derive n-blocks and mask-blocks from the ASCII.
        let mut n_blocks: Vec<(u32, u32)> = Vec::new();
        let mut m_blocks: Vec<(u32, u32)> = Vec::new();
        {
            let mut n_start: Option<u32> = None;
            let mut m_start: Option<u32> = None;
            for (i, &b) in bases.iter().enumerate() {
                let up = b.to_ascii_uppercase();
                let is_n = up == b'N';
                let is_mask = b.is_ascii_lowercase() && !is_n;
                if is_n {
                    if n_start.is_none() {
                        n_start = Some(i as u32);
                    }
                } else if let Some(s) = n_start.take() {
                    n_blocks.push((s, i as u32 - s));
                }
                if is_mask {
                    if m_start.is_none() {
                        m_start = Some(i as u32);
                    }
                } else if let Some(s) = m_start.take() {
                    m_blocks.push((s, i as u32 - s));
                }
            }
            if let Some(s) = n_start {
                n_blocks.push((s, bases.len() as u32 - s));
            }
            if let Some(s) = m_start {
                m_blocks.push((s, bases.len() as u32 - s));
            }
        }

        out.extend_from_slice(&(bases.len() as u32).to_le_bytes());
        out.extend_from_slice(&(n_blocks.len() as u32).to_le_bytes());
        for (s, _) in &n_blocks {
            out.extend_from_slice(&s.to_le_bytes());
        }
        for (_, sz) in &n_blocks {
            out.extend_from_slice(&sz.to_le_bytes());
        }
        out.extend_from_slice(&(m_blocks.len() as u32).to_le_bytes());
        for (s, _) in &m_blocks {
            out.extend_from_slice(&s.to_le_bytes());
        }
        for (_, sz) in &m_blocks {
            out.extend_from_slice(&sz.to_le_bytes());
        }
        out.extend_from_slice(&0u32.to_le_bytes()); // reserved

        // Packed DNA. 2bit codes: T=00, C=01, A=10, G=11. Arbitrary for
        // N (we use A).
        let packed_len = (bases.len() + 3) / 4;
        let mut packed = vec![0u8; packed_len];
        for (i, &b) in bases.iter().enumerate() {
            let code = match b.to_ascii_uppercase() {
                b'T' => 0b00u8,
                b'C' => 0b01u8,
                b'A' | b'N' => 0b10u8,
                b'G' => 0b11u8,
                _ => 0b10u8,
            };
            let shift = 6 - (i & 3) * 2;
            packed[i >> 2] |= code << shift;
        }
        out.extend_from_slice(&packed);

        let _ = seq_start;
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_magic_is_detected() {
        let bogus = vec![0u8; 32];
        assert!(matches!(parse_2bit(&bogus), Err(Bit2Error::BadMagic { .. })));
    }

    #[test]
    fn big_endian_round_trip() {
        // Build a 2bit with BE byte order by hand.
        let seqs: &[(&str, &[u8])] = &[("chr1", b"ACGTACGT")];
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&TWOBIT_MAGIC_LE.to_be_bytes());
        bytes.extend_from_slice(&0u32.to_be_bytes()); // version
        bytes.extend_from_slice(&1u32.to_be_bytes()); // seq_count
        bytes.extend_from_slice(&0u32.to_be_bytes()); // reserved
        bytes.push(4u8); // name_size
        bytes.extend_from_slice(b"chr1");
        let offset_patch = bytes.len();
        bytes.extend_from_slice(&0u32.to_be_bytes()); // placeholder offset
        let seq_start = bytes.len();
        bytes[offset_patch..offset_patch + 4]
            .copy_from_slice(&(seq_start as u32).to_be_bytes());
        bytes.extend_from_slice(&(seqs[0].1.len() as u32).to_be_bytes());
        bytes.extend_from_slice(&0u32.to_be_bytes()); // n_count
        bytes.extend_from_slice(&0u32.to_be_bytes()); // m_count
        bytes.extend_from_slice(&0u32.to_be_bytes()); // reserved
        // Packed DNA: ACGT = 10 01 11 00 = 0b10011100 = 0x9C
        // Second ACGT same.
        bytes.extend_from_slice(&[0x9C, 0x9C]);
        let parsed = parse_2bit(&bytes).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "chr1");
        assert_eq!(parsed[0].seq.to_ascii(), b"ACGTACGT");
    }

    #[test]
    fn roundtrip_single_sequence() {
        let seqs = &[("chr1", b"ACGTACGTACGT" as &[u8])];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].name, "chr1");
        assert_eq!(parsed[0].seq.to_ascii(), b"ACGTACGTACGT");
    }

    #[test]
    fn roundtrip_multi_sequence() {
        let seqs: &[(&str, &[u8])] = &[
            ("chrA", b"ACGT"),
            ("chrB", b"TTTT"),
            ("chrC", b"GATTACACATG"),
        ];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed[0].name, "chrA");
        assert_eq!(parsed[0].seq.to_ascii(), b"ACGT");
        assert_eq!(parsed[1].name, "chrB");
        assert_eq!(parsed[1].seq.to_ascii(), b"TTTT");
        assert_eq!(parsed[2].name, "chrC");
        assert_eq!(parsed[2].seq.to_ascii(), b"GATTACACATG");
    }

    #[test]
    fn n_blocks_round_trip() {
        let seqs: &[(&str, &[u8])] = &[("chr1", b"ACGTNNNNACGT")];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        assert_eq!(parsed[0].seq.to_ascii(), b"ACGTNNNNACGT");
        for i in 4..8 {
            assert!(!parsed[0].seq.is_valid(i), "pos {i} should be N");
        }
    }

    #[test]
    fn mask_blocks_round_trip() {
        // Uppercase ACGT + one lowercase run acgt + uppercase ACGT.
        let seqs: &[(&str, &[u8])] = &[("chr1", b"ACGTacgtACGT")];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        // to_ascii uppercases; use is_masked for the semantic check.
        assert_eq!(parsed[0].seq.to_ascii(), b"ACGTACGTACGT");
        for i in 4..8 {
            assert!(parsed[0].seq.is_masked(i), "pos {i} should be masked");
        }
        for i in 0..4 {
            assert!(!parsed[0].seq.is_masked(i));
        }
        for i in 8..12 {
            assert!(!parsed[0].seq.is_masked(i));
        }
    }

    #[test]
    fn overlapping_n_and_mask() {
        // N-block takes precedence over mask-block in the decoded output
        // (a position can't be both N and lowercase ACGT).
        let seqs: &[(&str, &[u8])] = &[("chr1", b"ACntACGT")];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        // positions 2,3 are 'nt' in source — 'n' is N (invalid), 't' is mask.
        assert!(!parsed[0].seq.is_valid(2));
        assert!(parsed[0].seq.is_valid(3) && parsed[0].seq.is_masked(3));
    }

    #[test]
    fn non_multiple_of_four_length_parses() {
        let seqs: &[(&str, &[u8])] = &[
            ("chr1", b"A"),      // 1 base
            ("chr2", b"AC"),     // 2 bases
            ("chr3", b"ACG"),    // 3 bases
            ("chr4", b"ACGT"),   // 4 bases
            ("chr5", b"ACGTA"),  // 5 bases
        ];
        let bytes = build_2bit(seqs);
        let parsed = parse_2bit(&bytes).unwrap();
        for (expected, p) in seqs.iter().zip(parsed.iter()) {
            assert_eq!(p.seq.to_ascii(), expected.1);
        }
    }

    #[test]
    fn truncated_file_errors_cleanly() {
        let seqs: &[(&str, &[u8])] = &[("chr1", b"ACGT")];
        let mut bytes = build_2bit(seqs);
        bytes.truncate(bytes.len() - 2);
        assert!(matches!(parse_2bit(&bytes), Err(Bit2Error::Truncated { .. })));
    }
}
