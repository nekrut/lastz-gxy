//! SAM (Sequence Alignment / Map) writer.
//!
//! Output follows the SAMv1 spec (samtools.github.io/hts-specs/SAMv1.pdf).
//! Each `Record` emits one line:
//!
//! ```text
//! QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL [TAGS]
//! ```
//!
//! `@SQ` header lines are emitted once per unique target name seen. `@PG`
//! identifies lastz-gxy as the aligner.

use std::collections::BTreeMap;
use std::io::{self, Write};

use super::{Record, Strand};
use crate::edit_script::{EditOp, EditScript};

pub struct SamWriter<W: Write> {
    inner: W,
    header_written: bool,
    /// Accumulates `(target_name → target_len)` seen so far so we can emit
    /// an `@SQ` header on first use.
    sq_lines: BTreeMap<String, u32>,
}

impl<W: Write> SamWriter<W> {
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            header_written: false,
            sq_lines: BTreeMap::new(),
        }
    }

    /// Emit headers for every target that any emitted record references.
    /// Safe to call before `write_record`; otherwise run `finish_headers`
    /// before the first record.
    pub fn declare_target(&mut self, name: &str, len: u32) {
        self.sq_lines.insert(name.to_string(), len);
    }

    fn write_header(&mut self) -> io::Result<()> {
        if self.header_written {
            return Ok(());
        }
        writeln!(self.inner, "@HD\tVN:1.6\tSO:unsorted")?;
        for (name, len) in &self.sq_lines {
            writeln!(self.inner, "@SQ\tSN:{name}\tLN:{len}")?;
        }
        writeln!(
            self.inner,
            "@PG\tID:lastz-gxy\tPN:lastz-gxy\tVN:{}",
            env!("CARGO_PKG_VERSION")
        )?;
        self.header_written = true;
        Ok(())
    }

    pub fn write_record(&mut self, rec: &Record) -> io::Result<()> {
        // Declare the target on first sight so @SQ is complete even if the
        // caller didn't pre-declare.
        self.sq_lines
            .entry(rec.target_name.clone())
            .or_insert(rec.target_len);
        self.write_header()?;

        let flag = match rec.query_strand {
            Strand::Plus => 0u16,
            Strand::Minus => 16u16, // 0x10 = reverse complement
        };

        // SEQ: the query bases used for this alignment block. lastz emits
        // them in aligned orientation (matching the CIGAR); for minus-strand
        // records the caller has already reverse-complemented the query,
        // which matches SAM semantics (SEQ is always on the reference strand).
        let seq = if rec.query_bases.is_empty() {
            b"*".to_vec()
        } else {
            rec.query_bases.clone()
        };
        let cigar = soft_clip(&rec.script, rec.q_start, rec.query_len, rec.q_span);
        let nm = nm_from(rec);

        writeln!(
            self.inner,
            "{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t*\tNM:i:{nm}\tAS:i:{as_}",
            qname = rec.query_name,
            flag = flag,
            rname = rec.target_name,
            pos = rec.t_start + 1, // SAM is 1-based
            mapq = 60,
            cigar = cigar,
            seq = std::str::from_utf8(&seq).unwrap_or("*"),
            nm = nm,
            as_ = rec.score,
        )
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Build the CIGAR string the way SAM expects, adding soft-clips for query
/// bases outside the aligned span. The query is clipped on the left by
/// `q_start` bases and on the right by `q_len - q_start - q_span`.
fn soft_clip(script: &EditScript, q_start: u32, q_len: u32, q_span: u32) -> String {
    let mut out = String::new();
    let left_clip = q_start;
    let right_clip = q_len.saturating_sub(q_start + q_span);
    if left_clip > 0 {
        out.push_str(&format!("{left_clip}S"));
    }
    for (op, n) in script.runs() {
        out.push_str(&format!("{n}{c}", c = op.cigar_char()));
    }
    if right_clip > 0 {
        out.push_str(&format!("{right_clip}S"));
    }
    out
}

/// Count mismatches + indels (CIGAR ops) to populate the `NM:i:` tag.
fn nm_from(rec: &Record) -> u32 {
    let mut nm = 0u32;
    let mut ti = 0usize;
    let mut qi = 0usize;
    for &(op, n) in rec.script.runs() {
        match op {
            EditOp::Match => {
                for _ in 0..n {
                    if rec.target_bases.get(ti) != rec.query_bases.get(qi) {
                        nm += 1;
                    }
                    ti += 1;
                    qi += 1;
                }
            }
            EditOp::InsertQuery => {
                nm += n;
                qi += n as usize;
            }
            EditOp::DeleteQuery => {
                nm += n;
                ti += n as usize;
            }
        }
    }
    nm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edit_script::EditScript;

    fn match_script(n: u32) -> EditScript {
        let mut s = EditScript::new();
        s.push(EditOp::Match, n);
        s
    }

    fn rec(strand: Strand, script: EditScript, target_bases: &[u8], query_bases: &[u8]) -> Record {
        let t_span = script.target_span();
        let q_span = script.query_span();
        Record {
            target_name: "chrT".into(),
            target_len: 100,
            query_name: "read1".into(),
            query_len: 30,
            query_strand: strand,
            t_start: 10,
            q_start: 5,
            t_span,
            q_span,
            score: 400,
            script,
            target_bases: target_bases.to_vec(),
            query_bases: query_bases.to_vec(),
        }
    }

    #[test]
    fn emits_header_once_and_records() {
        let mut out = Vec::new();
        {
            let mut w = SamWriter::new(&mut out);
            let r = rec(Strand::Plus, match_script(4), b"ACGT", b"ACGT");
            w.write_record(&r).unwrap();
            w.write_record(&r).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text.matches("@HD").count(), 1);
        assert_eq!(text.matches("@SQ\tSN:chrT").count(), 1);
        assert_eq!(text.matches("@PG").count(), 1);
        // Two data lines after the headers.
        let data_lines = text.lines().filter(|l| !l.starts_with('@')).count();
        assert_eq!(data_lines, 2);
    }

    #[test]
    fn soft_clip_bookends_query_outside_alignment() {
        let mut out = Vec::new();
        {
            let mut w = SamWriter::new(&mut out);
            let r = rec(Strand::Plus, match_script(4), b"ACGT", b"ACGT");
            // q_start=5, q_span=4, q_len=30 → left clip 5, right clip 21.
            w.write_record(&r).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let line = text.lines().find(|l| !l.starts_with('@')).unwrap();
        let cigar = line.split('\t').nth(5).unwrap();
        assert_eq!(cigar, "5S4M21S");
    }

    #[test]
    fn minus_strand_sets_flag_16() {
        let mut out = Vec::new();
        {
            let mut w = SamWriter::new(&mut out);
            let r = rec(Strand::Minus, match_script(4), b"ACGT", b"ACGT");
            w.write_record(&r).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let line = text.lines().find(|l| !l.starts_with('@')).unwrap();
        let flag = line.split('\t').nth(1).unwrap();
        assert_eq!(flag, "16");
    }

    #[test]
    fn nm_counts_mismatches_and_indels() {
        let mut script = EditScript::new();
        script.push(EditOp::Match, 2);
        script.push(EditOp::InsertQuery, 1);
        script.push(EditOp::Match, 2);
        let r = rec(Strand::Plus, script, b"ACGT", b"ACCGA"); // AC==AC, insert C, GA vs GT
        let mut out = Vec::new();
        {
            let mut w = SamWriter::new(&mut out);
            w.write_record(&r).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("NM:i:2"), "{text}"); // 1 insert + 1 mismatch
    }

    #[test]
    fn pos_is_one_based() {
        let mut out = Vec::new();
        {
            let mut w = SamWriter::new(&mut out);
            let r = rec(Strand::Plus, match_script(4), b"ACGT", b"ACGT");
            w.write_record(&r).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let line = text.lines().find(|l| !l.starts_with('@')).unwrap();
        let pos = line.split('\t').nth(3).unwrap();
        assert_eq!(pos, "11"); // t_start=10 → 1-based POS=11
    }
}
