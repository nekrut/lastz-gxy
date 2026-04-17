//! Streaming MAF (Multiple Alignment Format) writer.
//!
//! Output follows the MAF 1 spec as used by upstream lastz (`docs/maf_format.html`).
//! Phase 1 emits only ungapped `s` records; there are no gap characters since
//! HSPs are by definition gap-free.
//!
//! Coordinates:
//! - Plus strand: `s_start` is 0-based position on the `+` strand.
//! - Minus strand: `s_start` is 0-based position on the `-` strand
//!   (MAF spec §coordinate), i.e. `seq_len - q_end` of the forward coords.

use std::io::{self, Write};

use super::Record;

/// A writer that emits one MAF block per `write_record` call. Holds no state
/// between records; alignments appear in the order they are written.
pub struct MafWriter<W: Write> {
    inner: W,
    header_written: bool,
    scoring_desc: Option<String>,
}

impl<W: Write> MafWriter<W> {
    pub fn new(inner: W) -> Self {
        Self {
            inner,
            header_written: false,
            scoring_desc: None,
        }
    }

    /// Attach a scoring-matrix description for the `##maf` header line.
    pub fn with_scoring_desc(mut self, desc: impl Into<String>) -> Self {
        self.scoring_desc = Some(desc.into());
        self
    }

    fn write_header(&mut self) -> io::Result<()> {
        if self.header_written {
            return Ok(());
        }
        writeln!(self.inner, "##maf version=1 program=lastz-gxy")?;
        if let Some(desc) = &self.scoring_desc {
            writeln!(self.inner, "# scoring={desc}")?;
        }
        self.header_written = true;
        Ok(())
    }

    pub fn write_record(&mut self, rec: &Record) -> io::Result<()> {
        self.write_header()?;

        let t_name = &rec.target_name;
        let q_name = &rec.query_name;
        let len = rec.hsp.length;

        // For the minus strand, the MAF spec wants coordinates expressed on
        // the reverse-complement strand. Our `hsp.q_start` is already in the
        // coordinate frame of the sequence we extracted the HSP against
        // (which is the rc sequence when strand is Minus), so no additional
        // translation is needed here — the caller has already flipped.
        let q_strand_char = rec.query_strand.as_char();

        writeln!(self.inner, "a score={}", rec.hsp.score)?;
        writeln!(
            self.inner,
            "s {name:<20} {start:>10} {size:>6} {strand} {src_size:>10} {text}",
            name = t_name,
            start = rec.hsp.t_start,
            size = len,
            strand = '+',
            src_size = rec.target_len,
            text = std::str::from_utf8(&rec.target_bases).unwrap_or("?"),
        )?;
        writeln!(
            self.inner,
            "s {name:<20} {start:>10} {size:>6} {strand} {src_size:>10} {text}",
            name = q_name,
            start = rec.hsp.q_start,
            size = len,
            strand = q_strand_char,
            src_size = rec.query_len,
            text = std::str::from_utf8(&rec.query_bases).unwrap_or("?"),
        )?;
        writeln!(self.inner)?;
        Ok(())
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }

    pub fn into_inner(self) -> W {
        self.inner
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hsp::Hsp;
    use crate::output::Strand;

    #[test]
    fn writes_header_once() {
        let mut out = Vec::new();
        {
            let mut w = MafWriter::new(&mut out);
            let rec = Record {
                target_name: "chrT".into(),
                target_len: 100,
                query_name: "chrQ".into(),
                query_len: 50,
                query_strand: Strand::Plus,
                hsp: Hsp { t_start: 10, q_start: 5, length: 4, score: 400 },
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert_eq!(
            text.matches("##maf").count(),
            1,
            "header should appear once:\n{text}"
        );
        assert_eq!(text.matches("a score=400").count(), 2);
    }

    #[test]
    fn minus_strand_marked_on_query_line() {
        let mut out = Vec::new();
        {
            let mut w = MafWriter::new(&mut out);
            let rec = Record {
                target_name: "t".into(),
                target_len: 10,
                query_name: "q".into(),
                query_len: 10,
                query_strand: Strand::Minus,
                hsp: Hsp { t_start: 0, q_start: 0, length: 4, score: 100 },
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        // Target line always +; query line picks up the strand.
        let s_lines: Vec<&&str> = lines.iter().filter(|l| l.starts_with("s ")).collect();
        assert_eq!(s_lines.len(), 2);
        assert!(s_lines[0].contains(" + "), "{}", s_lines[0]);
        assert!(s_lines[1].contains(" - "), "{}", s_lines[1]);
    }
}
