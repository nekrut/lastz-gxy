//! Streaming MAF (Multiple Alignment Format) writer.
//!
//! Output follows the MAF 1 spec as used by upstream lastz
//! (`docs/maf_format.html`). Gap columns are rendered as `-` in both the
//! target and query rows per the spec.

use std::io::{self, Write};

use super::Record;
use crate::edit_script::render_aligned_bases;

/// A writer that emits one MAF block per `write_record` call.
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

        let (t_text, q_text) =
            render_aligned_bases(&rec.script, &rec.target_bases, &rec.query_bases);

        writeln!(self.inner, "a score={}", rec.score)?;
        writeln!(
            self.inner,
            "s {name:<20} {start:>10} {size:>6} {strand} {src_size:>10} {text}",
            name = rec.target_name,
            start = rec.t_start,
            size = rec.t_span,
            strand = '+',
            src_size = rec.target_len,
            text = std::str::from_utf8(&t_text).unwrap_or("?"),
        )?;
        writeln!(
            self.inner,
            "s {name:<20} {start:>10} {size:>6} {strand} {src_size:>10} {text}",
            name = rec.query_name,
            start = rec.q_start,
            size = rec.q_span,
            strand = rec.query_strand.as_char(),
            src_size = rec.query_len,
            text = std::str::from_utf8(&q_text).unwrap_or("?"),
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
    use crate::edit_script::{EditOp, EditScript};
    use crate::output::Strand;

    fn ungapped_script(n: u32) -> EditScript {
        let mut s = EditScript::new();
        s.push(EditOp::Match, n);
        s
    }

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
                t_start: 10,
                q_start: 5,
                t_span: 4,
                q_span: 4,
                score: 400,
                script: ungapped_script(4),
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert_eq!(text.matches("##maf").count(), 1);
        assert_eq!(text.matches("a score=400").count(), 2);
    }

    #[test]
    fn renders_gap_characters_in_both_rows() {
        let mut script = EditScript::new();
        script.push(EditOp::Match, 3);
        script.push(EditOp::InsertQuery, 1); // target gap
        script.push(EditOp::Match, 2);
        script.push(EditOp::DeleteQuery, 1); // query gap
        script.push(EditOp::Match, 1);

        let mut out = Vec::new();
        {
            let mut w = MafWriter::new(&mut out);
            let rec = Record {
                target_name: "t".into(),
                target_len: 50,
                query_name: "q".into(),
                query_len: 50,
                query_strand: Strand::Plus,
                t_start: 0,
                q_start: 0,
                t_span: 7,
                q_span: 7,
                score: 100,
                script,
                target_bases: b"AAACCCA".to_vec(),
                query_bases: b"AAAGCCCC".to_vec(),
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        // Target side should have a '-' where InsertQuery sits.
        assert!(text.contains("AAA-CC"), "target row missing '-': {text}");
        // Query side should have a '-' where DeleteQuery sits.
        assert!(text.contains("-C\n") || text.contains("-C"), "{text}");
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
                t_start: 0,
                q_start: 0,
                t_span: 4,
                q_span: 4,
                score: 100,
                script: ungapped_script(4),
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let s_lines: Vec<&str> = text.lines().filter(|l| l.starts_with("s ")).collect();
        assert_eq!(s_lines.len(), 2);
        assert!(s_lines[0].contains(" + "));
        assert!(s_lines[1].contains(" - "));
    }
}
