//! PAF (Pairwise Alignment Format) writer.
//!
//! PAF is the tab-separated format introduced by minimap2 and widely used in
//! modern WGA pipelines. Each record is one line with 12 mandatory columns
//! plus optional `tag:type:value` entries. We emit `cg:Z:<CIGAR>`, `NM:i`,
//! `AS:i`, and `tp:A:P` (primary alignment), which covers the tags most
//! downstream tools consume.

use std::io::{self, Write};

use super::{Record, Strand};
use crate::edit_script::EditOp;

pub struct PafWriter<W: Write> {
    inner: W,
}

impl<W: Write> PafWriter<W> {
    pub fn new(inner: W) -> Self {
        Self { inner }
    }

    pub fn write_record(&mut self, rec: &Record) -> io::Result<()> {
        let strand_char = match rec.query_strand {
            Strand::Plus => '+',
            Strand::Minus => '-',
        };

        let (matches, mismatches, indels) = score_script(rec);
        let aligned_cols = matches + mismatches + indels;

        writeln!(
            self.inner,
            "{qname}\t{qlen}\t{qstart}\t{qend}\t{strand}\t{tname}\t{tlen}\t{tstart}\t{tend}\t{nmatch}\t{alen}\t{mapq}\tNM:i:{nm}\tAS:i:{as_}\ttp:A:P\tcg:Z:{cigar}",
            qname = rec.query_name,
            qlen = rec.query_len,
            qstart = rec.q_start,
            qend = rec.q_end(),
            strand = strand_char,
            tname = rec.target_name,
            tlen = rec.target_len,
            tstart = rec.t_start,
            tend = rec.t_end(),
            nmatch = matches,
            alen = aligned_cols,
            mapq = 60,
            nm = mismatches + indels,
            as_ = rec.score,
            cigar = rec.script.to_cigar(),
        )
    }

    pub fn flush(&mut self) -> io::Result<()> {
        self.inner.flush()
    }
}

/// Count (matches, mismatches, indel columns) in the aligned region.
fn score_script(rec: &Record) -> (u32, u32, u32) {
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut indels = 0u32;
    let mut ti = 0usize;
    let mut qi = 0usize;
    for &(op, n) in rec.script.runs() {
        match op {
            EditOp::Match => {
                for _ in 0..n {
                    if rec.target_bases.get(ti) == rec.query_bases.get(qi) {
                        matches += 1;
                    } else {
                        mismatches += 1;
                    }
                    ti += 1;
                    qi += 1;
                }
            }
            EditOp::InsertQuery => {
                indels += n;
                qi += n as usize;
            }
            EditOp::DeleteQuery => {
                indels += n;
                ti += n as usize;
            }
        }
    }
    (matches, mismatches, indels)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::edit_script::{EditOp, EditScript};

    fn match_script(n: u32) -> EditScript {
        let mut s = EditScript::new();
        s.push(EditOp::Match, n);
        s
    }

    #[test]
    fn emits_twelve_mandatory_columns_plus_tags() {
        let mut out = Vec::new();
        {
            let mut w = PafWriter::new(&mut out);
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
                script: match_script(4),
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let line = text.lines().next().unwrap();
        let cols: Vec<&str> = line.split('\t').collect();
        assert!(cols.len() >= 12, "got {} columns: {}", cols.len(), line);
        assert_eq!(cols[0], "chrQ");
        assert_eq!(cols[4], "+");
        assert_eq!(cols[5], "chrT");
        assert_eq!(cols[7], "10");
        assert_eq!(cols[8], "14");
        assert!(line.contains("cg:Z:4M"));
        assert!(line.contains("NM:i:0"));
        assert!(line.contains("AS:i:400"));
    }

    #[test]
    fn counts_mismatches_and_indels() {
        let mut script = EditScript::new();
        script.push(EditOp::Match, 2); // both matches
        script.push(EditOp::InsertQuery, 1); // 1 indel col
        script.push(EditOp::Match, 2); // one match, one mismatch (we set query)
        let mut out = Vec::new();
        {
            let mut w = PafWriter::new(&mut out);
            let rec = Record {
                target_name: "t".into(),
                target_len: 10,
                query_name: "q".into(),
                query_len: 10,
                query_strand: Strand::Plus,
                t_start: 0,
                q_start: 0,
                t_span: 4,
                q_span: 5,
                score: 10,
                script,
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACCGA".to_vec(), // AC=match, then insert 'C', then GA vs GT: G=match, A!=T
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        assert!(text.contains("NM:i:2"), "{text}"); // 1 mismatch + 1 indel
        assert!(text.contains("cg:Z:2M1I2M"));
    }

    #[test]
    fn minus_strand_renders_dash() {
        let mut out = Vec::new();
        {
            let mut w = PafWriter::new(&mut out);
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
                score: 10,
                script: match_script(4),
                target_bases: b"ACGT".to_vec(),
                query_bases: b"ACGT".to_vec(),
            };
            w.write_record(&rec).unwrap();
        }
        let text = String::from_utf8(out).unwrap();
        let line = text.lines().next().unwrap();
        assert_eq!(line.split('\t').nth(4).unwrap(), "-");
    }
}
