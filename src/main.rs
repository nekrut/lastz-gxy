use std::fs::File;
use std::io::{stdout, BufWriter, Write};

use clap::Parser;

use lastz_gxy::bit2::load_2bit;
use lastz_gxy::cli::{build_config, Cli, Format};
use lastz_gxy::driver::run;
use lastz_gxy::output::{maf::MafWriter, paf::PafWriter, sam::SamWriter};
use lastz_gxy::sequences::{load_fasta, Sequence};

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    let config = build_config(&cli)?;
    let mut targets = load_sequences(&cli.target)?;
    let mut queries = load_sequences(&cli.query)?;

    if !config.respect_masking {
        for s in targets.iter_mut().chain(queries.iter_mut()) {
            s.seq.clear_masks();
        }
    }

    let records = run(&targets, &queries, &config);

    let out: Box<dyn Write> = if cli.out.as_os_str() == "-" {
        Box::new(BufWriter::new(stdout().lock()))
    } else {
        Box::new(BufWriter::new(File::create(&cli.out)?))
    };

    match cli.format {
        Format::Maf => {
            let mut w = MafWriter::new(out).with_scoring_desc("HOXD70");
            for rec in &records {
                w.write_record(rec)?;
            }
            w.flush()?;
        }
        Format::Paf => {
            let mut w = PafWriter::new(out);
            for rec in &records {
                w.write_record(rec)?;
            }
            w.flush()?;
        }
        Format::Sam => {
            let mut w = SamWriter::new(out);
            // Declare all targets up front so the header lists every chrom.
            for t in &targets {
                w.declare_target(&t.name, t.seq.len() as u32);
            }
            for rec in &records {
                w.write_record(rec)?;
            }
            w.flush()?;
        }
    }

    eprintln!(
        "lastz-gxy: {} target(s), {} query(s), {} alignment(s)",
        targets.len(),
        queries.len(),
        records.len()
    );
    Ok(())
}

/// Dispatch a sequence-file path to the right parser based on its
/// extension. `.2bit` routes to the UCSC 2bit parser; anything else is
/// treated as FASTA.
fn load_sequences(path: &std::path::Path) -> anyhow::Result<Vec<Sequence>> {
    let is_2bit = path
        .extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("2bit"));
    if is_2bit {
        Ok(load_2bit(path)?)
    } else {
        Ok(load_fasta(path)?)
    }
}
