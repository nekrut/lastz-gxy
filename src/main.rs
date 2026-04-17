use std::fs::File;
use std::io::{stdout, BufWriter, Write};

use clap::Parser;

use lastz_gxy::cli::{build_config, Cli, Format};
use lastz_gxy::driver::run;
use lastz_gxy::output::maf::MafWriter;
use lastz_gxy::sequences::load_fasta;

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok(); // If already initialised (e.g. in tests), proceed.
    }

    let config = build_config(&cli)?;
    let targets = load_fasta(&cli.target)?;
    let queries = load_fasta(&cli.query)?;

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
    }

    eprintln!(
        "lastz-gxy: {} target(s), {} query(s), {} alignment(s)",
        targets.len(),
        queries.len(),
        records.len()
    );
    Ok(())
}
