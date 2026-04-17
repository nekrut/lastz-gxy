//! CLI wiring.
//!
//! The v1 flag set is a subset of upstream lastz. Unknown / unsupported
//! flags fail fast with a pointer to PLAN.md §6 rather than silently being
//! ignored (a correctness trap upstream has had in practice).

use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use crate::driver::{Config, StrandSpec};
use crate::hsp::HspParams;
use crate::scoring::ScoringMatrix;
use crate::seeds::SeedPattern;

#[derive(Debug, Parser)]
#[command(
    name = "lastz-gxy",
    version,
    about = "Multi-core Rust reimplementation of lastz (ungapped HSP MVP)"
)]
pub struct Cli {
    /// Target FASTA (the reference side; gets indexed).
    #[arg(value_name = "TARGET")]
    pub target: PathBuf,

    /// Query FASTA (scanned against the indexed target).
    #[arg(value_name = "QUERY")]
    pub query: PathBuf,

    /// Output format. v1 Phase 1 ships MAF only; `axt`, `sam`, `paf`
    /// land with Phase 2 gapped extension.
    #[arg(long, value_enum, default_value_t = Format::Maf)]
    pub format: Format,

    /// Seed pattern.
    ///
    /// Accepts `12of19`, `19of20`, `match<k>` (e.g. `match12`), or a literal
    /// `0`/`1` spec like `11110111011110110111`.
    #[arg(long, default_value = "12of19")]
    pub seed: String,

    /// Stride between successive query / reference seed windows. `1` hits
    /// every position; larger values speed up at the cost of sensitivity.
    #[arg(long, default_value_t = 1)]
    pub step: usize,

    /// X-drop threshold for ungapped extension.
    #[arg(long, default_value_t = 910)]
    pub xdrop: i32,

    /// Minimum ungapped HSP score to report.
    #[arg(long, default_value_t = 3000)]
    pub hspthresh: i32,

    /// Which strand(s) of the query to align.
    #[arg(long, value_enum, default_value_t = StrandArg::Both)]
    pub strand: StrandArg,

    /// Optional lastz `--scores`-format matrix file. Default is HOXD70.
    #[arg(long)]
    pub scores: Option<PathBuf>,

    /// Drop seed words that occur more than this many times in the reference
    /// (maps to upstream `--maxwordcount`). `0` disables the filter.
    #[arg(long, default_value_t = 0)]
    pub max_word_count: u32,

    /// Number of worker threads. `0` (default) uses `rayon`'s default, which
    /// matches the number of logical CPUs.
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Output file. `-` means stdout.
    #[arg(long, short = 'o', default_value = "-")]
    pub out: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Format {
    Maf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StrandArg {
    Plus,
    Minus,
    Both,
}

impl From<StrandArg> for StrandSpec {
    fn from(a: StrandArg) -> Self {
        match a {
            StrandArg::Plus => StrandSpec::Plus,
            StrandArg::Minus => StrandSpec::Minus,
            StrandArg::Both => StrandSpec::Both,
        }
    }
}

/// Resolve the CLI into a runtime `Config`.
pub fn build_config(cli: &Cli) -> anyhow::Result<Config> {
    let pattern = parse_seed_spec(&cli.seed)?;
    let matrix = if let Some(path) = &cli.scores {
        let text = std::fs::read_to_string(path)?;
        ScoringMatrix::from_lastz_text(&text)?
    } else {
        ScoringMatrix::hoxd70()
    };

    Ok(Config {
        pattern,
        matrix,
        step: cli.step.max(1),
        hsp: HspParams { x_drop: cli.xdrop, hsp_threshold: cli.hspthresh },
        strand: cli.strand.into(),
        max_word_count: cli.max_word_count,
    })
}

fn parse_seed_spec(spec: &str) -> anyhow::Result<SeedPattern> {
    match spec {
        "12of19" => Ok(SeedPattern::twelve_of_nineteen()),
        "19of20" => Ok(SeedPattern::nineteen_of_twenty()),
        other if other.starts_with("match") => {
            let k: usize = other[5..].parse()?;
            Ok(SeedPattern::solid(k))
        }
        other if other.chars().all(|c| c == '0' || c == '1') => {
            Ok(SeedPattern::parse("custom", other)?)
        }
        other => anyhow::bail!(
            "unrecognised seed spec '{other}'; use '12of19', '19of20', 'match<k>', or a 0/1 pattern"
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_named_seeds() {
        assert_eq!(parse_seed_spec("12of19").unwrap().weight(), 12);
        let p = parse_seed_spec("19of20").unwrap();
        assert_eq!(p.weight(), 19);
        assert_eq!(p.len(), 20);
        assert_eq!(parse_seed_spec("match8").unwrap().weight(), 8);
    }

    #[test]
    fn parses_literal_pattern() {
        let p = parse_seed_spec("1101").unwrap();
        assert_eq!(p.weight(), 3);
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn rejects_unknown_seed_spec() {
        assert!(parse_seed_spec("deadbeef").is_err());
    }
}
