//! CLI wiring.
//!
//! The v1 flag set is a subset of upstream lastz. Unknown / unsupported
//! flags fail fast with a pointer to PLAN.md §6 rather than silently being
//! ignored (a correctness trap upstream has had in practice).

use std::path::PathBuf;

use clap::{Parser, ValueEnum};

use crate::driver::{Config, StrandSpec};
use crate::gapped_extend::GappedParams;
use crate::hsp::HspParams;
use crate::scoring::ScoringMatrix;
use crate::seeds::SeedPattern;
use crate::tweener::TweenerConfig;

#[derive(Debug, Parser)]
#[command(
    name = "lastz-gxy",
    version,
    about = "Multi-core Rust reimplementation of lastz"
)]
pub struct Cli {
    /// Target FASTA (the reference side; gets indexed).
    #[arg(value_name = "TARGET")]
    pub target: PathBuf,

    /// Query FASTA (scanned against the indexed target).
    #[arg(value_name = "QUERY")]
    pub query: PathBuf,

    /// Output format.
    #[arg(long, value_enum, default_value_t = Format::Maf)]
    pub format: Format,

    /// Seed pattern.
    ///
    /// Accepts `12of19`, `19of20`, `match<k>` (e.g. `match12`), or a literal
    /// `0`/`1` spec like `11110111011110110111`.
    #[arg(long, default_value = "12of19")]
    pub seed: String,

    /// Stride between successive query / reference seed windows.
    #[arg(long, default_value_t = 1)]
    pub step: usize,

    /// X-drop threshold for ungapped HSP extension.
    #[arg(long, default_value_t = 910)]
    pub xdrop: i32,

    /// Minimum ungapped HSP score to keep for gapped extension.
    #[arg(long, default_value_t = 3_000)]
    pub hspthresh: i32,

    /// Y-drop threshold for gapped affine extension.
    #[arg(long, default_value_t = 9_400)]
    pub ydrop: i32,

    /// Minimum gapped alignment score to report.
    #[arg(long, default_value_t = 3_000)]
    pub gappedthresh: i32,

    /// Skip gapped extension; output ungapped HSPs only.
    #[arg(long, default_value_t = false)]
    pub nogapped: bool,

    /// Filter HSPs through lastz-style chaining before extension.
    #[arg(long, default_value_t = false)]
    pub chain: bool,

    /// Run inter-alignment interpolation (upstream `tweener.c`) after
    /// chaining. Scans the gap between adjacent chain members with a
    /// denser seed pattern to recover borderline alignments.
    #[arg(long, default_value_t = false)]
    pub inner: bool,

    /// Denser seed pattern used by `--inner`. Same syntax as `--seed`.
    #[arg(long, default_value = "match12")]
    pub inner_seed: String,

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

    /// Ignore soft-masking (lowercase `acgt` in the FASTA) when picking
    /// seed positions. Matches upstream lastz's `--nomasking`. The default
    /// behaviour — honoring soft-masks — matches upstream's default.
    #[arg(long, default_value_t = false)]
    pub nomasking: bool,

    /// Number of transition substitutions allowed per seed. Matches
    /// upstream: 0 = `--notransition`, 1 = `--transition` (the default),
    /// 2 = `--transition=2`.
    #[arg(long, default_value_t = 1)]
    pub transition: u8,

    /// Shorthand for `--transition=0`. Takes precedence over `--transition`
    /// when both are given.
    #[arg(long, default_value_t = false)]
    pub notransition: bool,

    /// Anchor window width (columns) used to pick the gapped-extension
    /// start position inside each HSP.
    #[arg(long, default_value_t = 31)]
    pub anchor_window: u32,

    /// Worker thread count. `0` uses `rayon`'s default (logical CPU count).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    /// Output file. `-` means stdout.
    #[arg(long, short = 'o', default_value = "-")]
    pub out: PathBuf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Format {
    Maf,
    Paf,
    Sam,
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

pub fn build_config(cli: &Cli) -> anyhow::Result<Config> {
    let pattern = parse_seed_spec(&cli.seed)?;
    let matrix = if let Some(path) = &cli.scores {
        let text = std::fs::read_to_string(path)?;
        ScoringMatrix::from_lastz_text(&text)?
    } else {
        ScoringMatrix::hoxd70()
    };

    let tweener = if cli.inner {
        Some(TweenerConfig {
            pattern: parse_seed_spec(&cli.inner_seed)?,
            ..TweenerConfig::default()
        })
    } else {
        None
    };

    let transitions = if cli.notransition { 0 } else { cli.transition.min(2) };

    Ok(Config {
        pattern,
        matrix,
        step: cli.step.max(1),
        hsp: HspParams { x_drop: cli.xdrop, hsp_threshold: cli.hspthresh },
        gapped: GappedParams { y_drop: cli.ydrop, gapped_threshold: cli.gappedthresh },
        strand: cli.strand.into(),
        max_word_count: cli.max_word_count,
        respect_masking: !cli.nomasking,
        gapped_enabled: !cli.nogapped,
        chain_enabled: cli.chain,
        anchor_window: cli.anchor_window.max(1),
        transitions,
        tweener,
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
