//! `gxy-compare` — compare two MAF files and report parity metrics.
//!
//! Used to drive the PLAN.md §5 release gate: run upstream lastz and
//! lastz-gxy on the same fixture, feed both MAFs to `gxy-compare`, and
//! confirm Jaccard ≥ 0.99, median score Δ = 0, max |score Δ| ≤ 1, and
//! aligned-bp delta ≤ 0.1 %.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;

use lastz_gxy::parity::{compare, load_maf};

#[derive(Parser, Debug)]
#[command(name = "gxy-compare", about = "Diff two MAF files as a parity gate")]
struct Cli {
    /// MAF file treated as the baseline (typically upstream lastz).
    #[arg(value_name = "BASELINE_MAF")]
    baseline: PathBuf,

    /// MAF file under test (typically lastz-gxy).
    #[arg(value_name = "TEST_MAF")]
    test: PathBuf,

    /// Enforce the PLAN.md release gate: exit code 0 only on pass.
    #[arg(long, default_value_t = false)]
    enforce: bool,

    /// Maximum number of per-side divergent block signatures to list.
    #[arg(long, default_value_t = 10)]
    diff_limit: usize,
}

fn main() -> ExitCode {
    let cli = Cli::parse();

    let baseline = match load_maf(&cli.baseline) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("failed to read {}: {e}", cli.baseline.display());
            return ExitCode::from(2);
        }
    };
    let test = match load_maf(&cli.test) {
        Ok(b) => b,
        Err(e) => {
            eprintln!("failed to read {}: {e}", cli.test.display());
            return ExitCode::from(2);
        }
    };

    let report = compare(&baseline, &test);

    println!("baseline: {} blocks ({})", report.left_blocks, cli.baseline.display());
    println!("test:     {} blocks ({})", report.right_blocks, cli.test.display());
    println!("shared:   {}", report.intersection);
    println!("union:    {}", report.union_size);
    println!("Jaccard:  {:.6}", report.jaccard);
    println!();
    println!("score delta (test - baseline, over shared blocks):");
    println!("  median:  {}", report.score_delta_median);
    println!("  max|Δ|:  {}", report.score_delta_max_abs);
    println!("  mean:    {:.3}", report.score_delta_mean);
    println!();
    println!("aligned bp (non-gap columns):");
    println!("  baseline: {}", report.left_aligned_bp);
    println!("  test:     {}", report.right_aligned_bp);
    println!("  Δ:        {:+.3}%", report.aligned_bp_delta_pct);
    println!();
    println!("mean identity:");
    println!("  baseline: {:.4}", report.left_identity_mean);
    println!("  test:     {:.4}", report.right_identity_mean);

    if !report.left_only.is_empty() || !report.right_only.is_empty() {
        println!();
        println!("divergent blocks (up to {}):", cli.diff_limit);
        for s in report.left_only.iter().take(cli.diff_limit) {
            println!(
                "  ONLY_BASELINE {} {} {} {}..{} {}..{}",
                s.target_name, s.query_name, s.strand, s.t_start, s.t_end, s.q_start, s.q_end
            );
        }
        for s in report.right_only.iter().take(cli.diff_limit) {
            println!(
                "  ONLY_TEST     {} {} {} {}..{} {}..{}",
                s.target_name, s.query_name, s.strand, s.t_start, s.t_end, s.q_start, s.q_end
            );
        }
    }

    println!();
    let pass = report.passes_release_gate();
    println!(
        "release gate (Jaccard≥0.99, median Δ=0, max|Δ|≤1, |bp Δ|≤0.1%): {}",
        if pass { "PASS" } else { "FAIL" }
    );

    if cli.enforce && !pass {
        ExitCode::from(1)
    } else {
        ExitCode::SUCCESS
    }
}
