#!/usr/bin/env bash
# Run gxy-compare across every fixture pair in parity/corpus/, producing
# a one-line-per-pair summary table. Useful for spotting regressions
# after any change that might affect parity.
#
# Usage:  parity/scripts/matrix.sh [extra-lastz-flags...]
# Env:    same as compare.sh (LASTZ_UPSTREAM_BIN, LASTZ_GXY_BIN, etc.)

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CORPUS="$ROOT/parity/corpus"
COMPARE="$ROOT/parity/scripts/compare.sh"

EXTRA_FLAGS=("$@")

# (target, query, description) triples — one line per run. Keep names short
# so the output table stays readable.
PAIRS=(
    "pseudocat.fa|pseudocat.fa|cat self-alignment"
    "pseudopig1.fa|pseudopig1.fa|pig1 self-alignment"
    "pseudocat.fa|pseudopig1.fa|cat vs pig1 (single-contig cross)"
    "pseudocat.fa|pseudopig.fa|cat vs pig (multi-contig cross)"
    "sars_cov2.fa|sars_cov2.fa|sars-cov-2 self-alignment"
    "sars_cov2.fa|sars_cov1.fa|sars-cov-2 vs sars-cov-1 (real virus pair)"
)

printf '%-40s %6s %6s %6s %8s %8s %8s %6s\n' \
    "fixture" "base" "gxy" "shared" "RECALL" "PREC" "bpΔ" "gate"
printf '%-40s %6s %6s %6s %8s %8s %8s %6s\n' \
    "---" "---" "---" "---" "---" "---" "---" "---"

for row in "${PAIRS[@]}"; do
    IFS='|' read -r target query label <<< "$row"
    t_path="$CORPUS/$target"
    q_path="$CORPUS/$query"
    report=$("$COMPARE" "$t_path" "$q_path" "${EXTRA_FLAGS[@]}" 2>&1)
    base=$(echo "$report" | grep -E '^baseline:' | awk '{print $2}')
    gxy=$(echo "$report" | grep -E '^test:' | awk '{print $2}')
    shared=$(echo "$report" | grep -E '^shared:' | awk '{print $2}')
    recall=$(echo "$report" | grep -E '^RECALL:' | awk '{print $2}')
    precision=$(echo "$report" | grep -E '^PRECISION:' | awk '{print $2}')
    bpdelta=$(echo "$report" | grep -E '^\s*Δ:' | awk '{print $2}')
    gate=$(echo "$report" | grep -E '^release gate' | sed 's/.*: //')
    printf '%-40s %6s %6s %6s %8s %8s %8s %6s\n' \
        "$label" "$base" "$gxy" "$shared" "$recall" "$precision" "$bpdelta" "$gate"
done
