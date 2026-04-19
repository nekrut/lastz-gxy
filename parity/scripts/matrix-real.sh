#!/usr/bin/env bash
# Run the parity matrix over the real-world fixtures fetched by
# fetch-genomes.sh. Kept separate from matrix.sh because these fixtures
# are optional (they're not committed to the repo), can be large, and
# chr21 specifically takes hours to align.
#
# Usage:  parity/scripts/matrix-real.sh [extra-lastz-flags...]
# Env:    REAL_CHR21=1   include the chr21 pair (default: skipped)

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REAL_DIR="$ROOT/parity/corpus/real"
COMPARE="$ROOT/parity/scripts/compare.sh"

if [[ ! -d "$REAL_DIR" ]]; then
    echo "no real fixtures found. Fetch them first:" >&2
    echo "    parity/scripts/fetch-genomes.sh" >&2
    exit 2
fi

EXTRA_FLAGS=("$@")

PAIRS=(
    "hg38_chrM.fa|panTro_chrM.fa|human chrM vs chimp chrM (real mitochondrial pair)"
    "hg38_chrM.fa|hg38_chrM.fa|human chrM self-alignment"
)
if [[ "${REAL_CHR21:-0}" = "1" ]]; then
    PAIRS+=(
        "hg38_chr21.fa|panTro6_chr21.fa|human chr21 vs chimp chr21 (real mammalian pair)"
    )
fi

printf '%-60s %6s %6s %6s %8s %8s %8s %6s\n' \
    "fixture" "base" "gxy" "shared" "RECALL" "PREC" "bpΔ" "gate"
printf '%-60s %6s %6s %6s %8s %8s %8s %6s\n' \
    "---" "---" "---" "---" "---" "---" "---" "---"

for row in "${PAIRS[@]}"; do
    IFS='|' read -r target query label <<< "$row"
    t_path="$REAL_DIR/$target"
    q_path="$REAL_DIR/$query"
    if [[ ! -s "$t_path" || ! -s "$q_path" ]]; then
        printf '%-60s %-40s\n' "$label" "MISSING — run fetch-genomes.sh"
        continue
    fi
    report=$("$COMPARE" "$t_path" "$q_path" "${EXTRA_FLAGS[@]}" 2>&1)
    base=$(echo "$report" | grep -E '^baseline:' | awk '{print $2}')
    gxy=$(echo "$report" | grep -E '^test:' | awk '{print $2}')
    shared=$(echo "$report" | grep -E '^shared:' | awk '{print $2}')
    recall=$(echo "$report" | grep -E '^RECALL:' | awk '{print $2}')
    precision=$(echo "$report" | grep -E '^PRECISION:' | awk '{print $2}')
    bpdelta=$(echo "$report" | grep -E '^\s*Δ:' | awk '{print $2}')
    gate=$(echo "$report" | grep -E '^release gate' | sed 's/.*: //')
    printf '%-60s %6s %6s %6s %8s %8s %8s %6s\n' \
        "$label" "$base" "$gxy" "$shared" "$recall" "$precision" "$bpdelta" "$gate"
done
