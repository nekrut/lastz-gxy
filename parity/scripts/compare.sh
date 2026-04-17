#!/usr/bin/env bash
# Run upstream lastz and lastz-gxy on a FASTA pair, produce both MAFs, and
# diff them with `gxy-compare`. Returns 0 on parity-gate pass, 1 on fail,
# 2 on setup errors. Intended for CI and for the manual dev loop.
#
# Usage:  parity/scripts/compare.sh <target.fa> <query.fa> [lastz-flags...]
# Env:
#   LASTZ_UPSTREAM_BIN  — override upstream binary (default parity/upstream/bin/lastz)
#   LASTZ_GXY_BIN       — override test binary (default target/release/lastz-gxy)
#   GXY_COMPARE_BIN     — override compare binary (default target/release/gxy-compare)
#   ENFORCE_GATE        — set to 1 to make non-passing runs exit 1

set -euo pipefail

if [[ $# -lt 2 ]]; then
    echo "usage: $0 <target.fa> <query.fa> [lastz-flags...]" >&2
    exit 2
fi

TARGET="$1"
QUERY="$2"
shift 2
EXTRA_FLAGS=("$@")

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="$ROOT/parity/out"
mkdir -p "$OUT_DIR"

UPSTREAM_BIN="${LASTZ_UPSTREAM_BIN:-$ROOT/parity/upstream/bin/lastz}"
GXY_BIN="${LASTZ_GXY_BIN:-$ROOT/target/release/lastz-gxy}"
COMPARE_BIN="${GXY_COMPARE_BIN:-$ROOT/target/release/gxy-compare}"

if [[ ! -x "$UPSTREAM_BIN" ]]; then
    echo "upstream binary not found: $UPSTREAM_BIN" >&2
    echo "build it first:  parity/scripts/build-upstream.sh" >&2
    exit 2
fi
if [[ ! -x "$GXY_BIN" ]]; then
    echo "lastz-gxy binary not found: $GXY_BIN" >&2
    echo "build it first:  cargo build --release" >&2
    exit 2
fi
if [[ ! -x "$COMPARE_BIN" ]]; then
    echo "gxy-compare binary not found: $COMPARE_BIN" >&2
    echo "build it first:  cargo build --release --bin gxy-compare" >&2
    exit 2
fi

tag="$(basename "$TARGET" .fa).$(basename "$QUERY" .fa)"
UP_MAF="$OUT_DIR/$tag.upstream.maf"
GXY_MAF="$OUT_DIR/$tag.gxy.maf"

echo "==> running upstream lastz"
"$UPSTREAM_BIN" "$TARGET" "$QUERY" --format=maf "${EXTRA_FLAGS[@]}" > "$UP_MAF"

echo "==> running lastz-gxy"
"$GXY_BIN" "$TARGET" "$QUERY" --format maf "${EXTRA_FLAGS[@]}" > "$GXY_MAF"

echo "==> comparing"
if [[ "${ENFORCE_GATE:-0}" = "1" ]]; then
    "$COMPARE_BIN" --enforce "$UP_MAF" "$GXY_MAF"
else
    "$COMPARE_BIN" "$UP_MAF" "$GXY_MAF"
fi
