#!/usr/bin/env bash
# Clone upstream lastz at a pinned tag and build it. The resulting binary
# lands in parity/upstream/bin/lastz and serves as the ground truth for
# PLAN.md §5's release-gate comparisons. Idempotent: re-running is a no-op
# if the pinned commit is already built.

set -euo pipefail

LASTZ_TAG="${LASTZ_TAG:-1.04.52}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
UPSTREAM_DIR="$ROOT/parity/upstream"
CHECKOUT="$UPSTREAM_DIR/lastz-$LASTZ_TAG"
BIN_DIR="$UPSTREAM_DIR/bin"
BIN="$BIN_DIR/lastz"

if [[ -x "$BIN" ]]; then
    # Sanity-check the existing binary matches the requested tag.
    if "$BIN" --version 2>&1 | grep -q "$LASTZ_TAG"; then
        echo "upstream lastz $LASTZ_TAG already built at $BIN"
        exit 0
    fi
    echo "existing $BIN does not match tag $LASTZ_TAG; rebuilding"
    rm -f "$BIN"
fi

mkdir -p "$UPSTREAM_DIR" "$BIN_DIR"

if [[ ! -d "$CHECKOUT/.git" ]]; then
    echo "cloning lastz at tag $LASTZ_TAG"
    git clone --depth 1 --branch "$LASTZ_TAG" https://github.com/lastz/lastz.git "$CHECKOUT"
fi

pushd "$CHECKOUT" > /dev/null
echo "building lastz (this may take a minute)"
make -C src -j"$(nproc 2>/dev/null || echo 2)" lastz
popd > /dev/null

cp "$CHECKOUT/src/lastz" "$BIN"
echo "wrote $BIN"
"$BIN" --version | head -1
