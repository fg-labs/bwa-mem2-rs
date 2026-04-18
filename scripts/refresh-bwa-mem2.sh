#!/usr/bin/env bash
# Refresh the vendored bwa-mem2 snapshot to a given commit on fg-labs/bwa-mem2.
# Default branch vendored is `fg-main` (our integration branch).
#
# Usage: scripts/refresh-bwa-mem2.sh <commit-hash> [source-path]
#   - If source-path is provided, the vendor is copied from a local working tree
#     (useful when the commit is not yet pushed). Otherwise, a fresh clone is made.
set -euo pipefail
HASH="${1:?usage: refresh-bwa-mem2.sh <commit-hash> [source-path]}"
SRC="${2:-}"
REPO_URL="https://github.com/fg-labs/bwa-mem2.git"
DST="bwa-mem2-sys/vendor/bwa-mem2"

rm -rf "$DST"
mkdir -p "$(dirname "$DST")"

if [ -n "$SRC" ]; then
    # Copy from a local working tree at the given hash (submodules initialized).
    actual=$(git -C "$SRC" rev-parse HEAD)
    if [ "$actual" != "$HASH" ] && [ "${actual:0:${#HASH}}" != "$HASH" ]; then
        echo "ERROR: $SRC HEAD is $actual, not $HASH" >&2
        exit 1
    fi
    # rsync preserves submodule contents; exclude vcs, build artifacts, large tests.
    mkdir -p "$DST"
    rsync -a \
        --exclude='.git' --exclude='.github' \
        --exclude='*.o' --exclude='*.a' --exclude='obj/' --exclude='Debug/' --exclude='objtest/' \
        --exclude='/bwa-mem2' \
        --exclude='/bwa-mem2.sse41' --exclude='/bwa-mem2.sse42' \
        --exclude='/bwa-mem2.avx' --exclude='/bwa-mem2.avx2' --exclude='/bwa-mem2.avx512bw' \
        --exclude='/bwa-mem2.arm64' --exclude='/bwa-mem2.native' \
        "$SRC/" "$DST/"
else
    git clone --recurse-submodules "$REPO_URL" "$DST"
    git -C "$DST" checkout "$HASH"
    git -C "$DST" submodule update --init --recursive
    find "$DST" -name '.git' -print0 | xargs -0 rm -rf
fi

# Record the commit; used by build.rs as a sanity marker.
echo "$HASH" > bwa-mem2-sys/vendor/COMMIT

# Verify MATE_SORT=0 default (see spec "Pairing without patching upstream").
if ! grep -E '^CPPFLAGS\+?=.*-DMATE_SORT=0' "$DST/Makefile" >/dev/null; then
    echo "ERROR: vendored Makefile does not have -DMATE_SORT=0; shim semantics would diverge." >&2
    exit 1
fi

echo "Vendored $HASH at $DST"
