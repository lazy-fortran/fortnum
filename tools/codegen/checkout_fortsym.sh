#!/usr/bin/env bash
# Check out the one fortsym revision that every codegen lock is satisfied by.
#
# Usage: checkout_fortsym.sh FORTSYM_CHECKOUT
#
# fpm builds every generator in tools/codegen together against one fortsym
# tree, but generators are contracted to different revisions through the
# tools/codegen/fortsym*.lock files. The locks must lie on one line of
# history; the build uses the newest of them. Provenance stays per file:
# check_documentation.py asserts each generated kernel records the lock it is
# contracted to, not the revision it was built against.
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 FORTSYM_CHECKOUT" >&2
    exit 2
fi
fortsym=$1
codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

build_revision=""
build_lock=""
for lock in "$codegen_dir"/fortsym*.lock; do
    revision=$(tr -d '[:space:]' < "$lock")
    # A clone only carries objects reachable from a branch or tag. A locked
    # revision whose branch was squash-merged and deleted is absent, so fetch
    # it by SHA before deciding anything.
    if ! git -C "$fortsym" cat-file -e "$revision^{commit}" 2>/dev/null; then
        if ! git -C "$fortsym" fetch --depth 1 origin "$revision" \
            2>/dev/null; then
            echo "::error::fortsym revision $revision ($(basename "$lock"))" \
                "is unreachable. Its branch was probably squash-merged and" \
                "deleted; bump the lock to a revision on fortsym main."
            exit 1
        fi
    fi
    if [ -z "$build_revision" ]; then
        build_revision=$revision
        build_lock=$(basename "$lock")
    elif git -C "$fortsym" merge-base --is-ancestor "$build_revision" \
        "$revision"; then
        build_revision=$revision
        build_lock=$(basename "$lock")
    elif ! git -C "$fortsym" merge-base --is-ancestor "$revision" \
        "$build_revision"; then
        echo "::error::$(basename "$lock") ($revision) and $build_lock" \
            "($build_revision) have diverged; no single checkout satisfies" \
            "both"
        exit 1
    fi
done

if [ -z "$build_revision" ]; then
    echo "::error::no fortsym lock found in $codegen_dir" >&2
    exit 1
fi
echo "Building generators against fortsym $build_revision ($build_lock)"
git -C "$fortsym" checkout --quiet "$build_revision"
