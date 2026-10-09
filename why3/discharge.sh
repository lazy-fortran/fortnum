#!/usr/bin/env bash
# Discharge the Why3 contracts on the DOP853 numerical API.
#
# Requires the Why3 proof platform (why3 CLI) and at least one supported SMT
# prover (alt-ergo, z3, cvc5). The check is intentionally not part of the
# default CTest suite: Why3 is an optional, heavyweight external toolchain.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

if ! command -v why3 >/dev/null 2>&1; then
    echo "why3 CLI not found; install Why3 and an SMT prover to discharge." >&2
    exit 2
fi

provers=""
for p in alt-ergo z3 cvc5 cvc4; do
    if command -v "$p" >/dev/null 2>&1; then
        provers="$provers $p"
    fi
done
if [ -z "$provers" ]; then
    echo "no SMT prover found; install alt-ergo, z3, or cvc5." >&2
    exit 2
fi

failed=0
for p in $provers; do
    echo "discharging with $p ..."
    if ! why3 prove dop853_step.mlw -P "$p"; then
        echo "why3 ($p) did not discharge every goal" >&2
        failed=1
    fi
done
exit "$failed"
