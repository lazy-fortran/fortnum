#!/usr/bin/env bash
# Discharge the Why3 obligations on the exact derived DOP853 tableau.
#
# Requires the Why3 proof platform (why3 CLI) and at least one supported SMT
# prover (z3, cvc5, alt-ergo). The check is not part of the default CTest
# suite: Why3 is an optional, heavyweight external toolchain.
set -euo pipefail

why3_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$why3_dir/.." && pwd)
temporary_dir=$(mktemp -d)
trap 'rm -rf -- "$temporary_dir"' EXIT

# The committed theory must be the one the committed tableau generates.
python3 "$why3_dir/gen_dop853_model.py" \
    "$repository_dir/src/generated/fortnum_dop853_tableau.f90" \
    "$temporary_dir/dop853_tableau.mlw"
if ! cmp -- "$why3_dir/dop853_tableau.mlw" "$temporary_dir/dop853_tableau.mlw"
then
    echo "why3/dop853_tableau.mlw is stale; regenerate it with" \
        "why3/gen_dop853_model.py" >&2
    exit 1
fi

if ! command -v why3 >/dev/null 2>&1; then
    echo "why3 CLI not found; install Why3 and an SMT prover to discharge." >&2
    exit 2
fi

provers=""
for p in z3 cvc5 alt-ergo; do
    if command -v "$p" >/dev/null 2>&1; then
        provers="$provers $p"
    fi
done
if [ -z "$provers" ]; then
    echo "no SMT prover found; install z3, cvc5, or alt-ergo." >&2
    exit 2
fi

failed=0
for p in $provers; do
    echo "discharging with $p ..."
    # why3 prove exits non-zero unless every goal is Valid.
    if ! why3 prove -P "$p" -t 60 "$why3_dir/dop853_tableau.mlw"; then
        echo "why3 ($p) did not discharge every goal" >&2
        failed=1
    fi
done
exit "$failed"
