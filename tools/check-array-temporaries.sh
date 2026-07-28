#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/.." && pwd)
baseline="$repo_root/tools/array_temporary_baseline.csv"
work_dir=$(mktemp -d)
trap 'rm -rf "$work_dir"' EXIT

set +e
(cd "$repo_root" && fo lint) >"$work_dir/lint.out" 2>&1
lint_status=$?
set -e
if (( lint_status > 1 )); then
    cat "$work_dir/lint.out"
    exit "$lint_status"
fi

rg 'Warning: Creating array temporary' "$work_dir/lint.out" |
    awk -v prefix="$repo_root/" '
        {
            marker = match($0, /:[0-9]+:[0-9]+:/)
            if (marker == 0) next
            path = substr($0, 1, marker - 1)
            sub("^" prefix, "", path)
            count[path]++
        }
        END {
            for (path in count) print path "," count[path]
        }
    ' |
    sort >"$work_dir/actual.csv"

if ! diff -u "$baseline" "$work_dir/actual.csv"; then
    echo "array-temporary baseline drifted; inspect performance before updating" >&2
    exit 1
fi

warning_count=$(awk -F, '{total += $2} END {print total + 0}' "$baseline")
echo "PASS: $warning_count array-temporary warnings match the baseline"
