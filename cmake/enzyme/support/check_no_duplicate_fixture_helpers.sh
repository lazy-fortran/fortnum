#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$script_dir/../../.." && pwd)
target_dir=${1:-"$repository_dir/cmake/enzyme/hybrid"}

if [[ ! -d "$target_dir" ]]; then
    echo "fixture directory does not exist: $target_dir" >&2
    exit 2
fi

forbidden='__enzyme_(fwddiff|autodiff)|subroutine[[:space:]]+sort_values|call[[:space:]]+system_clock|bind\(c,[[:space:]]*name="fortnum_peak_rss_bytes"'
if grep -R -n -E --include='*.f90' "$forbidden" "$target_dir"; then
    echo "hybrid Enzyme fixtures must use generated wrappers and shared support" >&2
    exit 1
fi
