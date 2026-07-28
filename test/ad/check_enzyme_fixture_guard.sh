#!/usr/bin/env bash
set -euo pipefail

repository_dir=$1
guard="$repository_dir/cmake/enzyme/support/check_no_duplicate_fixture_helpers.sh"
temporary_dir=$(mktemp -d)
trap 'rm -rf -- "$temporary_dir"' EXIT

"$guard" "$repository_dir/cmake/enzyme/hybrid"

printf '%s\n' 'program clean_fixture' 'end program clean_fixture' \
    >"$temporary_dir/fixture.f90"
"$guard" "$temporary_dir"

for forbidden in \
    '__enzyme_fwddiff' \
    'subroutine sort_values(values)' \
    'call system_clock(start)' \
    'bind(c, name="fortnum_peak_rss_bytes")'; do
    printf '%s\n' "$forbidden" >"$temporary_dir/fixture.f90"
    if "$guard" "$temporary_dir" >/dev/null 2>&1; then
        echo "guard accepted forbidden fixture helper: $forbidden" >&2
        exit 1
    fi
done
