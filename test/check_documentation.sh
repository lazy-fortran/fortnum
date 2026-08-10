#!/usr/bin/env bash
set -euo pipefail

source_root=${1:?repository root is required}
temp_dir=$(mktemp -d)
trap 'rm -rf "${temp_dir}"' EXIT
base="${temp_dir}/base"

mkdir -p "${base}/benchmark/report" "${base}/tools/codegen/app"
cp "${source_root}/README.md" "${source_root}/CONTRIBUTING.md" \
    "${source_root}/ROADMAP.md" "${source_root}/AGENTS.md" \
    "${source_root}/LICENSE" "${base}/"
cp -a "${source_root}/docs" "${source_root}/src" "${source_root}/scripts" \
    "${source_root}/.github" "${source_root}/include" "${base}/"
cp "${source_root}/benchmark/README.md" "${base}/benchmark/"
cp -a "${source_root}/benchmark/reference" "${base}/benchmark/"
cp -a "${source_root}/benchmark/report/data" "${base}/benchmark/report/"
cp "${source_root}/tools/codegen/fortsym.lock" "${base}/tools/codegen/"
cp "${source_root}/tools/codegen/fortsym-rk54.lock" "${base}/tools/codegen/"
cp "${source_root}/tools/codegen/app/gen_enzyme_scalar_wrappers.f90" \
    "${base}/tools/codegen/app/"

checker="${base}/scripts/check_documentation.py"
python3 "${checker}" "${base}"

expect_failure() {
    local name=$1
    local case_root=$2
    if python3 "${case_root}/scripts/check_documentation.py" \
            "${case_root}" >/dev/null 2>&1; then
        echo "documentation checker accepted ${name}" >&2
        exit 1
    fi
}

case_root="${temp_dir}/broken-link"
cp -a "${base}" "${case_root}"
printf '\n[missing](absent.md)\n' >> "${case_root}/README.md"
expect_failure "a broken link" "${case_root}"

case_root="${temp_dir}/legacy-term"
cp -a "${base}" "${case_root}"
printf '\nanalytic_rule\n' >> "${case_root}/docs/design/ad.md"
expect_failure "legacy derivative terminology" "${case_root}"

case_root="${temp_dir}/hard-coded-revision"
cp -a "${base}" "${case_root}"
printf '\n0123456789abcdef0123456789abcdef01234567\n' \
    >> "${case_root}/docs/design/ad.md"
expect_failure "a hard-coded dependency revision" "${case_root}"

case_root="${temp_dir}/stale-claim"
cp -a "${base}" "${case_root}"
printf '\nNo derivative code ships now.\n' >> "${case_root}/docs/api.md"
expect_failure "a stale future-tense claim" "${case_root}"

case_root="${temp_dir}/api-module"
cp -a "${base}" "${case_root}"
sed -i 's/`fortnum_fft`/`transform_module`/g' "${case_root}/docs/api.md"
expect_failure "an undocumented production module" "${case_root}"

case_root="${temp_dir}/report-count"
cp -a "${base}" "${case_root}"
sed -i 's/contains [0-9][0-9]* derivative/contains 0 derivative/' \
    "${case_root}/docs/design/differentiation_report.md"
expect_failure "stale aggregate benchmark statistics" "${case_root}"

case_root="${temp_dir}/generator-lock"
cp -a "${base}" "${case_root}"
sed -i 's/Generator revision: fortsym@[0-9a-f]*/Generator revision: fortsym@0000000000000000000000000000000000000000/' \
    "${case_root}/src/generated/fortnum_dawson_outer_kernel.f90"
expect_failure "a generated revision mismatch" "${case_root}"

case_root="${temp_dir}/docs-map"
cp -a "${base}" "${case_root}"
sed -i '/performance_optimal_differentiation/d' \
    "${case_root}/docs/README.md"
expect_failure "an unindexed maintained document" "${case_root}"
