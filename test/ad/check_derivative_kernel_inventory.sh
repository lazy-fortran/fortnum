#!/usr/bin/env bash
set -euo pipefail

repo_root=${1:?repository root is required}
checker="${repo_root}/scripts/check_derivative_kernel_inventory.py"
inventory="${repo_root}/docs/design/derivative_kernel_inventory.csv"
temp_dir=$(mktemp -d)
trap 'rm -rf "${temp_dir}"' EXIT

python3 "${checker}" --root "${repo_root}"

cp "${inventory}" "${temp_dir}/missing.csv"
sed -i '/fortnum_dawson_identity_jvp_kernel/d' "${temp_dir}/missing.csv"
if python3 "${checker}" --root "${repo_root}" \
        --inventory "${temp_dir}/missing.csv" >/dev/null 2>&1; then
    echo "checker accepted an unclassified production kernel" >&2
    exit 1
fi

cp "${inventory}" "${temp_dir}/class.csv"
sed -i '2s/fortsym-generated/unclassified/' "${temp_dir}/class.csv"
if python3 "${checker}" --root "${repo_root}" \
        --inventory "${temp_dir}/class.csv" >/dev/null 2>&1; then
    echo "checker accepted an invalid implementation class" >&2
    exit 1
fi

cp "${inventory}" "${temp_dir}/stale.csv"
sed -i 's/fft_c2c_jvp/absent_jvp/' "${temp_dir}/stale.csv"
if python3 "${checker}" --root "${repo_root}" \
        --inventory "${temp_dir}/stale.csv" >/dev/null 2>&1; then
    echo "checker accepted a stale source symbol" >&2
    exit 1
fi

cp "${inventory}" "${temp_dir}/wrapper.csv"
sed -i 's/;fortnum_enzyme_square.f90//' "${temp_dir}/wrapper.csv"
if python3 "${checker}" --root "${repo_root}" \
        --inventory "${temp_dir}/wrapper.csv" >/dev/null 2>&1; then
    echo "checker accepted an unclassified generated hybrid boundary" >&2
    exit 1
fi
