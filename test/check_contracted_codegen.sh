#!/usr/bin/env bash
set -euo pipefail

source_root=${1:?repository root is required}
checker="${source_root}/scripts/check_contracted_codegen.py"
python3 "${checker}" "${source_root}"

temp_dir=$(mktemp -d)
trap 'rm -rf "${temp_dir}"' EXIT
mkdir -p "${temp_dir}/tools/codegen/app" "${temp_dir}/src/generated" \
    "${temp_dir}/scripts"
cp "${checker}" "${temp_dir}/scripts/"
cp "${source_root}"/tools/codegen/app/*.f90 "${temp_dir}/tools/codegen/app/"
cp "${source_root}"/src/generated/*_jvp_kernel.f90 \
    "${source_root}"/src/generated/*_vjp_kernel.f90 \
    "${temp_dir}/src/generated/"

printf '\ncall jacobian(x)\n' >> \
    "${temp_dir}/tools/codegen/app/gen_erf_products.f90"
if python3 "${temp_dir}/scripts/check_contracted_codegen.py" \
        "${temp_dir}" >/dev/null 2>&1; then
    echo "contracted-codegen checker accepted Jacobian materialization" >&2
    exit 1
fi
sed -i '$d' "${temp_dir}/tools/codegen/app/gen_erf_products.f90"
sed -i 's/Generator: gen_erf_products/Generator: absent_generator/' \
    "${temp_dir}/src/generated/fortnum_erf_jvp_kernel.f90"
if python3 "${temp_dir}/scripts/check_contracted_codegen.py" \
        "${temp_dir}" >/dev/null 2>&1; then
    echo "contracted-codegen checker accepted an absent generator" >&2
    exit 1
fi
sed -i 's/Generator: absent_generator/Generator: gen_erf_products/' \
    "${temp_dir}/src/generated/fortnum_erf_jvp_kernel.f90"
sed -i 's/vjp(/not_vjp(/g' \
    "${temp_dir}/tools/codegen/app/gen_erf_products.f90"
if python3 "${temp_dir}/scripts/check_contracted_codegen.py" \
        "${temp_dir}" >/dev/null 2>&1; then
    echo "contracted-codegen checker accepted a non-contracted VJP" >&2
    exit 1
fi
