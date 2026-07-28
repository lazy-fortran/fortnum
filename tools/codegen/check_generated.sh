#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
temporary_dir=$(mktemp -d)
trap 'rm -rf -- "$temporary_dir"' EXIT

(
    cd "$codegen_dir"
    fo build
    for generator in gen_dawson_outer gen_determinant_products \
        gen_inverse_products gen_multi_input_scalar gen_implicit_root_residual \
        gen_lagrange4_interpolation gen_dawson_jvp_variants \
        gen_erf_products; do
        FORTNUM_CODEGEN_OUTPUT_DIR="$temporary_dir" fo exec --no-build "$generator"
    done
)

for generated in "$repository_dir"/src/generated/fortnum_*_kernel.f90; do
    cmp -- "$generated" "$temporary_dir/$(basename "$generated")"
done

echo "generated kernels match committed sources"
