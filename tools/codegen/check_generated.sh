#!/usr/bin/env bash
set -euo pipefail

codegen_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_dir=$(cd "$codegen_dir/../.." && pwd)
temporary_dir=$(mktemp -d)
trap 'rm -rf -- "$temporary_dir"' EXIT

python3 "$repository_dir/scripts/check_contracted_codegen.py" "$repository_dir"

(
    cd "$codegen_dir"
    fo build
    for generator in gen_special_outer gen_determinant_products \
        gen_inverse_products gen_multi_input_scalar gen_implicit_root_residual \
        gen_lagrange4_interpolation gen_dawson_jvp_variants \
        gen_erf_products gen_stable_sqrt_difference \
        gen_special_region_candidates gen_jacobi_recurrence \
        gen_rk54_device gen_rk54_cpu_tableau; do
        FORTNUM_CODEGEN_OUTPUT_DIR="$temporary_dir" fo exec --no-build "$generator"
    done
)

for generated in "$temporary_dir"/fortnum_*.f90; do
    cmp -- "$repository_dir/src/generated/$(basename "$generated")" "$generated"
done

echo "generated kernels match committed sources"
