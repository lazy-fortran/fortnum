#!/usr/bin/env bash
# Builds fortnum's BLAS and LAPACK for WebAssembly and runs their tests under
# node, the way a browser runs them:
#   tools/wasm/test_wasm.sh [BUILD_DIR] [TEST ...]
# Needs LFortran built with WITH_TARGET_WASM (LFORTRAN, default lfortran), an
# emscripten SDK (EMSDK_PATH) and node. Without TEST names it runs every BLAS
# and LAPACK test.
set -euo pipefail
top=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out=$(mkdir -p "${1:-$top/build/wasm}" && cd "${1:-$top/build/wasm}" && pwd)
shift || true
emcc=${EMCC:-${EMSDK_PATH:?set EMSDK_PATH to the emscripten SDK}/upstream/emscripten/emcc}
lf=${LFORTRAN:-lfortran}
node=${NODE:-node}
runtime=$(dirname "$(command -v "$lf")")/../lib/lfortran_runtime_wasm_emcc.o
src=$top/src
flags=(--target=wasm32-unknown-emscripten --separate-compilation --realloc-lhs-arrays -O3
    --legacy-array-sections --cpp -I"$src/blas" -I"$src/lapack")
modules=(fortnum_kinds blas/fortnum_blas blas/fortnum_blas12 blas/fortnum_blas3
    lapack/fortnum_lapack_lu lapack/fortnum_lapack_chol lapack/fortnum_lapack_ldlt
    lapack/fortnum_lapack_qr lapack/fortnum_lapack_refine lapack/fortnum_lapack_syev
    lapack/fortnum_lapack_geev lapack/fortnum_lapack_svd lapack77/fortnum_lapack77)
tests=("$@")
if [ ${#tests[@]} -eq 0 ]; then
    tests=(test_fortnum_blas test_fortnum_blas12 test_fortnum_blas3 test_fortnum_lapack_lu
        test_fortnum_lapack_chol test_fortnum_lapack_ldlt test_fortnum_lapack_qr
        test_fortnum_lapack_refine test_fortnum_lapack_syev test_fortnum_lapack_geev
        test_fortnum_lapack_svd test_fortnum_lapack77)
fi
cd "$out"
objects=()
for m in "${modules[@]}"; do
    f=$src/$m.f90
    [ -f "$f" ] || f=$src/$m.F90
    "$lf" -c "${flags[@]}" "$f" -o "$(basename "$m").o"
    objects+=("$(basename "$m").o")
done
# the SIMD128 tile kernels and their dispatch
for f in "$src"/blas/fortnum_blas_cpu.c "$src"/generated/fortnum_dtile_wasm_*.c; do
    "$emcc" -O2 -msimd128 -c "$f" -o "$(basename "$f" .c).o"
    objects+=("$(basename "$f" .c).o")
done
failed=0
for t in "${tests[@]}"; do
    "$lf" -c "${flags[@]}" "$top/test/linalg/$t.f90" -o "$t.o"
    "$emcc" -O2 -sALLOW_MEMORY_GROWTH -sSTACK_SIZE=8mb -sENVIRONMENT=node -o "$t.js" \
        "$t.o" "${objects[@]}" "$runtime"
    if "$node" "$t.js"; then
        echo "PASS $t"
    else
        echo "FAIL $t"
        failed=1
    fi
done
exit $failed
