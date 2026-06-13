# Self-contained benchmark build for fortnum.
#
# Standalone by design so it never touches src/CMakeLists.txt or
# test/CMakeLists.txt.  Configure on its own:
#
#   cmake -S benchmark -B build-bench -G Ninja
#   cmake --build build-bench -j
#   ./build-bench/bin/bench_main
#
# It compiles the small set of fortnum core sources it depends on directly,
# rather than linking an installed library, to keep the build hermetic.

cmake_minimum_required(VERSION 3.22)

project(fortnum_benchmark LANGUAGES Fortran)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Core fortnum sources the benchmark exercises.  Add more here only as the
# seeded kernels grow to need them; the test/src builds are untouched.
set(FORTNUM_CORE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)

add_executable(bench_main
    ${FORTNUM_CORE_DIR}/fortnum_version.f90
    fortnum_benchmark.f90
    bench_main.f90
)

# Derivative-overhead benchmark (#39). Seeded with a placeholder kernel; #40
# adds real JVP/VJP/grad/HVP entries as the derivative products land.
add_executable(bench_ad_overheads
    ${FORTNUM_CORE_DIR}/fortnum_version.f90
    fortnum_benchmark.f90
    bench_ad_overheads.f90
)
