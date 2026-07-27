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

# Shared core compiled once into a static library. Building it per-executable
# made two targets emit the same .mod into the shared module directory, which
# ninja rejects ("multiple rules generate include/fortnum_version.mod").
add_library(bench_core STATIC
    ${FORTNUM_CORE_DIR}/ad/fortnum_build_selection.f90
    ${FORTNUM_CORE_DIR}/fortnum_kinds.f90
    ${FORTNUM_CORE_DIR}/fortnum_status.f90
    ${FORTNUM_CORE_DIR}/fortnum_version.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_kernel.f90
    ${FORTNUM_CORE_DIR}/linalg/fortnum_linalg.f90
    ${FORTNUM_CORE_DIR}/roots/fortnum_multiroot.f90
    fortnum_benchmark.f90
)

add_executable(bench_main bench_main.f90)
target_link_libraries(bench_main PRIVATE bench_core)

# Derivative-overhead benchmark (#39). Seeded with a placeholder kernel; #40
# adds real JVP/VJP/grad/HVP entries as the derivative products land.
add_executable(bench_ad_overheads bench_ad_overheads.f90)
target_link_libraries(bench_ad_overheads PRIVATE bench_core)

add_executable(bench_linear_solve_jvp bench_linear_solve_jvp.f90)
target_link_libraries(bench_linear_solve_jvp PRIVATE bench_core)

add_executable(bench_linear_solve_vjp bench_linear_solve_vjp.f90)
target_link_libraries(bench_linear_solve_vjp PRIVATE bench_core)

add_executable(bench_multiroot_preconditioner bench_multiroot_preconditioner.f90)
target_link_libraries(bench_multiroot_preconditioner PRIVATE bench_core)

add_executable(bench_multiroot_condition bench_multiroot_condition.f90)
target_link_libraries(bench_multiroot_condition PRIVATE bench_core)

add_executable(bench_static_registry bench_static_registry.f90)
target_link_libraries(bench_static_registry PRIVATE bench_core)
