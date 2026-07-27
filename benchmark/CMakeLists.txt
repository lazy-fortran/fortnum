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

project(fortnum_benchmark LANGUAGES Fortran C)

enable_testing()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake")
include(FortnumGpuBackend)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)

# Core fortnum sources the benchmark exercises.  Add more here only as the
# seeded kernels grow to need them; the test/src builds are untouched.
set(FORTNUM_CORE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../src)
set(FORTNUM_MULTI_INPUT_SOURCES
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p2_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p2_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p4_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p4_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p8_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p8_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p16_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_multi_input_p16_vjp_kernel.f90
)

# Shared core compiled once into a static library. Building it per-executable
# made two targets emit the same .mod into the shared module directory, which
# ninja rejects ("multiple rules generate include/fortnum_version.mod").
add_library(bench_core STATIC
    ${FORTNUM_CORE_DIR}/ad/fortnum_build_selection.f90
    ${FORTNUM_CORE_DIR}/fortnum_kinds.f90
    ${FORTNUM_CORE_DIR}/fortnum_status.f90
    ${FORTNUM_CORE_DIR}/fortnum_version.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_value_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_value_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_dawson_outer_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_det2_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_det2_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_det3_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_det3_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_inv2_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_inv2_vjp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_inv3_jvp_kernel.f90
    ${FORTNUM_CORE_DIR}/generated/fortnum_inv3_vjp_kernel.f90
    ${FORTNUM_MULTI_INPUT_SOURCES}
    ${FORTNUM_CORE_DIR}/quadrature/fortnum_integrate_gk.f90
    ${FORTNUM_CORE_DIR}/quadrature/fortnum_integrate.f90
    ${FORTNUM_CORE_DIR}/linalg/fortnum_linalg.f90
    ${FORTNUM_CORE_DIR}/interp/fortnum_polynomial.f90
    ${FORTNUM_CORE_DIR}/interp/fortnum_bspline.f90
    ${FORTNUM_CORE_DIR}/interp/fortnum_interp.f90
    ${FORTNUM_CORE_DIR}/roots/fortnum_fixed_point.f90
    ${FORTNUM_CORE_DIR}/roots/fortnum_multiroot.f90
    ${FORTNUM_CORE_DIR}/roots/fortnum_roots.f90
    gpu/fortnum_gpu_batch_wrappers.f90
    gpu/fortnum_gpu_multi_input_wrappers.f90
    gpu/test_gpu_multi_input_support.f90
    fortnum_benchmark.f90
    fortnum_benchmark_memory.f90
    peak_rss.c
)
target_link_libraries(bench_core PUBLIC fortnum_gpu_backend)

add_library(bench_multi_input_cpu_core STATIC
    ${FORTNUM_CORE_DIR}/fortnum_kinds.f90
    ${FORTNUM_MULTI_INPUT_SOURCES}
    gpu/fortnum_gpu_multi_input_wrappers.f90
    fortnum_benchmark_memory.f90
    peak_rss.c
)
set_target_properties(bench_multi_input_cpu_core PROPERTIES
    Fortran_MODULE_DIRECTORY
        ${PROJECT_BINARY_DIR}/include-multi-input-cpu)
target_include_directories(bench_multi_input_cpu_core PUBLIC
    ${PROJECT_BINARY_DIR}/include-multi-input-cpu)

add_executable(bench_main bench_main.f90)
target_link_libraries(bench_main PRIVATE bench_core)

# Derivative-overhead benchmark (#39). Seeded with a placeholder kernel; #40
# adds real JVP/VJP/grad/HVP entries as the derivative products land.
add_executable(bench_ad_overheads bench_ad_overheads.f90)
target_link_libraries(bench_ad_overheads PRIVATE bench_core)

add_executable(bench_dawson_generated_family
    bench_dawson_generated_family.f90)
target_link_libraries(bench_dawson_generated_family PRIVATE bench_core)

add_executable(bench_gpu_batch_wrapper
    gpu/bench_gpu_batch_wrapper.f90)
target_link_libraries(bench_gpu_batch_wrapper PRIVATE bench_core)

if(FORTNUM_GPU_BACKEND STREQUAL "OPENACC")
    add_executable(test_openacc_dawson_offload
        gpu/test_openacc_dawson_offload.f90)
    target_link_libraries(test_openacc_dawson_offload PRIVATE bench_core)
    add_test(NAME openacc_dawson_offload
        COMMAND test_openacc_dawson_offload)

    add_executable(test_openacc_multi_input_offload
        gpu/test_openacc_multi_input_offload.f90)
    target_link_libraries(test_openacc_multi_input_offload PRIVATE bench_core)
    add_test(NAME openacc_multi_input_offload
        COMMAND test_openacc_multi_input_offload)

    add_executable(bench_openacc_dawson_offload
        gpu/bench_openacc_dawson_offload.f90)
    target_link_libraries(bench_openacc_dawson_offload PRIVATE bench_core)
elseif(FORTNUM_GPU_BACKEND STREQUAL "OPENMP")
    add_executable(test_openmp_dawson_offload
        gpu/test_openmp_dawson_offload.f90)
    target_link_libraries(test_openmp_dawson_offload PRIVATE bench_core)
    add_test(NAME openmp_dawson_offload
        COMMAND test_openmp_dawson_offload)

    add_executable(test_openmp_multi_input_offload
        gpu/test_openmp_multi_input_offload.f90)
    target_link_libraries(test_openmp_multi_input_offload PRIVATE bench_core)
    add_test(NAME openmp_multi_input_offload
        COMMAND test_openmp_multi_input_offload)

    add_executable(bench_openmp_dawson_offload
        gpu/bench_openmp_dawson_offload.f90)
    target_link_libraries(bench_openmp_dawson_offload PRIVATE bench_core)
endif()

if(NOT FORTNUM_GPU_BACKEND STREQUAL "NONE")
    add_executable(bench_dawson_gpu_fusion
        gpu/bench_dawson_gpu_fusion.f90)
    target_link_libraries(bench_dawson_gpu_fusion PRIVATE bench_core)

    add_executable(bench_multi_input_gpu
        gpu/bench_multi_input_scaling.f90)
    target_link_libraries(bench_multi_input_gpu PRIVATE bench_core)

    add_executable(bench_multi_input_layout_gpu
        gpu/bench_multi_input_layout.f90)
    target_link_libraries(bench_multi_input_layout_gpu PRIVATE bench_core)
endif()

add_executable(bench_multi_input_cpu
    gpu/bench_multi_input_scaling.f90)
target_link_libraries(bench_multi_input_cpu PRIVATE
    bench_multi_input_cpu_core)

add_executable(bench_multi_input_layout_cpu
    gpu/bench_multi_input_layout.f90)
target_link_libraries(bench_multi_input_layout_cpu PRIVATE
    bench_multi_input_cpu_core)

add_executable(bench_linear_solve_jvp bench_linear_solve_jvp.f90)
target_link_libraries(bench_linear_solve_jvp PRIVATE bench_core)

add_executable(bench_linear_solve_vjp bench_linear_solve_vjp.f90)
target_link_libraries(bench_linear_solve_vjp PRIVATE bench_core)

add_executable(bench_lagrange_nodes bench_lagrange_nodes.f90)
target_link_libraries(bench_lagrange_nodes PRIVATE bench_core)

add_executable(bench_bspline_knots bench_bspline_knots.f90)
target_link_libraries(bench_bspline_knots PRIVATE bench_core)

add_executable(bench_bspline_fit bench_bspline_fit.f90)
target_link_libraries(bench_bspline_fit PRIVATE bench_core)

add_executable(bench_interp_cell_status bench_interp_cell_status.f90)
target_link_libraries(bench_interp_cell_status PRIVATE bench_core)

add_executable(bench_bspline_span_status bench_bspline_span_status.f90)
target_link_libraries(bench_bspline_span_status PRIVATE bench_core)

add_executable(bench_lagrange_combined bench_lagrange_combined.f90)
target_link_libraries(bench_lagrange_combined PRIVATE bench_core)

add_executable(bench_bspline_combined bench_bspline_combined.f90)
target_link_libraries(bench_bspline_combined PRIVATE bench_core)

add_executable(bench_multiroot_preconditioner bench_multiroot_preconditioner.f90)
target_link_libraries(bench_multiroot_preconditioner PRIVATE bench_core)

add_executable(bench_multiroot_jacobian_reuse bench_multiroot_jacobian_reuse.f90)
target_link_libraries(bench_multiroot_jacobian_reuse PRIVATE bench_core)

add_executable(bench_multiroot_jvp_factorization
    bench_multiroot_jvp_factorization.f90)
target_link_libraries(bench_multiroot_jvp_factorization PRIVATE bench_core)

add_executable(bench_multiroot_vjp_factorization
    bench_multiroot_vjp_factorization.f90)
target_link_libraries(bench_multiroot_vjp_factorization PRIVATE bench_core)

add_executable(bench_multiroot_implicit_jvp_reliability
    bench_multiroot_implicit_jvp_reliability.f90)
target_link_libraries(bench_multiroot_implicit_jvp_reliability PRIVATE bench_core)

add_executable(bench_multiroot_implicit_vjp_reliability
    bench_multiroot_implicit_vjp_reliability.f90)
target_link_libraries(bench_multiroot_implicit_vjp_reliability PRIVATE bench_core)

add_executable(bench_multiroot_condition bench_multiroot_condition.f90)
target_link_libraries(bench_multiroot_condition PRIVATE bench_core)

add_executable(bench_integrate_fixed_jvp bench_integrate_fixed_jvp.f90)
target_link_libraries(bench_integrate_fixed_jvp PRIVATE bench_core)

add_executable(bench_integrate_frozen_jvp bench_integrate_frozen_jvp.f90)
target_link_libraries(bench_integrate_frozen_jvp PRIVATE bench_core)

add_executable(bench_integrate_moving_lower_jvp
    bench_integrate_moving_lower_jvp.f90)
target_link_libraries(bench_integrate_moving_lower_jvp PRIVATE bench_core)

add_executable(bench_integrate_moving_upper_jvp
    bench_integrate_moving_upper_jvp.f90)
target_link_libraries(bench_integrate_moving_upper_jvp PRIVATE bench_core)

add_executable(bench_scalar_root_tangent bench_scalar_root_tangent.f90)
target_link_libraries(bench_scalar_root_tangent PRIVATE bench_core)

add_executable(bench_scalar_root_adjoint bench_scalar_root_adjoint.f90)
target_link_libraries(bench_scalar_root_adjoint PRIVATE bench_core)

add_executable(bench_vector_root_tangent bench_vector_root_tangent.f90)
target_link_libraries(bench_vector_root_tangent PRIVATE bench_core)

add_executable(bench_vector_root_adjoint bench_vector_root_adjoint.f90)
target_link_libraries(bench_vector_root_adjoint PRIVATE bench_core)

add_executable(bench_fixed_point_tangent bench_fixed_point_tangent.f90)
target_link_libraries(bench_fixed_point_tangent PRIVATE bench_core)

add_executable(bench_fixed_point_adjoint bench_fixed_point_adjoint.f90)
target_link_libraries(bench_fixed_point_adjoint PRIVATE bench_core)

add_executable(bench_static_registry bench_static_registry.f90)
target_link_libraries(bench_static_registry PRIVATE bench_core)

add_executable(test_peak_rss test_peak_rss.f90)
target_link_libraries(test_peak_rss PRIVATE bench_core)
add_test(NAME benchmark_peak_rss COMMAND test_peak_rss)
