# AD derivative-test harness (#39) plus the optional Flang+Enzyme path (#38).
# The harness self-tests run with plain gfortran. The Enzyme smoke tests build
# only when FORTNUM_ENABLE_ENZYME=ON and the pass plugin resolves; otherwise
# they are skipped so the suite stays green without Enzyme installed.
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(FortnumEnzyme)

# --- #39 derivative-test harness ---
add_library(fortnum_ad_test_utils STATIC fortnum_ad_test_utils.f90)
target_link_libraries(fortnum_ad_test_utils PUBLIC fortnum)
set_source_files_properties(fortnum_ad_test_utils.f90 PROPERTIES Fortran_PREPROCESS ON)
if(FORTNUM_ENABLE_ENZYME)
    target_compile_definitions(fortnum_ad_test_utils PRIVATE FORTNUM_ENABLE_ENZYME)
endif()

add_executable(test_finite_difference test_finite_difference.f90)
target_link_libraries(test_finite_difference PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_finite_difference COMMAND test_finite_difference)

add_executable(test_complex_step test_complex_step.f90)
target_link_libraries(test_complex_step PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_complex_step COMMAND test_complex_step)

add_executable(test_dot_products test_dot_products.f90)
target_link_libraries(test_dot_products PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_dot_products COMMAND test_dot_products)

# --- #40 module derivative sweep: analytic JVP/VJP/grad, verified by the
#     harness (finite difference + dot-product adjoint identity) ---
add_executable(test_special_ad test_special_ad.f90)
target_link_libraries(test_special_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_special COMMAND test_special_ad)

add_executable(test_complex_bessel_ad test_complex_bessel_ad.f90)
target_link_libraries(test_complex_bessel_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_complex_bessel COMMAND test_complex_bessel_ad)

add_executable(test_hyperg_ad test_hyperg_ad.f90)
target_link_libraries(test_hyperg_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_hyperg COMMAND test_hyperg_ad)

add_executable(test_erf_cbind_ad test_erf_cbind_ad.f90)
target_link_libraries(test_erf_cbind_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_erf_cbind COMMAND test_erf_cbind_ad)

add_executable(test_bspline_ad test_bspline_ad.f90)
target_link_libraries(test_bspline_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_bspline COMMAND test_bspline_ad)

add_executable(test_fft_ad test_fft_ad.f90)
target_link_libraries(test_fft_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_fft COMMAND test_fft_ad)

add_executable(test_quadrature_ad test_quadrature_ad.f90)
target_link_libraries(test_quadrature_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_quadrature COMMAND test_quadrature_ad)

add_executable(test_integrate_ad test_integrate_ad.f90)
target_link_libraries(test_integrate_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_integrate COMMAND test_integrate_ad)

add_executable(test_ode_ad test_ode_ad.f90)
target_link_libraries(test_ode_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_ode COMMAND test_ode_ad)

add_executable(test_roots_ad test_roots_ad.f90)
target_link_libraries(test_roots_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_roots COMMAND test_roots_ad)

add_executable(test_fixed_point_ad test_fixed_point_ad.f90)
target_link_libraries(test_fixed_point_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_fixed_point COMMAND test_fixed_point_ad)

add_executable(test_interp_ad test_interp_ad.f90)
target_link_libraries(test_interp_ad PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_interp COMMAND test_interp_ad)

# --- #41 optimizer-facing API ---
add_executable(test_optimizer_api test_optimizer_api.f90)
target_link_libraries(test_optimizer_api PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_optimizer_api COMMAND test_optimizer_api)

add_executable(test_derivative_registry test_derivative_registry.f90)
target_link_libraries(test_derivative_registry PRIVATE fortnum)
add_test(NAME ad_derivative_registry COMMAND test_derivative_registry)

# --- #38 Enzyme smoke tests (compiler-AD path; skipped when unavailable) ---
if(FORTNUM_ENABLE_ENZYME)
    include(FortnumAddEnzymeTest)
    fortnum_add_enzyme_test(
        NAME enzyme_square
        SOURCES ${CMAKE_SOURCE_DIR}/cmake/enzyme/smoke/square.f90
        MODE reverse
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_dot
        SOURCES ${CMAKE_SOURCE_DIR}/cmake/enzyme/smoke/dot_product.f90
        MODE forward
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_array_sum
        SOURCES ${CMAKE_SOURCE_DIR}/cmake/enzyme/smoke/array_sum.f90
        MODE reverse
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_dawson_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/generated/fortnum_dawson_outer_kernel.f90
            ${FORTNUM_BUILD_SELECTION_SOURCE}
            ${CMAKE_SOURCE_DIR}/src/special/fortnum_special_dawson.f90
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/dawson_rule.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/dawson_outer.f90
        FORWARD_RULES
            fortnum_dawson_kernel=fortnum_dawson_kernel_derivative
        MODE forward
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_scalar_root_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/fortnum_status.f90
            ${CMAKE_SOURCE_DIR}/src/roots/fortnum_roots.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/scalar_root_jvp.f90
        MODE forward
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_scalar_root_vjp_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/fortnum_status.f90
            ${CMAKE_SOURCE_DIR}/src/roots/fortnum_roots.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/scalar_root_vjp.f90
        MODE reverse
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_vector_root_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/fortnum_status.f90
            ${CMAKE_SOURCE_DIR}/src/linalg/fortnum_linalg.f90
            ${CMAKE_SOURCE_DIR}/src/roots/fortnum_multiroot.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/vector_root_jvp.f90
        MODE forward
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_vector_root_vjp_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/fortnum_status.f90
            ${CMAKE_SOURCE_DIR}/src/linalg/fortnum_linalg.f90
            ${CMAKE_SOURCE_DIR}/src/roots/fortnum_multiroot.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/vector_root_vjp.f90
        MODE reverse
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_fixed_quadrature_jvp_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/quadrature/fortnum_quadrature.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/fixed_quadrature_jvp.f90
        MODE forward
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_fixed_quadrature_vjp_hybrid
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/quadrature/fortnum_quadrature.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/fixed_quadrature_vjp.f90
        MODE reverse
        OPT_LEVEL O2)
    fortnum_add_enzyme_test(
        NAME enzyme_adaptive_frozen_trace_jvp
        SOURCES
            ${CMAKE_SOURCE_DIR}/src/fortnum_kinds.f90
            ${CMAKE_SOURCE_DIR}/src/fortnum_status.f90
            ${CMAKE_SOURCE_DIR}/src/quadrature/fortnum_integrate_gk.f90
            ${CMAKE_SOURCE_DIR}/src/quadrature/fortnum_integrate.f90
            ${CMAKE_SOURCE_DIR}/benchmark/peak_rss.c
            ${CMAKE_SOURCE_DIR}/cmake/enzyme/hybrid/adaptive_frozen_trace_jvp.f90
        MODE forward
        OPT_LEVEL O2)
endif()
