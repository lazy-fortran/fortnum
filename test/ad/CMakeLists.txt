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

add_executable(test_enzyme_vs_analytic test_enzyme_vs_analytic.f90)
target_link_libraries(test_enzyme_vs_analytic PRIVATE fortnum_ad_test_utils)
add_test(NAME ad_enzyme_vs_analytic COMMAND test_enzyme_vs_analytic)

# --- #38 Enzyme smoke tests (compiler-AD path; skipped when unavailable) ---
if(FORTNUM_ENABLE_ENZYME)
    include(FortnumAddEnzymeTest)
    fortnum_add_enzyme_test(NAME enzyme_square    SOURCES smoke/square.f90      MODE reverse OPT_LEVEL O2)
    fortnum_add_enzyme_test(NAME enzyme_dot       SOURCES smoke/dot_product.f90 MODE forward OPT_LEVEL O2)
    fortnum_add_enzyme_test(NAME enzyme_array_sum SOURCES smoke/array_sum.f90   MODE reverse OPT_LEVEL O2)
endif()
