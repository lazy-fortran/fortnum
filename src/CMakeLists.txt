# fortnum core library. Downstream agents append one module source per line
# under the "module sources" block below, kept alphabetical and one-per-line
# to minimise merge conflicts.
add_library(fortnum STATIC)

target_sources(fortnum PRIVATE
    # module sources (alphabetical, one per line; append here)
    ad/fortnum_active_vector.f90
    ad/fortnum_ad_interfaces.f90
    bindings/fortnum_capi.f90
    bindings/fortnum_capi_bspline.f90
    fft/fortnum_fft.f90
    fortnum_kinds.f90
    fortnum_status.f90
    fortnum_version.f90
    interp/fortnum_bspline.f90
    interp/fortnum_interp.f90
    interp/fortnum_polynomial.f90
    ode/fortnum_ode.f90
    ode/fortnum_ode_cash_karp.f90
    ode/fortnum_ode_dop853.f90
    ode/fortnum_ode_events.f90
    ode/fortnum_ode_wrapper.f90
    quadrature/fortnum_integrate.f90
    quadrature/fortnum_integrate_gk.f90
    quadrature/fortnum_levin.f90
    quadrature/fortnum_quadrature.f90
    rng/fortnum_rng.f90
    roots/fortnum_multiroot.f90
    roots/fortnum_roots.f90
    roots/fortnum_roots_complex.f90
    special/fortnum_special.f90
    special/fortnum_special_bessel.f90
    special/fortnum_special_complex_bessel.f90
    special/fortnum_special_dawson.f90
    special/fortnum_special_erf_cbind.f90
    special/fortnum_special_gamma.f90
    special/fortnum_special_hypergeometric_1f1.f90
    testing/fortnum_oracle.f90
)

# Consumers find .mod files in the module output directory.
target_include_directories(fortnum PUBLIC ${CMAKE_Fortran_MODULE_DIRECTORY})

# roots/fortnum_roots_complex.f90 calls LAPACK ZGGEV; link LAPACK (with BLAS)
# into the library so every consumer resolves the symbol.
target_link_libraries(fortnum PUBLIC LAPACK::LAPACK)

# Install the hand-written C ABI header so C/C++ consumers find it. Guard on a
# standalone build: as a FetchContent subproject CMAKE_SOURCE_DIR points at the
# consumer, and the consumer's install tree should not carry fortnum's header.
if(FORTNUM_IS_STANDALONE)
    install(FILES ${PROJECT_SOURCE_DIR}/include/fortnum.h DESTINATION include)
endif()
