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
    bindings/fortnum_capi_rk8pd.f90
    fft/fortnum_fft.f90
    fortnum_kinds.f90
    fortnum_status.f90
    fortnum_version.f90
    interp/fortnum_bspline.f90
    interp/fortnum_interp.f90
    interp/fortnum_polynomial.f90
    linalg/fortnum_linalg.f90
    ode/fortnum_ode.f90
    ode/fortnum_ode_cash_karp.f90
    ode/fortnum_ode_ddeabm.f90
    ode/fortnum_ode_dop853.f90
    ode/fortnum_ode_events.f90
    ode/fortnum_ode_rk8pd.f90
    ode/fortnum_ode_vode.f90
    ode/fortnum_ode_wrapper.f90
    quadrature/fortnum_cquad.f90
    quadrature/fortnum_integrate.f90
    quadrature/fortnum_integrate_gk.f90
    quadrature/fortnum_levin.f90
    quadrature/fortnum_quadrature.f90
    rng/fortnum_rng.f90
    roots/fortnum_multiroot.f90
    roots/fortnum_multiroot_rc.f90
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

# Adaptive quadrature, root polishing, and the NaN-sample handling depend on
# IEEE semantics. When fortnum is embedded via add_subdirectory/FetchContent it
# inherits the parent's directory-scope flags; a parent building with
# -ffast-math (finite-math-only, reassociation) silently deletes the NaN
# checks and shifts rule sums. Re-disable fast-math last so it wins on the
# command line, whatever the parent sets.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(fortnum PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-fno-fast-math>
        $<$<COMPILE_LANGUAGE:Fortran>:-ffp-contract=off>)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(fortnum PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:-fp-model=precise>)
endif()

# roots/fortnum_roots_complex.f90 calls LAPACK ZGGEV; link LAPACK (with BLAS)
# into the library so every consumer resolves the symbol.
target_link_libraries(fortnum PUBLIC LAPACK::LAPACK)

# Install the hand-written C ABI header so C/C++ consumers find it.
install(FILES ${PROJECT_SOURCE_DIR}/include/fortnum.h DESTINATION include)
