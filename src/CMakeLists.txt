# fortnum core library. Downstream agents append one module source per line
# under the "module sources" block below, kept alphabetical and one-per-line
# to minimise merge conflicts.
add_library(fortnum STATIC)

target_sources(fortnum PRIVATE
    ${FORTNUM_BUILD_SELECTION_SOURCE}
    # module sources (alphabetical, one per line; append here)
    ad/fortnum_active_vector.f90
    ad/fortnum_ad_interfaces.f90
    ad/fortnum_derivative_registry.f90
    bindings/fortnum_capi.f90
    bindings/fortnum_capi_bspline.f90
    bindings/fortnum_capi_rk8pd.f90
    fft/fortnum_fft.f90
    fft/fortnum_fft8_kernel.f90
    fortnum_ad_backend.f90
    fortnum_kinds.f90
    fortnum_status.f90
    fortnum_version.f90
    generated/fortad/dawson_outer_value.f90
    generated/fortad/det2.f90
    generated/fortad/det3.f90
    generated/fortad/erf.f90
    generated/fortad/erfc.f90
    generated/fortad/fortnum_dawson_outer_value_jvp_fortad.f90
    generated/fortad/fortnum_dawson_outer_value_vjp_fortad.f90
    generated/fortad/fortnum_det2_jvp_fortad.f90
    generated/fortad/fortnum_det2_vjp_fortad.f90
    generated/fortad/fortnum_det3_jvp_fortad.f90
    generated/fortad/fortnum_det3_vjp_fortad.f90
    generated/fortad/fortnum_dot_sin_primal.f90
    generated/fortad/fortnum_erf_jvp_fortad.f90
    generated/fortad/fortnum_erf_vjp_fortad.f90
    generated/fortad/fortnum_erfc_jvp_fortad.f90
    generated/fortad/fortnum_erfc_vjp_fortad.f90
    generated/fortad/fortnum_hypergeom_2f1_term_jvp_fortad.f90
    generated/fortad/fortnum_hypergeom_2f1_term_vjp_fortad.f90
    generated/fortad/fortnum_implicit_root_residual_jvp_fortad.f90
    generated/fortad/fortnum_implicit_root_residual_vjp_fortad.f90
    generated/fortad/fortnum_inv2_jvp_fortad.f90
    generated/fortad/fortnum_lagrange4_jvp_fortad.f90
    generated/fortad/fortnum_lagrange4_vjp_fortad.f90
    generated/fortad/fortnum_legendre_recurrence_derivative_jvp_fortad.f90
    generated/fortad/fortnum_legendre_recurrence_derivative_vjp_fortad.f90
    generated/fortad/fortnum_legendre_recurrence_jvp_fortad.f90
    generated/fortad/fortnum_legendre_recurrence_vjp_fortad.f90
    generated/fortad/fortnum_multi_input_p16_jvp_fortad.f90
    generated/fortad/fortnum_multi_input_p16_vjp_fortad.f90
    generated/fortad/fortnum_multi_input_p2_jvp_fortad.f90
    generated/fortad/fortnum_multi_input_p2_vjp_fortad.f90
    generated/fortad/fortnum_multi_input_p4_jvp_fortad.f90
    generated/fortad/fortnum_multi_input_p4_vjp_fortad.f90
    generated/fortad/fortnum_multi_input_p8_jvp_fortad.f90
    generated/fortad/fortnum_multi_input_p8_vjp_fortad.f90
    generated/fortad/fortnum_scalar_root_residual_jvp_fortad.f90
    generated/fortad/fortnum_scalar_root_residual_vjp_fortad.f90
    generated/fortad/fortnum_scaled_jacobi_recurrence_jvp_fortad.f90
    generated/fortad/fortnum_scaled_jacobi_recurrence_vjp_fortad.f90
    generated/fortad/fortnum_sqrt1pm1_raw_jvp_fortad.f90
    generated/fortad/fortnum_sqrt1pm1_raw_vjp_fortad.f90
    generated/fortad/fortnum_sqrt1pm1_stable_jvp_fortad.f90
    generated/fortad/fortnum_sqrt1pm1_stable_vjp_fortad.f90
    generated/fortad/fortnum_toroidal_order_jvp_fortad.f90
    generated/fortad/fortnum_toroidal_order_vjp_fortad.f90
    generated/fortad/hypergeom_2f1_term.f90
    generated/fortad/implicit_root_residual.f90
    generated/fortad/inv2.f90
    generated/fortad/lagrange4.f90
    generated/fortad/legendre_recurrence.f90
    generated/fortad/multi_input_p16.f90
    generated/fortad/multi_input_p2.f90
    generated/fortad/multi_input_p4.f90
    generated/fortad/multi_input_p8.f90
    generated/fortad/scalar_root_residual.f90
    generated/fortad/scaled_jacobi_recurrence.f90
    generated/fortad/sqrt1pm1_raw.f90
    generated/fortad/sqrt1pm1_stable.f90
    generated/fortad/toroidal_order.f90
    generated/fortnum_bessel_outer_jvp_kernel.f90
    generated/fortnum_dawson_outer_kernel.f90
    generated/fortnum_dawson_outer_jvp_kernel.f90
    generated/fortnum_dawson_outer_value_kernel.f90
    generated/fortnum_dawson_outer_value_vjp_kernel.f90
    generated/fortnum_dawson_outer_vjp_kernel.f90
    generated/fortnum_det2_jvp_kernel.f90
    generated/fortnum_det2_vjp_kernel.f90
    generated/fortnum_det3_jvp_kernel.f90
    generated/fortnum_det3_vjp_kernel.f90
    generated/fortnum_erf_jvp_kernel.f90
    generated/fortnum_erf_vjp_kernel.f90
    generated/fortnum_erfc_jvp_kernel.f90
    generated/fortnum_erfc_vjp_kernel.f90
    generated/fortnum_hyperg_asymptotic_outer_jvp_kernel.f90
    generated/fortnum_hyperg_asymptotic_outer_vjp_kernel.f90
    generated/fortnum_hypergeom_2f1_term_kernel.f90
    generated/fortnum_inv2_jvp_kernel.f90
    generated/fortnum_inv2_vjp_kernel.f90
    generated/fortnum_inv3_jvp_kernel.f90
    generated/fortnum_inv3_vjp_kernel.f90
    generated/fortnum_implicit_root_residual_kernel.f90
    generated/fortnum_legendre_recurrence_kernel.f90
    generated/fortnum_jacobi_recurrence_kernel.f90
    generated/fortnum_scaled_jacobi_recurrence_kernel.f90
    generated/fortnum_scaled_jacobi_gradient_recurrence_kernel.f90
    generated/fortnum_scaled_jacobi_hessian_recurrence_kernel.f90
    generated/fortnum_scalar_root_residual_jvp_kernel.f90
    generated/fortnum_scalar_root_residual_vjp_kernel.f90
    generated/fortnum_sqrt1pm1_raw_kernel.f90
    generated/fortnum_sqrt1pm1_stable_kernel.f90
    generated/fortnum_toroidal_order_kernel.f90
    generated/fortnum_vector_root_residual_jacobian_kernel.f90
    generated/fortnum_vector_root_residual_jvp_kernel.f90
    generated/fortnum_vector_root_residual_vjp_kernel.f90
    interp/fortnum_bspline.f90
    interp/fortnum_bspline_lsq.f90
    interp/fortnum_interp.f90
    interp/fortnum_polynomial.f90
    linalg/fortnum_cholesky.f90
    linalg/fortnum_krylov.f90
    linalg/fortnum_linalg.f90
    linalg/fortnum_tensor_product.f90
    linalg/fortnum_toeplitz.f90
    ode/fortnum_ode.f90
    ode/fortnum_ode_cash_karp.f90
    ode/fortnum_ode_ddeabm.f90
    ode/fortnum_ode_dop853.f90
    ode/fortnum_ode_tdrk.f90
    ode/fortnum_ode_extrapolation.f90
    ode/fortnum_ode_gauss_radau.f90
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
    roots/fortnum_fixed_point.f90
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
    special/fortnum_special_legendre.f90
    special/fortnum_special_jacobi.f90
    special/fortnum_special_spherical.f90
    special/fortnum_special_toroidal.f90
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
target_link_libraries(fortnum PUBLIC LAPACK::LAPACK fortnum_gpu_backend)

# Install the hand-written C ABI header so C/C++ consumers find it.
install(FILES ${PROJECT_SOURCE_DIR}/include/fortnum.h DESTINATION include)
