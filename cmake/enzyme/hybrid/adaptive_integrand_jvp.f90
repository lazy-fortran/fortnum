module adaptive_integrand_autodiff
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_adaptive_integrand, only: &
        fortnum_enzyme_adaptive_integrand_jvp
    use fortnum_generated_enzyme_singular_integrand, only: &
        fortnum_enzyme_singular_integrand_jvp
    implicit none
    private

    public :: integrand_jvp, singular_integrand_jvp

contains

    pure function integrand(x, p) result(value) &
            bind(c, name="fortnum_adaptive_trace_integrand")
        real(c_double), value :: x, p
        real(c_double) :: value

        value = exp(p*x)
    end function integrand

    pure function singular_integrand(x, p) result(value) &
            bind(c, name="fortnum_singular_trace_integrand")
        real(c_double), value :: x, p
        real(c_double) :: value

        value = exp(p)/sqrt(x)
    end function singular_integrand

    function integrand_jvp(x, p) result(value)
        real(c_double), intent(in) :: x, p
        real(c_double) :: value

        value = fortnum_enzyme_adaptive_integrand_jvp( &
            x, 0.0_c_double, &
            p, 1.0_c_double)
    end function integrand_jvp

    function singular_integrand_jvp(x, p) result(value)
        real(c_double), intent(in) :: x, p
        real(c_double) :: value

        value = fortnum_enzyme_singular_integrand_jvp( &
            x, 0.0_c_double, &
            p, 1.0_c_double)
    end function singular_integrand_jvp

end module adaptive_integrand_autodiff
