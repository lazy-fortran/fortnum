module adaptive_integrand_autodiff
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    public :: integrand_jvp

    interface
        function enzyme_fwddiff(f, x, dx, p, dp_seed) result(derivative) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x, dx, p, dp_seed
            real(c_double) :: derivative
        end function enzyme_fwddiff
    end interface

contains

    pure function integrand(x, p) result(value) &
            bind(c, name="fortnum_adaptive_trace_integrand")
        real(c_double), value :: x, p
        real(c_double) :: value

        value = exp(p*x)
    end function integrand

    function integrand_jvp(x, p) result(value)
        real(c_double), intent(in) :: x, p
        real(c_double) :: value

        value = enzyme_fwddiff(c_funloc(integrand), x, 0.0_c_double, &
            p, 1.0_c_double)
    end function integrand_jvp

end module adaptive_integrand_autodiff
