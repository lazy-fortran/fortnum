module fixed_quadrature_full_jvp_autodiff
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    public :: full_quadrature_jvp

    interface
        function quadrature_kernel(p1, p2, p3, p4) result(value) &
                bind(c, name="fortnum_fixed_quadrature_jvp_kernel")
            import :: c_double
            real(c_double), value :: p1, p2, p3, p4
            real(c_double) :: value
        end function quadrature_kernel

        function enzyme_fwddiff(f, p1, dp1, p2, dp2, p3, dp3, p4, dp4) &
                result(derivative) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: p1, dp1, p2, dp2
            real(c_double), value :: p3, dp3, p4, dp4
            real(c_double) :: derivative
        end function enzyme_fwddiff
    end interface

contains

    function full_quadrature_jvp(parameters, direction) result(derivative)
        real(c_double), intent(in) :: parameters(4), direction(4)
        real(c_double) :: derivative

        derivative = enzyme_fwddiff(c_funloc(quadrature_kernel), &
            parameters(1), direction(1), parameters(2), direction(2), &
            parameters(3), direction(3), parameters(4), direction(4))
    end function full_quadrature_jvp

end module fixed_quadrature_full_jvp_autodiff
