module fixed_quadrature_full_vjp_autodiff
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    type, bind(c) :: quadrature_gradient_t
        real(c_double) :: values(4)
    end type quadrature_gradient_t

    public :: full_quadrature_vjp

    interface
        function quadrature_kernel(p1, p2, p3, p4) result(value) &
                bind(c, name="fortnum_fixed_quadrature_vjp_kernel")
            import :: c_double
            real(c_double), value :: p1, p2, p3, p4
            real(c_double) :: value
        end function quadrature_kernel

        function enzyme_autodiff(f, p1, p2, p3, p4) result(gradient) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr, quadrature_gradient_t
            type(c_funptr), value :: f
            real(c_double), value :: p1, p2, p3, p4
            type(quadrature_gradient_t) :: gradient
        end function enzyme_autodiff
    end interface

contains

    function full_quadrature_vjp(parameters, cotangent) result(vjp)
        real(c_double), intent(in) :: parameters(4), cotangent
        real(c_double) :: vjp(4)
        type(quadrature_gradient_t) :: gradient

        gradient = enzyme_autodiff(c_funloc(quadrature_kernel), &
            parameters(1), parameters(2), parameters(3), parameters(4))
        vjp = cotangent*gradient%values
    end function full_quadrature_vjp

end module fixed_quadrature_full_vjp_autodiff
