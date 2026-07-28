module forward_over_reverse_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    implicit none
    private
    public :: least_squares_objective, reverse_gradient

    interface
        function enzyme_autodiff(f, x, x_bar, target, target_bar, n) result(value) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x(*), target(*), target_bar(*)
            real(c_double), intent(inout) :: x_bar(*)
            integer(c_int), value :: n
            real(c_double) :: value
        end function enzyme_autodiff
    end interface

contains

    pure function least_squares_objective(x, target, n) result(value) &
            bind(c, name="fortnum_hvp_least_squares")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), target(n)
        real(c_double) :: value, residual
        integer :: i

        value = 0.0_c_double
        do i = 1, n
            residual = x(i)*x(i) - target(i)
            value = value + 0.5_c_double*residual*residual
        end do
    end function least_squares_objective

    subroutine reverse_gradient(x, gradient, target, n) &
            bind(c, name="fortnum_hvp_reverse_gradient")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), target(n)
        real(c_double), intent(out) :: gradient(n)
        real(c_double) :: target_bar(n), value

        gradient = 0.0_c_double
        target_bar = 0.0_c_double
        value = enzyme_autodiff(c_funloc(least_squares_objective), x, gradient, &
            target, target_bar, n)
    end subroutine reverse_gradient

end module forward_over_reverse_kernel

program test_forward_over_reverse_hvp
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    use forward_over_reverse_kernel, only: reverse_gradient
    use fortnum_hvp_fixture_support, only: hvp_size, run_hvp_fixture
    implicit none

    interface
        subroutine enzyme_fwddiff(f, x, direction, gradient, product, target, &
                target_direction, n) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x(*), direction(*), target(*)
            real(c_double), intent(in) :: target_direction(*)
            real(c_double), intent(out) :: gradient(*), product(*)
            integer(c_int), value :: n
        end subroutine enzyme_fwddiff
    end interface

    call run_hvp_fixture("forward_over_reverse", evaluate_hvp)

contains

    subroutine evaluate_hvp(x, direction, target, product)
        real(c_double), intent(in) :: x(hvp_size), direction(hvp_size)
        real(c_double), intent(in) :: target(hvp_size)
        real(c_double), intent(out) :: product(hvp_size)
        real(c_double) :: target_direction(hvp_size), gradient(hvp_size)

        target_direction = 0.0_c_double
        call enzyme_fwddiff(c_funloc(reverse_gradient), x, direction, gradient, &
            product, target, target_direction, int(hvp_size, c_int))
    end subroutine evaluate_hvp

end program test_forward_over_reverse_hvp
