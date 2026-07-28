module reverse_over_forward_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    implicit none
    private
    public :: directional_derivative

    interface
        function enzyme_fwddiff(f, x, direction, target, target_direction, n) &
                result(value) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x(*), direction(*), target(*)
            real(c_double), intent(in) :: target_direction(*)
            integer(c_int), value :: n
            real(c_double) :: value
        end function enzyme_fwddiff
    end interface

contains

    pure function least_squares_objective(x, target, n) result(value) &
            bind(c, name="fortnum_rof_least_squares")
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

    function directional_derivative(x, direction, target, target_direction, n) &
            result(value) bind(c, name="fortnum_hvp_directional_derivative")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), direction(n), target(n)
        real(c_double), intent(in) :: target_direction(n)
        real(c_double) :: value

        value = enzyme_fwddiff(c_funloc(least_squares_objective), x, direction, &
            target, target_direction, n)
    end function directional_derivative

end module reverse_over_forward_kernel

program test_reverse_over_forward_hvp
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    use reverse_over_forward_kernel, only: directional_derivative
    use fortnum_hvp_fixture_support, only: hvp_size, run_hvp_fixture
    implicit none

    interface
        function enzyme_autodiff(f, x, x_bar, direction, direction_bar, target, &
                target_bar, target_direction, target_direction_bar, n) &
                result(value) bind(c, name="__enzyme_autodiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x(*), direction(*), target(*)
            real(c_double), intent(in) :: target_direction(*)
            real(c_double), intent(inout) :: x_bar(*), direction_bar(*)
            real(c_double), intent(inout) :: target_bar(*)
            real(c_double), intent(inout) :: target_direction_bar(*)
            integer(c_int), value :: n
            real(c_double) :: value
        end function enzyme_autodiff
    end interface

    call run_hvp_fixture("reverse_over_forward", evaluate_hvp)

contains

    subroutine evaluate_hvp(x, direction, target, product)
        real(c_double), intent(in) :: x(hvp_size), direction(hvp_size)
        real(c_double), intent(in) :: target(hvp_size)
        real(c_double), intent(out) :: product(hvp_size)
        real(c_double) :: direction_bar(hvp_size), target_bar(hvp_size)
        real(c_double) :: target_direction(hvp_size)
        real(c_double) :: target_direction_bar(hvp_size), value

        product = 0.0_c_double
        direction_bar = 0.0_c_double
        target_bar = 0.0_c_double
        target_direction = 0.0_c_double
        target_direction_bar = 0.0_c_double
        value = enzyme_autodiff(c_funloc(directional_derivative), x, product, &
            direction, direction_bar, target, target_bar, target_direction, &
            target_direction_bar, int(hvp_size, c_int))
    end subroutine evaluate_hvp

end program test_reverse_over_forward_hvp
