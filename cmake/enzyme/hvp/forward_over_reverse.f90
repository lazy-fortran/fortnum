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

    integer(c_int), parameter :: n = 4
    real(c_double) :: x(n), direction(n), target(n), target_direction(n)
    real(c_double) :: gradient(n), product(n), expected(n)

    x = [0.3_c_double, -0.7_c_double, 1.1_c_double, 0.5_c_double]
    direction = [0.2_c_double, 0.4_c_double, -0.3_c_double, 0.8_c_double]
    target = [0.1_c_double, -0.2_c_double, 0.7_c_double, 0.4_c_double]
    target_direction = 0.0_c_double
    call enzyme_fwddiff(c_funloc(reverse_gradient), x, direction, gradient, &
        product, target, target_direction, n)
    expected = (6.0_c_double*x*x - 2.0_c_double*target)*direction
    if (maxval(abs(product - expected)) > 2.0e-12_c_double) then
        error stop "forward-over-reverse HVP mismatch"
    end if
    write (*, "(a)") "PASS forward-over-reverse HVP"
end program test_forward_over_reverse_hvp
