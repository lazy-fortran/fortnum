program test_fortnum_special_second_derivatives
    ! Independent ODE and centered-difference checks for second derivatives.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_legendre, only: &
        legendre_p, legendre_p_derivative, legendre_p_second_derivative, &
        legendre_q, legendre_q_derivative, legendre_q_second_derivative
    use fortnum_special_toroidal, only: &
        toroidal_p, toroidal_p_derivative, toroidal_p_second_derivative, &
        toroidal_q, toroidal_q_derivative, toroidal_q_second_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: x_ferrers = 0.37_dp
    real(dp), parameter :: x_toroidal = 2.2_dp
    real(dp), parameter :: finite_difference_step = 2.0e-4_dp
    real(dp), parameter :: ode_tolerance = 2.0e-11_dp
    real(dp), parameter :: finite_difference_tolerance = 3.0e-7_dp

    nfail = 0
    call check_close("P_2 second derivative", &
                     legendre_p_second_derivative(2, 0, x_ferrers), 3.0_dp, 2.0e-13_dp)
    call check_close("P_3^2 second derivative", &
                    legendre_p_second_derivative(3, 2, x_ferrers), -90.0_dp*x_ferrers, &
                     2.0e-12_dp)
    call check_close("Q_0 second derivative", &
                     legendre_q_second_derivative(0, 2.0_dp), 4.0_dp/9.0_dp, 2.0e-13_dp)
    call check_ferrers_ode(4, 1, x_ferrers)
    call check_legendre_q_ode(4, 2.1_dp)
    call check_ferrers_fd(4, 1, x_ferrers)
    call check_legendre_q_fd(4, 2.1_dp)

    call check_toroidal_ode(.true., 3, 1, x_toroidal)
    call check_toroidal_ode(.false., 3, 1, x_toroidal)
    call check_toroidal_fd(.true., 3, 1, x_toroidal)
    call check_toroidal_fd(.false., 3, 1, x_toroidal)

    call check_true("Ferrers second derivative endpoint is undefined", &
                    ieee_is_nan(legendre_p_second_derivative(2, 0, 1.0_dp)))
    call check_true("ordinary Q second derivative domain", &
                    ieee_is_nan(legendre_q_second_derivative(2, 1.0_dp)))
    call check_true("toroidal P second derivative domain", &
                    ieee_is_nan(toroidal_p_second_derivative(2, 0, 1.0_dp)))
    call check_true("toroidal Q second derivative domain", &
                    ieee_is_nan(toroidal_q_second_derivative(2, 0, 1.0_dp)))

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_close(label, got, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, tolerance
        if (abs(got - expected) > tolerance*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check_close

    subroutine check_true(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition
        if (.not. condition) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: "//label
        end if
    end subroutine check_true

    subroutine check_ferrers_ode(degree, order, x)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: x
        real(dp) :: value, first, second, residual, denominator
        value = legendre_p(degree, order, x)
        first = legendre_p_derivative(degree, order, x)
        second = legendre_p_second_derivative(degree, order, x)
        denominator = 1.0_dp - x*x
        residual = denominator*second - 2.0_dp*x*first + &
                   (real(degree*(degree + 1), dp) - &
                    real(order*order, dp)/denominator)*value
        call check_close("associated Legendre ODE residual", residual, 0.0_dp, &
                         ode_tolerance)
    end subroutine check_ferrers_ode

    subroutine check_legendre_q_ode(degree, x)
        integer, intent(in) :: degree
        real(dp), intent(in) :: x
        real(dp) :: value, first, second, residual
        value = legendre_q(degree, x)
        first = legendre_q_derivative(degree, x)
        second = legendre_q_second_derivative(degree, x)
        residual = (1.0_dp - x*x)*second - 2.0_dp*x*first + &
                   real(degree*(degree + 1), dp)*value
        call check_close("ordinary Legendre Q ODE residual", residual, 0.0_dp, &
                         ode_tolerance)
    end subroutine check_legendre_q_ode

    subroutine check_ferrers_fd(degree, order, x)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: x
        real(dp) :: finite_difference
        finite_difference = (legendre_p_derivative(degree, order, x + &
                    finite_difference_step) - legendre_p_derivative(degree, order, x - &
                                finite_difference_step))/(2.0_dp*finite_difference_step)
        call check_close("associated Legendre second derivative finite difference", &
                    legendre_p_second_derivative(degree, order, x), finite_difference, &
                         finite_difference_tolerance)
    end subroutine check_ferrers_fd

    subroutine check_legendre_q_fd(degree, x)
        integer, intent(in) :: degree
        real(dp), intent(in) :: x
        real(dp) :: finite_difference
        finite_difference = (legendre_q_derivative(degree, x + &
                           finite_difference_step) - legendre_q_derivative(degree, x - &
                                finite_difference_step))/(2.0_dp*finite_difference_step)
        call check_close("ordinary Legendre Q second derivative finite difference", &
                         legendre_q_second_derivative(degree, x), finite_difference, &
                         finite_difference_tolerance)
    end subroutine check_legendre_q_fd

    subroutine check_toroidal_ode(first_kind, degree_index, order, x)
        logical, intent(in) :: first_kind
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: x
        real(dp) :: value, first, second, residual, degree, denominator
        degree = real(degree_index, dp) - 0.5_dp
        if (first_kind) then
            value = toroidal_p(degree_index, order, x)
            first = toroidal_p_derivative(degree_index, order, x)
            second = toroidal_p_second_derivative(degree_index, order, x)
        else
            value = toroidal_q(degree_index, order, x)
            first = toroidal_q_derivative(degree_index, order, x)
            second = toroidal_q_second_derivative(degree_index, order, x)
        end if
        denominator = 1.0_dp - x*x
        residual = denominator*second - 2.0_dp*x*first + &
                   (degree*(degree + 1.0_dp) - &
                    real(order*order, dp)/denominator)*value
        call check_close("toroidal associated ODE residual", residual, 0.0_dp, &
                         ode_tolerance)
    end subroutine check_toroidal_ode

    subroutine check_toroidal_fd(first_kind, degree_index, order, x)
        logical, intent(in) :: first_kind
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: x
        real(dp) :: finite_difference, value
        if (first_kind) then
            finite_difference = (toroidal_p_derivative(degree_index, order, x + &
              finite_difference_step) - toroidal_p_derivative(degree_index, order, x - &
                                finite_difference_step))/(2.0_dp*finite_difference_step)
            value = toroidal_p_second_derivative(degree_index, order, x)
        else
            finite_difference = (toroidal_q_derivative(degree_index, order, x + &
              finite_difference_step) - toroidal_q_derivative(degree_index, order, x - &
                                finite_difference_step))/(2.0_dp*finite_difference_step)
            value = toroidal_q_second_derivative(degree_index, order, x)
        end if
        call check_close("toroidal second derivative finite difference", value, &
                         finite_difference, finite_difference_tolerance)
    end subroutine check_toroidal_fd

end program test_fortnum_special_second_derivatives
