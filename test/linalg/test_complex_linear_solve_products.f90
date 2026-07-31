program test_complex_linear_solve_products
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_linalg, only: dense_solve, linear_solve_complex_jvp, &
        linear_solve_complex_vjp, LINALG_OK
    implicit none

    integer, parameter :: n = 3
    real(dp), parameter :: step = 2.0e-7_dp
    complex(dp) :: a(n, n), a_bar(n, n), a_dot(n, n)
    complex(dp) :: b(n), b_bar(n), b_dot(n), x(n), x_bar(n), x_dot(n)
    complex(dp) :: x_minus(n), x_plus(n)
    real(dp) :: lhs, rhs
    integer :: info, info_minus, info_plus

    a = reshape([ &
        cmplx(3.0_dp, 0.2_dp, dp), cmplx(-0.4_dp, 0.1_dp, dp), &
        cmplx(0.2_dp, -0.3_dp, dp), cmplx(0.5_dp, -0.2_dp, dp), &
        cmplx(2.7_dp, 0.4_dp, dp), cmplx(-0.3_dp, 0.25_dp, dp), &
        cmplx(-0.15_dp, 0.35_dp, dp), cmplx(0.45_dp, -0.1_dp, dp), &
        cmplx(2.4_dp, -0.3_dp, dp)], [n, n])
    a_dot = reshape([ &
        cmplx(0.02_dp, -0.01_dp, dp), cmplx(-0.03_dp, 0.015_dp, dp), &
        cmplx(0.01_dp, 0.025_dp, dp), cmplx(-0.015_dp, 0.02_dp, dp), &
        cmplx(0.025_dp, -0.005_dp, dp), cmplx(0.03_dp, 0.01_dp, dp), &
        cmplx(-0.02_dp, 0.005_dp, dp), cmplx(0.01_dp, -0.03_dp, dp), &
        cmplx(0.015_dp, 0.02_dp, dp)], [n, n])
    b = [cmplx(0.7_dp, -0.2_dp, dp), cmplx(-0.4_dp, 0.5_dp, dp), &
        cmplx(0.3_dp, 0.15_dp, dp)]
    b_dot = [cmplx(-0.02_dp, 0.03_dp, dp), &
        cmplx(0.04_dp, -0.01_dp, dp), cmplx(0.015_dp, 0.025_dp, dp)]
    x_bar = [cmplx(0.2_dp, -0.1_dp, dp), cmplx(-0.3_dp, 0.25_dp, dp), &
        cmplx(0.15_dp, 0.35_dp, dp)]

    call dense_solve(a, b, x, info)
    call linear_solve_complex_jvp(n, a, x, a_dot, b_dot, x_dot, info)
    call dense_solve(a + step*a_dot, b + step*b_dot, x_plus, info_plus)
    call dense_solve(a - step*a_dot, b - step*b_dot, x_minus, info_minus)
    call require( &
        info == LINALG_OK .and. info_plus == LINALG_OK .and. &
        info_minus == LINALG_OK, "complex solve JVP succeeds")
    call require(maxval(abs( &
        x_dot - (x_plus - x_minus)/(2.0_dp*step))) < 2.0e-10_dp, &
        "complex solve JVP matches central differences")

    call linear_solve_complex_vjp( &
        n, a, x, x_bar, a_bar, b_bar, info)
    lhs = real(sum(conjg(x_bar)*x_dot), dp)
    rhs = real(sum(conjg(a_bar)*a_dot) + sum(conjg(b_bar)*b_dot), dp)
    call require(info == LINALG_OK, "complex solve VJP succeeds")
    call require( &
        abs(lhs - rhs) < 3.0e-13_dp*max(1.0_dp, abs(lhs), abs(rhs)), &
        "complex solve products obey the real-complex adjoint identity")

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write(error_unit, "(a,a)") "FAIL: ", message
            error stop 1
        end if
    end subroutine require

end program test_complex_linear_solve_products
