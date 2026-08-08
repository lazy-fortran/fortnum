program test_complex_linear_solve_products
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_linalg, only: dense_solve, linear_solve_complex_jvp, &
        linear_solve_complex_vjp, LINALG_OK
    implicit none

    integer, parameter :: n = 3
    !! Base step for the finite-difference oracle, with the error budget
    !! written down rather than assumed.
    !!
    !! A central difference carries truncation error growing as `step**2` and
    !! cancellation error growing as `eps/step`. The previous value, 2e-7, put
    !! the cancellation term at about `2.2e-16 * 0.3 / 2e-7 = 3e-10` -- larger
    !! than the 2e-10 tolerance it was checked against. The test was therefore
    !! inside its own noise floor and passed on a particular rounding pattern:
    !! relinking BLAS changed the rounding and it failed, which is how this was
    !! found.
    !!
    !! A larger step with Richardson extrapolation fixes both terms at once.
    !! At 1e-4 the cancellation term falls to about 1e-12, and extrapolating
    !! two central differences cancels the leading truncation term so what
    !! remains scales as `step**4`, around 1e-16. The tolerance below now sits
    !! orders of magnitude above the error rather than beneath it.
    real(dp), parameter :: step = 1.0e-4_dp
    complex(dp) :: a(n, n), a_bar(n, n), a_dot(n, n)
    complex(dp) :: b(n), b_bar(n), b_dot(n), x(n), x_bar(n), x_dot(n)
    complex(dp) :: x_minus(n), x_plus(n)
    complex(dp) :: x_fine_minus(n), x_fine_plus(n)
    complex(dp) :: coarse(n), fine(n), extrapolated(n)
    integer :: info_fine_minus, info_fine_plus
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
    call dense_solve(a + 0.5_dp*step*a_dot, b + 0.5_dp*step*b_dot, &
        x_fine_plus, info_fine_plus)
    call dense_solve(a - 0.5_dp*step*a_dot, b - 0.5_dp*step*b_dot, &
        x_fine_minus, info_fine_minus)
    call require( &
        info == LINALG_OK .and. info_plus == LINALG_OK .and. &
        info_minus == LINALG_OK .and. info_fine_plus == LINALG_OK .and. &
        info_fine_minus == LINALG_OK, "complex solve JVP succeeds")

    ! Central differences are second order, so the extrapolation weights are
    ! (4*fine - coarse)/3 and the leading truncation term cancels.
    coarse = (x_plus - x_minus)/(2.0_dp*step)
    fine = (x_fine_plus - x_fine_minus)/step
    extrapolated = (4.0_dp*fine - coarse)/3.0_dp
    call require(maxval(abs(x_dot - extrapolated)) < 1.0e-9_dp, &
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
