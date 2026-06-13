program test_complex_step
    ! Proves the complex-step helper on analytic scalar kernels with known
    ! derivatives: sin (d=cos), exp (d=exp), and a polynomial. Complex step has
    ! no subtractive cancellation, so agreement is to a few ulp even with a
    ! tiny step.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_complex_step, complex_step_deriv
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846_dp
    integer :: nfail
    nfail = 0

    call test_sin(nfail)
    call test_exp(nfail)
    call test_poly(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    subroutine test_sin(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x
        x = 0.7_dp
        if (.not. check_complex_step("sin", f_sin, x, cos(x), &
                tol=1.0e-14_dp)) nfail = nfail + 1
        x = PI/3.0_dp
        if (.not. check_complex_step("sin_pi3", f_sin, x, cos(x), &
                tol=1.0e-14_dp)) nfail = nfail + 1
    end subroutine test_sin

    subroutine test_exp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x
        x = 1.3_dp
        if (.not. check_complex_step("exp", f_exp, x, exp(x), &
                tol=1.0e-14_dp)) nfail = nfail + 1
    end subroutine test_exp

    ! p(x) = 3 x^3 - 2 x + 5, p'(x) = 9 x^2 - 2.
    subroutine test_poly(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x, dwant
        x = -1.25_dp
        dwant = 9.0_dp*x*x - 2.0_dp
        if (.not. check_complex_step("poly", f_poly, x, dwant, &
                tol=1.0e-14_dp)) nfail = nfail + 1
    end subroutine test_poly

    ! ----------------------------------------------------- reference kernels

    pure function f_sin(z) result(w)
        complex(dp), intent(in) :: z
        complex(dp)             :: w
        w = sin(z)
    end function f_sin

    pure function f_exp(z) result(w)
        complex(dp), intent(in) :: z
        complex(dp)             :: w
        w = exp(z)
    end function f_exp

    pure function f_poly(z) result(w)
        complex(dp), intent(in) :: z
        complex(dp)             :: w
        w = 3.0_dp*z**3 - 2.0_dp*z + 5.0_dp
    end function f_poly

end program test_complex_step
