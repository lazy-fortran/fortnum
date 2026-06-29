program test_fortnum_integrate_gk
    ! Behavioural tests for integrate_gk.
    ! Each test integrates a function with a known analytic value and asserts
    ! the error lies within the requested tolerance.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_integrate_gk, only: integrate_gk
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer :: nfail

    nfail = 0
    call check_exp(nfail)
    call check_poly(nfail)
    call check_sin(nfail)
    call check_all_keys(nfail)
    call check_bad_input(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "integrate_gk: all tests passed"
    stop 0

contains

    ! int exp(x) dx from 0 to 1 = e - 1
    subroutine check_exp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: res, err
        integer  :: ierr
        call integrate_gk(f_exp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
            res, err, ierr)
        if (ierr /= 0 .or. abs(res - (exp(1.0_dp) - 1.0_dp)) > 1.0e-10_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_exp: res=", res, " ierr=", ierr
            nfail = nfail + 1
        end if
    end subroutine check_exp

    ! int (3x^2+2x+1) dx from 0 to 2 = [x^3+x^2+x]_0^2 = 14
    subroutine check_poly(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: res, err
        integer  :: ierr
        call integrate_gk(f_poly, 0.0_dp, 2.0_dp, 0.0_dp, 1.0e-12_dp, &
            res, err, ierr)
        if (ierr /= 0 .or. abs(res - 14.0_dp) > 1.0e-10_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_poly: res=", res, " ierr=", ierr
            nfail = nfail + 1
        end if
    end subroutine check_poly

    ! int sin(x) dx from 0 to pi = 2
    subroutine check_sin(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: res, err
        integer  :: ierr
        call integrate_gk(f_sin, 0.0_dp, PI, 0.0_dp, 1.0e-12_dp, &
            res, err, ierr)
        if (ierr /= 0 .or. abs(res - 2.0_dp) > 1.0e-10_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_sin: res=", res, " ierr=", ierr
            nfail = nfail + 1
        end if
    end subroutine check_sin

    ! All four GK keys agree on exp to tight tolerance.
    subroutine check_all_keys(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: keys(4) = [15, 21, 31, 61]
        real(dp), parameter :: exact = exp(1.0_dp) - 1.0_dp
        real(dp) :: res, err
        integer  :: ierr, k
        do k = 1, size(keys)
            call integrate_gk(f_exp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
                res, err, ierr, key=keys(k))
            if (ierr /= 0 .or. abs(res - exact) > 1.0e-10_dp) then
                write (error_unit, "(a,i0,a,es14.6,a,i0)") &
                    "FAIL check_all_keys key=", keys(k), &
                    " res=", res, " ierr=", ierr
                nfail = nfail + 1
            end if
        end do
    end subroutine check_all_keys

    ! Invalid key must return ierr == 6 without crashing.
    subroutine check_bad_input(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: res, err
        integer  :: ierr
        call integrate_gk(f_exp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, &
            res, err, ierr, key=99)
        if (ierr /= 6) then
            write (error_unit, "(a,i0)") &
                "FAIL check_bad_input: expected ierr=6, got ", ierr
            nfail = nfail + 1
        end if
    end subroutine check_bad_input

    function f_exp(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        fx = exp(x)
    end function f_exp

    function f_poly(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        fx = 3.0_dp*x**2 + 2.0_dp*x + 1.0_dp
    end function f_poly

    function f_sin(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        fx = sin(x)
    end function f_sin

end program test_fortnum_integrate_gk
