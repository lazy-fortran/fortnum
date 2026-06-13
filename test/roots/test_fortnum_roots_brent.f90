program test_fortnum_roots_brent
    ! Behavioral tests for root_brent.
    !
    ! Covers: IQI path on a smooth cubic; secant path (near-flat function);
    ! wide bracket forcing many bisection steps; sign-check domain error;
    ! convergence-error path (max_iter=1); exact zero at endpoint.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots, only: root_fn_t, root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none

    integer :: nfail
    nfail = 0

    call test_brent_sqrt2(nfail)
    call test_brent_cosxx(nfail)
    call test_brent_cubic(nfail)
    call test_brent_lambertw(nfail)
    call test_brent_wide_bracket(nfail)
    call test_brent_exact_endpoint(nfail)
    call test_brent_no_bracket(nfail)
    call test_brent_max_iter(nfail)
    call test_brent_ftol(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! ------------------------------------------------------------------ helpers

    subroutine check(label, got, expected, atol, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, expected, atol
        integer,      intent(inout) :: nfail
        real(dp) :: err
        err = abs(got - expected)
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es12.4,a,es24.16,a,es24.16)") &
                "FAIL [", label, "] abserr=", err, " got=", got, &
                " expected=", expected
        end if
    end subroutine check

    subroutine check_code(label, got_code, expected_code, nfail)
        character(*), intent(in)    :: label
        integer,      intent(in)    :: got_code, expected_code
        integer,      intent(inout) :: nfail
        if (got_code /= expected_code) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a,i0)") &
                "FAIL [", label, "] status=", got_code, &
                " expected=", expected_code
        end if
    end subroutine check_code

    ! ------------------------------------------------------------------ f defs

    pure function f_sqrt2(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 2.0_dp
    end function f_sqrt2

    pure function f_cosxx(x) result(y)
        ! cos(x) - x: near-flat around the Dottie number; stresses secant path.
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = cos(x) - x
    end function f_cosxx

    pure function f_cubic(x) result(y)
        ! x^3 - x - 2; root ~ 1.5213797068
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x*x - x - 2.0_dp
    end function f_cubic

    pure function f_lambertw(x) result(y)
        ! x*exp(x) - 1; root = W(1) ~ 0.5671
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * exp(x) - 1.0_dp
    end function f_lambertw

    pure function f_sqrt3(x) result(y)
        ! x^2 - 3 in [1, 3]: wide bracket tests bisection fallback.
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 3.0_dp
    end function f_sqrt3

    pure function f_x2m4(x) result(y)
        ! x^2 - 4: f(2.0) = 0 exactly in IEEE-754.
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 4.0_dp
    end function f_x2m4

    ! ------------------------------------------------------------------ tests

    subroutine test_brent_sqrt2(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.4142135623730951455_dp

        call root_brent(f_sqrt2, 1.0_dp, 2.0_dp, x, s)
        call check_code("brent_sqrt2_status", s%code, FORTNUM_OK, nfail)
        call check("brent_sqrt2", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_brent_sqrt2

    subroutine test_brent_cosxx(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 0.73908513321516067229_dp

        call root_brent(f_cosxx, 0.5_dp, 1.0_dp, x, s)
        call check_code("brent_cosxx_status", s%code, FORTNUM_OK, nfail)
        call check("brent_cosxx", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_brent_cosxx

    subroutine test_brent_cubic(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.5213797068045675775_dp

        call root_brent(f_cubic, 1.0_dp, 2.0_dp, x, s)
        call check_code("brent_cubic_status", s%code, FORTNUM_OK, nfail)
        call check("brent_cubic", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_brent_cubic

    subroutine test_brent_lambertw(nfail)
        ! W(1): root of x*exp(x) - 1 in [0, 1]; IQI path.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 0.56714329040978387300_dp

        call root_brent(f_lambertw, 0.0_dp, 1.0_dp, x, s)
        call check_code("brent_lambertw_status", s%code, FORTNUM_OK, nfail)
        call check("brent_lambertw", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_brent_lambertw

    subroutine test_brent_wide_bracket(nfail)
        ! sqrt(3) in [1, 3]: wide bracket forces early bisection steps.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.7320508075688772936_dp

        call root_brent(f_sqrt3, 1.0_dp, 3.0_dp, x, s)
        call check_code("brent_sqrt3_status", s%code, FORTNUM_OK, nfail)
        call check("brent_sqrt3", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_brent_wide_bracket

    subroutine test_brent_exact_endpoint(nfail)
        ! f(a) = 0 exactly (IEEE-754): must return immediately with x = a.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_brent(f_x2m4, 2.0_dp, 3.0_dp, x, s)
        call check_code("brent_endpoint_status", s%code, FORTNUM_OK, nfail)
        call check("brent_endpoint", x, 2.0_dp, 0.0_dp, nfail)
    end subroutine test_brent_exact_endpoint

    subroutine test_brent_no_bracket(nfail)
        ! [1, 1.3]: f(1)=-1, f(1.3)<0 -- same sign, domain error.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_brent(f_sqrt2, 1.0_dp, 1.3_dp, x, s)
        call check_code("brent_no_bracket", s%code, FORTNUM_DOMAIN_ERROR, nfail)
    end subroutine test_brent_no_bracket

    subroutine test_brent_max_iter(nfail)
        ! max_iter=1, xtol=0, ftol=0: cannot converge; convergence error.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_brent(f_sqrt2, 1.0_dp, 2.0_dp, x, s, &
            xtol=0.0_dp, ftol=0.0_dp, max_iter=1)
        call check_code("brent_max_iter", s%code, FORTNUM_CONVERGENCE_ERROR, nfail)
    end subroutine test_brent_max_iter

    subroutine test_brent_ftol(nfail)
        ! ftol large enough that the first midpoint satisfies |f| <= ftol.
        ! The root of x^2-2 is ~1.414; f(1.5)=0.25, f(1.6)=0.56.
        ! ftol=0.3 means f=0.25 at midpoint 1.5 triggers early exit.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_brent(f_sqrt2, 1.0_dp, 2.0_dp, x, s, ftol=0.3_dp)
        call check_code("brent_ftol_status", s%code, FORTNUM_OK, nfail)
        ! Just verify the result is inside [1, 2].
        if (x < 1.0_dp .or. x > 2.0_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,es24.16)") "FAIL [brent_ftol] x out of bracket: ", x
        end if
    end subroutine test_brent_ftol

end program test_fortnum_roots_brent
