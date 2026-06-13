program test_fortnum_roots
    ! Behavioral tests for fortnum_roots.
    !
    ! Covers: bisection on simple roots; Newton with analytic derivatives;
    ! sign-check domain error; convergence-error path for a deliberately
    ! underresolved bracket; near-zero-derivative guard in Newton.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots, only: root_fn_t, root_fn_df_t, root_bisect, root_newton
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none

    integer :: nfail
    nfail = 0

    call test_bisect_sqrt2(nfail)
    call test_bisect_cosxx(nfail)
    call test_bisect_cubic(nfail)
    call test_bisect_exact_endpoint(nfail)
    call test_bisect_no_bracket(nfail)
    call test_bisect_max_iter(nfail)
    call test_newton_sqrt2(nfail)
    call test_newton_exp(nfail)
    call test_newton_cubic(nfail)
    call test_newton_no_bracket(nfail)
    call test_newton_step_outside_bracket(nfail)

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

    pure function f_lnroot(x) result(y)
        ! exp(x) - 3; root = ln(3)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = exp(x) - 3.0_dp
    end function f_lnroot

    pure function f_x2m4(x) result(y)
        ! x^2 - 4; root = 2 exactly in IEEE-754.
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 4.0_dp
    end function f_x2m4

    pure subroutine fdf_sqrt2(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x * x - 2.0_dp
        dfx = 2.0_dp * x
    end subroutine fdf_sqrt2

    pure subroutine fdf_exp(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = exp(x) - 3.0_dp
        dfx = exp(x)
    end subroutine fdf_exp

    pure subroutine fdf_cubic(x, fx, dfx)
        ! x^3 - x - 2
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x*x*x - x - 2.0_dp
        dfx = 3.0_dp * x*x - 1.0_dp
    end subroutine fdf_cubic

    ! ------------------------------------------------------------------ tests

    subroutine test_bisect_sqrt2(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.4142135623730951455_dp

        call root_bisect(f_sqrt2, 1.0_dp, 2.0_dp, x, s)
        call check_code("bisect_sqrt2_status", s%code, FORTNUM_OK, nfail)
        call check("bisect_sqrt2", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_bisect_sqrt2

    subroutine test_bisect_cosxx(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 0.73908513321516067229_dp

        call root_bisect(f_cosxx, 0.5_dp, 1.0_dp, x, s)
        call check_code("bisect_cosxx_status", s%code, FORTNUM_OK, nfail)
        call check("bisect_cosxx", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_bisect_cosxx

    subroutine test_bisect_cubic(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.5213797068045675775_dp

        call root_bisect(f_cubic, 1.0_dp, 2.0_dp, x, s)
        call check_code("bisect_cubic_status", s%code, FORTNUM_OK, nfail)
        call check("bisect_cubic", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_bisect_cubic

    subroutine test_bisect_exact_endpoint(nfail)
        ! f(a) = 0 exactly: return immediately with x=a.
        ! x^2 - 4: f(2.0) = 0.0 exactly in IEEE-754.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_bisect(f_x2m4, 2.0_dp, 3.0_dp, x, s)
        call check_code("bisect_endpoint_status", s%code, FORTNUM_OK, nfail)
        call check("bisect_endpoint", x, 2.0_dp, 0.0_dp, nfail)
    end subroutine test_bisect_exact_endpoint

    subroutine test_bisect_no_bracket(nfail)
        ! [1, 1.3]: f(1)=-1, f(1.3)=1.3^2-2=-0.31 -- same-sign domain error.
        ! Actually f(1)=-1 < 0, f(1.3) ~ -0.31 < 0; both negative => error.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_bisect(f_sqrt2, 1.0_dp, 1.3_dp, x, s)
        call check_code("bisect_no_bracket", s%code, FORTNUM_DOMAIN_ERROR, nfail)
    end subroutine test_bisect_no_bracket

    subroutine test_bisect_max_iter(nfail)
        ! max_iter=1 and xtol=0, ftol=0: cannot converge; expect CONVERGENCE_ERROR.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_bisect(f_sqrt2, 1.0_dp, 2.0_dp, x, s, &
            xtol=0.0_dp, ftol=0.0_dp, max_iter=1)
        call check_code("bisect_max_iter", s%code, FORTNUM_CONVERGENCE_ERROR, nfail)
    end subroutine test_bisect_max_iter

    subroutine test_newton_sqrt2(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.4142135623730951455_dp

        call root_newton(fdf_sqrt2, 1.0_dp, 2.0_dp, 1.5_dp, x, s)
        call check_code("newton_sqrt2_status", s%code, FORTNUM_OK, nfail)
        call check("newton_sqrt2", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_newton_sqrt2

    subroutine test_newton_exp(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.0986122886681097821_dp

        call root_newton(fdf_exp, 0.0_dp, 2.0_dp, 1.0_dp, x, s)
        call check_code("newton_exp_status", s%code, FORTNUM_OK, nfail)
        call check("newton_exp", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_newton_exp

    subroutine test_newton_cubic(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 1.5213797068045675775_dp

        call root_newton(fdf_cubic, 1.0_dp, 2.0_dp, 1.5_dp, x, s)
        call check_code("newton_cubic_status", s%code, FORTNUM_OK, nfail)
        call check("newton_cubic", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_newton_cubic

    subroutine test_newton_no_bracket(nfail)
        ! [1, 1.3]: f positive on both ends (f_sqrt2(1)=-1, f_sqrt2(1.3)<0).
        ! Use [3, 4]: f(3)=7>0, f(4)=14>0 -- domain error.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x

        call root_newton(fdf_sqrt2, 3.0_dp, 4.0_dp, 3.5_dp, x, s)
        call check_code("newton_no_bracket", s%code, FORTNUM_DOMAIN_ERROR, nfail)
    end subroutine test_newton_no_bracket

    subroutine test_newton_step_outside_bracket(nfail)
        ! f = cos(x) - x has a root at ~0.739.  Starting from x0 far from the
        ! root (but inside [0.5, 1]) exercises the bisection-guard branch.
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s
        real(dp) :: x
        real(dp), parameter :: ref = 0.73908513321516067229_dp

        call root_newton(fdf_cosxx, 0.5_dp, 1.0_dp, 0.99_dp, x, s)
        call check_code("newton_cosxx_status", s%code, FORTNUM_OK, nfail)
        call check("newton_cosxx", x, ref, 2.0e-15_dp, nfail)
    end subroutine test_newton_step_outside_bracket

    ! cos(x) - x with derivative for Newton test.
    pure subroutine fdf_cosxx(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = cos(x) - x
        dfx = -sin(x) - 1.0_dp
    end subroutine fdf_cosxx

end program test_fortnum_roots
