program test_fortnum_cquad
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_cquad, only: integrate_cquad
    use fortnum_integrate, only: integrate_integrand_t
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    real(dp), parameter :: PI = acos(-1.0_dp)
    integer :: nfail

    nfail = 0

    call check("x^4 on [0,1]", f_x4, 0.0_dp, 1.0_dp, 0.2_dp, 1.0e-13_dp)
    call check("exp on [-1,1]", f_exp, -1.0_dp, 1.0_dp, exp(1.0_dp) - exp(-1.0_dp), &
        1.0e-13_dp)
    call check("1/(1+x) on [0,1]", f_recip, 0.0_dp, 1.0_dp, log(2.0_dp), 1.0e-13_dp)
    call check("exp(-x) on [0,2]", f_expm, 0.0_dp, 2.0_dp, 1.0_dp - exp(-2.0_dp), &
        1.0e-13_dp)
    call check("sin on [0,pi]", f_sin, 0.0_dp, PI, 2.0_dp, 1.0e-12_dp)
    call check("sqrt on [0,1]", f_sqrt, 0.0_dp, 1.0_dp, 2.0_dp/3.0_dp, 1.0e-10_dp)
    ! Collision-operator character: Maxwellian-weighted moment.
    ! integral_0^L x^2 exp(-x^2) dx = sqrt(pi)/4 erf(L) - (L/2) exp(-L^2).
    call check("x^2 exp(-x^2) on [0,5]", f_maxw, 0.0_dp, 5.0_dp, &
        0.25_dp*sqrt(PI)*erf(5.0_dp) - 2.5_dp*exp(-25.0_dp), 1.0e-13_dp)
    ! Estimator-degeneration regression: T40 + c*T44 with c chosen so the
    ! degree-32 and degree-16 rule sums coincide on the top interval (both
    ! aliases cancel). A rule-difference error estimate reads zero and accepts
    ! the first panel with a true error near 1e-2; the coefficient-difference
    ! estimate keeps bisecting to the true value.
    call check("T40 + c*T44 aliasing pair", f_alias, -1.0_dp, 1.0_dp, &
        alias_truth(), 1.0e-12_dp)
    ! Non-finite sample regression: x*log(x) evaluates to 0*(-inf) = NaN at
    ! the endpoint node x = 0. The NaN sample must be zeroed (as in GSL CQUAD)
    ! instead of poisoning the whole result. integral_0^1 x log x dx = -1/4.
    call check("x*log(x) on [0,1] (NaN at 0)", f_xlogx, 0.0_dp, 1.0_dp, &
        -0.25_dp, 1.0e-10_dp)

    if (nfail == 0) then
        print *, "All tests passed!"
    else
        print *, "FAIL: ", nfail, " case(s)"
        error stop 1
    end if

contains

    subroutine check(name, f, a, b, ref, tol)
        character(*), intent(in) :: name
        procedure(integrate_integrand_t) :: f
        real(dp),     intent(in) :: a, b, ref, tol
        real(dp) :: val, err
        type(fortnum_status_t) :: status
        real(dp) :: d

        call integrate_cquad(f, a, b, val, status, epsabs=1.0e-13_dp, &
            epsrel=1.0e-13_dp, abserr=err)
        d = abs(val - ref)
        if (status%code /= FORTNUM_OK .or. d > tol) then
            nfail = nfail + 1
            write (*, '(A,A,A,ES12.4,A,ES12.4,A,ES12.4)') "FAIL ", name, &
                "  err=", d, " tol=", tol, " abserr=", err
        else
            write (*, '(A,A,A,ES12.4)') "ok   ", name, "  err=", d
        end if
    end subroutine check

    function f_x4(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = x**4
    end function f_x4

    function f_exp(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = exp(x)
    end function f_exp

    function f_recip(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = 1.0_dp/(1.0_dp + x)
    end function f_recip

    function f_expm(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = exp(-x)
    end function f_expm

    function f_sin(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = sin(x)
    end function f_sin

    function f_sqrt(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = sqrt(x)
    end function f_sqrt

    function f_maxw(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = x*x*exp(-x*x)
    end function f_maxw

    ! integral of T_k over [-1,1]: 0 for odd k, 2/(1-k^2) for even k.
    pure function cheb_int(k) result(v)
        integer, intent(in) :: k
        real(dp) :: v
        if (mod(k, 2) == 1) then
            v = 0.0_dp
        else
            v = 2.0_dp/(1.0_dp - real(k, dp)**2)
        end if
    end function cheb_int

    ! On the 33-node grid T40 aliases to T24 and T44 to T20; on the 17-node
    ! subset they alias to T8 and T12. c equates the two rule sums.
    pure function alias_coeff() result(c)
        real(dp) :: c
        c = (cheb_int(24) - cheb_int(8))/(cheb_int(12) - cheb_int(20))
    end function alias_coeff

    pure function alias_truth() result(v)
        real(dp) :: v
        v = cheb_int(40) + alias_coeff()*cheb_int(44)
    end function alias_truth

    function f_alias(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y, theta
        theta = acos(max(-1.0_dp, min(1.0_dp, x)))
        y = cos(40.0_dp*theta) + alias_coeff()*cos(44.0_dp*theta)
    end function f_alias

    function f_xlogx(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = x*log(x)
    end function f_xlogx

end program test_fortnum_cquad
