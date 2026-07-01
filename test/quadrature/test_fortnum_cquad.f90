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

end program test_fortnum_cquad
