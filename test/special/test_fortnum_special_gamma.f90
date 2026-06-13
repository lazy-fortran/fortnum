program test_fortnum_special_gamma
    ! Behavioral tests for fortnum_special_gamma.
    ! Checks boundary values, the series/contfrac branch transition, and
    ! the identity gamma_lower(a,x) = gamma_reg_p(a,x) * Gamma(a).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_gamma, only: gamma_lower, gamma_reg_p
    implicit none

    integer  :: nfail
    real(dp) :: tol

    nfail = 0
    tol   = 1.0e-12_dp

    ! x == 0 -> result is 0 regardless of a.
    call check("gamma_lower(1.0, 0.0) == 0", &
               gamma_lower(1.0_dp, 0.0_dp), 0.0_dp, tol)
    call check("gamma_reg_p(2.0, 0.0) == 0", &
               gamma_reg_p(2.0_dp, 0.0_dp), 0.0_dp, tol)

    ! gamma_lower(1, x) = 1 - exp(-x)  (DLMF 8.6.1).
    call check("gamma_lower(1,0.5) = 1-exp(-0.5)", &
               gamma_lower(1.0_dp, 0.5_dp), 1.0_dp - exp(-0.5_dp), tol)
    call check("gamma_lower(1,2.0) = 1-exp(-2)", &
               gamma_lower(1.0_dp, 2.0_dp), 1.0_dp - exp(-2.0_dp), tol)

    ! gamma_lower(a, big x) -> Gamma(a) as x -> inf; use x = 50 as proxy.
    call check("gamma_lower(2, 50) ~ Gamma(2) = 1", &
               gamma_lower(2.0_dp, 50.0_dp), gamma(2.0_dp), tol)

    ! Identity: gamma_lower(a,x) == gamma_reg_p(a,x) * Gamma(a).
    ! Test across both branches (series: x < a+1; contfrac: x >= a+1).
    call check("identity series branch a=3 x=1", &
               gamma_lower(3.0_dp, 1.0_dp), &
               gamma_reg_p(3.0_dp, 1.0_dp) * gamma(3.0_dp), tol)
    call check("identity contfrac branch a=1 x=5", &
               gamma_lower(1.0_dp, 5.0_dp), &
               gamma_reg_p(1.0_dp, 5.0_dp) * gamma(1.0_dp), tol)

    ! gamma_reg_p(a, inf proxy) -> 1.
    call check("gamma_reg_p(1, 50) ~ 1", &
               gamma_reg_p(1.0_dp, 50.0_dp), 1.0_dp, tol)

    ! Non-integer a: compare against Fortran intrinsic gamma (sanity cross-check).
    ! gamma_lower(0.5, inf) = Gamma(0.5) = sqrt(pi).
    call check("gamma_lower(0.5, 40) ~ sqrt(pi)", &
               gamma_lower(0.5_dp, 40.0_dp), sqrt(acos(-1.0_dp)), tol)

    if (nfail > 0) then
        write(error_unit, "(i0,a)") nfail, " test(s) FAILED"
        stop 1
    end if
    write(*, "(a)") "PASS"

contains

    subroutine check(label, got, expected, t)
        character(*), intent(in) :: label
        real(dp),     intent(in) :: got, expected, t
        real(dp) :: err
        err = abs(got - expected)
        if (.not. (err <= t + t * abs(expected))) then
            nfail = nfail + 1
            write(error_unit, "(a,a,a,es12.5,a,es12.5,a,es12.5)") &
                "FAIL: ", label, "  got=", got, "  exp=", expected, "  err=", err
        end if
    end subroutine check

end program test_fortnum_special_gamma
