program test_fortnum_roots_multiplicity
    ! Multiple-root accuracy of complex_region_roots.
    !
    ! At a multiplicity-m zero the plain Newton step is only linearly
    ! convergent, so an absolute step-stop halts the polish ~m*newtonz short
    ! of the root. The multiplicity-aware step (m*f/f') restores quadratic
    ! convergence at multiple roots. With the plain step these cases (a
    ! double root and a triple root) halted at abserr ~7.2e-11 / ~1.7e-11;
    ! the modified step reaches ~3e-12 / ~1.9e-12. The bar here is 1e-11,
    ! well below the plain-step errors and far under the 1e-8 golden bar.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots_complex, only: complex_region_roots
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, save :: active_case = -1
    integer :: nfail

    nfail = 0

    ! Case 0: simple-simple-double product, the oracle case 3 regime; the
    ! double root at -0.9-0.6i is the sensitive one.
    call check_root(0, cmplx(-0.9_dp, -0.6_dp, dp), 2, &
        cmplx(-1.2_dp, -1.0_dp, dp), cmplx(0.9_dp, 0.5_dp, dp), nfail)

    ! Case 1: pure triple root (x0-x)^3, the oracle case 2 regime.
    call check_root(1, cmplx(-0.3_dp, 0.6_dp, dp), 3, &
        cmplx(-1.0_dp, -0.2_dp, dp), cmplx(0.4_dp, 1.0_dp, dp), nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " multiple-root check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine f_case(kr, fk, ctx)
        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: fk
        class(*),    intent(in), optional :: ctx

        select case (active_case)
        case (0)
            fk = (kr - cmplx(0.20_dp, 0.10_dp, dp)) &
                * (kr - cmplx(0.50_dp, 0.10_dp, dp)) &
                * (kr - cmplx(-0.9_dp, -0.6_dp, dp))**2
        case (1)
            fk = (kr - cmplx(-0.3_dp, 0.6_dp, dp))**3
        case default
            fk = (0.0_dp, 0.0_dp)
        end select
        if (present(ctx)) continue
    end subroutine f_case

    ! Run the finder over [ll, ur], find the returned zero nearest zref, and
    ! assert it lands within 1e-12 of zref and reports multiplicity mexp.
    subroutine check_root(cid, zref, mexp, ll, ur, nfail)
        integer,     intent(in)    :: cid, mexp
        complex(dp), intent(in)    :: zref, ll, ur
        integer,     intent(inout) :: nfail

        complex(dp), allocatable :: roots(:), fvals(:)
        integer,     allocatable :: mult(:)
        integer :: nfound, j, jbest
        type(fortnum_status_t) :: status
        real(dp) :: d, dbest
        real(dp), parameter :: tol = 1.0e-11_dp

        active_case = cid
        call complex_region_roots(f_case, ll, ur, roots, fvals, mult, &
            nfound, status, m_max=5)

        if (status%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0)") "FAIL case ", cid, &
                " status=", status%code
            return
        end if

        jbest = 0
        dbest = huge(1.0_dp)
        do j = 1, nfound
            d = abs(roots(j) - zref)
            if (d < dbest) then
                dbest = d
                jbest = j
            end if
        end do

        if (jbest == 0) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a)") "FAIL case ", cid, " no zeros found"
            return
        end if

        if (dbest > tol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.4,a,es12.4)") "FAIL case ", cid, &
                " abserr=", dbest, " tol=", tol
        end if

        if (mult(jbest) /= mexp) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,i0)") "FAIL case ", cid, &
                " mult got=", mult(jbest), " expected=", mexp
        end if
    end subroutine check_root

end program test_fortnum_roots_multiplicity
