program test_fortnum_region_roots_converge
    ! Convergence contract of complex_region_roots: a returned zero must be a
    ! real zero of f.
    !
    ! newton_polish refines each ZGGEV moment estimate. When a cell holds a
    ! pair of simple zeros too close to separate from the contour moments, the
    ! Hankel pencil collapses to a single double estimate at their midpoint.
    ! That midpoint is not a zero: |f| there is ~(sep/2)^2 * scale. The polish
    ! must not report it as a converged root. The acceptance guard therefore
    ! has to verify |f(root)| <= newtonf, not merely that the Newton step fell
    ! below the |dz| tolerance (the step-exit) or that |f| stayed under the
    ! loose 1e-8 floor (the max-iteration exit).
    !
    ! Case A places two simple zeros 3.16e-4 apart near the right portion of a
    ! cell with m_max=3, so the trio resolves in one cell without bisection and
    ! the pair collapses to a midpoint estimate. Contract: every returned root
    ! is a true zero, |f| <= 1e-10 and within 1e-10 of an exact construction
    ! zero. Case B runs the identical function with m_max=1, forcing bisection
    ! to isolate each zero, and requires all three found to the same bar; it
    ! guards that the tightened guard does not over-reject genuine zeros.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots_complex, only: complex_region_roots
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    real(dp), parameter :: sep = 3.1622776601683794e-4_dp
    real(dp), parameter :: cx  = -0.77_dp
    real(dp), parameter :: fbar = 1.0e-10_dp
    real(dp), parameter :: zbar = 1.0e-10_dp

    complex(dp), save :: zexact(3)
    integer :: nfail

    zexact(1) = cmplx(cx - 0.5_dp*sep, 0.13_dp, dp)
    zexact(2) = cmplx(cx + 0.5_dp*sep, 0.13_dp, dp)
    zexact(3) = cmplx(-0.83_dp, -0.41_dp, dp)

    nfail = 0
    call check_no_false_root(3, nfail)
    call check_all_found(1, nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " convergence contract check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine f_trio(kr, fk, ctx)
        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: fk
        class(*),    intent(in), optional :: ctx
        fk = (kr - zexact(1)) * (kr - zexact(2)) * (kr - zexact(3))
        if (present(ctx)) continue
    end subroutine f_trio

    ! Every returned root must be a genuine zero: |f| <= fbar and within zbar
    ! of an exact construction zero. A false midpoint root (|f| ~ 1e-8) fails.
    subroutine check_no_false_root(mmax, nfail)
        integer, intent(in)    :: mmax
        integer, intent(inout) :: nfail

        complex(dp), allocatable :: roots(:), fvals(:)
        integer,     allocatable :: mult(:)
        integer :: nfound, k
        type(fortnum_status_t) :: status
        real(dp) :: dmin

        call complex_region_roots(f_trio, cmplx(-1.0_dp, -1.0_dp, dp), &
            cmplx(1.0_dp, 1.0_dp, dp), roots, fvals, mult, nfound, status, &
            m_max=mmax)

        do k = 1, nfound
            if (abs(fvals(k)) > fbar) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,a,2es20.11,a,es12.4,a,es10.2)") &
                    "FAIL A root ", k, " (", real(roots(k)), aimag(roots(k)), &
                    ") |f|=", abs(fvals(k)), " > ", fbar
            end if
            dmin = nearest_exact(roots(k))
            if (dmin > zbar) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,a,2es20.11,a,es12.4,a,es10.2)") &
                    "FAIL A root ", k, " (", real(roots(k)), aimag(roots(k)), &
                    ") dist-to-exact=", dmin, " > ", zbar
            end if
        end do
    end subroutine check_no_false_root

    ! With bisection isolating every zero the finder must return all three to
    ! the same accuracy bar (guards against over-rejection by the |f| guard).
    subroutine check_all_found(mmax, nfail)
        integer, intent(in)    :: mmax
        integer, intent(inout) :: nfail

        complex(dp), allocatable :: roots(:), fvals(:)
        integer,     allocatable :: mult(:)
        integer :: nfound, j, k, hit
        type(fortnum_status_t) :: status
        real(dp) :: dbest

        call complex_region_roots(f_trio, cmplx(-1.0_dp, -1.0_dp, dp), &
            cmplx(1.0_dp, 1.0_dp, dp), roots, fvals, mult, nfound, status, &
            m_max=mmax)

        if (status%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL B status=", status%code
            return
        end if

        do j = 1, 3
            hit = 0
            do k = 1, nfound
                dbest = abs(roots(k) - zexact(j))
                if (dbest <= zbar .and. abs(fvals(k)) <= fbar) hit = k
            end do
            if (hit == 0) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,a,2es20.11,a)") &
                    "FAIL B exact zero ", j, " (", real(zexact(j)), &
                    aimag(zexact(j)), ") not found to bar"
            end if
        end do
    end subroutine check_all_found

    real(dp) function nearest_exact(z) result(d)
        complex(dp), intent(in) :: z
        integer :: j
        real(dp) :: dj
        d = huge(1.0_dp)
        do j = 1, 3
            dj = abs(z - zexact(j))
            if (dj < d) d = dj
        end do
    end function nearest_exact

end program test_fortnum_region_roots_converge
