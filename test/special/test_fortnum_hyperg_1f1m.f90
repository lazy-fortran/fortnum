program test_fortnum_hyperg_1f1m
    ! Focused test for the modified-form 1F1 entry point hyperg_1f1m_a1, which
    ! the KiLCA FLR reduction consumes as F11m where
    !   M(1,b,z) = 1 + z/b + z^2/(b(b+1)) * (1 + F11m).
    !
    ! Reconstructing F11m from M cancels two ~1 quantities and divides by z^2,
    ! so at the small-z consumer regime (z ~ 1e-2, b = 1 + t1 - i*x2) the
    ! M-then-cancel route loses ~8 digits and breaks the 1e-8 golden bar.
    ! hyperg_1f1m_a1 returns F11m = M(1,b+2,z) - 1 directly, with no
    ! cancellation.  Reference values from mpmath.hyp1f1 (dps=50) at the worst
    ! consumer grid points (x1 in {0.1,0.3,2}, x2 in {0,1,2,4}); the relative
    ! error must stay under the golden bar.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1m_a1
    implicit none

    real(dp), parameter :: rel_tol = 1.0e-8_dp
    integer,  parameter :: ncase = 8

    ! b_re, b_im, z_re, f11m_re, f11m_im (z_im = 0 over the consumer grid).
    ! The flre conductivity sweep reaches large z = x1^2 with b = 1 + z - i*x2,
    ! so |b| > |z| there: the last three cases (z = 64, 64, 400) guard the
    ! series/asymptotic crossover that the |z|-only gate routed wrongly.
    real(dp), parameter :: bre(ncase) = &
        [1.01_dp, 1.01_dp, 1.01_dp, 1.09_dp, 5.0_dp, 65.0_dp, 65.0_dp, 401.0_dp]
    real(dp), parameter :: bim(ncase) = &
        [-4.0_dp, -1.0_dp, 0.0_dp, -2.0_dp, -2.0_dp, 0.0_dp, -8.0_dp, 0.0_dp]
    real(dp), parameter :: zre(ncase) = &
        [0.01_dp, 0.01_dp, 0.01_dp, 0.09_dp, 4.0_dp, 64.0_dp, 64.0_dp, 400.0_dp]
    real(dp), parameter :: fre(ncase) = [ &
        0.00120061968696303871_dp, 0.00299847135353516077_dp, &
        0.00333056063606658202_dp, 0.0207773222910076061_dp, &
        0.879941883520409987_dp, 7.78607612691535186_dp, &
        5.25741317247839837_dp, 22.583254728782295_dp]
    real(dp), parameter :: fim(ncase) = [ &
        0.00159965944697732065_dp, 0.000998121726977174834_dp, &
        0.0_dp, 0.0137081430292269098_dp, 0.420815058286077092_dp, &
        0.0_dp, 4.12522932270957618_dp, 0.0_dp]

    integer :: i, nfail
    complex(dp) :: b, z, ref, got
    real(dp) :: err, scale, worst
    type(fortnum_status_t) :: status

    nfail = 0
    worst = 0.0_dp

    do i = 1, ncase
        b   = cmplx(bre(i), bim(i), dp)
        z   = cmplx(zre(i), 0.0_dp, dp)
        ref = cmplx(fre(i), fim(i), dp)
        call hyperg_1f1m_a1(b, z, got, status)
        if (.not. status_ok(status)) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,a)") "FAIL case ", i, &
                " status: ", trim(status%msg)
            cycle
        end if
        scale = max(abs(ref), 1.0_dp)
        err   = abs(got - ref) / scale
        worst = max(worst, err)
        if (err > rel_tol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.5,a,2es24.16,a,2es24.16)") &
                "FAIL case ", i, " relerr=", err, &
                " ref=", real(ref, dp), aimag(ref), &
                " got=", real(got, dp), aimag(got)
        end if
    end do

    write (*, "(a,es13.6)") "worst rel err: ", worst
    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, " failures in ", ncase, &
            " cases"
        stop 1
    end if
    write (*, "(a,i0,a)") "PASS: ", ncase, " modified-form cases within tolerance"
    stop 0
end program test_fortnum_hyperg_1f1m
