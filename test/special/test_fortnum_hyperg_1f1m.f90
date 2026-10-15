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
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_CONVERGENCE_ERROR
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1m_a1, hyperg_1f1
    implicit none

    real(dp), parameter :: rel_tol = 1.0e-8_dp
    integer,  parameter :: ncase = 12

    ! b_re, b_im, z_re, f11m_re, f11m_im (z_im = 0 over the consumer grid).
    ! The flre conductivity sweep reaches large z = x1^2 with b = 1 + z - i*x2,
    ! so |b| > |z| there: cases 6-8 (z = 64, 64, 400) guard the
    ! series/asymptotic crossover that the |z|-only gate routed wrongly.
    real(dp), parameter :: bre(ncase) = &
        [1.01_dp, 1.01_dp, 1.01_dp, 1.09_dp, 5.0_dp, 65.0_dp, 65.0_dp, &
        401.0_dp, 90001.0_dp, 9000001.0_dp, 36000001.0_dp, 81000001.0_dp]
    real(dp), parameter :: bim(ncase) = &
        [-4.0_dp, -1.0_dp, 0.0_dp, -2.0_dp, -2.0_dp, 0.0_dp, -8.0_dp, &
        0.0_dp, 30.0_dp, 300.0_dp, 600.0_dp, 900.0_dp]
    real(dp), parameter :: zre(ncase) = &
        [0.01_dp, 0.01_dp, 0.01_dp, 0.09_dp, 4.0_dp, 64.0_dp, 64.0_dp, &
        400.0_dp, 90000.0_dp, 9000000.0_dp, 36000000.0_dp, 81000000.0_dp]
    real(dp), parameter :: fre(ncase) = [ &
        0.00120061968696303871_dp, 0.00299847135353516077_dp, &
        0.00333056063606658202_dp, 0.0207773222910076061_dp, &
        0.879941883520409987_dp, 7.78607612691535186_dp, &
        5.25741317247839837_dp, 22.583254728782295_dp, &
        371.4882576664374786629473417784377835263_dp, &
        3738.547494503537187431787551031365396016_dp, &
        7479.736478815122997723588229889770878149_dp, &
        11220.92588749229619473858523142111182721_dp]
    real(dp), parameter :: fim(ncase) = [ &
        0.00159965944697732065_dp, 0.000998121726977174834_dp, &
        0.0_dp, 0.0137081430292269098_dp, 0.420815058286077092_dp, &
        0.0_dp, 4.12522932270957618_dp, 0.0_dp, &
        -29.58998894861820577541750682023263132703_dp, &
        -298.6902003639461459460768555566784092585_dp, &
        -597.6921090942329854733859910692603395003_dp, &
        -896.6940767612271714482627615862643013485_dp]

    ! Last four references independently integrate Euler's representation at
    ! mpmath.mp.dps=60, rather than using a hypergeometric series:
    ! (b+1)*quad(lambda u: exp(z*(-expm1(-u))-(b+1)*u), [0,1,inf])-1.
    ! Here b=1+x*x+0.1i*x, z=x*x for x=300,3000,6000,9000. The final three
    ! cases exceed the former 5000-term limit; require accurate SUCCESS.

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
        if (i >= 9 .and. err > 1.0e-11_dp) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.5)") &
                "FAIL large-order integral case ", i, " relerr=", err
        end if
    end do

    ! A much wider near-transition series must still fail explicitly when
    ! the bounded work budget is exhausted; increasing it is not convergence.
    call hyperg_1f1m_a1(cmplx(1.0e18_dp, 0.0_dp, dp), &
        cmplx(1.0e18_dp, 0.0_dp, dp), got, status)
    if (status%code /= FORTNUM_CONVERGENCE_ERROR) then
        nfail = nfail + 1
        write (error_unit, "(a)") "FAIL: exhausted series budget reported success"
    end if

    ! Preserve terminating-polynomial behavior when Re(b) is negative and
    ! the sufficient infinite-tail ratio bound cannot yet be used.
    b = cmplx(-1000000000.5_dp, 0.0_dp, dp)
    z = cmplx(1.0_dp, 0.0_dp, dp)
    call hyperg_1f1(cmplx(-2.0_dp, 0.0_dp, dp), b, z, got, status)
    if (.not. status_ok(status)) nfail = nfail + 1
    ref = 1.0_dp - 2.0_dp * z / b + z**2 / (b * (b + 1.0_dp))
    if (abs(got - ref) > epsilon(1.0_dp)) nfail = nfail + 1
    call hyperg_1f1(cmplx(1.5_dp, 0.0_dp, dp), b, &
        cmplx(0.0_dp, 0.0_dp, dp), got, status)
    if (.not. status_ok(status)) nfail = nfail + 1
    if (abs(got - 1.0_dp) > epsilon(1.0_dp)) nfail = nfail + 1

    write (*, "(a,es13.6)") "worst rel err: ", worst
    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, " failures in ", ncase, &
            " cases"
        stop 1
    end if
    write (*, "(a,i0,a)") "PASS: ", ncase, " modified-form cases within tolerance"
    stop 0
end program test_fortnum_hyperg_1f1m
