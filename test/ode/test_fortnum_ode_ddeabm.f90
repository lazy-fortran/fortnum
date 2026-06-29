program test_fortnum_ode_ddeabm
    ! Behavioural tests for the clean-room variable-order Adams PECE integrator
    ! (fortnum_ode_ddeabm), SLATEC DDEABM equivalent.
    !   - exponential decay y' = -y over a long interval to < tol
    !   - a small linear 2x2 system against its analytic solution
    !   - restart continuation: integrate grid point to grid point and match a
    !     single-shot integration (the KAMEL grid-march pattern)
    !   - tstop hard-stop: never step past tstop, return exactly at tout/tstop
    !
    ! Accuracy is the integrator's own (compared to analytic solutions); the
    ! direct SLATEC DDEABM cross-check lives in a separate driver.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_ode_ddeabm, only: ddeabm_state_t, ddeabm_init, &
        ddeabm_integrate_to
    implicit none

    integer :: nfail

    nfail = 0
    call check_exp_decay(nfail)
    call check_linear_system(nfail)
    call check_restart_continuation(nfail)
    call check_tstop_hard_stop(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_ddeabm: all tests passed"
    stop 0

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)
    end subroutine rhs_decay

    ! y1' = -2 y1 + y2, y2' = y1 - 2 y2; eigenvalues -1, -3.
    subroutine rhs_lin(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -2.0_dp * y(1) + y(2)
        dydt(2) =  y(1) - 2.0_dp * y(2)
    end subroutine rhs_lin

    ! y' = -y, y(0) = 1: assert y(tout) - exp(-tout) is below the requested
    ! tolerance over a long interval.
    subroutine check_exp_decay(nfail)
        integer, intent(inout) :: nfail
        type(ddeabm_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), tout, rtol, expected, err

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        tout = 20.0_dp
        call ddeabm_init(st, 1, 0.0_dp, [1.0_dp])
        call ddeabm_integrate_to(rhs_decay, st, tout, rtol, atol, yout, status)

        expected = exp(-tout)
        err = abs(yout(1) - expected)
        write (*, "(a,es12.4,a,es12.4,a,i0)") &
            "exp decay: y=", yout(1), " err=", err, " nstep=", st%nsteps
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  exp_decay: status not OK: "// &
                trim(status%msg)
            nfail = nfail + 1
        end if
        if (err > 10.0_dp * rtol * expected + atol(1)) then
            write (error_unit, "(a)") "  exp_decay: accuracy below tolerance"
            nfail = nfail + 1
        end if
    end subroutine check_exp_decay

    ! Linear system with y(0) = (1,0): analytic solution
    ! y1 = (e^-t + e^-3t)/2, y2 = (e^-t - e^-3t)/2.
    subroutine check_linear_system(nfail)
        integer, intent(inout) :: nfail
        type(ddeabm_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(2), rtol, tout, y1e, y2e, em1, em3

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        tout = 5.0_dp
        call ddeabm_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call ddeabm_integrate_to(rhs_lin, st, tout, rtol, atol, yout, status)
        em1 = exp(-tout)
        em3 = exp(-3.0_dp * tout)
        y1e = 0.5_dp * (em1 + em3)
        y2e = 0.5_dp * (em1 - em3)
        write (*, "(a,es12.4,a,es12.4)") &
            "linear: err1=", abs(yout(1) - y1e), " err2=", abs(yout(2) - y2e)
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  linear: status not OK"
            nfail = nfail + 1
        end if
        if (max(abs(yout(1) - y1e), abs(yout(2) - y2e)) > 1.0e-8_dp) then
            write (error_unit, "(a)") "  linear: accuracy below tolerance"
            nfail = nfail + 1
        end if
    end subroutine check_linear_system

    ! Restart continuation: integrating to tend in one shot must match a series
    ! of contiguous calls that carry state forward (KAMEL grid-to-grid usage).
    subroutine check_restart_continuation(nfail)
        integer, intent(inout) :: nfail
        type(ddeabm_state_t) :: st1, st2
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yfull(:), ystep(:)
        real(dp) :: atol(1), rtol, tend
        integer  :: i

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        tend = 10.0_dp

        call ddeabm_init(st1, 1, 0.0_dp, [1.0_dp])
        call ddeabm_integrate_to(rhs_decay, st1, tend, rtol, atol, yfull, &
                                 status)

        call ddeabm_init(st2, 1, 0.0_dp, [1.0_dp])
        do i = 1, 20
            call ddeabm_integrate_to(rhs_decay, st2, 0.5_dp * i, rtol, atol, &
                                     ystep, status)
        end do
        write (*, "(a,es12.4)") &
            "restart: |one-shot - stepped| = ", abs(yfull(1) - ystep(1))
        if (abs(yfull(1) - ystep(1)) > 1.0e-7_dp) then
            write (error_unit, "(a)") "  restart: continuation mismatch"
            nfail = nfail + 1
        end if
    end subroutine check_restart_continuation

    ! tstop hard-stop (DDEABM INFO(4)=1): set tstop = tend; march the grid up
    ! to and including tend. The integrator must never step past tstop and must
    ! return exactly at each tout (the last one equal to tstop).
    subroutine check_tstop_hard_stop(nfail)
        integer, intent(inout) :: nfail
        type(ddeabm_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), rtol, tend, tout, expected, err, errmax
        integer  :: i, ngrid

        rtol = 1.0e-11_dp
        atol = 1.0e-13_dp
        tend = 4.0_dp
        ngrid = 8
        errmax = 0.0_dp

        call ddeabm_init(st, 1, 0.0_dp, [1.0_dp])
        do i = 1, ngrid
            tout = tend * real(i, dp) / real(ngrid, dp)
            call ddeabm_integrate_to(rhs_decay, st, tout, rtol, atol, yout, &
                                     status, tstop=tend)
            if (.not. status_ok(status)) then
                write (error_unit, "(a,i0,a)") "  tstop: status not OK at i=", &
                    i, " "//trim(status%msg)
                nfail = nfail + 1
                return
            end if
            ! Reported point must be exactly tout, and the internal mesh top x
            ! must never have moved past tstop.
            if (abs(st%t - tout) > 1.0e-13_dp) then
                write (error_unit, "(a)") "  tstop: did not return at tout"
                nfail = nfail + 1
            end if
            if ((st%x - tend) * st%delsgn > 1.0e-12_dp) then
                write (error_unit, "(a,es20.12)") &
                    "  tstop: internal mesh stepped past tstop, x=", st%x
                nfail = nfail + 1
            end if
            expected = exp(-tout)
            err = abs(yout(1) - expected)
            errmax = max(errmax, err)
        end do
        write (*, "(a,es20.12,a,es12.4)") &
            "tstop: final x=", st%x, " max err vs analytic=", errmax
        if (errmax > 1.0e-8_dp) then
            write (error_unit, "(a)") "  tstop: accuracy below tolerance"
            nfail = nfail + 1
        end if
    end subroutine check_tstop_hard_stop

end program test_fortnum_ode_ddeabm
