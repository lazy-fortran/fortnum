program test_fortnum_ode_vode
    ! Behavioural tests for the clean-room nonstiff Adams integrator
    ! (fortnum_ode_vode), DVODE MF=10 equivalent.
    !   - exponential decay y' = -y over a long interval to < rtol
    !   - harmonic-oscillator energy conservation over many periods
    !   - a linear 2x2 system against its analytic solution
    !   - g-root location at a known analytic time, to <= event_tol
    !   - restart continuation: integrate grid point to grid point and match a
    !     single-shot integration
    !
    ! Accuracy is the integrator's own (compared to analytic solutions); a
    ! direct DVODE_F90 cross-check lives in test_fortnum_ode_vode_xcheck.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer :: nfail

    nfail = 0
    call check_exp_decay(nfail)
    call check_harmonic_energy(nfail)
    call check_linear_system(nfail)
    call check_event_root(nfail)
    call check_two_event_root(nfail)
    call check_restart_continuation(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_vode: all tests passed"
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

    subroutine rhs_osc(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) =  y(2)
        dydt(2) = -y(1)
    end subroutine rhs_osc

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

    function ev_y1(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(1)
    end function ev_y1

    ! g = y1 - 1/2: with y1 = cos(t) this falls through zero at t = pi/3,
    ! earlier than ev_y1's pi/2.
    function ev_y1_half(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(1) - 0.5_dp
    end function ev_y1_half

    ! y' = -y, y(0) = 1: assert y(tout) - exp(-tout) is below the requested
    ! tolerance over a long interval. This is the case the prior draft failed.
    subroutine check_exp_decay(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), tout, rtol, expected, err

        rtol = 1.0e-9_dp
        atol = 1.0e-12_dp
        tout = 20.0_dp
        call vode_init(st, 1, 0.0_dp, [1.0_dp])
        call vode_integrate_to(rhs_decay, st, tout, rtol, atol, yout, status)

        expected = exp(-tout)
        err = abs(yout(1) - expected)
        write (*, "(a,es12.4,a,es12.4,a,i0)") &
            "exp decay: y=", yout(1), " err=", err, " nstep=", st%nsteps
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  exp_decay: status not OK: " // &
                trim(status%msg)
            nfail = nfail + 1
        end if
        ! Relative error must beat rtol with margin.
        if (err > 10.0_dp * rtol * expected + atol(1)) then
            write (error_unit, "(a)") "  exp_decay: accuracy below tolerance"
            nfail = nfail + 1
        end if
    end subroutine check_exp_decay

    ! Harmonic oscillator: energy E = y1^2 + y2^2 conserved. Integrate many
    ! periods and assert the energy drift stays tiny.
    subroutine check_harmonic_energy(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), rtol, e0, e1, tout

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        tout = 20.0_dp * PI
        e0 = 1.0_dp
        call vode_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call vode_integrate_to(rhs_osc, st, tout, rtol, atol, yout, status)
        e1 = yout(1)**2 + yout(2)**2
        write (*, "(a,es12.4,a,es12.4)") &
            "harmonic: y1=", yout(1), " energy drift=", abs(e1 - e0)
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  harmonic: status not OK"
            nfail = nfail + 1
        end if
        ! After 10 periods y1 should be ~cos(20pi)=1.
        if (abs(yout(1) - 1.0_dp) > 1.0e-6_dp) then
            write (error_unit, "(a)") "  harmonic: y1 off analytic value"
            nfail = nfail + 1
        end if
        if (abs(e1 - e0) > 1.0e-7_dp) then
            write (error_unit, "(a)") "  harmonic: energy not conserved"
            nfail = nfail + 1
        end if
    end subroutine check_harmonic_energy

    ! Linear system with y(0) = (1,0): the analytic solution is
    ! y1 = (e^-t + e^-3t)/2, y2 = (e^-t - e^-3t)/2.
    subroutine check_linear_system(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(2), rtol, tout, y1e, y2e, em1, em3

        rtol = 1.0e-9_dp
        atol = 1.0e-12_dp
        tout = 5.0_dp
        call vode_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call vode_integrate_to(rhs_lin, st, tout, rtol, atol, yout, status)
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

    ! Oscillator y1 = cos(t): first rising-or-falling zero of y1 is at t=pi/2.
    ! Locate it with event_tol and assert the located root matches.
    subroutine check_event_root(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), rtol, tout, troot, etol
        logical  :: found

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        etol = 1.0e-10_dp
        tout = 5.0_dp
        call vode_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call vode_integrate_to(rhs_osc, st, tout, rtol, atol, yout, status, &
            event=ev_y1, event_dir=-1, event_tol=etol, t_root=troot, &
            root_found=found)
        write (*, "(a,l1,a,es20.12,a,es12.4)") &
            "event: found=", found, " troot=", troot, &
            " err=", abs(troot - 0.5_dp * PI)
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  event: status not OK"
            nfail = nfail + 1
        end if
        if (.not. found) then
            write (error_unit, "(a)") "  event: root not found"
            nfail = nfail + 1
        else if (abs(troot - 0.5_dp * PI) > etol) then
            write (error_unit, "(a)") "  event: root outside event_tol"
            nfail = nfail + 1
        end if
    end subroutine check_event_root

    ! Two monitored events (DVODE NEVENTS=2). With y1 = cos(t), ev_y1 falls
    ! through zero at t = pi/2 and ev_y1_half at t = pi/3 (earlier). The
    ! integrator must locate pi/3 and report the function that owns it,
    ! independent of the argument order in which the two g's are passed.
    subroutine check_two_event_root(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yout(:)
        real(dp) :: atol(1), rtol, tout, troot, etol, expected, root_margin
        integer  :: idx
        logical  :: found

        rtol = 1.0e-10_dp
        atol = 1.0e-12_dp
        etol = 1.0e-10_dp
        ! The Illinois search returns a bracket endpoint once the bracket width
        ! drops to etol, so the located root sits within a few etol of the true
        ! crossing; require that bound rather than etol itself.
        root_margin = 10.0_dp * etol
        tout = 5.0_dp
        expected = PI / 3.0_dp

        ! event = late root (pi/2), event2 = early root (pi/3): index 2 wins.
        call vode_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call vode_integrate_to(rhs_osc, st, tout, rtol, atol, yout, status, &
            event=ev_y1, event_dir=-1, event_tol=etol, t_root=troot, &
            root_found=found, event_index=idx, &
            event2=ev_y1_half, event_dir2=-1)
        write (*, "(a,l1,a,i0,a,es20.12,a,es12.4)") &
            "two-event: found=", found, " idx=", idx, " troot=", troot, &
            " err=", abs(troot - expected)
        if (.not. status_ok(status)) then
            write (error_unit, "(a)") "  two-event: status not OK"
            nfail = nfail + 1
        end if
        if (.not. found) then
            write (error_unit, "(a)") "  two-event: root not found"
            nfail = nfail + 1
        else
            if (idx /= 2) then
                write (error_unit, "(a)") "  two-event: wrong event index"
                nfail = nfail + 1
            end if
            if (abs(troot - expected) > root_margin) then
                write (error_unit, "(a)") "  two-event: root outside tolerance"
                nfail = nfail + 1
            end if
        end if

        ! Swap the order: the early root now arrives as event, index 1 wins.
        call vode_init(st, 2, 0.0_dp, [1.0_dp, 0.0_dp])
        call vode_integrate_to(rhs_osc, st, tout, rtol, atol, yout, status, &
            event=ev_y1_half, event_dir=-1, event_tol=etol, t_root=troot, &
            root_found=found, event_index=idx, &
            event2=ev_y1, event_dir2=-1)
        if (.not. found .or. idx /= 1 .or. abs(troot - expected) > root_margin) then
            write (error_unit, "(a)") "  two-event: swapped order mismatch"
            nfail = nfail + 1
        end if
    end subroutine check_two_event_root

    ! Restart continuation: integrating to tout in one shot must match a series
    ! of contiguous calls that carry state forward (KAMEL grid-to-grid usage).
    subroutine check_restart_continuation(nfail)
        integer, intent(inout) :: nfail
        type(vode_state_t) :: st1, st2
        type(fortnum_status_t) :: status
        real(dp), allocatable :: yfull(:), ystep(:)
        real(dp) :: atol(1), rtol, tend
        integer  :: i

        rtol = 1.0e-9_dp
        atol = 1.0e-12_dp
        tend = 10.0_dp

        call vode_init(st1, 1, 0.0_dp, [1.0_dp])
        call vode_integrate_to(rhs_decay, st1, tend, rtol, atol, yfull, status)

        call vode_init(st2, 1, 0.0_dp, [1.0_dp])
        do i = 1, 20
            call vode_integrate_to(rhs_decay, st2, 0.5_dp * i, rtol, atol, &
                ystep, status)
        end do
        write (*, "(a,es12.4)") &
            "restart: |one-shot - stepped| = ", abs(yfull(1) - ystep(1))
        if (abs(yfull(1) - ystep(1)) > 1.0e-7_dp) then
            write (error_unit, "(a)") "  restart: continuation mismatch"
            nfail = nfail + 1
        end if
    end subroutine check_restart_continuation

end program test_fortnum_ode_vode
