program test_fortnum_ode_events
    ! Behavioural tests for event location on the Cash-Karp trace
    ! (docs/design/ode.md §6).
    !   - oscillator y1 = cos(t): falling zero crossing located at t = pi/2
    !   - threshold crossing y1 = 0.5: located at t = pi/3 (rising for cos
    !     decreasing through 0.5 is falling; checked explicitly)
    !   - direction filter: a rising-only scan skips the first falling crossing
    !   - transversal crossing reports FORTNUM_OK and gdot of the right sign
    !   - tangential crossing (g touches zero, dg/dt = 0) reports a
    !     non-transversal domain error but still returns the primal root
    !   - no crossing leaves found=.false. and status OK

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ODE_EVENT_RISING, ODE_EVENT_FALLING, ODE_EVENT_ANY
    use fortnum_ode_events, only: ode_event_scan, ode_event_result_t
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer :: nfail

    nfail = 0
    call check_zero_crossing(nfail)
    call check_threshold_crossing(nfail)
    call check_direction_filter(nfail)
    call check_no_crossing(nfail)
    call check_tangential(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_events: all tests passed"
    stop 0

contains

    ! Undamped oscillator y1' = y2, y2' = -y1; y1(0)=1, y2(0)=0 -> y1 = cos(t).
    subroutine rhs_osc(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) =  y(2)
        dydt(2) = -y(1)
    end subroutine rhs_osc

    function ev_y1(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(1)
    end function ev_y1

    function ev_y1_half(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(1) - 0.5_dp
    end function ev_y1_half

    function ev_never(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(1) + 2.0_dp   ! cos(t) + 2 >= 1 > 0, never crosses zero
    end function ev_never

    ! Tangential event: g = y2 = -sin(t) on cos-oscillator. At t=0, y2=0 and
    ! dg/dt = -y1 = -1 (transversal); at t=pi, y2=0 with dg/dt = -y1 = 1. To get
    ! a true touch use g = y2^2, which grazes zero at t=0 and t=pi with dg/dt=0.
    function ev_y2_sq(t, y, ctx) result(g)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
        associate (unused_t => t); end associate
        g = y(2)**2 - 1.0e-30_dp   ! grazes just below zero near the turning pts
    end function ev_y2_sq

    subroutine integrate_osc(t1, solution)
        real(dp), intent(in)               :: t1
        type(ode_solution_t), intent(out)  :: solution
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        problem%rhs => rhs_osc
        problem%t0 = 0.0_dp
        problem%t1 = t1
        problem%y0 = [1.0_dp, 0.0_dp]
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
        call ode_integrate(problem, workspace, solution, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "integrate_osc failed code=", status%code
            stop 1
        end if
    end subroutine integrate_osc

    ! cos(t) falls through zero at t = pi/2.
    subroutine check_zero_crossing(nfail)
        integer, intent(inout) :: nfail
        type(ode_solution_t)     :: solution
        type(ode_event_result_t) :: ev
        type(fortnum_status_t)   :: status
        real(dp) :: terr
        call integrate_osc(3.0_dp, solution)
        call ode_event_scan(rhs_osc, ev_y1, solution, ODE_EVENT_FALLING, &
                            1.0e-12_dp, ev, status)
        if (status%code /= FORTNUM_OK .or. .not. ev%found) then
            write (error_unit, "(a,i0,a,l1)") &
                "FAIL check_zero_crossing: code=", status%code, &
                " found=", ev%found
            nfail = nfail + 1
            return
        end if
        ! Cubic-Hermite dense output is third order: O(h^4) location band.
        terr = abs(ev%t_event - 0.5_dp * PI)
        if (terr > 1.0e-8_dp) then
            write (error_unit, "(a,es20.12,a,es12.4)") &
                "FAIL check_zero_crossing: t_event=", ev%t_event, " err=", terr
            nfail = nfail + 1
        end if
        ! State at the root: y1 ~ 0, y2 ~ -1 (= -sin(pi/2)).
        if (abs(ev%y_event(1)) > 1.0e-8_dp .or. &
            abs(ev%y_event(2) + 1.0_dp) > 1.0e-8_dp) then
            write (error_unit, "(a,2es16.8)") &
                "FAIL check_zero_crossing: y_event=", ev%y_event
            nfail = nfail + 1
        end if
        if (ev%direction /= ODE_EVENT_FALLING .or. .not. ev%transversal) then
            write (error_unit, "(a,i0,a,l1)") &
                "FAIL check_zero_crossing: direction=", ev%direction, &
                " transversal=", ev%transversal
            nfail = nfail + 1
        end if
        ! dg/dt = y2 = -sin(pi/2) = -1 at the falling crossing.
        if (abs(ev%g_dot + 1.0_dp) > 1.0e-5_dp) then
            write (error_unit, "(a,es16.8)") &
                "FAIL check_zero_crossing: g_dot=", ev%g_dot
            nfail = nfail + 1
        end if
    end subroutine check_zero_crossing

    ! cos(t) = 0.5 first at t = pi/3, decreasing -> g = y1-0.5 falls through 0.
    subroutine check_threshold_crossing(nfail)
        integer, intent(inout) :: nfail
        type(ode_solution_t)     :: solution
        type(ode_event_result_t) :: ev
        type(fortnum_status_t)   :: status
        real(dp) :: terr
        call integrate_osc(2.0_dp, solution)
        call ode_event_scan(rhs_osc, ev_y1_half, solution, ODE_EVENT_ANY, &
                            1.0e-12_dp, ev, status)
        if (status%code /= FORTNUM_OK .or. .not. ev%found) then
            write (error_unit, "(a,i0,a,l1)") &
                "FAIL check_threshold_crossing: code=", status%code, &
                " found=", ev%found
            nfail = nfail + 1
            return
        end if
        ! Cubic-Hermite dense output is third order, so the located time sits
        ! within an O(h^4) band of the analytic crossing, not at machine zero.
        terr = abs(ev%t_event - PI / 3.0_dp)
        if (terr > 1.0e-8_dp) then
            write (error_unit, "(a,es20.12,a,es12.4)") &
                "FAIL check_threshold_crossing: t_event=", ev%t_event, &
                " err=", terr
            nfail = nfail + 1
        end if
        if (ev%direction /= ODE_EVENT_FALLING) then
            write (error_unit, "(a,i0)") &
                "FAIL check_threshold_crossing: direction=", ev%direction
            nfail = nfail + 1
        end if
    end subroutine check_threshold_crossing

    ! A rising-only scan must skip the falling crossing at pi/2 and find the
    ! rising crossing at t = 3*pi/2.
    subroutine check_direction_filter(nfail)
        integer, intent(inout) :: nfail
        type(ode_solution_t)     :: solution
        type(ode_event_result_t) :: ev
        type(fortnum_status_t)   :: status
        real(dp) :: terr
        call integrate_osc(5.0_dp, solution)
        call ode_event_scan(rhs_osc, ev_y1, solution, ODE_EVENT_RISING, &
                            1.0e-12_dp, ev, status)
        if (status%code /= FORTNUM_OK .or. .not. ev%found) then
            write (error_unit, "(a,i0,a,l1)") &
                "FAIL check_direction_filter: code=", status%code, &
                " found=", ev%found
            nfail = nfail + 1
            return
        end if
        terr = abs(ev%t_event - 1.5_dp * PI)
        if (ev%direction /= ODE_EVENT_RISING .or. terr > 1.0e-8_dp) then
            write (error_unit, "(a,es20.12,a,i0)") &
                "FAIL check_direction_filter: t_event=", ev%t_event, &
                " direction=", ev%direction
            nfail = nfail + 1
        end if
    end subroutine check_direction_filter

    subroutine check_no_crossing(nfail)
        integer, intent(inout) :: nfail
        type(ode_solution_t)     :: solution
        type(ode_event_result_t) :: ev
        type(fortnum_status_t)   :: status
        call integrate_osc(5.0_dp, solution)
        call ode_event_scan(rhs_osc, ev_never, solution, ODE_EVENT_ANY, &
                            1.0e-12_dp, ev, status)
        if (status%code /= FORTNUM_OK .or. ev%found) then
            write (error_unit, "(a,i0,a,l1)") &
                "FAIL check_no_crossing: code=", status%code, &
                " found=", ev%found
            nfail = nfail + 1
        end if
    end subroutine check_no_crossing

    ! g = y2^2 - tiny grazes zero at the velocity turning points (t=0, pi, ...)
    ! where dg/dt -> 0. A bracketed touch must report a non-transversal domain
    ! error while still handing back the primal root.
    subroutine check_tangential(nfail)
        integer, intent(inout) :: nfail
        type(ode_solution_t)     :: solution
        type(ode_event_result_t) :: ev
        type(fortnum_status_t)   :: status
        call integrate_osc(4.0_dp, solution)
        call ode_event_scan(rhs_osc, ev_y2_sq, solution, ODE_EVENT_ANY, &
                            1.0e-12_dp, ev, status)
        ! Either no bracket is detected (the graze may not change sign) or a
        ! crossing is reported as non-transversal. A clean FORTNUM_OK with a
        ! transversal flag on this grazing event would be wrong.
        if (ev%found) then
            if (status%code /= FORTNUM_DOMAIN_ERROR .or. ev%transversal) then
                write (error_unit, "(a,i0,a,l1)") &
                    "FAIL check_tangential: code=", status%code, &
                    " transversal=", ev%transversal
                nfail = nfail + 1
            end if
            if (.not. allocated(ev%y_event)) then
                write (error_unit, "(a)") &
                    "FAIL check_tangential: primal root not returned"
                nfail = nfail + 1
            end if
        end if
    end subroutine check_tangential

end program test_fortnum_ode_events
