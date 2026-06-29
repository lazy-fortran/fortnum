module fortnum_ode_events
    ! Event location for the adaptive Cash-Karp integrator (ADR docs/design/ode.md
    ! §6). Given an event function g(t, y, ctx) and the accepted-step trace from
    ! ode_integrate, find the first time g crosses zero in the requested
    ! direction, locate it to event_tol, and report the root time, the state
    ! there, and the bracketing step index. No module-level state: the caller
    ! owns the trace, the rhs, and the event function; this module reads them.
    !
    ! Derivative policy: trace_rule (ad.md §1, §4; ode.md §6).
    !   The event time t_event and event state y_event are differentiable only
    !   at a smooth, transversal crossing. Transversal means the total time
    !   derivative dg/dt = grad_y g . f is nonzero at the root, so the implicit
    !   function theorem turns dg = 0 into dt_event. A tangential crossing
    !   (|dg/dt| <= event_tol) or a non-smooth g at the root has no such
    !   derivative; ode_event_scan then sets FORTNUM_DOMAIN_ERROR naming the
    !   non-transversal event, while the primal root it returns stays valid. The
    !   #40 event sensitivity checks transversal_ok before propagating.
    !   Active: y0, ctx parameters, the frozen trace. Inactive: event_direction,
    !   event_tol, the bracket index, transversal_ok, status.
    !
    ! Dense location. Between two accepted mesh points the integrator advanced
    ! one Cash-Karp step. ode_integrate records the endpoint states t(i), y(:,i)
    ! and t(i+1), y(:,i+1). The step's continuous extension is a cubic Hermite
    ! polynomial matched to both endpoint states and both endpoint derivatives
    ! f(t,y); this is the standard third-order continuous output for an embedded
    ! RK pair (Hairer II.6) and is enough to root-find g to event_tol. The dense
    ! state feeds g, and Brent's method (Numerical Recipes, ch. 9) brackets the
    ! sign change. Endpoint derivatives come from two extra rhs calls per
    ! candidate step; the count is reported through nfev so a caller can audit
    ! the event-detection cost separately from the primal stepping.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_rhs_t, ode_event_t, &
        ode_solution_t, &
        ODE_EVENT_RISING, ODE_EVENT_FALLING, ODE_EVENT_ANY
    implicit none
    private

    public :: ode_event_scan
    public :: ode_event_result_t

    integer, parameter :: BRENT_MAX_ITER = 200

    ! Outcome of an event scan over the recorded trace. Holds enough primal
    ! metadata for a later sensitivity (#40): the root time and state, the
    ! bracketing step index into solution%t/%y, the realised crossing
    ! direction, the total time derivative dg/dt at the root, and whether the
    ! crossing was transversal (the precondition for differentiating t_event).
    type :: ode_event_result_t
        logical               :: found = .false.
        real(dp)              :: t_event = 0.0_dp
        real(dp), allocatable :: y_event(:)
        integer               :: bracket_step = 0 ! root lies in [t(i), t(i+1)]
        integer               :: direction = ODE_EVENT_ANY
        real(dp)              :: g_dot = 0.0_dp ! dg/dt at the root
        logical               :: transversal = .false.
        integer               :: nfev = 0 ! rhs calls spent locating
    end type ode_event_result_t

contains

    ! Scan the accepted-step trace in solution for the first event consistent
    ! with direction, then locate it to event_tol. The caller supplies the same
    ! rhs and event function used to produce the trace and the optional ctx.
    !
    ! direction is one of ODE_EVENT_RISING / _FALLING / _ANY. event_tol is the
    ! absolute tolerance on the located time. On a clean transversal crossing
    ! status is FORTNUM_OK; on a tangential or non-smooth crossing the root is
    ! still returned but status is FORTNUM_DOMAIN_ERROR (ode.md §6). When no
    ! crossing is bracketed result%found stays .false. and status is
    ! FORTNUM_OK.
    subroutine ode_event_scan(rhs, event, solution, direction, event_tol, &
            result, status, ctx)
        procedure(ode_rhs_t)                 :: rhs
        procedure(ode_event_t)               :: event
        type(ode_solution_t),  intent(in)    :: solution
        integer,               intent(in)    :: direction
        real(dp),              intent(in)    :: event_tol
        type(ode_event_result_t), intent(out) :: result
        type(fortnum_status_t), intent(out)  :: status
        class(*), intent(in), optional       :: ctx

        integer  :: neq, npts, i, cross
        real(dp) :: ga, gb, tr, gr
        real(dp), allocatable :: ya(:), yb(:), fa(:), fb(:), yr(:), fr(:)

        call status_set(status, FORTNUM_OK, "")
        call clear_result(result)

        if (.not. allocated(solution%t) .or. .not. allocated(solution%y)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_event_scan: empty trace")
            return
        end if
        npts = size(solution%t)
        neq  = size(solution%y, 1)
        if (npts < 2 .or. neq < 1) return

        allocate(ya(neq), yb(neq), fa(neq), fb(neq), yr(neq), fr(neq))

        ! g and f at the left end of the first interval.
        ya = solution%y(:, 1)
        ga = event(solution%t(1), ya, ctx)
        call rhs(solution%t(1), ya, fa, ctx)
        result%nfev = result%nfev + 1

        do i = 1, npts - 1
            yb = solution%y(:, i + 1)
            gb = event(solution%t(i + 1), yb, ctx)
            call rhs(solution%t(i + 1), yb, fb, ctx)
            result%nfev = result%nfev + 1

            cross = crossing_direction(ga, gb)
            if (cross /= 0 .and. direction_admits(direction, cross)) then
                call locate_in_step(event, &
                    solution%t(i), solution%t(i + 1), ya, yb, fa, fb, &
                    event_tol, tr, yr, gr, ctx)
                call total_derivative(rhs, event, tr, yr, &
                    event_tol, result%g_dot, &
                    result%transversal, ctx, result%nfev, fr)

                result%found        = .true.
                result%t_event      = tr
                result%y_event      = yr
                result%bracket_step = i
                result%direction    = cross

                if (.not. result%transversal) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        "ode_event_scan: non-transversal crossing, " // &
                        "event time not differentiable")
                end if
                return
            end if

            ! Advance the left end to this point for the next interval.
            ya = yb
            ga = gb
            fa = fb
        end do
    end subroutine ode_event_scan

    ! --- internals ---

    subroutine clear_result(result)
        type(ode_event_result_t), intent(out) :: result
        result%found = .false.
        result%t_event = 0.0_dp
        result%bracket_step = 0
        result%direction = ODE_EVENT_ANY
        result%g_dot = 0.0_dp
        result%transversal = .false.
        result%nfev = 0
    end subroutine clear_result

    ! Sign-change classifier on a step. Returns +1 for a g - to + crossing
    ! (rising), -1 for + to -, 0 for no bracketed sign change. A zero exactly at
    ! the left endpoint is consumed by the previous interval; an exact zero at
    ! the right endpoint is treated as a crossing in the sign direction of ga.
    integer function crossing_direction(ga, gb) result(cross)
        real(dp), intent(in) :: ga, gb
        cross = 0
        if (ga < 0.0_dp .and. gb >= 0.0_dp) then
            cross = ODE_EVENT_RISING
        else if (ga > 0.0_dp .and. gb <= 0.0_dp) then
            cross = ODE_EVENT_FALLING
        end if
    end function crossing_direction

    logical function direction_admits(direction, cross) result(ok)
        integer, intent(in) :: direction, cross
        ok = (direction == ODE_EVENT_ANY) .or. (direction == cross)
    end function direction_admits

    ! Cubic-Hermite dense state on one step [ta, tb], matched to the endpoint
    ! states and endpoint derivatives (Hairer II.6 continuous output).
    pure subroutine hermite_state(ta, tb, ya, yb, fa, fb, t, yt)
        real(dp), intent(in)  :: ta, tb
        real(dp), intent(in)  :: ya(:), yb(:), fa(:), fb(:)
        real(dp), intent(in)  :: t
        real(dp), intent(out) :: yt(:)
        real(dp) :: h, s, h00, h10, h01, h11
        h = tb - ta
        s = (t - ta) / h
        h00 =  2.0_dp * s**3 - 3.0_dp * s**2 + 1.0_dp
        h10 =  s**3 - 2.0_dp * s**2 + s
        h01 = -2.0_dp * s**3 + 3.0_dp * s**2
        h11 =  s**3 - s**2
        yt = h00 * ya + h10 * h * fa + h01 * yb + h11 * h * fb
    end subroutine hermite_state

    ! Locate the root of g along the dense state on [ta, tb] by Brent's method.
    ! Brackets are guaranteed by the caller's sign change. tr is the root time,
    ! yr the dense state there, gr the residual g at the root. No rhs call is
    ! spent: the Hermite interpolant reuses the endpoint derivatives the scan
    ! already evaluated, so location is pure root-finding on g.
    subroutine locate_in_step(event, ta, tb, ya, yb, fa, fb, &
            event_tol, tr, yr, gr, ctx)
        procedure(ode_event_t)            :: event
        real(dp), intent(in)              :: ta, tb
        real(dp), intent(in)              :: ya(:), yb(:), fa(:), fb(:)
        real(dp), intent(in)              :: event_tol
        real(dp), intent(out)             :: tr
        real(dp), intent(out)             :: yr(:)
        real(dp), intent(out)             :: gr
        class(*), intent(in), optional    :: ctx

        real(dp) :: a, b, c, fpa, fpb, fpc, d, e
        real(dp) :: tol1, xm, p, q, r, s, eps
        integer  :: iter
        logical  :: only_two

        eps = epsilon(1.0_dp)
        a = ta
        b = tb
        fpa = g_at(event, ta, tb, ya, yb, fa, fb, a, yr, ctx)
        fpb = g_at(event, ta, tb, ya, yb, fa, fb, b, yr, ctx)
        c = b
        fpc = fpb
        only_two = .true. ! a and c coincide until the first contraction
        d = b - a
        e = d

        do iter = 1, BRENT_MAX_ITER
            if ((fpb > 0.0_dp .and. fpc > 0.0_dp) .or. &
                (fpb < 0.0_dp .and. fpc < 0.0_dp)) then
                c = a; fpc = fpa; d = b - a; e = d
                only_two = .true.
            end if
            if (abs(fpc) < abs(fpb)) then
                a = b; b = c; c = a
                fpa = fpb; fpb = fpc; fpc = fpa
                only_two = .true.
            end if
            tol1 = 2.0_dp * eps * abs(b) + 0.5_dp * event_tol
            xm = 0.5_dp * (c - b)
            ! Converged on a tight bracket or an exact residual zero.
            if (abs(xm) <= tol1 .or. .not. (abs(fpb) > 0.0_dp)) exit
            if (abs(e) >= tol1 .and. abs(fpa) > abs(fpb)) then
                s = fpb / fpa
                if (only_two) then
                    p = 2.0_dp * xm * s
                    q = 1.0_dp - s
                else
                    q = fpa / fpc
                    r = fpb / fpc
                    p = s * (2.0_dp * xm * q * (q - r) - (b - a) * (r - 1.0_dp))
                    q = (q - 1.0_dp) * (r - 1.0_dp) * (s - 1.0_dp)
                end if
                if (p > 0.0_dp) q = -q
                p = abs(p)
                if (2.0_dp * p < min(3.0_dp * xm * q - abs(tol1 * q), &
                    abs(e * q))) then
                    e = d
                    d = p / q
                else
                    d = xm; e = d
                end if
            else
                d = xm; e = d
            end if
            a = b; fpa = fpb
            only_two = .false. ! a now differs from c: three distinct points
            if (abs(d) > tol1) then
                b = b + d
            else
                b = b + sign(tol1, xm)
            end if
            fpb = g_at(event, ta, tb, ya, yb, fa, fb, b, yr, ctx)
        end do

        tr = b
        call hermite_state(ta, tb, ya, yb, fa, fb, tr, yr)
        gr = event(tr, yr, ctx)
    end subroutine locate_in_step

    ! g evaluated on the dense Hermite state at time t within [ta, tb].
    real(dp) function g_at(event, ta, tb, ya, yb, fa, fb, t, ytmp, ctx) &
            result(g)
        procedure(ode_event_t)         :: event
        real(dp), intent(in)           :: ta, tb
        real(dp), intent(in)           :: ya(:), yb(:), fa(:), fb(:)
        real(dp), intent(in)           :: t
        real(dp), intent(inout)        :: ytmp(:)
        class(*), intent(in), optional :: ctx
        call hermite_state(ta, tb, ya, yb, fa, fb, t, ytmp)
        g = event(t, ytmp, ctx)
    end function g_at

    ! Total time derivative dg/dt = grad_y g . f at the root, by a central
    ! difference of g along the trajectory. The state derivative there is f
    ! from one rhs call; g is differentiated by a small time perturbation that
    ! advances y by f. Transversality holds when |dg/dt| exceeds event_tol.
    subroutine total_derivative(rhs, event, tr, yr, event_tol, gdot, &
            transversal, ctx, nfev, fr)
        procedure(ode_rhs_t)             :: rhs
        procedure(ode_event_t)           :: event
        real(dp), intent(in)             :: tr
        real(dp), intent(in)             :: yr(:)
        real(dp), intent(in)             :: event_tol
        real(dp), intent(out)            :: gdot
        logical,  intent(out)            :: transversal
        class(*), intent(in), optional   :: ctx
        integer,  intent(inout)          :: nfev
        real(dp), intent(inout)          :: fr(:)
        real(dp) :: dt, gp, gm
        real(dp), allocatable :: yp(:), ym(:)
        integer :: n
        n = size(yr)
        allocate(yp(n), ym(n))
        call rhs(tr, yr, fr, ctx)
        nfev = nfev + 1
        ! Step along the trajectory: y(t +/- dt) ~= yr +/- dt * f(tr, yr).
        dt = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(tr))
        yp = yr + dt * fr
        ym = yr - dt * fr
        gp = event(tr + dt, yp, ctx)
        gm = event(tr - dt, ym, ctx)
        gdot = (gp - gm) / (2.0_dp * dt)
        transversal = abs(gdot) > event_tol
    end subroutine total_derivative

end module fortnum_ode_events
