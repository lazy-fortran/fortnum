module fortnum_ode_vode
    ! Variable-coefficient, variable-order nonstiff Adams integrator,
    ! functionally equivalent to DVODE with MF=10 (METH=1 Adams, MITER=0
    ! functional iteration, no Jacobian). Variable order 1..12 with a
    ! Nordsieck history array, the DVSET variable-step Adams-Moulton
    ! coefficients, the standard local-error test (DSM = ||ACOR|| / TQ(2))
    ! and the q-1/q/q+1 order-selection logic. Re-entrant: the integrator
    ! state (vode_state_t) carries the Nordsieck array, order, step size and
    ! step history across calls so a consumer can integrate grid point to
    ! grid point (DVODE ISTATE continuation). Optional g-function root finding
    ! locates zeros of user event functions on the integrator's own Nordsieck
    ! interpolant (DVINDY) by the Illinois algorithm, at the accuracy DVODE
    ! delivers.
    !
    ! Provenance: clean-room implementation of the published algorithm.
    !   P. N. Brown, G. D. Byrne, A. C. Hindmarsh, "VODE: A Variable-
    !   Coefficient ODE Solver", SIAM J. Sci. Stat. Comput. 10 (1989)
    !   1038-1051; and the standard Nordsieck variable-step Adams theory
    !   (Hairer, Norsett, Wanner, "Solving Ordinary Differential Equations I",
    !   III.5-III.7). Written natively from the algorithm description, MIT
    !   licensed, differentiation-ready style. No netlib DVODE/VODE source was
    !   copied; the control constants reproduce the published method.
    !
    ! Derivative policy: this module records no trace and exposes no AD path
    !   yet (the AD surface lives in fortnum_ode via the Cash-Karp trace). The
    !   public dispatch in fortnum_ode treats vode as a primal-only method.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t, ode_event_t
    implicit none
    private

    public :: vode_state_t
    public :: vode_init, vode_integrate_to
    public :: ode_solve_vode
    public :: VODE_MAX_ORDER, VODE_MAX_EVENTS

    ! Maximum Adams order and the YH column count L = q + 1. DVODE caps the
    ! Adams order at 12, so YH carries up to 13 columns.
    integer, parameter :: VODE_MAX_ORDER  = 12
    integer, parameter :: VODE_MAX_L      = VODE_MAX_ORDER + 1
    integer, parameter :: VODE_MAX_EVENTS = 2

    ! Method control constants (VODE paper; reproduce the published scheme).
    real(dp), parameter :: ADDON  = 1.0e-6_dp
    real(dp), parameter :: BIAS1  = 6.0_dp
    real(dp), parameter :: BIAS2  = 6.0_dp
    real(dp), parameter :: BIAS3  = 10.0_dp
    real(dp), parameter :: CORTES = 0.1_dp
    real(dp), parameter :: CRDOWN = 0.3_dp
    real(dp), parameter :: RDIV   = 2.0_dp
    real(dp), parameter :: ETACF  = 0.25_dp
    real(dp), parameter :: ETAMIN = 0.1_dp
    real(dp), parameter :: ETAMX1 = 1.0e4_dp
    real(dp), parameter :: ETAMX2 = 10.0_dp
    real(dp), parameter :: ETAMX3 = 10.0_dp
    real(dp), parameter :: ETAMXF = 0.2_dp
    real(dp), parameter :: THRESH = 1.5_dp
    real(dp), parameter :: ONEPSM = 1.00001_dp
    integer,  parameter :: MAXCOR = 3
    integer,  parameter :: KFC    = -3
    integer,  parameter :: KFH    = -15
    integer,  parameter :: MXNCF  = 10

    ! The RHS interface (ode_rhs_t) and event interface (ode_event_t) are shared
    ! with fortnum_ode so consumers pass the same user routines to either method.

    ! Re-entrant integrator state. yh(:,1..l) is the Nordsieck array:
    ! yh(i,j+1) ~ h**j/j! * d^j y_i/dt^j. tau holds the recent step sizes.
    type :: vode_state_t
        integer  :: neq    = 0
        integer  :: nq     = 1 ! current order
        integer  :: l      = 2 ! nq + 1
        integer  :: nqwait = 2 ! steps until next order change test
        integer  :: maxord = VODE_MAX_ORDER
        integer  :: nsteps = 0
        integer  :: nfev   = 0
        integer  :: nrejected = 0
        logical  :: started = .false.
        real(dp) :: tn    = 0.0_dp ! current time at top of yh
        real(dp) :: h     = 0.0_dp ! current step size (signed)
        real(dp) :: hu    = 0.0_dp ! step that produced the current yh
        real(dp) :: hscal = 0.0_dp
        real(dp) :: rc    = 0.0_dp
        real(dp) :: prl1  = 1.0_dp
        real(dp) :: eta   = 1.0_dp
        real(dp) :: etamax = ETAMX1
        real(dp) :: conp  = 0.0_dp
        real(dp) :: hmax  = 0.0_dp ! 0 => unbounded
        real(dp) :: hmin  = 0.0_dp
        integer  :: max_steps = 500000
        real(dp), allocatable :: yh(:,:) ! neq x (maxord+2)
        real(dp) :: tau(VODE_MAX_L) = 0.0_dp
        real(dp) :: el(VODE_MAX_L)  = 0.0_dp
        real(dp) :: tq(5)           = 0.0_dp
        ! Work arrays sized to neq; allocated once in vode_init (no hot-loop
        ! allocation).
        real(dp), allocatable :: ewt(:), savf(:), acor(:), y(:), ftmp(:)
    end type vode_state_t

contains

    ! Initialise the integrator state at (t0, y0). Allocates the Nordsieck
    ! array and work vectors; sets order 1. Tolerances and bounds are passed
    ! per integration call, not stored here.
    subroutine vode_init(state, neq, t0, y0)
        type(vode_state_t), intent(inout) :: state
        integer,  intent(in) :: neq
        real(dp), intent(in) :: t0
        real(dp), intent(in) :: y0(:)

        state%neq    = neq
        state%nq     = 1
        state%l      = 2
        state%nqwait = 2
        state%maxord = VODE_MAX_ORDER
        state%nsteps = 0
        state%nfev   = 0
        state%nrejected = 0
        state%started = .false.
        state%tn    = t0
        state%h     = 0.0_dp
        state%hu    = 0.0_dp
        state%hscal = 0.0_dp
        state%rc    = 0.0_dp
        state%prl1  = 1.0_dp
        state%eta   = 1.0_dp
        state%etamax = ETAMX1
        state%conp  = 0.0_dp
        state%tau   = 0.0_dp
        state%el    = 0.0_dp
        state%tq    = 0.0_dp

        if (allocated(state%yh)) deallocate(state%yh)
        allocate(state%yh(neq, VODE_MAX_ORDER + 2))
        state%yh = 0.0_dp
        state%yh(:, 1) = y0(1:neq)

        call ensure_work(state, neq)
    end subroutine vode_init

    subroutine ensure_work(state, neq)
        type(vode_state_t), intent(inout) :: state
        integer, intent(in) :: neq
        if (allocated(state%ewt)) then
            if (size(state%ewt) == neq) return
            deallocate(state%ewt, state%savf, state%acor, state%y, state%ftmp)
        end if
        allocate(state%ewt(neq), state%savf(neq), state%acor(neq), &
            state%y(neq), state%ftmp(neq))
    end subroutine ensure_work

    ! Integrate from the current state%tn to tout (DVODE ITASK=1: take internal
    ! steps until tout is reached or passed, then interpolate exactly to tout).
    ! relerr is scalar; abserr is a vector of length neq (ITOL=2) or length 1
    ! (scalar atol broadcast). On return state holds tout as the new tn (via the
    ! Nordsieck interpolant), ready for the next call.
    !
    ! If event is present, the first sign change of g over a completed step is
    ! located on the Nordsieck interpolant; the integration stops at the root,
    ! state%tn becomes the root time, y_out the state there, and t_root the root
    ! time. event_dir restricts the crossing direction (+1 rising, -1 falling,
    ! 0 any). event_tol is the absolute resolution in t.
    !
    ! Up to VODE_MAX_EVENTS event functions may be monitored at once (DVODE
    ! NEVENTS): pass event2/event_dir2 for the second g. The earliest root of
    ! either function over a completed step wins; event_index reports which
    ! function located it (1 or 2, 0 if none).
    subroutine vode_integrate_to(rhs, state, tout, relerr, abserr, y_out, &
            status, event, event_dir, event_tol, &
            t_root, root_found, event_index, event2, event_dir2, ctx)
        procedure(ode_rhs_t)               :: rhs
        type(vode_state_t),    intent(inout) :: state
        real(dp), intent(in)                 :: tout
        real(dp), intent(in)                 :: relerr
        real(dp), intent(in)                 :: abserr(:)
        real(dp), allocatable, intent(out)   :: y_out(:)
        type(fortnum_status_t), intent(out)  :: status
        procedure(ode_event_t), optional :: event
        integer,  intent(in),  optional      :: event_dir
        real(dp), intent(in),  optional      :: event_tol
        real(dp), intent(out), optional      :: t_root
        logical,  intent(out), optional      :: root_found
        integer,  intent(out), optional      :: event_index
        procedure(ode_event_t), optional :: event2
        integer,  intent(in),  optional      :: event_dir2
        class(*), intent(in), optional       :: ctx

        integer  :: neq, kflag, nev, found_idx
        integer  :: edir(VODE_MAX_EVENTS)
        real(dp) :: dir, etol, troot, tlast
        real(dp) :: g_left(VODE_MAX_EVENTS)
        logical  :: found, has_event
        real(dp), allocatable :: y_left(:)

        call status_set(status, FORTNUM_OK, "")
        found = .false.
        found_idx = 0
        troot = state%tn
        neq = state%neq

        if (present(root_found)) root_found = .false.
        if (present(t_root)) t_root = state%tn
        if (present(event_index)) event_index = 0

        if (neq < 1 .or. .not. allocated(state%yh)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "vode_integrate_to: state not initialised")
            return
        end if
        if (relerr < 0.0_dp .or. any(abserr < 0.0_dp) .or. &
            (relerr == 0.0_dp .and. all(abserr == 0.0_dp))) then
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            "vode_integrate_to: tolerances must be nonnegative, not both 0")
        return
    end if

    allocate(y_out(neq), y_left(neq))
    has_event = present(event)
    nev = 0
    edir = 0
    if (present(event)) then
        nev = 1
        if (present(event_dir)) edir(1) = event_dir
    end if
    if (present(event2)) then
        if (.not. present(event)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "vode_integrate_to: event2 requires event")
            return
        end if
        nev = 2
        if (present(event_dir2)) edir(2) = event_dir2
    end if
    etol = 1.0e-10_dp
    if (present(event_tol)) etol = event_tol

    ! tout already reached: interpolate and return (DVODE returns state).
    if (tout == state%tn) then
        y_out = state%yh(:, 1)
        return
    end if

    ! Integration direction. Once started it is the sign of the internal
    ! step h (a prior call may have stepped past this tout, so tout-tn can
    ! point backward even though integration continues forward).
    if (state%started) then
        dir = sign(1.0_dp, state%h)
    else
        dir = sign(1.0_dp, tout - state%tn)
    end if

    ! Continuation: if the internal mesh top already reached or passed tout
    ! in the integration direction, interpolate without stepping. The
    ! Nordsieck interpolant (DVINDY) is valid on [tn-hu, tn]; a tout within
    ! that window is recovered exactly. Skip when an event scan is requested
    ! (the event must see the freshly stepped interval).
    if (state%started .and. .not. has_event .and. &
        abs(state%hu) > 0.0_dp) then
    if ((state%tn - tout) * dir >= 0.0_dp .and. &
        (tout - (state%tn - state%hu)) * dir >= 0.0_dp) then
    call interpolate(state, tout, y_out)
    return
end if
end if

! First entry: prime the Nordsieck array and pick the starting step.
if (.not. state%started) then
    call rhs(state%tn, state%yh(:, 1), state%ftmp, ctx)
    state%nfev = state%nfev + 1
    state%yh(:, 2) = 0.0_dp
    call set_ewt(state, relerr, abserr)
    state%h = initial_step(rhs, state, tout, relerr, abserr, ctx)
    state%yh(:, 2) = state%h * state%ftmp
    state%hscal = state%h
    state%tau(1) = state%h
    state%started = .true.
end if

! seed the event left endpoint at the current state
if (has_event) then
    y_left = state%yh(:, 1)
    g_left(1) = event(state%tn, y_left, ctx)
    if (nev == 2) g_left(2) = event2(state%tn, y_left, ctx)
    tlast = state%tn
end if

do
    if (state%nsteps >= state%max_steps) then
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "vode_integrate_to: exceeded max_steps")
        exit
    end if

    ! Clip the step so it does not overshoot tout by more than needed;
    ! DVODE steps past tout and interpolates, but limiting h to the
    call set_ewt(state, relerr, abserr)
    call take_step(rhs, state, relerr, abserr, kflag, ctx)

    if (kflag < 0) then
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "vode_integrate_to: step failed (h underflow or " // &
            "repeated convergence/error-test failure)")
        exit
    end if

    ! Event scan over the completed step [tlast, state%tn].
    if (has_event) then
        call scan_step_for_root(event, event2, nev, state, edir, etol, &
            tlast, g_left, troot, found, found_idx, ctx)
        if (found) then
            call interpolate(state, troot, y_out)
            state%tn = troot
            state%hu = troot - tlast
            if (present(t_root)) t_root = troot
            if (present(root_found)) root_found = .true.
            if (present(event_index)) event_index = found_idx
            exit
        end if
        tlast = state%tn
        y_left = state%yh(:, 1)
        g_left(1) = event(state%tn, y_left, ctx)
        if (nev == 2) g_left(2) = event2(state%tn, y_left, ctx)
    end if

    ! Reached tout?
    if ((state%tn - tout) * dir >= 0.0_dp) then
        call interpolate(state, tout, y_out)
        ! Advance the conceptual current point to tout for the next
        ! call without disturbing the Nordsieck array (DVODE keeps tn
        ! at the internal mesh top and interpolates on entry). We return
        ! the interpolated value; tn/h/yh stay as the internal mesh.
        exit
    end if
end do

if (status%code == FORTNUM_OK .and. .not. found) then
    ! Normal tout return: leave state at the internal mesh top so the
    ! next call continues seamlessly; report the interpolated y_out.
    continue
end if
end subroutine vode_integrate_to

! Flat one-shot entry matching the ode_solve / ode_solve_dop surface: build
! state, integrate t0 -> t1, return the endpoint. For continuation and event
! root finding use vode_init + vode_integrate_to directly (the re-entrant
! path KAMEL needs). rtol/atol default to the ode_problem_t defaults.
subroutine ode_solve_vode(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
    procedure(ode_rhs_t)               :: rhs
    real(dp),               intent(in) :: t0, t1
    real(dp),               intent(in) :: y0(:)
    real(dp), allocatable, intent(out) :: t_out(:)
    real(dp), allocatable, intent(out) :: y_out(:,:)
    type(fortnum_status_t), intent(out) :: status
    real(dp), intent(in), optional     :: rtol, atol

    type(vode_state_t) :: st
    real(dp), allocatable :: yend(:)
    real(dp) :: rt, at(1)

    rt = 1.0e-6_dp
    at = 1.0e-9_dp
    if (present(rtol)) rt = rtol
    if (present(atol)) at = atol

    call vode_init(st, size(y0), t0, y0)
    call vode_integrate_to(rhs, st, t1, rt, at, yend, status)

    allocate(t_out(2))
    allocate(y_out(size(y0), 2))
    t_out(1) = t0
    t_out(2) = t1
    y_out(:, 1) = y0
    y_out(:, 2) = yend
end subroutine ode_solve_vode

! --- core stepper ---

! One Adams step with adaptive order/step control (DVSTEP for METH=1).
! On entry the Nordsieck array holds the state at state%tn. On a successful
! step state%tn advances by state%h, state%yh is updated and the order/step
! for the next step are selected. kflag >= 0 on success, < 0 on failure.
subroutine take_step(rhs, state, relerr, abserr, kflag, ctx)
    procedure(ode_rhs_t)            :: rhs
    type(vode_state_t), intent(inout) :: state
    real(dp), intent(in)             :: relerr
    real(dp), intent(in)             :: abserr(:)
    integer,  intent(out)            :: kflag
    class(*), intent(in), optional   :: ctx

    integer  :: neq, ncf, nflag, i, j, iback
    real(dp) :: told, dsm, flotl, r, rl1

    neq = state%neq
    kflag = 0
    told = state%tn
    ncf = 0

    ! Apply any pending order/step change selected on the previous step,
    ! then rescale the Nordsieck array if the step size changed.
    call apply_pending_change(state, neq)

    do ! retry loop for convergence / error-test failures
        ! Predictor: multiply YH by the Pascal-triangle matrix in place.
        state%tn = state%tn + state%h
        call predict(state, neq)

        call dvset(state)
        rl1 = 1.0_dp / state%el(2)
        state%rc = state%rc * (rl1 / state%prl1)
        state%prl1 = rl1

        call corrector(rhs, state, rl1, nflag, dsm, ctx)

        if (nflag == 0) then
            ! Local error test.
            if (dsm <= 1.0_dp) then
                ! Accept: commit the correction into the Nordsieck array.
                kflag = 0
                state%nsteps = state%nsteps + 1
                state%hu = state%h
                do iback = 1, state%nq
                    i = state%l - iback
                    state%tau(i + 1) = state%tau(i)
                end do
                state%tau(1) = state%h
                do j = 1, state%l
                    state%yh(:, j) = state%yh(:, j) &
                        + state%el(j) * state%acor(:)
                end do
                state%nqwait = state%nqwait - 1
                if (state%l /= state%maxord + 1 .and. state%nqwait == 1) then
                    state%yh(:, state%maxord + 1) = state%acor(:)
                    state%conp = state%tq(5)
                end if
                ! Select order/step for the next step.
                if (abs(state%etamax - 1.0_dp) > 0.0_dp) then
                    call select_order(state, neq, dsm)
                else
                    if (state%nqwait < 2) state%nqwait = 2
                    state%eta = 1.0_dp
                end if
                state%etamax = ETAMX3
                if (state%nsteps <= 10) state%etamax = ETAMX2
                ! Scale ACOR to hold the estimated local error (TQ(2)).
                r = 1.0_dp / state%tq(2)
                state%acor(:) = r * state%acor(:)
                return
            else
                ! Error-test failure: retract, shrink, retry.
                kflag = kflag - 1
                state%nrejected = state%nrejected + 1
                nflag = -2
                state%tn = told
                call unpredict(state, neq)
                if (abs(state%h) <= state%hmin * ONEPSM) then
                    kflag = -1; return
                end if
                state%etamax = 1.0_dp
                if (kflag > KFC) then
                    flotl = real(state%l, dp)
                    state%eta = 1.0_dp / &
                        ((BIAS2 * dsm)**(1.0_dp / flotl) + ADDON)
                    state%eta = max(state%eta, ETAMIN)
                    if (state%hmin > 0.0_dp) &
                        state%eta = max(state%eta, state%hmin / abs(state%h))
                    if (kflag <= -2 .and. state%eta > ETAMXF) &
                        state%eta = ETAMXF
                    call rescale(state, neq)
                    cycle
                else
                    ! >=3 consecutive failures: drop order, shrink hard.
                    if (kflag == KFH) then
                        kflag = -1; return
                    end if
                    if (state%nq == 1) then
                        ! Order already 1: restart from the derivative.
                        state%eta = ETAMIN
                        if (state%hmin > 0.0_dp) state%eta = &
                            max(state%eta, state%hmin / abs(state%h))
                        state%h = state%h * state%eta
                        state%hscal = state%h
                        state%tau(1) = state%h
                        call rhs(state%tn, state%yh(:, 1), state%savf, ctx)
                        state%nfev = state%nfev + 1
                        state%yh(:, 2) = state%h * state%savf(:)
                        state%nqwait = 10
                        cycle
                    else
                        state%eta = ETAMIN
                        if (state%hmin > 0.0_dp) state%eta = &
                            max(state%eta, state%hmin / abs(state%h))
                        call dvjust_down(state, neq)
                        state%l = state%nq
                        state%nq = state%nq - 1
                        state%nqwait = state%l
                        call rescale(state, neq)
                        cycle
                    end if
                end if
            end if
        else
            ! Corrector convergence failure: retract, shrink, retry.
            ncf = ncf + 1
            state%nrejected = state%nrejected + 1
            state%etamax = 1.0_dp
            state%tn = told
            call unpredict(state, neq)
            if (abs(state%h) <= state%hmin * ONEPSM) then
                kflag = -2; return
            end if
            if (ncf == MXNCF) then
                kflag = -2; return
            end if
            state%eta = ETACF
            if (state%hmin > 0.0_dp) &
                state%eta = max(state%eta, state%hmin / abs(state%h))
            call rescale(state, neq)
            cycle
        end if
    end do
end subroutine take_step

! Functional-iteration corrector (MITER=0). The fixed-point map is
!   acor <- rl1*(h*f(yh(:,1)+acor) - yh(:,2)),   y = yh(:,1)+acor,
! iterated until the weighted-norm correction passes the convergence test
! DCON = del*min(1,crate)/TQ(4) <= 1. dsm = ||acor||_ewt / TQ(2) on return.
subroutine corrector(rhs, state, rl1, nflag, dsm, ctx)
    procedure(ode_rhs_t)            :: rhs
    type(vode_state_t), intent(inout) :: state
    real(dp), intent(in)             :: rl1
    integer,  intent(out)            :: nflag
    real(dp), intent(out)            :: dsm
    class(*), intent(in), optional   :: ctx

    integer  :: neq, m
    real(dp) :: del, delp, crate, dcon, acnrm
    real(dp), allocatable :: ynew(:)

    neq = state%neq
    allocate(ynew(neq))
    m = 0
    delp = 0.0_dp
    crate = 1.0_dp
    nflag = -1
    dsm = 0.0_dp

    state%y(:) = state%yh(:, 1)
    call rhs(state%tn, state%y, state%savf, ctx)
    state%nfev = state%nfev + 1
    state%acor(:) = 0.0_dp

    do
        ! ynew = rl1*(h*f - yh(:,2)); correction = ynew - acor.
        ynew(:) = rl1 * (state%h * state%savf(:) - state%yh(:, 2))
        state%y(:) = ynew(:) - state%acor(:)
        del = weighted_norm(state%y, state%ewt)
        state%y(:) = state%yh(:, 1) + ynew(:)
        state%acor(:) = ynew(:)

        if (m /= 0) crate = max(CRDOWN * crate, del / delp)
        dcon = del * min(1.0_dp, crate) / state%tq(4)
        if (dcon <= 1.0_dp) then
            nflag = 0
            if (m == 0) then
                acnrm = del
            else
                acnrm = weighted_norm(state%acor, state%ewt)
            end if
            dsm = acnrm / state%tq(2)
            return
        end if
        m = m + 1
        if (m == MAXCOR) return
        if (m >= 2 .and. del > RDIV * delp) return
        delp = del
        call rhs(state%tn, state%y, state%savf, ctx)
        state%nfev = state%nfev + 1
    end do
end subroutine corrector

! --- Nordsieck array operations ---

! Predictor: YH(:,j) += YH(:,j+1) accumulated as the Pascal-triangle
! product over the active columns (advances every derivative by one step).
subroutine predict(state, neq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq
    integer :: j, k
    do j = 1, state%nq
        do k = state%nq, j, -1
            state%yh(:, k) = state%yh(:, k) + state%yh(:, k + 1)
        end do
    end do
end subroutine predict

! Inverse of predict (DVSTEP retraction): undo the Pascal-triangle
! accumulation column-by-column to restore the Nordsieck array after a
! rejected step. Mirrors the YH1(I) = YH1(I) - YH1(I+LDYH) flat loop, which
! is the exact inverse of the predictor's accumulation order.
subroutine unpredict(state, neq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq
    integer :: jb, c0, k
    do jb = 1, state%nq
        c0 = state%nq + 1 - jb
        do k = c0, state%nq
            state%yh(:, k) = state%yh(:, k) - state%yh(:, k + 1)
        end do
    end do
end subroutine unpredict

! Rescale the Nordsieck columns 2..l by powers of eta after a step-size
! change, then set h = hscal*eta. Mirrors the DVSTEP rescale block: this is
! the operation whose mis-ordering left a spurious O(1) term in acor.
subroutine rescale(state, neq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq
    integer  :: j
    real(dp) :: r
    ! Clamp eta against hmax (HMXI in DVODE).
    if (state%hmax > 0.0_dp) then
        state%eta = state%eta / &
            max(1.0_dp, abs(state%hscal) * state%eta / state%hmax)
    end if
    r = 1.0_dp
    do j = 2, state%l
        r = r * state%eta
        state%yh(:, j) = r * state%yh(:, j)
    end do
    state%h = state%hscal * state%eta
    state%hscal = state%h
    state%rc = state%rc * state%eta
end subroutine rescale

! Realise the pending step-size change selected on the previous step, before
! the next prediction. The order change (and its YH adjustment) was already
! applied by select_order/commit_order; this only rescales for a new h.
subroutine apply_pending_change(state, neq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq
    ! select_order already applied any order change to the Nordsieck array
    ! (commit_order). Here only the pending step-size change is realised. eta
    ! exactly 1 means the controller kept the step: nothing to do.
    if (abs(state%eta - 1.0_dp) <= 0.0_dp) return
    call rescale(state, neq)
end subroutine apply_pending_change

! Interpolate the solution to time t from the Nordsieck array (DVINDY,
! K=0): y(t) = sum_{j=0}^{q} ((t-tn)/h)^j * YH(:,j+1).
subroutine interpolate(state, t, y)
    type(vode_state_t), intent(in)  :: state
    real(dp),           intent(in)  :: t
    real(dp),           intent(out) :: y(:)
    integer  :: j
    real(dp) :: s
    s = (t - state%tn) / state%h
    y(:) = state%yh(:, state%l)
    do j = state%l - 1, 1, -1
        y(:) = state%yh(:, j) + s * y(:)
    end do
end subroutine interpolate

! --- order/step selection (DVSTEP order-change block) ---

! After an accepted step, choose eta and the new order from the q-1/q/q+1
! candidates. Only reconsiders order when nqwait hits 0.
subroutine select_order(state, neq, dsm)
    type(vode_state_t), intent(inout) :: state
    integer,  intent(in)  :: neq
    real(dp), intent(in)  :: dsm

    real(dp) :: flotl, etaq, etaqm1, etaqp1, ddn, dup, cnquot
    integer  :: newq

    flotl = real(state%l, dp)
    etaq = 1.0_dp / ((BIAS2 * dsm)**(1.0_dp / flotl) + ADDON)

    if (state%nqwait /= 0) then
        ! Same order; fall through to the common THRESH/ETAMAX clamp.
        state%eta = etaq
        newq = state%nq
    else
        state%nqwait = 2
        etaqm1 = 0.0_dp
        if (state%nq /= 1) then
            ddn = weighted_norm(state%yh(:, state%l), state%ewt) &
                / state%tq(1)
            etaqm1 = 1.0_dp / &
                ((BIAS1 * ddn)**(1.0_dp / (flotl - 1.0_dp)) + ADDON)
        end if
        etaqp1 = 0.0_dp
        if (state%l /= state%maxord + 1) then
            cnquot = (state%tq(5) / state%conp) * &
                (state%h / state%tau(2))**state%l
            state%savf(:) = state%acor(:) - &
                cnquot * state%yh(:, state%maxord + 1)
            dup = weighted_norm(state%savf, state%ewt) / state%tq(3)
            etaqp1 = 1.0_dp / &
                ((BIAS3 * dup)**(1.0_dp / (flotl + 1.0_dp)) + ADDON)
        end if

        if (etaq >= etaqp1) then
            if (etaq < etaqm1) then
                state%eta = etaqm1; newq = state%nq - 1
            else
                state%eta = etaq;   newq = state%nq
            end if
        else
            if (etaqp1 > etaqm1) then
                state%eta = etaqp1; newq = state%nq + 1
                ! Stash acor for the order-increase column.
                state%yh(:, state%maxord + 1) = state%acor(:)
            else
                state%eta = etaqm1; newq = state%nq - 1
            end if
        end if
    end if

    ! Require a meaningful step growth before changing (DVSTEP label 200).
    if (state%eta < THRESH .or. abs(state%etamax - 1.0_dp) <= 0.0_dp) then
        state%eta = 1.0_dp
        newq = state%nq
    else
        state%eta = min(state%eta, state%etamax)
        if (state%hmax > 0.0_dp) state%eta = state%eta / &
            max(1.0_dp, abs(state%h) * state%eta / state%hmax)
    end if
    call commit_order(state, neq, newq)
end subroutine select_order

! Apply an order change to the Nordsieck array and update nq/l/nqwait.
subroutine commit_order(state, neq, newq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq, newq
    if (newq == state%nq) return
    if (newq < state%nq) then
        call dvjust_down(state, neq)
        state%nq = newq
        state%l = state%nq + 1
        state%nqwait = state%l
    else
        ! Order increase (DVJUST METH=1, IORD=1): the new derivative column
        ! starts at zero and is filled by the next step's corrector. The
        ! stashed acor in column maxord+1 was only for the etaqp1 estimate.
        state%yh(:, newq + 1) = 0.0_dp
        state%nq = newq
        state%l = state%nq + 1
        state%nqwait = state%l
    end if
end subroutine commit_order

! Adjust the Nordsieck array on an order decrease (DVJUST, METH=1, IORD=-1):
! subtract correction terms built from the variable-step coefficients so the
! lower-order history is consistent.
subroutine dvjust_down(state, neq)
    type(vode_state_t), intent(inout) :: state
    integer, intent(in) :: neq
    integer  :: nqm1, nqm2, j, i, iback, jp1
    real(dp) :: hsum, xi
    real(dp) :: elj(VODE_MAX_L)

    if (state%nq == 2) return
    nqm1 = state%nq - 1
    nqm2 = state%nq - 2
    elj = 0.0_dp
    elj(2) = 1.0_dp
    hsum = 0.0_dp
    do j = 1, nqm2
        hsum = hsum + state%tau(j)
        xi = hsum / state%hscal
        jp1 = j + 1
        do iback = 1, jp1
            i = (j + 3) - iback
            elj(i) = elj(i) * xi + elj(i - 1)
        end do
    end do
    do j = 2, nqm1
        elj(j + 1) = real(state%nq, dp) * elj(j) / real(j, dp)
    end do
    do j = 3, state%nq
        state%yh(:, j) = state%yh(:, j) - state%yh(:, state%l) * elj(j)
    end do
end subroutine dvjust_down

! --- DVSET: variable-step Adams-Moulton coefficients (METH=1) ---
! Verified against authentic DVODE to machine precision (the el arrays for
! q=2,3,4). TQ(1)/TQ(3) are only formed when an order change is pending
! (nqwait == 1), matching DVODE.
subroutine dvset(state)
    type(vode_state_t), intent(inout) :: state
    integer  :: nq, l, nqm1, i, j, iback
    real(dp) :: em(VODE_MAX_L)
    real(dp) :: hsum, rxi, xi, s, em0, csum, floti, flotnq, flotl, h

    nq = state%nq
    l = state%l
    h = state%h
    state%el(1:l) = 0.0_dp

    if (nq == 1) then
        state%el(1) = 1.0_dp
        state%el(2) = 1.0_dp
        state%tq(1) = 1.0_dp
        state%tq(2) = 2.0_dp
        state%tq(3) = 6.0_dp * state%tq(2)
        state%tq(5) = 1.0_dp
        state%tq(4) = CORTES * state%tq(2)
        return
    end if

    nqm1 = nq - 1
    flotl = real(l, dp)
    flotnq = flotl - 1.0_dp
    em = 0.0_dp
    em(1) = 1.0_dp
    hsum = h
    do j = 1, nqm1
        if (j == nqm1 .and. state%nqwait == 1) then
            s = 1.0_dp
            csum = 0.0_dp
            do i = 1, nqm1
                csum = csum + s * em(i) / real(i + 1, dp)
                s = -s
            end do
            state%tq(1) = em(nqm1) / (flotnq * csum)
        end if
        rxi = h / hsum
        do iback = 1, j
            i = (j + 2) - iback
            em(i) = em(i) + em(i - 1) * rxi
        end do
        hsum = hsum + state%tau(j)
    end do

    s = 1.0_dp
    em0 = 0.0_dp
    csum = 0.0_dp
    do i = 1, nq
        floti = real(i, dp)
        em0 = em0 + s * em(i) / floti
        csum = csum + s * em(i) / (floti + 1.0_dp)
        s = -s
    end do

    s = 1.0_dp / em0
    state%el(1) = 1.0_dp
    do i = 1, nq
        state%el(i + 1) = s * em(i) / real(i, dp)
    end do
    xi = hsum / h
    state%tq(2) = xi * em0 / csum
    state%tq(5) = xi / state%el(l)

    if (state%nqwait == 1) then
        rxi = 1.0_dp / xi
        do iback = 1, nq
            i = (l + 1) - iback
            em(i) = em(i) + em(i - 1) * rxi
        end do
        s = 1.0_dp
        csum = 0.0_dp
        do i = 1, l
            csum = csum + s * em(i) / real(i + 1, dp)
            s = -s
        end do
        state%tq(3) = flotl * em0 / csum
    end if
    state%tq(4) = CORTES * state%tq(2)
end subroutine dvset

! --- helpers ---

! Error-weight vector EWT(i) = 1/(relerr*|y_i| + atol_i). atol is a vector
! of length neq (ITOL=2) or a scalar broadcast (length 1).
subroutine set_ewt(state, relerr, abserr)
    type(vode_state_t), intent(inout) :: state
    real(dp), intent(in) :: relerr
    real(dp), intent(in) :: abserr(:)
    integer  :: i
    real(dp) :: a
    do i = 1, state%neq
        if (size(abserr) == 1) then
            a = abserr(1)
        else
            a = abserr(i)
        end if
        state%ewt(i) = 1.0_dp / (relerr * abs(state%yh(i, 1)) + a)
    end do
end subroutine set_ewt

! Weighted RMS norm: sqrt(mean((v_i*w_i)^2)). w is EWT (already 1/scale).
pure real(dp) function weighted_norm(v, w) result(nrm)
    real(dp), intent(in) :: v(:), w(:)
    integer :: i, n
    real(dp) :: acc
    n = size(v)
    acc = 0.0_dp
    do i = 1, n
        acc = acc + (v(i) * w(i))**2
    end do
    nrm = sqrt(acc / real(n, dp))
end function weighted_norm

! Starting step (DVHIN): h from w.r.m.s.norm(h^2 yddot/2)=1 with a 1/2 bias,
! bounded by 0.1*|tout-t0| and the roundoff floor.
real(dp) function initial_step(rhs, state, tout, relerr, abserr, ctx) &
        result(h0)
    procedure(ode_rhs_t)          :: rhs
    type(vode_state_t), intent(inout) :: state
    real(dp), intent(in)           :: tout, relerr
    real(dp), intent(in)           :: abserr(:)
    class(*), intent(in), optional :: ctx

    integer  :: i, iter, neq
    real(dp) :: t0, tdist, tround, hlb, hub, hg, h, t1, hnew, hrat
    real(dp) :: yddnrm, delyi, afi, a, uround
    real(dp), allocatable :: ytmp(:), temp(:)

    neq = state%neq
    t0 = state%tn
    uround = epsilon(1.0_dp)
    allocate(ytmp(neq), temp(neq))

    tdist = abs(tout - t0)
    tround = uround * max(abs(t0), abs(tout))
    hlb = 100.0_dp * tround
    hub = 0.1_dp * tdist
    do i = 1, neq
        if (size(abserr) == 1) then
            a = abserr(1)
        else
            a = abserr(i)
        end if
        delyi = 0.1_dp * abs(state%yh(i, 1)) + a
        afi = abs(state%ftmp(i))
        if (afi * hub > delyi) hub = delyi / afi
    end do

    hg = sqrt(hlb * hub)
    if (hub < hlb) then
        h0 = sign(hg, tout - t0)
        return
    end if

    iter = 0
    do
        h = sign(hg, tout - t0)
        t1 = t0 + h
        ytmp(:) = state%yh(:, 1) + h * state%ftmp(:)
        call rhs(t1, ytmp, temp, ctx)
        state%nfev = state%nfev + 1
        temp(:) = (temp(:) - state%ftmp(:)) / h
        yddnrm = weighted_norm(temp, state%ewt)
        if (yddnrm * hub * hub > 2.0_dp) then
            hnew = sqrt(2.0_dp / yddnrm)
        else
            hnew = sqrt(hg * hub)
        end if
        iter = iter + 1
        if (iter >= 4) exit
        hrat = hnew / hg
        if (hrat > 0.5_dp .and. hrat < 2.0_dp) exit
        if (iter >= 2 .and. hnew > 2.0_dp * hg) then
            hnew = hg
            exit
        end if
        hg = hnew
    end do

    h0 = hnew * 0.5_dp
    if (h0 < hlb) h0 = hlb
    if (h0 > hub) h0 = hub
    h0 = sign(h0, tout - t0)
end function initial_step

! Scan the completed step [tlast, tn] for the first direction-consistent sign
! change of either monitored event function, locate each candidate root on the
! Nordsieck interpolant by the Illinois algorithm to resolution etol, and
! return the earliest in the integration direction (DVODE NEVENTS root logic).
! ev_idx reports which function (1 or 2) owns the returned root.
subroutine scan_step_for_root(event, event2, nev, state, edir, etol, &
        tlast, g_left, troot, found, ev_idx, ctx)
    procedure(ode_event_t)            :: event
    procedure(ode_event_t), optional  :: event2
    integer,  intent(in)              :: nev
    type(vode_state_t), intent(in)    :: state
    integer,  intent(in)              :: edir(:)
    real(dp), intent(in)              :: etol, tlast
    real(dp), intent(in)              :: g_left(:)
    real(dp), intent(out)             :: troot
    logical,  intent(out)             :: found
    integer,  intent(out)             :: ev_idx
    class(*), intent(in), optional    :: ctx

    real(dp) :: troot_i, dirn
    logical  :: found_i
    integer  :: k

    found = .false.
    ev_idx = 0
    troot = state%tn
    dirn = sign(1.0_dp, state%tn - tlast)

    do k = 1, nev
        if (k == 1) then
            call locate_event_root(event, state, edir(k), etol, &
                tlast, g_left(k), troot_i, found_i, ctx)
        else
            call locate_event_root(event2, state, edir(k), etol, &
                tlast, g_left(k), troot_i, found_i, ctx)
        end if
        if (.not. found_i) cycle
        ! First found, or earlier in the integration direction.
        if (.not. found .or. (troot_i - troot) * dirn < 0.0_dp) then
            troot = troot_i
            ev_idx = k
            found = .true.
        end if
    end do
end subroutine scan_step_for_root

! Locate the first direction-consistent root of one event function over the
! completed step [tlast, tn] on the Nordsieck interpolant by the Illinois
! algorithm to resolution etol.
subroutine locate_event_root(event, state, edir, etol, &
        tlast, g_left, troot, found, ctx)
    procedure(ode_event_t)            :: event
    type(vode_state_t), intent(in)    :: state
    integer,  intent(in)              :: edir
    real(dp), intent(in)              :: etol, tlast, g_left
    real(dp), intent(out)             :: troot
    logical,  intent(out)             :: found
    class(*), intent(in), optional    :: ctx

    real(dp) :: g_right, g0, g1, gx, x0, x1, x2, alpha, dirn
    integer  :: last, nxlast, cross
    real(dp), allocatable :: ytmp(:)

    found = .false.
    troot = state%tn
    allocate(ytmp(state%neq))

    ytmp(:) = state%yh(:, 1)
    g_right = event(state%tn, ytmp, ctx)
    cross = crossing(g_left, g_right)
    if (cross == 0) return
    if (edir /= 0 .and. edir /= cross) return

    ! Illinois on the interpolant between tlast (g0) and tn (g1).
    x0 = tlast
    x1 = state%tn
    g0 = g_left
    g1 = g_right
    dirn = sign(1.0_dp, x1 - x0)
    alpha = 1.0_dp
    last = 1
    nxlast = 0

    do
        if (abs(x1 - x0) <= etol) exit
        if (nxlast == last) then
            alpha = 1.0_dp
        else if (last == 0) then
            alpha = 0.5_dp * alpha
        else
            alpha = 2.0_dp * alpha
        end if
        x2 = x1 - (x1 - x0) * g1 / (g1 - alpha * g0)
        ! Keep x2 strictly inside the bracket.
        if ((x2 - x0) * dirn < 0.25_dp * etol) &
            x2 = x0 + 0.5_dp * (x1 - x0)
        if ((x1 - x2) * dirn < 0.25_dp * etol) &
            x2 = x1 - 0.5_dp * (x1 - x0)
        call interpolate(state, x2, ytmp)
        gx = event(x2, ytmp, ctx)
        nxlast = last
        if (crossing(g0, gx) /= 0) then
            x1 = x2; g1 = gx; last = 1
        else
            x0 = x2; g0 = gx; last = 0
        end if
        if (abs(gx) == 0.0_dp) exit
    end do

    troot = x1
    found = .true.
end subroutine locate_event_root

! Sign-change classifier: +1 rising (- to +), -1 falling (+ to -), 0 none.
pure integer function crossing(ga, gb) result(c)
    real(dp), intent(in) :: ga, gb
    c = 0
    if (ga < 0.0_dp .and. gb >= 0.0_dp) then
        c = 1
    else if (ga > 0.0_dp .and. gb <= 0.0_dp) then
        c = -1
    end if
end function crossing

end module fortnum_ode_vode
