module fortnum_ode_rk8pd
    ! Clean-room Prince-Dormand RK8(7)13M stepper with adaptive step
    ! control, bit-faithful to the documented rk8pd stepper driven by
    ! the y_new error-per-step control / an accepted-step evolve advance. KiLCA's background
    ! equilibrium ODE (calc_back.cpp) was tuned against that exact documented behavior,
    ! and Hairer's DOP853 (a different 8(7) error norm) drifts from the recorded
    ! golden near a resonant grid point. This module reproduces the documented behavior
    ! so the background profile matches the golden continuously.
    !
    ! Method and coefficients: Prince and Dormand, "High order embedded
    !   Runge-Kutta formulae", J. Comput. Appl. Math. 7 (1981) 67-75, the
    !   RK8(7)13M pair. The thirteen nodes, the coupling matrix, the eighth-order
    !   propagation weights b, and the seventh-order embedded weights bhat are the
    !   published exact rationals. The error coefficients ec = b - bhat drive the
    !   embedded estimate yerr = h * sum_i ec_i k_i (the the reference rk8pd error vector).
    !   The propagation weights match fortnum_ode_dop853's B* constants; the two
    !   modules differ only in the error estimator and the step controller.
    !
    ! Step control: the reference "standard" controller (the standard error-per-step controller
    !   with a_y = 1, a_dydt = 0, i.e. control_y_new). Per component the scaled
    !   tolerance is D0 = eps_abs + eps_rel*|y|; the failure ratio is
    !   rmax = max_i |yerr_i| / D0_i. the standard step controller uses three bands
    !   keyed on the step's reported order (rk8pd order = 8):
    !     rmax > 1.1  : decrease, h *= max(0.2, S / rmax^(1/order)), step rejected;
    !     rmax < 0.5  : increase, h *= min(5.0, S / rmax^(1/(order+1)));
    !     otherwise   : leave h unchanged.
    !   S = 0.9. This is the documented the reference formula, written clean-room
    !   source.
    !
    ! Re-entrant evolve: rk8pd_state_t carries the current step size and the
    !   first-stage derivative across output points so a caller can step from one
    !   output abscissa to the next while keeping the adaptive schedule, mirroring
    !   a loop over an accepted-step evolve advance. No module-level state; the state and
    !   the stage scratch are caller-owned.
    !
    ! Derivative policy: trace_rule (ad.md sec 1, 4; ode.md sec 1, 4). The
    !   adaptive schedule is data-dependent; a sensitivity differentiates the
    !   frozen accepted-step mesh. This stepper records nothing globally and only
    !   writes its output arguments, so the fortnum_ode trace_rule products attach
    !   the same way they do for the Cash-Karp and DOP853 paths.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t
    implicit none
    private

    public :: rk8pd_step, rk8pd_state_t
    public :: rk8pd_evolve_init, rk8pd_evolve_apply

    ! the standard controller defaults (a_y = 1, a_dydt = 0): the safety factor,
    ! the accept band, and the growth/shrink clamps.
    real(dp), parameter :: SAFETY  = 0.9_dp
    real(dp), parameter :: DEC_BAND = 1.1_dp ! rmax above this => decrease+reject
    real(dp), parameter :: INC_BAND = 0.5_dp ! rmax below this => increase
    real(dp), parameter :: GROW_MAX = 5.0_dp
    real(dp), parameter :: SHRINK_MIN = 0.2_dp
    ! the standard step controller keys its exponents on the step's reported order.
    ! the documented rk8pd stepper reports order 8: decrease uses 1/order, increase uses
    ! 1/(order+1).
    integer, parameter :: STEP_ORDER = 8

    ! Stage nodes c2..c13 (c1 = 0 is implicit). RK8(7)13M, exact rationals.
    real(dp), parameter :: C2  = 1.0_dp / 18.0_dp
    real(dp), parameter :: C3  = 1.0_dp / 12.0_dp
    real(dp), parameter :: C4  = 1.0_dp / 8.0_dp
    real(dp), parameter :: C5  = 5.0_dp / 16.0_dp
    real(dp), parameter :: C6  = 3.0_dp / 8.0_dp
    real(dp), parameter :: C7  = 59.0_dp / 400.0_dp
    real(dp), parameter :: C8  = 93.0_dp / 200.0_dp
    real(dp), parameter :: C9  = 5490023248.0_dp / 9719169821.0_dp
    real(dp), parameter :: C10 = 13.0_dp / 20.0_dp
    real(dp), parameter :: C11 = 1201146811.0_dp / 1299019798.0_dp
    real(dp), parameter :: C12 = 1.0_dp
    real(dp), parameter :: C13 = 1.0_dp

    ! Coupling coefficients a(i,j), j < i. RK8(7)13M, exact rationals.
    real(dp), parameter :: A2_1 = 1.0_dp / 18.0_dp
    real(dp), parameter :: A3_1 = 1.0_dp / 48.0_dp
    real(dp), parameter :: A3_2 = 1.0_dp / 16.0_dp
    real(dp), parameter :: A4_1 = 1.0_dp / 32.0_dp
    real(dp), parameter :: A4_3 = 3.0_dp / 32.0_dp
    real(dp), parameter :: A5_1 = 5.0_dp / 16.0_dp
    real(dp), parameter :: A5_3 = -75.0_dp / 64.0_dp
    real(dp), parameter :: A5_4 = 75.0_dp / 64.0_dp
    real(dp), parameter :: A6_1 = 3.0_dp / 80.0_dp
    real(dp), parameter :: A6_4 = 3.0_dp / 16.0_dp
    real(dp), parameter :: A6_5 = 3.0_dp / 20.0_dp
    real(dp), parameter :: A7_1 = 29443841.0_dp / 614563906.0_dp
    real(dp), parameter :: A7_4 = 77736538.0_dp / 692538347.0_dp
    real(dp), parameter :: A7_5 = -28693883.0_dp / 1125000000.0_dp
    real(dp), parameter :: A7_6 = 23124283.0_dp / 1800000000.0_dp
    real(dp), parameter :: A8_1 = 16016141.0_dp / 946692911.0_dp
    real(dp), parameter :: A8_4 = 61564180.0_dp / 158732637.0_dp
    real(dp), parameter :: A8_5 = 22789713.0_dp / 633445777.0_dp
    real(dp), parameter :: A8_6 = 545815736.0_dp / 2771057229.0_dp
    real(dp), parameter :: A8_7 = -180193667.0_dp / 1043307555.0_dp
    real(dp), parameter :: A9_1 = 39632708.0_dp / 573591083.0_dp
    real(dp), parameter :: A9_4 = -433636366.0_dp / 683701615.0_dp
    real(dp), parameter :: A9_5 = -421739975.0_dp / 2616292301.0_dp
    real(dp), parameter :: A9_6 = 100302831.0_dp / 723423059.0_dp
    real(dp), parameter :: A9_7 = 790204164.0_dp / 839813087.0_dp
    real(dp), parameter :: A9_8 = 800635310.0_dp / 3783071287.0_dp
    real(dp), parameter :: A10_1 = 246121993.0_dp / 1340847787.0_dp
    real(dp), parameter :: A10_4 = -37695042795.0_dp / 15268766246.0_dp
    real(dp), parameter :: A10_5 = -309121744.0_dp / 1061227803.0_dp
    real(dp), parameter :: A10_6 = -12992083.0_dp / 490766935.0_dp
    real(dp), parameter :: A10_7 = 6005943493.0_dp / 2108947869.0_dp
    real(dp), parameter :: A10_8 = 393006217.0_dp / 1396673457.0_dp
    real(dp), parameter :: A10_9 = 123872331.0_dp / 1001029789.0_dp
    real(dp), parameter :: A11_1 = -1028468189.0_dp / 846180014.0_dp
    real(dp), parameter :: A11_4 = 8478235783.0_dp / 508512852.0_dp
    real(dp), parameter :: A11_5 = 1311729495.0_dp / 1432422823.0_dp
    real(dp), parameter :: A11_6 = -10304129995.0_dp / 1701304382.0_dp
    real(dp), parameter :: A11_7 = -48777925059.0_dp / 3047939560.0_dp
    real(dp), parameter :: A11_8 = 15336726248.0_dp / 1032824649.0_dp
    real(dp), parameter :: A11_9 = -45442868181.0_dp / 3398467696.0_dp
    real(dp), parameter :: A11_10 = 3065993473.0_dp / 597172653.0_dp
    real(dp), parameter :: A12_1 = 185892177.0_dp / 718116043.0_dp
    real(dp), parameter :: A12_4 = -3185094517.0_dp / 667107341.0_dp
    real(dp), parameter :: A12_5 = -477755414.0_dp / 1098053517.0_dp
    real(dp), parameter :: A12_6 = -703635378.0_dp / 230739211.0_dp
    real(dp), parameter :: A12_7 = 5731566787.0_dp / 1027545527.0_dp
    real(dp), parameter :: A12_8 = 5232866602.0_dp / 850066563.0_dp
    real(dp), parameter :: A12_9 = -4093664535.0_dp / 808688257.0_dp
    real(dp), parameter :: A12_10 = 3962137247.0_dp / 1805957418.0_dp
    real(dp), parameter :: A12_11 = 65686358.0_dp / 487910083.0_dp
    real(dp), parameter :: A13_1 = 403863854.0_dp / 491063109.0_dp
    real(dp), parameter :: A13_4 = -5068492393.0_dp / 434740067.0_dp
    real(dp), parameter :: A13_5 = -411421997.0_dp / 543043805.0_dp
    real(dp), parameter :: A13_6 = 652783627.0_dp / 914296604.0_dp
    real(dp), parameter :: A13_7 = 11173962825.0_dp / 925320556.0_dp
    real(dp), parameter :: A13_8 = -13158990841.0_dp / 6184727034.0_dp
    real(dp), parameter :: A13_9 = 3936647629.0_dp / 1978049680.0_dp
    real(dp), parameter :: A13_10 = -160528059.0_dp / 685178525.0_dp
    real(dp), parameter :: A13_11 = 248638103.0_dp / 1413531060.0_dp

    ! Eighth-order propagation weights b (stages 2..5 carry zero weight). These
    ! equal fortnum_ode_dop853's B* by construction (same RK8(7)13M tableau).
    real(dp), parameter :: B1  = 14005451.0_dp / 335480064.0_dp
    real(dp), parameter :: B6  = -59238493.0_dp / 1068277825.0_dp
    real(dp), parameter :: B7  = 181606767.0_dp / 758867731.0_dp
    real(dp), parameter :: B8  = 561292985.0_dp / 797845732.0_dp
    real(dp), parameter :: B9  = -1041891430.0_dp / 1371343529.0_dp
    real(dp), parameter :: B10 = 760417239.0_dp / 1151165299.0_dp
    real(dp), parameter :: B11 = 118820643.0_dp / 751138087.0_dp
    real(dp), parameter :: B12 = -528747749.0_dp / 2220607170.0_dp
    real(dp), parameter :: B13 = 1.0_dp / 4.0_dp

    ! Seventh-order embedded weights bhat (Prince-Dormand 1981; the reference rk8pd A[]).
    ! yerr is formed as h*(seventh-order sum - eighth-order sum) at runtime, the
    ! same separate-then-subtract order the reference uses; bhat13 = 0 (no k13 term).
    real(dp), parameter :: BH1  = 13451932.0_dp / 455176623.0_dp
    real(dp), parameter :: BH6  = -808719846.0_dp / 976000145.0_dp
    real(dp), parameter :: BH7  = 1757004468.0_dp / 5645159321.0_dp
    real(dp), parameter :: BH8  = 656045339.0_dp / 265891186.0_dp
    real(dp), parameter :: BH9  = -3867574721.0_dp / 1518517206.0_dp
    real(dp), parameter :: BH10 = 465885868.0_dp / 322736535.0_dp
    real(dp), parameter :: BH11 = 53011238.0_dp / 667516719.0_dp
    real(dp), parameter :: BH12 = 2.0_dp / 45.0_dp

    ! Re-entrant evolve state. Carries the adaptive step size and the cached
    ! first-stage derivative (rk8pd is not FSAL, but f(t,y) at the segment start
    ! is reused across the evolve loop) across output points, plus the stage
    ! scratch so the caller owns all storage. neq guards reallocation.
    type :: rk8pd_state_t
        integer  :: neq = 0
        real(dp) :: h = 0.0_dp ! current step size magnitude
        logical  :: have_dydt = .false.
        real(dp) :: t_dydt = 0.0_dp ! abscissa where dydt0 is valid
        real(dp), allocatable :: dydt0(:)
        real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), allocatable :: k7(:), k8(:), k9(:), k10(:), k11(:), k12(:)
        real(dp), allocatable :: k13(:)
        real(dp), allocatable :: ytmp(:), y8(:), yerr(:)
    end type rk8pd_state_t

contains

    ! Advance one RK8(7)13M step of size h from (t, y). rhs matches ode_rhs_t.
    ! k1..k13 are caller-owned stage-derivative slots (length neq); ytmp is
    ! caller-owned scratch. y8 receives the eighth-order solution, yerr the
    ! embedded error vector (already scaled by h). nfev counts stages evaluated.
    ! have_k1 reuses an externally supplied first-stage derivative; the stepper
    ! only sets that flag when k1 already holds f(t, y).
    subroutine rk8pd_step(rhs, t, y, h, have_k1, k1, k2, k3, k4, k5, k6, &
            k7, k8, k9, k10, k11, k12, k13, ytmp, y8, yerr, &
            nfev, ctx)
        procedure(ode_rhs_t)            :: rhs
        real(dp), intent(in)            :: t
        real(dp), intent(in)            :: y(:)
        real(dp), intent(in)            :: h
        logical,  intent(in)            :: have_k1
        real(dp), intent(inout)         :: k1(:)
        real(dp), intent(out)           :: k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), intent(out)           :: k7(:), k8(:), k9(:), k10(:)
        real(dp), intent(out)           :: k11(:), k12(:), k13(:)
        real(dp), intent(out)           :: ytmp(:)
        real(dp), intent(out)           :: y8(:)
        real(dp), intent(out)           :: yerr(:)
        integer,  intent(inout)         :: nfev
        class(*), intent(in), optional  :: ctx

        real(dp) :: ksum7(size(y)), ksum8(size(y))

        if (.not. have_k1) then
            call rhs(t, y, k1, ctx)
            nfev = nfev + 1
        end if

        ytmp = y + h * (A2_1 * k1)
        call rhs(t + C2 * h, ytmp, k2, ctx)

        ytmp = y + h * (A3_1 * k1 + A3_2 * k2)
        call rhs(t + C3 * h, ytmp, k3, ctx)

        ytmp = y + h * (A4_1 * k1 + A4_3 * k3)
        call rhs(t + C4 * h, ytmp, k4, ctx)

        ytmp = y + h * (A5_1 * k1 + A5_3 * k3 + A5_4 * k4)
        call rhs(t + C5 * h, ytmp, k5, ctx)

        ytmp = y + h * (A6_1 * k1 + A6_4 * k4 + A6_5 * k5)
        call rhs(t + C6 * h, ytmp, k6, ctx)

        ytmp = y + h * (A7_1 * k1 + A7_4 * k4 + A7_5 * k5 + A7_6 * k6)
        call rhs(t + C7 * h, ytmp, k7, ctx)

        ytmp = y + h * (A8_1 * k1 + A8_4 * k4 + A8_5 * k5 + A8_6 * k6 &
            + A8_7 * k7)
        call rhs(t + C8 * h, ytmp, k8, ctx)

        ytmp = y + h * (A9_1 * k1 + A9_4 * k4 + A9_5 * k5 + A9_6 * k6 &
            + A9_7 * k7 + A9_8 * k8)
        call rhs(t + C9 * h, ytmp, k9, ctx)

        ytmp = y + h * (A10_1 * k1 + A10_4 * k4 + A10_5 * k5 + A10_6 * k6 &
            + A10_7 * k7 + A10_8 * k8 + A10_9 * k9)
        call rhs(t + C10 * h, ytmp, k10, ctx)

        ytmp = y + h * (A11_1 * k1 + A11_4 * k4 + A11_5 * k5 + A11_6 * k6 &
            + A11_7 * k7 + A11_8 * k8 + A11_9 * k9 + A11_10 * k10)
        call rhs(t + C11 * h, ytmp, k11, ctx)

        ytmp = y + h * (A12_1 * k1 + A12_4 * k4 + A12_5 * k5 + A12_6 * k6 &
            + A12_7 * k7 + A12_8 * k8 + A12_9 * k9 + A12_10 * k10 &
            + A12_11 * k11)
        call rhs(t + C12 * h, ytmp, k12, ctx)

        ytmp = y + h * (A13_1 * k1 + A13_4 * k4 + A13_5 * k5 + A13_6 * k6 &
            + A13_7 * k7 + A13_8 * k8 + A13_9 * k9 + A13_10 * k10 &
            + A13_11 * k11)
        call rhs(t + C13 * h, ytmp, k13, ctx)

        nfev = nfev + 12

        ! Stages 2 and 3 carry zero propagation weight; k13 (the doubled final
        ! node) contributes only to the eighth-order sum. Form the eighth- and
        ! seventh-order weighted sums separately and take yerr from their
        ! difference, exactly as the reference rk8pd does (yerr = h*(ksum7 - ksum8)).
        ! Precombining the per-stage difference (b_i - bhat_i)*k_i instead rounds
        ! differently and drifts rmax, hence the step schedule, from the reference's.
        ksum8 = B1 * k1 + B6 * k6 + B7 * k7 + B8 * k8 + B9 * k9 &
            + B10 * k10 + B11 * k11 + B12 * k12 + B13 * k13
        ksum7 = BH1 * k1 + BH6 * k6 + BH7 * k7 + BH8 * k8 + BH9 * k9 &
            + BH10 * k10 + BH11 * k11 + BH12 * k12

        y8 = y + h * ksum8
        yerr = h * (ksum7 - ksum8)
    end subroutine rk8pd_step

    ! Initialise a re-entrant evolve state for neq equations with starting step
    ! magnitude h0 (must be positive). Clears the cached first-stage derivative.
    subroutine rk8pd_evolve_init(state, neq, h0, status)
        type(rk8pd_state_t),    intent(inout) :: state
        integer,                intent(in)    :: neq
        real(dp),               intent(in)    :: h0
        type(fortnum_status_t), intent(out)   :: status
        call status_set(status, FORTNUM_OK, "")
        if (neq < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rk8pd_evolve_init: neq < 1")
            return
        end if
        if (h0 <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rk8pd_evolve_init: h0 must be positive")
            return
        end if
        call ensure_state(state, neq)
        state%h = h0
        state%have_dydt = .false.
        state%t_dydt = 0.0_dp
    end subroutine rk8pd_evolve_init

    ! Evolve y from t to t1 with adaptive error-per-step control, mirroring a loop
    ! over an accepted-step evolve advance: the step size carried in state%h is reused
    ! and re-adapted continuously, and the final step is clipped to land exactly
    ! on t1 (the reference caps h at t1 - t without recording the clip as a permanent
    ! shrink). On return t == t1, y holds the solution there, state%h carries the
    ! step proposed for continuing past t1, and the cached first-stage derivative
    ! is valid at t1. nfev accumulates RHS evaluations.
    !
    ! single_step (default .false.): when .true. the routine returns after one
    ! accepted step exactly as a single an accepted-step evolve advance call does, so a
    ! caller that records the adaptive mesh (one row per accepted step) reproduces
    ! the adaptive evolve loop point for point. t then lands either on the accepted
    ! step endpoint or on t1 when that step was clipped to the interval end.
    subroutine rk8pd_evolve_apply(rhs, state, t, t1, y, eps_abs, eps_rel, &
            max_steps, nfev, status, ctx, single_step)
        procedure(ode_rhs_t)            :: rhs
        type(rk8pd_state_t),    intent(inout) :: state
        real(dp),               intent(inout) :: t
        real(dp),               intent(in)    :: t1
        real(dp),               intent(inout) :: y(:)
        real(dp),               intent(in)    :: eps_abs, eps_rel
        integer,                intent(in)    :: max_steps
        integer,                intent(inout) :: nfev
        type(fortnum_status_t), intent(out)   :: status
        class(*), intent(in), optional        :: ctx
        logical,  intent(in), optional        :: single_step

        integer  :: neq, nstep
        real(dp) :: dir, h, h_used, rmax, h_clip
        logical  :: final_clip, one_step

        call status_set(status, FORTNUM_OK, "")
        neq = size(y)
        if (state%neq /= neq .or. .not. allocated(state%k1)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rk8pd_evolve_apply: state not initialised for this neq")
            return
        end if
        if (eps_abs <= 0.0_dp .or. eps_rel < 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rk8pd_evolve_apply: invalid tolerances")
            return
        end if

        one_step = .false.
        if (present(single_step)) one_step = single_step

        dir = sign(1.0_dp, t1 - t)
        if (t1 == t) return

        nstep = 0
        do
            if (abs(t1 - t) <= 0.0_dp) exit
            if (nstep >= max_steps) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "rk8pd_evolve_apply: exceeded max_steps")
                return
            end if

            ! Cache f(t, y) for the step's first stage. the reference recomputes it at the
            ! start of each evolve_apply; here it is reused if already valid at t.
            if (.not. (state%have_dydt .and. state%t_dydt == t)) then
                call rhs(t, y, state%dydt0, ctx)
                nfev = nfev + 1
                state%have_dydt = .true.
                state%t_dydt = t
            end if

            ! Clip the trial step to the interval, mirroring the h0 cap at the
            ! remaining span. The clip is local: state%h keeps the unclipped
            ! proposal so the continuous schedule survives across output points.
            h = state%h
            h_clip = abs(t1 - t)
            final_clip = h >= h_clip
            if (final_clip) h = h_clip

            state%k1 = state%dydt0
            call rk8pd_step(rhs, t, y, dir * h, .true., &
                state%k1, state%k2, state%k3, state%k4, state%k5, state%k6, &
                state%k7, state%k8, state%k9, state%k10, state%k11, &
                state%k12, state%k13, state%ytmp, state%y8, state%yerr, &
                nfev, ctx)

            rmax = failure_ratio(y, state%y8, state%yerr, eps_abs, eps_rel)

            h_used = h
            if (rmax > DEC_BAND) then
                ! Reject: shrink and retry from the same t. The carried step
                ! state%h shrinks by the reference's decrease rule; the clip is not
                ! recorded so the continuous schedule is unaffected.
                state%h = adjust_step(h, rmax, decrease=.true.)
                state%have_dydt = .true. ! dydt0 still valid at t
                nstep = nstep + 1
                cycle
            end if

            ! Accept. Advance t and y; recompute f at the new point next loop.
            if (final_clip) then
                t = t1
            else
                t = t + dir * h_used
            end if
            y = state%y8
            state%have_dydt = .false.

            ! Carry the next step. the reference leaves h unchanged in the dead band
            ! (0.5 <= rmax <= 1.1) and grows it when rmax < 0.5. The clip on a
            ! final step is local, so grow/keep is applied to the unclipped
            ! carried value, not to the clipped h_used.
            if (rmax < INC_BAND) then
                state%h = adjust_step(state%h, rmax, decrease=.false.)
            end if

            nstep = nstep + 1
            if (final_clip) exit
            if (one_step) exit
        end do
    end subroutine rk8pd_evolve_apply

    ! --- internals ---

    ! the y_new error-per-step control failure ratio: max over components of |yerr_i| / D0_i,
    ! D0_i = eps_abs + eps_rel*|y_i|. y is the pre-step state; the reference scales by the
    ! input y, not the proposed y8.
    real(dp) function failure_ratio(y, y8, yerr, eps_abs, eps_rel) result(r)
        real(dp), intent(in) :: y(:), y8(:), yerr(:)
        real(dp), intent(in) :: eps_abs, eps_rel
        real(dp) :: d0, ratio
        integer  :: i, n
        n = size(y)
        r = 0.0_dp
        do i = 1, n
            d0 = eps_abs + eps_rel * abs(y(i))
            if (d0 > 0.0_dp) then
                ratio = abs(yerr(i)) / d0
                if (ratio > r) r = ratio
            end if
        end do
    end function failure_ratio

    ! the standard step controller step adjustment. Decrease (rejected step):
    !   h_new = h * max(0.2, S / rmax^(1/order)). Increase (accepted, rmax<0.5):
    !   h_new = h * min(5.0, S / rmax^(1/(order+1))), then floored at 1.0 so the
    !   safety factor S<1 never turns an accepted step into a shrink (the reference v1
    !   std_control_hadjust: "don't allow any decrease caused by S<1"). Omitting
    !   that floor lets rmax just under 0.5 nudge the carried step down by a hair,
    !   which drifts the step schedule from the reference on sensitive ODEs. order =
    !   STEP_ORDER.
    real(dp) function adjust_step(h, rmax, decrease) result(h_new)
        real(dp), intent(in) :: h, rmax
        logical,  intent(in) :: decrease
        real(dp) :: fac
        if (decrease) then
            fac = SAFETY / rmax**(1.0_dp / real(STEP_ORDER, dp))
            if (fac < SHRINK_MIN) fac = SHRINK_MIN
        else
            if (rmax <= 0.0_dp) then
                h_new = GROW_MAX * h
                return
            end if
            fac = SAFETY / rmax**(1.0_dp / real(STEP_ORDER + 1, dp))
            if (fac > GROW_MAX) fac = GROW_MAX
            if (fac < 1.0_dp) fac = 1.0_dp
        end if
        h_new = fac * h
    end function adjust_step

    ! Allocate the thirteen stage slots, the cached first-stage derivative, and
    ! the scratch states once when the size changes.
    subroutine ensure_state(state, neq)
        type(rk8pd_state_t), intent(inout) :: state
        integer,             intent(in)    :: neq
        if (state%neq == neq .and. allocated(state%k1)) return
        state%neq = neq
        call realloc(state%dydt0, neq)
        call realloc(state%k1, neq)
        call realloc(state%k2, neq)
        call realloc(state%k3, neq)
        call realloc(state%k4, neq)
        call realloc(state%k5, neq)
        call realloc(state%k6, neq)
        call realloc(state%k7, neq)
        call realloc(state%k8, neq)
        call realloc(state%k9, neq)
        call realloc(state%k10, neq)
        call realloc(state%k11, neq)
        call realloc(state%k12, neq)
        call realloc(state%k13, neq)
        call realloc(state%ytmp, neq)
        call realloc(state%y8, neq)
        call realloc(state%yerr, neq)
    end subroutine ensure_state

    subroutine realloc(a, n)
        real(dp), allocatable, intent(inout) :: a(:)
        integer,               intent(in)    :: n
        if (allocated(a)) then
            if (size(a) == n) return
            deallocate(a)
        end if
        allocate(a(n))
    end subroutine realloc

end module fortnum_ode_rk8pd
