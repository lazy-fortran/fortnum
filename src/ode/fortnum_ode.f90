module fortnum_ode
    ! Adaptive ODE integration: Cash-Karp RK5(4) with PI step-size control.
    ! API and numerics fixed by docs/design/ode.md (ADR #14, M2.1).
    !
    ! Derivative policy: trace_rule (ad.md §1, §4; ode.md §1, §4).
    !   Adaptive step control makes the step schedule data-dependent. A
    !   sensitivity differentiates the recorded accepted-step mesh with that
    !   schedule held fixed: the variational system rides the frozen trace
    !   (solution%t, solution%y, solution%h). The integrator records that mesh
    !   and owns no global state, so the #40 forward/reverse products attach
    !   additively. ode_integrate_jvp integrates the variational equation
    !   d/dt(dy/dp.v) = J_f(t,y)(dy/dp.v) + df/dp.vp over the frozen schedule;
    !   ode_integrate_vjp is its discrete adjoint, the exact transpose walked
    !   backward over the same trace (verified by the dot-product identity to
    !   machine precision). Both take the caller-supplied variational RHS
    !   (ode_var_rhs_t): the forward one applies J_f, the adjoint one applies
    !   J_f^T. The Jacobian-w.r.t.-y action is supplied this way; the
    !   PARAMETER-gradient quadrature df/dp^T lam (so the VJP returns the
    !   parameter adjoint, not only the y0 adjoint) is the documented next step
    !   and is not accumulated here. HVP is deferred: ode has no scalar primal
    !   output, so a Hessian-vector product needs a caller-defined scalar loss
    !   and a second tangent pass; left to the loss-aware layer.
    !   Active: y0 (via the JVP seed s0), ctx parameters (via var_rhs), the
    !   recorded trace. Inactive: rtol, atol, h0, hmin, hmax, max_steps,
    !   event_direction, and the counts/err/status.
    !
    ! Method: Cash and Karp, ACM TOMS 16 (1990) 201-222 (stepper in
    !   fortnum_ode_cash_karp). PI controller and starting-step estimate:
    !   Hairer, Norsett, Wanner, "Solving Ordinary Differential Equations I",
    !   2nd ed., II.4. No module-level save; the workspace is caller-owned.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode_cash_karp, only: cash_karp_step
    implicit none
    private

    public :: ode_rhs_t, ode_event_t, ode_var_rhs_t
    public :: ode_integrate, ode_solve
    public :: ode_integrate_jvp, ode_integrate_vjp

    integer, parameter, public :: ODE_EVENT_RISING  =  1
    integer, parameter, public :: ODE_EVENT_FALLING = -1
    integer, parameter, public :: ODE_EVENT_ANY     =  0

    ! PI controller constants (Hairer I, II.4). alpha/beta are the error
    ! exponents for an order-5 method; fac is the safety factor; the step
    ! growth ratio is clamped to [fac_min, fac_max].
    real(dp), parameter :: PI_ALPHA = 0.7_dp / 5.0_dp
    real(dp), parameter :: PI_BETA  = 0.4_dp / 5.0_dp
    real(dp), parameter :: SAFETY   = 0.9_dp
    real(dp), parameter :: FAC_MIN  = 0.2_dp
    real(dp), parameter :: FAC_MAX  = 5.0_dp
    real(dp), parameter :: TRACE_GROWTH = 2.0_dp

    ! Cash-Karp tableau constants for the augmented and adjoint steppers.
    ! Must match fortnum_ode_cash_karp exactly; that module keeps them private.
    real(dp), parameter :: CK_C2 = 0.2_dp,   CK_C3 = 0.3_dp
    real(dp), parameter :: CK_C4 = 0.6_dp,   CK_C5 = 1.0_dp,  CK_C6 = 0.875_dp
    real(dp), parameter :: CK_A21 = 0.2_dp
    real(dp), parameter :: CK_A31 = 3.0_dp/40.0_dp,  CK_A32 = 9.0_dp/40.0_dp
    real(dp), parameter :: CK_A41 = 0.3_dp,  CK_A42 = -0.9_dp, CK_A43 = 1.2_dp
    real(dp), parameter :: CK_A51 = -11.0_dp/54.0_dp, CK_A52 = 2.5_dp
    real(dp), parameter :: CK_A53 = -70.0_dp/27.0_dp, CK_A54 = 35.0_dp/27.0_dp
    real(dp), parameter :: CK_A61 = 1631.0_dp/55296.0_dp
    real(dp), parameter :: CK_A62 = 175.0_dp/512.0_dp
    real(dp), parameter :: CK_A63 = 575.0_dp/13824.0_dp
    real(dp), parameter :: CK_A64 = 44275.0_dp/110592.0_dp
    real(dp), parameter :: CK_A65 = 253.0_dp/4096.0_dp
    real(dp), parameter :: CK_B5_1 = 37.0_dp/378.0_dp
    real(dp), parameter :: CK_B5_3 = 250.0_dp/621.0_dp
    real(dp), parameter :: CK_B5_4 = 125.0_dp/594.0_dp
    real(dp), parameter :: CK_B5_6 = 512.0_dp/1771.0_dp

    abstract interface
        subroutine ode_rhs_t(t, y, dydt, ctx)
            import :: dp
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
        end subroutine ode_rhs_t
    end interface

    abstract interface
        ! Variational (tangent) RHS for forward sensitivity (#40, trace_rule).
        ! Returns the directional derivative of ode_rhs_t along the augmented
        ! tangent: dsdt = J_f(t,y) s + df/dp . vp, where J_f is the RHS Jacobian
        ! w.r.t. y and the df/dp.vp term carries the parameter-direction part of
        ! the seed. The caller supplies the closed form (or its own AD) so the
        ! tangent rides the same trace the primal recorded. y is the primal state
        ! on the frozen trace; s is the current sensitivity dy/dp.v.
        subroutine ode_var_rhs_t(t, y, s, dsdt, ctx)
            import :: dp
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(in)  :: s(:)
            real(dp), intent(out) :: dsdt(:)
            class(*), intent(in), optional :: ctx
        end subroutine ode_var_rhs_t
    end interface

    abstract interface
        function ode_event_t(t, y, ctx) result(g)
            import :: dp
            real(dp), intent(in) :: t
            real(dp), intent(in) :: y(:)
            class(*), intent(in), optional :: ctx
            real(dp) :: g
        end function ode_event_t
    end interface

    type, public :: ode_problem_t
        procedure(ode_rhs_t), pointer, nopass :: rhs => null()
        real(dp)              :: t0   = 0.0_dp
        real(dp)              :: t1   = 0.0_dp
        real(dp), allocatable :: y0(:)
        real(dp)              :: rtol = 1.0e-6_dp
        real(dp)              :: atol = 1.0e-9_dp
        real(dp)              :: h0   = 0.0_dp
        real(dp)              :: hmin = 0.0_dp
        real(dp)              :: hmax = 0.0_dp
        integer               :: max_steps = 100000
        procedure(ode_event_t), pointer, nopass :: event => null()
        integer               :: event_direction = ODE_EVENT_ANY
        logical               :: terminal_event  = .true.
        real(dp)              :: event_tol = 1.0e-12_dp
    end type ode_problem_t

    type, public :: ode_workspace_t
        integer               :: neq = 0
        real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), allocatable :: ytmp(:), y5(:), y4(:), yerr(:)
        ! Higher stages and error vectors used by the RK8(7)13M integrator
        ! (fortnum_ode_dop853). The Cash-Karp paths ignore these slots.
        real(dp), allocatable :: k7(:), k8(:), k9(:), k10(:), k11(:), k12(:)
        real(dp), allocatable :: y8(:), err5(:), err3(:)
    end type ode_workspace_t

    type, public :: ode_solution_t
        integer               :: nsteps     = 0
        integer               :: nrejected  = 0
        integer               :: nfev       = 0
        real(dp), allocatable :: t(:)
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: h(:)
        real(dp), allocatable :: err(:)
        logical               :: event_found = .false.
        real(dp)              :: t_event = 0.0_dp
        real(dp), allocatable :: y_event(:)
        type(fortnum_status_t) :: status
    end type ode_solution_t

contains

    ! Integrate problem%rhs from t0 to t1 with adaptive Cash-Karp RK5(4).
    ! Records the accepted-step mesh into solution. Event detection lands in
    ! #16; this primal ignores problem%event.
    subroutine ode_integrate(problem, workspace, solution, status)
        type(ode_problem_t),    intent(in)    :: problem
        type(ode_workspace_t),  intent(inout) :: workspace
        type(ode_solution_t),   intent(inout) :: solution
        type(fortnum_status_t), intent(out)   :: status

        integer  :: neq, nstep, cap
        real(dp) :: t, dir, span, h, hmax, hmin
        real(dp) :: err_norm, err_prev, fac, h_new
        logical  :: accepted, first_step, after_reject, final_step

        call status_set(status, FORTNUM_OK, "")
        solution%nsteps = 0
        solution%nrejected = 0
        solution%nfev = 0
        solution%event_found = .false.
        solution%t_event = 0.0_dp

        if (.not. validate_problem(problem, status)) then
            solution%status = status
            return
        end if

        neq = size(problem%y0)
        call ensure_workspace(workspace, neq)

        dir  = sign(1.0_dp, problem%t1 - problem%t0)
        span = abs(problem%t1 - problem%t0)
        if (problem%hmax > 0.0_dp) then
            hmax = min(problem%hmax, span)
        else
            hmax = span
        end if
        hmin = problem%hmin

        cap = 64
        call alloc_trace(solution, neq, cap)
        solution%t(1)   = problem%t0
        solution%y(:,1) = problem%y0
        nstep = 0

        ! t1 == t0: nothing to integrate, single recorded point.
        if (span <= 0.0_dp) then
            call trim_trace(solution, 0, neq)
            solution%status = status
            return
        end if

        h = initial_step(problem, workspace, hmax)
        err_prev = 1.0_dp
        first_step = .true.
        after_reject = .false.
        t = problem%t0

        do
            if (nstep >= problem%max_steps) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ode_integrate: exceeded max_steps")
                exit
            end if

            ! Clip the final step to land exactly on t1.
            final_step = abs(problem%t1 - t) <= h
            if (final_step) h = abs(problem%t1 - t)
            if (h <= 0.0_dp) exit

            call cash_karp_step(problem%rhs, t, &
                solution%y(:,nstep+1), dir * h, .false., &
                workspace%k1, workspace%k2, workspace%k3, workspace%k4, &
                workspace%k5, workspace%k6, workspace%ytmp, &
                workspace%y5, workspace%yerr, solution%nfev)

            err_norm = error_norm(solution%y(:,nstep+1), workspace%y5, &
                                  workspace%yerr, problem%rtol, problem%atol)
            accepted = err_norm <= 1.0_dp

            if (accepted) then
                ! Snap the clipped final step onto t1 exactly; avoids drift in
                ! the recorded endpoint from accumulated h additions.
                if (final_step) then
                    t = problem%t1
                else
                    t = t + dir * h
                end if
                nstep = nstep + 1
                if (nstep + 1 > cap) then
                    cap = max(cap + 1, nint(cap * TRACE_GROWTH))
                    call grow_trace(solution, neq, cap)
                end if
                solution%t(nstep+1) = t
                solution%y(:,nstep+1) = workspace%y5
                solution%h(nstep) = dir * h
                solution%err(nstep) = err_norm
                if (final_step) exit
            else
                solution%nrejected = solution%nrejected + 1
            end if

            ! PI control on accepted history; I-rule on first step or after a
            ! rejection (Hairer I, II.4).
            fac = control_factor(err_norm, err_prev, first_step .or. after_reject)
            h_new = h * fac
            if (h_new > hmax) h_new = hmax

            if (accepted) then
                err_prev = max(err_norm, 1.0e-10_dp)
                first_step = .false.
                after_reject = .false.
            else
                after_reject = .true.
            end if

            if (hmin > 0.0_dp .and. h_new < hmin .and. .not. accepted) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ode_integrate: step forced below hmin")
                exit
            end if
            h = h_new
        end do

        solution%nsteps = nstep
        call trim_trace(solution, nstep, neq)
        solution%status = status
    end subroutine ode_integrate

    ! Flat call: build a problem, integrate, hand back the trace. Events are
    ! not exposed (ode.md §5).
    subroutine ode_solve(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
        procedure(ode_rhs_t)               :: rhs
        real(dp),               intent(in) :: t0, t1
        real(dp),               intent(in) :: y0(:)
        real(dp), allocatable, intent(out) :: t_out(:)
        real(dp), allocatable, intent(out) :: y_out(:,:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional     :: rtol, atol

        type(ode_problem_t)   :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t)  :: solution
        integer :: npts

        problem%rhs => rhs
        problem%t0 = t0
        problem%t1 = t1
        problem%y0 = y0
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol

        call ode_integrate(problem, workspace, solution, status)

        npts = solution%nsteps + 1
        if (allocated(solution%t)) then
            t_out = solution%t(1:npts)
            y_out = solution%y(:, 1:npts)
        else
            allocate(t_out(0))
            allocate(y_out(size(y0), 0))
        end if
    end subroutine ode_solve

    ! Forward sensitivity (JVP), trace_rule (ad.md sec 1, 4). Integrates the
    ! variational equation d/dt(dy/dp.v) = J_f(t,y)(dy/dp.v) + df/dp.vp alongside
    ! the primal over the FROZEN accepted-step schedule recorded in solution
    ! (solution%t, solution%h). Returns s1 = dy(t1)/dp . v, the forward product
    ! J(y0,p->y(t1)) applied to the perturbation direction encoded by the seed
    ! s0 (the y0-part: s0 = dy0/dp.v) and by var_rhs (the parameter part).
    !
    ! Re-runs the primal stepper in lockstep with the tangent so the stage
    ! states match what the primal saw; the step sizes are taken verbatim from
    ! the recorded trace and never re-adapted. Caller must have called
    ! ode_integrate first to fill solution. Active: y0 (via s0), parameters (via
    ! var_rhs). Inactive: tolerances, step controls, status (ad.md sec 3).
    subroutine ode_integrate_jvp(problem, var_rhs, s0, solution, s1, status)
        type(ode_problem_t),    intent(in)  :: problem
        procedure(ode_var_rhs_t)            :: var_rhs
        real(dp),               intent(in)  :: s0(:)
        type(ode_solution_t),   intent(in)  :: solution
        real(dp), allocatable,  intent(out) :: s1(:)
        type(fortnum_status_t), intent(out) :: status

        integer  :: neq, i, nstep
        real(dp) :: t, h
        real(dp), allocatable :: y(:), s(:)

        call status_set(status, FORTNUM_OK, "")
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_jvp: rhs not associated")
            return
        end if
        if (.not. allocated(solution%t) .or. .not. allocated(solution%y)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_jvp: empty trace (run ode_integrate first)")
            return
        end if

        neq = size(solution%y, 1)
        if (size(s0) /= neq) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_jvp: s0 size mismatch")
            return
        end if

        allocate(s1(neq), y(neq), s(neq))
        s = s0

        nstep = solution%nsteps
        do i = 1, nstep
            t = solution%t(i)
            h = solution%h(i)
            y = solution%y(:, i)
            call augmented_step(problem%rhs, var_rhs, t, y, s, h, neq)
        end do

        s1 = s
    end subroutine ode_integrate_jvp

    ! Reverse sensitivity (VJP) over the frozen trace. Discrete adjoint: the
    ! forward product is the composition of per-step linear maps S = M_n...M_1,
    ! so J^T u = M_1^T...M_n^T u, walking the recorded trace backward. Each
    ! M_i^T is applied matrix-free by the variational RHS, transposed at the
    ! stage level. Returns jtu = dy(t1)/dy0 ^T . u in the y0-input space.
    !
    ! var_rhs_adj(t,y,lam,out) must apply J_f(t,y)^T lam (the transpose of the
    ! tangent RHS w.r.t. s). The parameter-gradient quadrature (df/dp^T lam) is
    ! NOT accumulated here: this returns the y0-adjoint only. The parameter
    ! adjoint is the documented next step (see header / report).
    subroutine ode_integrate_vjp(problem, var_rhs_adj, u, solution, jtu, status)
        type(ode_problem_t),    intent(in)  :: problem
        procedure(ode_var_rhs_t)            :: var_rhs_adj
        real(dp),               intent(in)  :: u(:)
        type(ode_solution_t),   intent(in)  :: solution
        real(dp), allocatable,  intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status

        integer  :: neq, i, nstep
        real(dp) :: t, h
        real(dp), allocatable :: y(:), lam(:)

        call status_set(status, FORTNUM_OK, "")
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_vjp: rhs not associated")
            return
        end if
        if (.not. allocated(solution%t) .or. .not. allocated(solution%y)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_vjp: empty trace (run ode_integrate first)")
            return
        end if

        neq = size(solution%y, 1)
        if (size(u) /= neq) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_vjp: u size mismatch")
            return
        end if

        allocate(jtu(neq), y(neq), lam(neq))
        lam = u

        nstep = solution%nsteps
        do i = nstep, 1, -1
            t = solution%t(i)
            h = solution%h(i)
            y = solution%y(:, i)
            call augmented_step_adj(problem%rhs, var_rhs_adj, t, y, lam, h, neq)
        end do

        jtu = lam
    end subroutine ode_integrate_vjp

    ! One Cash-Karp step on the augmented state (y, s): the primal stages drive
    ! the tangent stages so they share the same internal states. h is taken from
    ! the frozen trace (signed). On return s holds the propagated sensitivity;
    ! y is advanced too but discarded by the caller (the trace already has it).
    subroutine augmented_step(rhs, var_rhs, t, y, s, h, neq)
        procedure(ode_rhs_t)     :: rhs
        procedure(ode_var_rhs_t) :: var_rhs
        real(dp), intent(in)     :: t
        real(dp), intent(inout)  :: y(:)
        real(dp), intent(inout)  :: s(:)
        real(dp), intent(in)     :: h
        integer,  intent(in)     :: neq

        real(dp) :: yk1(neq), yk2(neq), yk3(neq), yk4(neq), yk5(neq), yk6(neq)
        real(dp) :: sk1(neq), sk2(neq), sk3(neq), sk4(neq), sk5(neq), sk6(neq)
        real(dp) :: yt(neq), st(neq)

        call ck_stage(rhs, var_rhs, t, y, s, yk1, sk1)

        yt = y + h * (CK_A21 * yk1)
        st = s + h * (CK_A21 * sk1)
        call ck_stage(rhs, var_rhs, t + CK_C2 * h, yt, st, yk2, sk2)

        yt = y + h * (CK_A31 * yk1 + CK_A32 * yk2)
        st = s + h * (CK_A31 * sk1 + CK_A32 * sk2)
        call ck_stage(rhs, var_rhs, t + CK_C3 * h, yt, st, yk3, sk3)

        yt = y + h * (CK_A41 * yk1 + CK_A42 * yk2 + CK_A43 * yk3)
        st = s + h * (CK_A41 * sk1 + CK_A42 * sk2 + CK_A43 * sk3)
        call ck_stage(rhs, var_rhs, t + CK_C4 * h, yt, st, yk4, sk4)

        yt = y + h * (CK_A51 * yk1 + CK_A52 * yk2 + CK_A53 * yk3 + CK_A54 * yk4)
        st = s + h * (CK_A51 * sk1 + CK_A52 * sk2 + CK_A53 * sk3 + CK_A54 * sk4)
        call ck_stage(rhs, var_rhs, t + CK_C5 * h, yt, st, yk5, sk5)

        yt = y + h * (CK_A61 * yk1 + CK_A62 * yk2 + CK_A63 * yk3 &
                      + CK_A64 * yk4 + CK_A65 * yk5)
        st = s + h * (CK_A61 * sk1 + CK_A62 * sk2 + CK_A63 * sk3 &
                      + CK_A64 * sk4 + CK_A65 * sk5)
        call ck_stage(rhs, var_rhs, t + CK_C6 * h, yt, st, yk6, sk6)

        y = y + h * (CK_B5_1 * yk1 + CK_B5_3 * yk3 + CK_B5_4 * yk4 + CK_B5_6 * yk6)
        s = s + h * (CK_B5_1 * sk1 + CK_B5_3 * sk3 + CK_B5_4 * sk4 + CK_B5_6 * sk6)
    end subroutine augmented_step

    ! Evaluate primal and tangent stage derivatives at (t, y, s).
    subroutine ck_stage(rhs, var_rhs, t, y, s, ydot, sdot)
        procedure(ode_rhs_t)     :: rhs
        procedure(ode_var_rhs_t) :: var_rhs
        real(dp), intent(in)  :: t, y(:), s(:)
        real(dp), intent(out) :: ydot(:), sdot(:)
        call rhs(t, y, ydot)
        call var_rhs(t, y, s, sdot)
    end subroutine ck_stage

    ! Discrete-adjoint step: the transpose of augmented_step's linear s-map.
    ! The forward s-update is s' = s + h sum b_i sk_i with sk_i a linear chain
    ! in the earlier sk_j; the adjoint propagates lam through the transposed
    ! chain. Implemented matrix-free: each stage's contribution to lam is
    ! var_rhs_adj(t_i, y_i, .)^T applied to the accumulated stage adjoint.
    subroutine augmented_step_adj(rhs, var_rhs_adj, t, y, lam, h, neq)
        procedure(ode_rhs_t)     :: rhs
        procedure(ode_var_rhs_t) :: var_rhs_adj
        real(dp), intent(in)     :: t
        real(dp), intent(in)     :: y(:)
        real(dp), intent(inout)  :: lam(:)
        real(dp), intent(in)     :: h
        integer,  intent(in)     :: neq

        real(dp) :: yk1(neq), yk2(neq), yk3(neq), yk4(neq), yk5(neq), yk6(neq)
        real(dp) :: yt(neq)
        real(dp) :: tn(6)
        real(dp) :: yn(neq, 6)          ! primal stage states
        real(dp) :: skbar(neq, 6)       ! cotangent of stage outputs sk_i
        real(dp) :: inbar(neq)          ! cotangent of stage input in_i
        real(dp) :: sbar(neq)
        real(dp) :: a(6, 6), b(6)
        integer  :: i, j

        ! Recompute the primal stage STATES so the transposed Jacobian g_i^T is
        ! evaluated at the same nodes the forward tangent used.
        call rhs(t, y, yk1)
        yn(:, 1) = y;  tn(1) = t
        yt = y + h * (CK_A21 * yk1)
        call rhs(t + CK_C2 * h, yt, yk2); yn(:, 2) = yt; tn(2) = t + CK_C2 * h
        yt = y + h * (CK_A31 * yk1 + CK_A32 * yk2)
        call rhs(t + CK_C3 * h, yt, yk3); yn(:, 3) = yt; tn(3) = t + CK_C3 * h
        yt = y + h * (CK_A41 * yk1 + CK_A42 * yk2 + CK_A43 * yk3)
        call rhs(t + CK_C4 * h, yt, yk4); yn(:, 4) = yt; tn(4) = t + CK_C4 * h
        yt = y + h * (CK_A51 * yk1 + CK_A52 * yk2 + CK_A53 * yk3 + CK_A54 * yk4)
        call rhs(t + CK_C5 * h, yt, yk5); yn(:, 5) = yt; tn(5) = t + CK_C5 * h
        yt = y + h * (CK_A61 * yk1 + CK_A62 * yk2 + CK_A63 * yk3 &
                      + CK_A64 * yk4 + CK_A65 * yk5)
        call rhs(t + CK_C6 * h, yt, yk6); yn(:, 6) = yt; tn(6) = t + CK_C6 * h

        ! Tableau coupling/weights as a matrix, zeros for unused entries.
        a = 0.0_dp; b = 0.0_dp
        a(2,1) = CK_A21
        a(3,1) = CK_A31; a(3,2) = CK_A32
        a(4,1) = CK_A41; a(4,2) = CK_A42; a(4,3) = CK_A43
        a(5,1) = CK_A51; a(5,2) = CK_A52; a(5,3) = CK_A53; a(5,4) = CK_A54
        a(6,1) = CK_A61; a(6,2) = CK_A62; a(6,3) = CK_A63; a(6,4) = CK_A64
        a(6,5) = CK_A65
        b(1) = CK_B5_1; b(3) = CK_B5_3; b(4) = CK_B5_4; b(6) = CK_B5_6

        ! Forward s-map:  in_i = s + h sum_{j<i} a_ij sk_j,  sk_i = g_i in_i,
        !                 s' = s + h sum_i b_i sk_i.
        ! Reverse (sbar = cotangent of s, seeded by lam = cotangent of s'):
        sbar = lam
        do i = 1, 6
            skbar(:, i) = h * b(i) * lam
        end do
        do i = 6, 1, -1
            call var_rhs_adj(tn(i), yn(:, i), skbar(:, i), inbar)  ! g_i^T skbar_i
            sbar = sbar + inbar
            do j = 1, i - 1
                skbar(:, j) = skbar(:, j) + h * a(i, j) * inbar
            end do
        end do
        lam = sbar
    end subroutine augmented_step_adj

    ! --- internals ---

    logical function validate_problem(problem, status) result(ok)
        type(ode_problem_t),    intent(in)  :: problem
        type(fortnum_status_t), intent(out) :: status
        ok = .false.
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate: rhs not associated")
            return
        end if
        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate: y0 not allocated")
            return
        end if
        if (size(problem%y0) < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate: neq < 1")
            return
        end if
        if (problem%rtol <= 0.0_dp .or. problem%atol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate: tolerances must be positive")
            return
        end if
        if (problem%max_steps < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate: max_steps < 1")
            return
        end if
        ok = .true.
    end function validate_problem

    ! Allocate stage and temporary arrays once when the size changes.
    subroutine ensure_workspace(workspace, neq)
        type(ode_workspace_t), intent(inout) :: workspace
        integer,               intent(in)    :: neq
        if (workspace%neq == neq .and. allocated(workspace%k1)) return
        workspace%neq = neq
        call realloc(workspace%k1, neq)
        call realloc(workspace%k2, neq)
        call realloc(workspace%k3, neq)
        call realloc(workspace%k4, neq)
        call realloc(workspace%k5, neq)
        call realloc(workspace%k6, neq)
        call realloc(workspace%ytmp, neq)
        call realloc(workspace%y5, neq)
        call realloc(workspace%y4, neq)
        call realloc(workspace%yerr, neq)
    end subroutine ensure_workspace

    subroutine realloc(a, n)
        real(dp), allocatable, intent(inout) :: a(:)
        integer,               intent(in)    :: n
        if (allocated(a)) then
            if (size(a) == n) return
            deallocate(a)
        end if
        allocate(a(n))
    end subroutine realloc

    ! RMS norm of the local error scaled by atol + rtol*max(|y|,|y5|).
    real(dp) function error_norm(y, y5, yerr, rtol, atol) result(en)
        real(dp), intent(in) :: y(:), y5(:), yerr(:)
        real(dp), intent(in) :: rtol, atol
        real(dp) :: acc, sc
        integer  :: i, n
        n = size(y)
        acc = 0.0_dp
        do i = 1, n
            sc = atol + rtol * max(abs(y(i)), abs(y5(i)))
            acc = acc + (yerr(i) / sc)**2
        end do
        en = sqrt(acc / real(n, dp))
    end function error_norm

    ! Step-growth factor. PI form on accepted history; pure I form when
    ! restarting (Hairer I, II.4). Clamped to [FAC_MIN, FAC_MAX].
    real(dp) function control_factor(err_norm, err_prev, restart) result(fac)
        real(dp), intent(in) :: err_norm, err_prev
        logical,  intent(in) :: restart
        real(dp) :: e
        e = max(err_norm, 1.0e-10_dp)
        if (restart) then
            fac = SAFETY * e**(-1.0_dp / 5.0_dp)
        else
            fac = SAFETY * e**(-PI_ALPHA) * err_prev**(PI_BETA)
        end if
        fac = max(FAC_MIN, min(FAC_MAX, fac))
    end function control_factor

    ! Starting step from the scaled RHS magnitude at t0 (Hairer I, II.4),
    ! falling back to a fraction of the span. Returns a positive magnitude.
    real(dp) function initial_step(problem, workspace, hmax) result(h)
        type(ode_problem_t),   intent(in)    :: problem
        type(ode_workspace_t), intent(inout) :: workspace
        real(dp),              intent(in)    :: hmax
        real(dp) :: d0, d1, sc
        integer  :: i, n
        if (problem%h0 > 0.0_dp) then
            h = min(problem%h0, hmax)
            return
        end if
        n = size(problem%y0)
        call problem%rhs(problem%t0, problem%y0, workspace%k1)
        d0 = 0.0_dp
        d1 = 0.0_dp
        do i = 1, n
            sc = problem%atol + problem%rtol * abs(problem%y0(i))
            d0 = d0 + (problem%y0(i) / sc)**2
            d1 = d1 + (workspace%k1(i) / sc)**2
        end do
        d0 = sqrt(d0 / real(n, dp))
        d1 = sqrt(d1 / real(n, dp))
        if (d1 <= 1.0e-10_dp) then
            h = 1.0e-6_dp
        else
            h = 0.01_dp * d0 / d1
        end if
        h = max(h, 1.0e-10_dp)
        h = min(h, hmax)
        if (h <= 0.0_dp) h = min(1.0e-4_dp, hmax)
    end function initial_step

    subroutine alloc_trace(solution, neq, cap)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: neq, cap
        if (allocated(solution%t)) deallocate(solution%t)
        if (allocated(solution%y)) deallocate(solution%y)
        if (allocated(solution%h)) deallocate(solution%h)
        if (allocated(solution%err)) deallocate(solution%err)
        allocate(solution%t(cap))
        allocate(solution%y(neq, cap))
        allocate(solution%h(cap))
        allocate(solution%err(cap))
    end subroutine alloc_trace

    ! Grow the recorded trace to capacity cap, preserving content.
    subroutine grow_trace(solution, neq, cap)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: neq, cap
        real(dp), allocatable :: tt(:), yy(:,:), hh(:), ee(:)
        integer :: old
        old = size(solution%t)
        allocate(tt(cap))
        allocate(yy(neq, cap))
        allocate(hh(cap))
        allocate(ee(cap))
        tt(1:old) = solution%t
        yy(:,1:old) = solution%y
        hh(1:old) = solution%h
        ee(1:old) = solution%err
        call move_alloc(tt, solution%t)
        call move_alloc(yy, solution%y)
        call move_alloc(hh, solution%h)
        call move_alloc(ee, solution%err)
    end subroutine grow_trace

    ! Trim trace arrays to the final recorded length: t/y carry nstep+1
    ! points, h/err carry nstep accepted steps.
    subroutine trim_trace(solution, nstep, neq)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: nstep, neq
        real(dp), allocatable :: tt(:), yy(:,:), hh(:), ee(:)
        integer :: npts
        npts = nstep + 1
        allocate(tt(npts))
        allocate(yy(neq, npts))
        tt = solution%t(1:npts)
        yy = solution%y(:, 1:npts)
        call move_alloc(tt, solution%t)
        call move_alloc(yy, solution%y)
        if (nstep > 0) then
            allocate(hh(nstep))
            allocate(ee(nstep))
            hh = solution%h(1:nstep)
            ee = solution%err(1:nstep)
            call move_alloc(hh, solution%h)
            call move_alloc(ee, solution%err)
        else
            if (allocated(solution%h)) deallocate(solution%h)
            if (allocated(solution%err)) deallocate(solution%err)
            allocate(solution%h(0))
            allocate(solution%err(0))
        end if
    end subroutine trim_trace

end module fortnum_ode
