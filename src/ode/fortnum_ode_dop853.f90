module fortnum_ode_dop853
    ! Adaptive Prince-Dormand RK8(7)-13M integrator (dop853 class). Twelve
    ! explicit stages give an eighth-order solution; two embedded error
    ! estimators (orders 5 and 3) drive the adaptive step, for the very tight
    ! tolerances KiLCA needs (Hairer, Norsett, Wanner; Prince-Dormand 1981).
    !
    ! Derivative policy: trace_rule (ad.md sec 1, 4; ode.md sec 1, 4).
    !   Same policy as the Cash-Karp integrator in fortnum_ode: the adaptive
    !   schedule is data-dependent, and a sensitivity differentiates the frozen
    !   accepted-step mesh (solution%t, solution%y, solution%h) with that
    !   schedule held fixed. The high-order forward/reverse products ride the
    !   recorded trace exactly as fortnum_ode's do; this module shares the
    !   ode_problem_t / ode_workspace_t / ode_solution_t carriers and the
    !   recorded mesh those products walk. Active: y0, ctx parameters. Inactive:
    !   rtol, atol, h0, hmin, hmax, max_steps.
    !
    ! Method and coefficients: Prince and Dormand, "High order embedded
    !   Runge-Kutta formulae", J. Comput. Appl. Math. 7 (1981) 67-75; Hairer,
    !   Norsett, Wanner, "Solving Ordinary Differential Equations I", 2nd ed.,
    !   II.5 (DOP853). The 13M node/coupling/weight constants and the order-5
    !   and order-3 error weights are the published RK8(7)13M values. PI step
    !   control and the starting-step estimate follow Hairer I, II.4.
    !
    ! No module-level state. The caller owns the workspace and the recorded
    ! trace; the step writes only its output arguments and stage slots.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    ! The nodes, coupling matrix and weights are derived from the reduced
    ! system, not declared here. See tools/codegen/app/gen_dop853_tableau.f90.
    use fortnum_dop853_tableau, only: dop853_c, dop853_b, dop853_a, &
        dop853_e5, dop853_e3
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_rhs_t
    implicit none
    private

    public :: dop853_step, ode_integrate_dop, ode_solve_dop

    ! PI controller for an order-8 method (Hairer I, II.4). The error estimate
    ! the step accepts is order 8 (err = err5^2 / sqrt(err5^2 + 0.01 err3^2)
    ! has the order-8 leading term), so the controller exponents use order 8.
    real(dp), parameter :: PI_ALPHA = 1.0_dp / 8.0_dp - 0.2_dp * 0.75_dp / 8.0_dp
    real(dp), parameter :: PI_BETA  = 0.2_dp * 0.75_dp / 8.0_dp
    real(dp), parameter :: SAFETY   = 0.9_dp
    real(dp), parameter :: FAC_MIN  = 0.333_dp
    real(dp), parameter :: FAC_MAX  = 6.0_dp
    real(dp), parameter :: TRACE_GROWTH = 2.0_dp

    ! Stage nodes c2..c12 (c1 = 0 is implicit).

    ! Coupling coefficients a(i,j), j < i.

    ! Eighth-order solution weights (stages 2..5 carry zero weight).

    ! Order-5 embedded error weights (err5 = sum E5_i k_i).

    ! Order-3 embedded error weights (err3 = sum E3_i k_i).

contains

    ! Advance one RK8(7)13M step of size h from (t, y). rhs matches ode_rhs_t.
    ! k1..k12 are caller-owned stage-derivative slots (length neq); ytmp is a
    ! caller-owned scratch state. y8 receives the eighth-order solution, err5
    ! and err3 the two embedded error vectors (already scaled by h). nfev is
    ! incremented by the stages evaluated.
    !
    ! have_k1 reuses an externally supplied first-stage derivative (the method
    ! is not FSAL, so the integrator only sets it when k1 already holds
    ! f(t, y)).
    subroutine dop853_step(rhs, t, y, h, have_k1, k1, k2, k3, k4, k5, k6, &
            k7, k8, k9, k10, k11, k12, ytmp, y8, err5, err3, &
            nfev, ctx)
        procedure(ode_rhs_t)            :: rhs
        real(dp), intent(in)            :: t
        real(dp), intent(in)            :: y(:)
        real(dp), intent(in)            :: h
        logical,  intent(in)            :: have_k1
        real(dp), intent(inout)         :: k1(:)
        real(dp), intent(out)           :: k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), intent(out)           :: k7(:), k8(:), k9(:), k10(:)
        real(dp), intent(out)           :: k11(:), k12(:)
        real(dp), intent(out)           :: ytmp(:)
        real(dp), intent(out)           :: y8(:)
        real(dp), intent(out)           :: err5(:), err3(:)
        integer,  intent(inout)         :: nfev
        class(*), intent(in), optional  :: ctx

        if (.not. have_k1) then
            call rhs(t, y, k1, ctx)
            nfev = nfev + 1
        end if

        ytmp = y + h * (dop853_a(2,1) * k1)
        call rhs(t + dop853_c(2) * h, ytmp, k2, ctx)

        ytmp = y + h * (dop853_a(3,1) * k1 + dop853_a(3,2) * k2)
        call rhs(t + dop853_c(3) * h, ytmp, k3, ctx)

        ytmp = y + h * (dop853_a(4,1) * k1 + dop853_a(4,3) * k3)
        call rhs(t + dop853_c(4) * h, ytmp, k4, ctx)

        ytmp = y + h * (dop853_a(5,1) * k1 + dop853_a(5,3) * k3 + dop853_a(5,4) * k4)
        call rhs(t + dop853_c(5) * h, ytmp, k5, ctx)

        ytmp = y + h * (dop853_a(6,1) * k1 + dop853_a(6,4) * k4 + dop853_a(6,5) * k5)
        call rhs(t + dop853_c(6) * h, ytmp, k6, ctx)

        ytmp = y + h * (dop853_a(7,1) * k1 + dop853_a(7,4) * k4 + dop853_a(7,5) * k5 + dop853_a(7,6) * k6)
        call rhs(t + dop853_c(7) * h, ytmp, k7, ctx)

        ytmp = y + h * (dop853_a(8,1) * k1 + dop853_a(8,4) * k4 + dop853_a(8,5) * k5 + dop853_a(8,6) * k6 &
            + dop853_a(8,7) * k7)
        call rhs(t + dop853_c(8) * h, ytmp, k8, ctx)

        ytmp = y + h * (dop853_a(9,1) * k1 + dop853_a(9,4) * k4 + dop853_a(9,5) * k5 + dop853_a(9,6) * k6 &
            + dop853_a(9,7) * k7 + dop853_a(9,8) * k8)
        call rhs(t + dop853_c(9) * h, ytmp, k9, ctx)

        ytmp = y + h * (dop853_a(10,1) * k1 + dop853_a(10,4) * k4 + dop853_a(10,5) * k5 + dop853_a(10,6) * k6 &
            + dop853_a(10,7) * k7 + dop853_a(10,8) * k8 + dop853_a(10,9) * k9)
        call rhs(t + dop853_c(10) * h, ytmp, k10, ctx)

        ytmp = y + h * (dop853_a(11,1) * k1 + dop853_a(11,4) * k4 + dop853_a(11,5) * k5 + dop853_a(11,6) * k6 &
            + dop853_a(11,7) * k7 + dop853_a(11,8) * k8 + dop853_a(11,9) * k9 + dop853_a(11,10) * k10)
        call rhs(t + dop853_c(11) * h, ytmp, k11, ctx)

        ytmp = y + h * (dop853_a(12,1) * k1 + dop853_a(12,4) * k4 + dop853_a(12,5) * k5 + dop853_a(12,6) * k6 &
            + dop853_a(12,7) * k7 + dop853_a(12,8) * k8 + dop853_a(12,9) * k9 + dop853_a(12,10) * k10 &
            + dop853_a(12,11) * k11)
        call rhs(t + dop853_c(12) * h, ytmp, k12, ctx)

        nfev = nfev + 11

        y8 = y + h * (dop853_b(1) * k1 + dop853_b(6) * k6 + dop853_b(7) * k7 + dop853_b(8) * k8 + dop853_b(9) * k9 &
            + dop853_b(10) * k10 + dop853_b(11) * k11 + dop853_b(12) * k12)

        err5 = h * (dop853_e5(1) * k1 + dop853_e5(6) * k6 + dop853_e5(7) * k7 + dop853_e5(8) * k8 + dop853_e5(9) * k9 &
            + dop853_e5(10) * k10 + dop853_e5(11) * k11 + dop853_e5(12) * k12)
        err3 = h * (dop853_e3(1) * k1 + dop853_e3(6) * k6 + dop853_e3(7) * k7 + dop853_e3(8) * k8 + dop853_e3(9) * k9 &
            + dop853_e3(10) * k10 + dop853_e3(11) * k11 + dop853_e3(12) * k12)
    end subroutine dop853_step

    ! Integrate problem%rhs from t0 to t1 with adaptive RK8(7)13M. Records the
    ! accepted-step mesh into solution, matching ode_integrate's layout so the
    ! fortnum_ode trace_rule products can walk it. Events are not exposed.
    subroutine ode_integrate_dop(problem, workspace, solution, status)
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
                    "ode_integrate_dop: exceeded max_steps")
                exit
            end if

            final_step = abs(problem%t1 - t) <= h
            if (final_step) h = abs(problem%t1 - t)
            if (h <= 0.0_dp) exit

            call dop853_step(problem%rhs, t, solution%y(:,nstep+1), dir * h, &
                .false., workspace%k1, workspace%k2, workspace%k3, &
                workspace%k4, workspace%k5, workspace%k6, workspace%k7, &
                workspace%k8, workspace%k9, workspace%k10, workspace%k11, &
                workspace%k12, workspace%ytmp, workspace%y8, &
                workspace%err5, workspace%err3, solution%nfev)

            err_norm = error_norm(solution%y(:,nstep+1), workspace%y8, &
                workspace%err5, workspace%err3, &
                problem%rtol, problem%atol)
            accepted = err_norm <= 1.0_dp

            if (accepted) then
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
                solution%y(:,nstep+1) = workspace%y8
                solution%h(nstep) = dir * h
                solution%err(nstep) = err_norm
                if (final_step) exit
            else
                solution%nrejected = solution%nrejected + 1
            end if

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
                    "ode_integrate_dop: step forced below hmin")
                exit
            end if
            h = h_new
        end do

        solution%nsteps = nstep
        call trim_trace(solution, nstep, neq)
        solution%status = status
    end subroutine ode_integrate_dop

    ! Flat call: build a problem, integrate with RK8(7)13M, hand back the trace.
    subroutine ode_solve_dop(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
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
        ! Explicit allocation, not allocation-on-assignment: consumers that
        ! compile with -fno-realloc-lhs (libneo does, and SIMPLE inherits it)
        ! turn "lhs = rhs" on an unallocated allocatable into a write through a
        ! null descriptor rather than an allocation, which segfaults.
        allocate (problem%y0(size(y0)))
        problem%y0 = y0
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol

        call ode_integrate_dop(problem, workspace, solution, status)

        npts = solution%nsteps + 1
        if (allocated(solution%t)) then
            allocate (t_out(npts))
            allocate (y_out(size(solution%y, 1), npts))
            t_out = solution%t(1:npts)
            y_out = solution%y(:, 1:npts)
        else
            allocate(t_out(0))
            allocate(y_out(size(y0), 0))
        end if
    end subroutine ode_solve_dop

    ! --- internals ---

    logical function validate_problem(problem, status) result(ok)
        type(ode_problem_t),    intent(in)  :: problem
        type(fortnum_status_t), intent(out) :: status
        ok = .false.
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: rhs not associated")
            return
        end if
        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: y0 not allocated")
            return
        end if
        if (size(problem%y0) < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: neq < 1")
            return
        end if
        if (problem%rtol <= 0.0_dp .or. problem%atol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: tolerances must be positive")
            return
        end if
        if (problem%max_steps < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: max_steps < 1")
            return
        end if
        ok = .true.
    end function validate_problem

    ! Allocate the twelve stage slots and temporaries when the size changes.
    ! Reuses ode_workspace_t: k1..k6 + ytmp are shared with Cash-Karp; the
    ! higher stages live in the workspace's spare slots.
    subroutine ensure_workspace(workspace, neq)
        type(ode_workspace_t), intent(inout) :: workspace
        integer,               intent(in)    :: neq
        if (workspace%neq == neq .and. allocated(workspace%k7)) return
        workspace%neq = neq
        call realloc(workspace%k1, neq)
        call realloc(workspace%k2, neq)
        call realloc(workspace%k3, neq)
        call realloc(workspace%k4, neq)
        call realloc(workspace%k5, neq)
        call realloc(workspace%k6, neq)
        call realloc(workspace%k7, neq)
        call realloc(workspace%k8, neq)
        call realloc(workspace%k9, neq)
        call realloc(workspace%k10, neq)
        call realloc(workspace%k11, neq)
        call realloc(workspace%k12, neq)
        call realloc(workspace%ytmp, neq)
        call realloc(workspace%y8, neq)
        call realloc(workspace%err5, neq)
        call realloc(workspace%err3, neq)
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

    ! Hairer-Norsett-Wanner DOP853 error norm: combine the order-5 and order-3
    ! estimates so the accepted error has the order-8 leading term. Both errN
    ! arrays already carry the h factor (see dop853_step).
    real(dp) function error_norm(y, y8, err5, err3, rtol, atol) result(en)
        real(dp), intent(in) :: y(:), y8(:), err5(:), err3(:)
        real(dp), intent(in) :: rtol, atol
        real(dp) :: sc, e5sq, e3sq, denom
        integer  :: i, n
        n = size(y)
        e5sq = 0.0_dp
        e3sq = 0.0_dp
        do i = 1, n
            sc = atol + rtol * max(abs(y(i)), abs(y8(i)))
            e5sq = e5sq + (err5(i) / sc)**2
            e3sq = e3sq + (err3(i) / sc)**2
        end do
        if (e5sq <= 0.0_dp .and. e3sq <= 0.0_dp) then
            en = 0.0_dp
            return
        end if
        denom = e5sq + 0.01_dp * e3sq
        en = sqrt(e5sq / (denom * real(n, dp))) * sqrt(e5sq)
    end function error_norm

    ! Step-growth factor. PI form on accepted history; pure I form when
    ! restarting (Hairer I, II.4). Clamped to [FAC_MIN, FAC_MAX].
    real(dp) function control_factor(err_norm, err_prev, restart) result(fac)
        real(dp), intent(in) :: err_norm, err_prev
        logical,  intent(in) :: restart
        real(dp) :: e
        e = max(err_norm, 1.0e-10_dp)
        if (restart) then
            fac = SAFETY * e**(-1.0_dp / 8.0_dp)
        else
            fac = SAFETY * e**(-PI_ALPHA) * err_prev**(PI_BETA)
        end if
        fac = max(FAC_MIN, min(FAC_MAX, fac))
    end function control_factor

    ! Starting step from the scaled RHS magnitude at t0 (Hairer I, II.4),
    ! tuned for an order-8 method via the 1/8 power. Returns a positive
    ! magnitude bounded by hmax.
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

end module fortnum_ode_dop853
