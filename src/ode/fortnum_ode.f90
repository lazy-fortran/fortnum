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
    !   additively (reserved names ode_integrate_jvp, ode_integrate_vjp).
    !   Active: y0, ctx parameters, the recorded trace. Inactive: rtol, atol,
    !   h0, hmin, hmax, max_steps, event_direction, and the counts/err/status.
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

    public :: ode_rhs_t, ode_event_t
    public :: ode_integrate, ode_solve

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
