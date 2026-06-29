module fortnum_ode_wrapper
    ! Convenience wrapper: evaluate the ODE solution at a caller-specified set
    ! of output times t_eval. The adaptive integrator (ode_integrate) runs from
    ! one output point to the next, reusing the workspace across segments.
    !
    ! Derivative policy: trace_rule (ad.md §1, §4; ode.md §1, §5).
    !   Each segment's trace follows the same trace_rule as ode_integrate; a
    !   sensitivity would concatenate the per-segment variational systems on the
    !   frozen per-segment meshes. Reserved derivative names:
    !   ode_at_jvp, ode_at_vjp (defined under issue #40; not implemented here).
    !
    ! Caller contract for t_eval:
    !   - Must be strictly monotone (increasing or decreasing).
    !   - t_eval(1) is the first output time; the integrator starts from y0 at
    !     t0 = t_eval(1). Equivalently the caller may set t0 = t_eval(1).
    !   - Events are not exposed; use ode_integrate directly for event detection.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, &
        ode_solution_t, ode_integrate
    implicit none
    private

    public :: ode_at

contains

    ! Integrate problem%rhs from t_eval(1) to t_eval(n), storing the solution
    ! at each requested output time in y_out(:, i) = y(t_eval(i)).
    !
    ! t_eval must be strictly monotone. Each segment [t_eval(i), t_eval(i+1)]
    ! is a separate ode_integrate call; the workspace is reused across segments
    ! so stage arrays are allocated once.
    !
    ! On entry problem%y0 is the initial condition at t_eval(1). On return
    ! y_out is allocated (neq, n) and y_out(:, 1) = problem%y0. If any segment
    ! fails, y_out is left in the partially filled state, status carries the
    ! first error, and the subroutine returns immediately.
    subroutine ode_at(problem, t_eval, workspace, y_out, status)
        type(ode_problem_t),    intent(in)    :: problem
        real(dp),               intent(in)    :: t_eval(:)
        type(ode_workspace_t),  intent(inout) :: workspace
        real(dp), allocatable,  intent(out)   :: y_out(:,:)
        type(fortnum_status_t), intent(out)   :: status

        type(ode_problem_t)   :: seg
        type(ode_solution_t)  :: sol
        integer :: n, neq, i
        real(dp) :: dir, step_dir

        call status_set(status, FORTNUM_OK, "")

        n = size(t_eval)
        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_at: y0 not allocated")
            allocate(y_out(0, 0))
            return
        end if
        neq = size(problem%y0)

        if (n < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_at: t_eval must have at least one point")
            allocate(y_out(neq, 0))
            return
        end if

        allocate(y_out(neq, n))
        y_out(:, 1) = problem%y0

        if (n == 1) return   ! single output point: nothing to integrate

        ! Direction must be consistent across all intervals.
        dir = sign(1.0_dp, t_eval(n) - t_eval(1))
        do i = 1, n - 1
            step_dir = sign(1.0_dp, t_eval(i+1) - t_eval(i))
            if (abs(t_eval(i+1) - t_eval(i)) == 0.0_dp .or. &
                step_dir /= dir) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "ode_at: t_eval must be strictly monotone")
                return
            end if
        end do

        ! Copy problem fields that stay constant across segments.
        seg%rhs  => problem%rhs
        seg%rtol = problem%rtol
        seg%atol = problem%atol
        seg%h0   = problem%h0
        seg%hmin = problem%hmin
        seg%hmax = problem%hmax
        seg%max_steps = problem%max_steps

        ! Integrate segment by segment. The endpoint of segment i becomes
        ! y0 for segment i+1 (the last column of y from the previous call).
        seg%y0 = problem%y0
        do i = 1, n - 1
            seg%t0 = t_eval(i)
            seg%t1 = t_eval(i+1)
            call ode_integrate(seg, workspace, sol, status)
            if (status%code /= FORTNUM_OK) return
            ! Last accepted point is the endpoint t_eval(i+1).
            y_out(:, i+1) = sol%y(:, sol%nsteps + 1)
            seg%y0 = y_out(:, i+1)
        end do
    end subroutine ode_at

end module fortnum_ode_wrapper
