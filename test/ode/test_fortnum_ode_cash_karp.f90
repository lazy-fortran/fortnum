program test_fortnum_ode_cash_karp
    ! Behavioural tests for the adaptive Cash-Karp RK5(4) integrator.
    !   - scalar decay y' = -y reaches exp(-t) within tolerance
    !   - harmonic oscillator conserves energy over many periods
    !   - the adaptive step count is finite and the trace is populated
    !   - the procedural wrapper returns the same endpoint
    !   - backward integration (t1 < t0) lands on t1
    !   - invalid input is rejected before stepping

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_solve
    implicit none

    integer :: nfail

    nfail = 0
    call check_decay(nfail)
    call check_oscillator_energy(nfail)
    call check_trace_finite(nfail)
    call check_solve_wrapper(nfail)
    call check_backward(nfail)
    call check_bad_input(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_cash_karp: all tests passed"
    stop 0

contains

    ! y' = -y, y(0) = 1 -> y(t) = exp(-t).
    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)
    end subroutine rhs_decay

    ! Harmonic oscillator y1' = y2, y2' = -y1; energy 0.5*(y1^2+y2^2) constant.
    subroutine rhs_osc(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) =  y(2)
        dydt(2) = -y(1)
    end subroutine rhs_osc

    subroutine check_decay(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(ode_solution_t)   :: solution
        type(fortnum_status_t) :: status
        real(dp) :: yend, exact, errabs
        problem%rhs => rhs_decay
        problem%t0 = 0.0_dp
        problem%t1 = 4.0_dp
        problem%y0 = [1.0_dp]
        problem%rtol = 1.0e-8_dp
        problem%atol = 1.0e-10_dp
        call ode_integrate(problem, workspace, solution, status)
        yend  = solution%y(1, solution%nsteps + 1)
        exact = exp(-problem%t1)
        errabs = abs(yend - exact)
        if (status%code /= FORTNUM_OK .or. errabs > 1.0e-7_dp) then
            write (error_unit, "(a,es14.6,a,es14.6,a,i0)") &
                "FAIL check_decay: yend=", yend, " err=", errabs, &
                " code=", status%code
            nfail = nfail + 1
        end if
        ! Endpoint must be exactly t1.
        if (abs(solution%t(solution%nsteps + 1) - problem%t1) > 1.0e-12_dp) then
            write (error_unit, "(a,es14.6)") &
                "FAIL check_decay: endpoint t=", solution%t(solution%nsteps + 1)
            nfail = nfail + 1
        end if
    end subroutine check_decay

    subroutine check_oscillator_energy(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(ode_solution_t)   :: solution
        type(fortnum_status_t) :: status
        real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
        real(dp) :: e0, ek, maxdrift
        integer  :: k
        problem%rhs => rhs_osc
        problem%t0 = 0.0_dp
        problem%t1 = 20.0_dp * PI
        problem%y0 = [1.0_dp, 0.0_dp]
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-11_dp
        call ode_integrate(problem, workspace, solution, status)
        e0 = 0.5_dp * (solution%y(1,1)**2 + solution%y(2,1)**2)
        maxdrift = 0.0_dp
        do k = 1, solution%nsteps + 1
            ek = 0.5_dp * (solution%y(1,k)**2 + solution%y(2,k)**2)
            maxdrift = max(maxdrift, abs(ek - e0))
        end do
        if (status%code /= FORTNUM_OK .or. maxdrift > 1.0e-5_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_oscillator_energy: maxdrift=", maxdrift, &
                " code=", status%code
            nfail = nfail + 1
        end if
        ! Returns to the start after 10 full periods within tolerance.
        if (abs(solution%y(1, solution%nsteps + 1) - 1.0_dp) > 1.0e-5_dp) then
            write (error_unit, "(a,es14.6)") &
                "FAIL check_oscillator_energy: y1_end=", &
                solution%y(1, solution%nsteps + 1)
            nfail = nfail + 1
        end if
    end subroutine check_oscillator_energy

    subroutine check_trace_finite(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(ode_solution_t)   :: solution
        type(fortnum_status_t) :: status
        problem%rhs => rhs_decay
        problem%t0 = 0.0_dp
        problem%t1 = 10.0_dp
        problem%y0 = [1.0_dp]
        call ode_integrate(problem, workspace, solution, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") &
                "FAIL check_trace_finite: code=", status%code
            nfail = nfail + 1
            return
        end if
        ! Trace must be populated, finite, and consistent in size.
        if (solution%nsteps < 1 .or. solution%nsteps >= problem%max_steps) then
            write (error_unit, "(a,i0)") &
                "FAIL check_trace_finite: nsteps=", solution%nsteps
            nfail = nfail + 1
        end if
        if (size(solution%t) /= solution%nsteps + 1 .or. &
            size(solution%y, 2) /= solution%nsteps + 1 .or. &
            size(solution%h) /= solution%nsteps .or. &
            size(solution%err) /= solution%nsteps) then
            write (error_unit, "(a)") &
                "FAIL check_trace_finite: trace array sizes inconsistent"
            nfail = nfail + 1
        end if
        if (solution%nfev <= 0) then
            write (error_unit, "(a,i0)") &
                "FAIL check_trace_finite: nfev=", solution%nfev
            nfail = nfail + 1
        end if
        ! Monotone advancing mesh from t0.
        if (abs(solution%t(1) - problem%t0) > 0.0_dp) then
            write (error_unit, "(a)") "FAIL check_trace_finite: t(1) /= t0"
            nfail = nfail + 1
        end if
    end subroutine check_trace_finite

    subroutine check_solve_wrapper(nfail)
        integer, intent(inout) :: nfail
        real(dp), allocatable  :: t_out(:), y_out(:,:)
        type(fortnum_status_t) :: status
        real(dp) :: yend, exact
        call ode_solve(rhs_decay, 0.0_dp, 3.0_dp, [2.0_dp], t_out, y_out, &
                       status, rtol=1.0e-9_dp, atol=1.0e-11_dp)
        if (status%code /= FORTNUM_OK .or. size(t_out) < 2) then
            write (error_unit, "(a,i0)") &
                "FAIL check_solve_wrapper: code=", status%code
            nfail = nfail + 1
            return
        end if
        yend  = y_out(1, size(t_out))
        exact = 2.0_dp * exp(-3.0_dp)
        if (abs(yend - exact) > 1.0e-7_dp) then
            write (error_unit, "(a,es14.6,a,es14.6)") &
                "FAIL check_solve_wrapper: yend=", yend, " exact=", exact
            nfail = nfail + 1
        end if
        if (abs(t_out(size(t_out)) - 3.0_dp) > 1.0e-12_dp) then
            write (error_unit, "(a,es14.6)") &
                "FAIL check_solve_wrapper: t_end=", t_out(size(t_out))
            nfail = nfail + 1
        end if
    end subroutine check_solve_wrapper

    ! Backward integration: t1 < t0 must land on t1 and stay accurate.
    subroutine check_backward(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(ode_solution_t)   :: solution
        type(fortnum_status_t) :: status
        real(dp) :: yend, exact
        problem%rhs => rhs_decay
        problem%t0 = 2.0_dp
        problem%t1 = 0.0_dp
        problem%y0 = [exp(-2.0_dp)]
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-11_dp
        call ode_integrate(problem, workspace, solution, status)
        yend  = solution%y(1, solution%nsteps + 1)
        exact = 1.0_dp
        if (status%code /= FORTNUM_OK .or. abs(yend - exact) > 1.0e-7_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_backward: yend=", yend, " code=", status%code
            nfail = nfail + 1
        end if
        if (abs(solution%t(solution%nsteps + 1) - problem%t1) > 1.0e-12_dp) then
            write (error_unit, "(a,es14.6)") &
                "FAIL check_backward: endpoint=", &
                solution%t(solution%nsteps + 1)
            nfail = nfail + 1
        end if
    end subroutine check_backward

    subroutine check_bad_input(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(ode_solution_t)   :: solution
        type(fortnum_status_t) :: status
        ! Null rhs must be rejected with a domain error.
        problem%t0 = 0.0_dp
        problem%t1 = 1.0_dp
        problem%y0 = [1.0_dp]
        call ode_integrate(problem, workspace, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_bad_input (null rhs): code=", status%code
            nfail = nfail + 1
        end if
        ! Non-positive tolerance must be rejected.
        problem%rhs => rhs_decay
        problem%rtol = 0.0_dp
        call ode_integrate(problem, workspace, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_bad_input (rtol=0): code=", status%code
            nfail = nfail + 1
        end if
    end subroutine check_bad_input

end program test_fortnum_ode_cash_karp
