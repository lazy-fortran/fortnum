program test_fortnum_ode_dop853
    ! Behavioural tests for the adaptive RK8(7)13M integrator
    ! (fortnum_ode_dop853).
    !   - scalar decay y' = -y reaches exp(-t) to near machine precision
    !   - harmonic oscillator conserves energy over many periods
    !   - the adaptive step count is finite and the trace is populated/exact-end
    !   - the procedural wrapper ode_solve_dop returns the same endpoint
    !   - backward integration lands on t1
    !   - invalid input is rejected before stepping
    !   - the observed convergence order exceeds Cash-Karp RK5(4): halving the
    !     fixed step shrinks the endpoint error far faster (order >= 7)

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t
    use fortnum_ode_dop853, only: ode_integrate_dop, ode_solve_dop, dop853_step
    use fortnum_ode_cash_karp, only: cash_karp_step
    implicit none

    integer :: nfail

    nfail = 0
    call check_decay(nfail)
    call check_oscillator_energy(nfail)
    call check_trace_finite(nfail)
    call check_solve_wrapper(nfail)
    call check_backward(nfail)
    call check_bad_input(nfail)
    call check_order_vs_cash_karp(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_dop853: all tests passed"
    stop 0

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        dydt(1) = -y(1)
    end subroutine rhs_decay

    subroutine rhs_osc(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
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
        problem%rtol = 1.0e-13_dp
        problem%atol = 1.0e-15_dp
        call ode_integrate_dop(problem, workspace, solution, status)
        yend  = solution%y(1, solution%nsteps + 1)
        exact = exp(-problem%t1)
        errabs = abs(yend - exact)
        ! Eighth-order at this tolerance: expect well under 1e-12.
        if (status%code /= FORTNUM_OK .or. errabs > 1.0e-12_dp) then
            write (error_unit, "(a,es14.6,a,es14.6,a,i0)") &
                "FAIL check_decay: yend=", yend, " err=", errabs, &
                " code=", status%code
            nfail = nfail + 1
        end if
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
        problem%rtol = 1.0e-12_dp
        problem%atol = 1.0e-14_dp
        call ode_integrate_dop(problem, workspace, solution, status)
        e0 = 0.5_dp * (solution%y(1,1)**2 + solution%y(2,1)**2)
        maxdrift = 0.0_dp
        do k = 1, solution%nsteps + 1
            ek = 0.5_dp * (solution%y(1,k)**2 + solution%y(2,k)**2)
            maxdrift = max(maxdrift, abs(ek - e0))
        end do
        if (status%code /= FORTNUM_OK .or. maxdrift > 1.0e-9_dp) then
            write (error_unit, "(a,es14.6,a,i0)") &
                "FAIL check_oscillator_energy: maxdrift=", maxdrift, &
                " code=", status%code
            nfail = nfail + 1
        end if
        if (abs(solution%y(1, solution%nsteps + 1) - 1.0_dp) > 1.0e-9_dp) then
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
        call ode_integrate_dop(problem, workspace, solution, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") &
                "FAIL check_trace_finite: code=", status%code
            nfail = nfail + 1
            return
        end if
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
    end subroutine check_trace_finite

    subroutine check_solve_wrapper(nfail)
        integer, intent(inout) :: nfail
        real(dp), allocatable  :: t_out(:), y_out(:,:)
        type(fortnum_status_t) :: status
        real(dp) :: yend, exact
        call ode_solve_dop(rhs_decay, 0.0_dp, 3.0_dp, [2.0_dp], t_out, y_out, &
                           status, rtol=1.0e-13_dp, atol=1.0e-15_dp)
        if (status%code /= FORTNUM_OK .or. size(t_out) < 2) then
            write (error_unit, "(a,i0)") &
                "FAIL check_solve_wrapper: code=", status%code
            nfail = nfail + 1
            return
        end if
        yend  = y_out(1, size(t_out))
        exact = 2.0_dp * exp(-3.0_dp)
        if (abs(yend - exact) > 1.0e-12_dp) then
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
        problem%rtol = 1.0e-12_dp
        problem%atol = 1.0e-14_dp
        call ode_integrate_dop(problem, workspace, solution, status)
        yend  = solution%y(1, solution%nsteps + 1)
        exact = 1.0_dp
        if (status%code /= FORTNUM_OK .or. abs(yend - exact) > 1.0e-11_dp) then
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
        problem%t0 = 0.0_dp
        problem%t1 = 1.0_dp
        problem%y0 = [1.0_dp]
        call ode_integrate_dop(problem, workspace, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_bad_input (null rhs): code=", status%code
            nfail = nfail + 1
        end if
        problem%rhs => rhs_decay
        problem%rtol = 0.0_dp
        call ode_integrate_dop(problem, workspace, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_bad_input (rtol=0): code=", status%code
            nfail = nfail + 1
        end if
    end subroutine check_bad_input

    ! Empirical order check: one fixed step of size h on y'=-y from y(0)=1 has
    ! local error ~ C h^(p+1). Halving h must shrink the error by ~2^(p+1).
    ! For RK8 (p=8) the ratio is ~512; for Cash-Karp RK5 (p=5) it is ~64. The
    ! DOP853 ratio must comfortably beat the Cash-Karp ratio, confirming the
    ! tableau is genuinely higher order. Steps are taken small enough that the
    ! leading-order term dominates but large enough to stay above rounding.
    subroutine check_order_vs_cash_karp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: y0(1), h1, h2
        real(dp) :: e_dop_1, e_dop_2, e_ck_1, e_ck_2
        real(dp) :: ratio_dop, ratio_ck
        y0 = [1.0_dp]
        h1 = 0.25_dp
        h2 = 0.125_dp
        e_dop_1 = dop_local_err(y0, h1)
        e_dop_2 = dop_local_err(y0, h2)
        e_ck_1  = ck_local_err(y0, h1)
        e_ck_2  = ck_local_err(y0, h2)
        ratio_dop = e_dop_1 / e_dop_2
        ratio_ck  = e_ck_1 / e_ck_2
        write (*, "(a,es12.4,a,es12.4)") &
            "dop853 local err h,h/2 =", e_dop_1, " ,", e_dop_2
        write (*, "(a,f10.2,a,f10.2)") &
            "halving ratio dop853 =", ratio_dop, "  cash-karp =", ratio_ck
        ! RK8 observed ratio must be at least 2^7 = 128 (one order of safety
        ! below the theoretical 2^9), and must exceed the Cash-Karp ratio.
        if (ratio_dop < 128.0_dp) then
            write (error_unit, "(a,f12.2)") &
                "FAIL order: dop853 halving ratio too low ", ratio_dop
            nfail = nfail + 1
        end if
        if (ratio_dop <= ratio_ck) then
            write (error_unit, "(a,f12.2,a,f12.2)") &
                "FAIL order: dop853 ratio ", ratio_dop, &
                " not above cash-karp ", ratio_ck
            nfail = nfail + 1
        end if
    end subroutine check_order_vs_cash_karp

    ! Local error of one DOP853 step of size h on y'=-y from y0 vs exp(-h)*y0.
    function dop_local_err(y0, h) result(e)
        real(dp), intent(in) :: y0(1), h
        real(dp) :: e
        real(dp) :: k1(1), k2(1), k3(1), k4(1), k5(1), k6(1)
        real(dp) :: k7(1), k8(1), k9(1), k10(1), k11(1), k12(1)
        real(dp) :: ytmp(1), y8(1), err5(1), err3(1)
        integer  :: nfev
        nfev = 0
        call dop853_step(rhs_decay, 0.0_dp, y0, h, .false., k1, k2, k3, k4, &
            k5, k6, k7, k8, k9, k10, k11, k12, ytmp, y8, err5, err3, nfev)
        e = abs(y8(1) - y0(1) * exp(-h))
    end function dop_local_err

    ! Local error of one Cash-Karp RK5 step of size h on y'=-y from y0.
    function ck_local_err(y0, h) result(e)
        real(dp), intent(in) :: y0(1), h
        real(dp) :: e
        real(dp) :: k1(1), k2(1), k3(1), k4(1), k5(1), k6(1)
        real(dp) :: ytmp(1), y5(1), yerr(1)
        integer  :: nfev
        nfev = 0
        call cash_karp_step(rhs_decay, 0.0_dp, y0, h, .false., k1, k2, k3, &
            k4, k5, k6, ytmp, y5, yerr, nfev)
        e = abs(y5(1) - y0(1) * exp(-h))
    end function ck_local_err

end program test_fortnum_ode_dop853
