module ode_rhs_forward_autodiff
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_ode_scalar_rhs, only: &
        fortnum_enzyme_ode_scalar_rhs_jvp
    implicit none
    private

    public :: rhs_jvp

contains

    pure function scalar_rhs(t, y, k) result(value) &
            bind(c, name="fortnum_ode_scalar_rhs")
        real(c_double), value :: t, y, k
        real(c_double) :: value

        value = -k*y + 0.0_c_double*t
    end function scalar_rhs

    function rhs_jvp(t, y, dy, k, dk) result(derivative)
        real(c_double), intent(in) :: t, y, dy, k, dk
        real(c_double) :: derivative

        derivative = fortnum_enzyme_ode_scalar_rhs_jvp(t, 0.0_c_double, &
            y, dy, k, dk)
    end function rhs_jvp

end module ode_rhs_forward_autodiff

module ode_forward_sensitivity_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    use ode_rhs_forward_autodiff, only: rhs_jvp
    implicit none
    private

    real(dp) :: active_k = 0.7_dp
    real(dp) :: active_direction = 1.0_dp

    public :: analytical_value_jvp, hybrid_value_jvp, diagnostic_value_jvp
    public :: exact_jvp

contains

    function analytical_value_jvp(k) result(derivative)
        real(dp), intent(in) :: k
        real(dp) :: derivative

        active_k = k
        active_direction = 1.0_dp
        derivative = frozen_sensitivity(analytical_var_rhs)
    end function analytical_value_jvp

    function hybrid_value_jvp(k) result(derivative)
        real(dp), intent(in) :: k
        real(dp) :: derivative

        active_k = k
        active_direction = 1.0_dp
        derivative = frozen_sensitivity(hybrid_var_rhs)
    end function hybrid_value_jvp

    function diagnostic_value_jvp(k) result(derivative)
        real(dp), intent(in) :: k
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: derivative, base_value

        base_value = primal_value(k)
        derivative = (primal_value(k + h) - primal_value(k - h))/(2.0_dp*h)
        if (base_value == huge(base_value)) error stop "unreachable ODE value"
    end function diagnostic_value_jvp

    function frozen_sensitivity(var_rhs) result(derivative)
        interface
            subroutine var_rhs(t, y, s, dsdt, ctx)
                import :: dp
                real(dp), intent(in) :: t, y(:), s(:)
                real(dp), intent(out) :: dsdt(:)
                class(*), intent(in), optional :: ctx
            end subroutine var_rhs
        end interface
        real(dp) :: derivative
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp), allocatable :: sensitivity(:)

        call configure_problem(problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "ODE primal trace failed"
        call ode_integrate_jvp(problem, var_rhs, [0.0_dp], solution, &
            sensitivity, status)
        if (.not. status_ok(status)) error stop "ODE forward sensitivity failed"
        derivative = sensitivity(1)
    end function frozen_sensitivity

    function primal_value(k) result(value)
        real(dp), intent(in) :: k
        real(dp) :: value
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        active_k = k
        call configure_problem(problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "ODE diagnostic primal failed"
        value = solution%y(1, solution%nsteps + 1)
    end function primal_value

    subroutine configure_problem(problem)
        type(ode_problem_t), intent(out) :: problem

        problem%rhs => primal_rhs
        problem%t0 = 0.0_dp
        problem%t1 = 2.0_dp
        problem%y0 = [1.3_dp]
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
    end subroutine configure_problem

    subroutine primal_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt(1) = -active_k*y(1) + 0.0_dp*t
    end subroutine primal_rhs

    subroutine analytical_var_rhs(t, y, s, dsdt, ctx)
        real(dp), intent(in) :: t, y(:), s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx

        dsdt(1) = -active_k*s(1) - active_direction*y(1) + 0.0_dp*t
    end subroutine analytical_var_rhs

    subroutine hybrid_var_rhs(t, y, s, dsdt, ctx)
        real(dp), intent(in) :: t, y(:), s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx

        dsdt(1) = rhs_jvp(t, y(1), s(1), active_k, active_direction)
    end subroutine hybrid_var_rhs

    pure function exact_jvp(k) result(derivative)
        real(dp), intent(in) :: k
        real(dp) :: derivative

        derivative = -2.0_dp*1.3_dp*exp(-2.0_dp*k)
    end function exact_jvp

end module ode_forward_sensitivity_kernel

program enzyme_ode_forward_sensitivity_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_result
    use ode_forward_sensitivity_kernel, only: analytical_value_jvp, &
        hybrid_value_jvp, diagnostic_value_jvp, exact_jvp
    implicit none

    character(32) :: action, candidate
    real(dp) :: samples(fixture_sample_count), sink
    integer :: candidate_kind, repetitions

    call get_command_argument(1, action)
    call get_command_argument(2, candidate)
    if (len_trim(action) == 0) then
        call get_environment_variable("FORTNUM_ODE_ACTION", action)
        call get_environment_variable("FORTNUM_ODE_CANDIDATE", candidate)
    end if

    if (trim(action) == "--benchmark") then
        call run_benchmark(trim(candidate))
    else if (trim(action) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        call validate_candidates()
    end if

contains

    subroutine validate_candidates()
        real(dp) :: values(3), reference

        reference = exact_jvp(0.7_dp)
        values(1) = analytical_value_jvp(0.7_dp)
        values(2) = hybrid_value_jvp(0.7_dp)
        values(3) = diagnostic_value_jvp(0.7_dp)
        if (maxval(abs(values(1:2) - reference)) > 2.0e-8_dp .or. &
            abs(values(3) - reference) > 2.0e-7_dp) then
            print *, "ODE hybrid forward-sensitivity mismatch", values - reference
            error stop 1
        end if
        print *, "PASS ODE hybrid forward sensitivity", values - reference
    end subroutine validate_candidates

    subroutine run_benchmark(name)
        character(*), intent(in) :: name
        real(dp) :: median, mad

        call validate_name(name)
        call select_candidate(name)
        repetitions = 2000
        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_result(name, repetitions, median, mad)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name

        call validate_name(name)
        call select_candidate(name)
        repetitions = 5000
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if ((name /= "analytical") .and. (name /= "hybrid") .and. &
            (name /= "diagnostic")) then
            error stop "ODE candidate must be analytical, hybrid, or diagnostic"
        end if
    end subroutine validate_name

    subroutine select_candidate(name)
        character(*), intent(in) :: name

        select case (name)
        case ("analytical")
            candidate_kind = 1
        case ("hybrid")
            candidate_kind = 2
        case default
            candidate_kind = 3
        end select
    end subroutine select_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        real(dp) :: k, local_sink
        integer :: iteration

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            k = 0.7_dp + 1.0e-6_dp*real(mod(iteration, 101), dp)
            select case (candidate_kind)
            case (1)
                local_sink = local_sink + analytical_value_jvp(k)
            case (2)
                local_sink = local_sink + hybrid_value_jvp(k)
            case default
                local_sink = local_sink + diagnostic_value_jvp(k)
            end select
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "ODE forward-sensitivity benchmark failed"
        end if
        sink = local_sink
    end function measure_candidate

end program enzyme_ode_forward_sensitivity_jvp
