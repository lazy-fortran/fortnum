module ode_rhs_forward_autodiff
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    public :: rhs_jvp

    interface
        function enzyme_fwddiff(f, t, dt, y, dy, k, dk) result(derivative) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: t, dt, y, dy, k, dk
            real(c_double) :: derivative
        end function enzyme_fwddiff
    end interface

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

        derivative = enzyme_fwddiff(c_funloc(scalar_rhs), t, 0.0_c_double, &
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
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use ode_forward_sensitivity_kernel, only: analytical_value_jvp, &
        hybrid_value_jvp, diagnostic_value_jvp, exact_jvp
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    character(32) :: action, candidate

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
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 2000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, reps, elapsed(sample))
        end do
        call report(name, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name
        integer(int64), parameter :: reps = 5000_int64
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if ((name /= "analytical") .and. (name /= "hybrid") .and. &
            (name /= "diagnostic")) then
            error stop "ODE candidate must be analytical, hybrid, or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        real(dp) :: k, sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            k = 0.7_dp + 1.0e-6_dp*real(mod(iteration, 101_int64), dp)
            select case (name)
            case ("analytical")
                sink = sink + analytical_value_jvp(k)
            case ("hybrid")
                sink = sink + hybrid_value_jvp(k)
            case ("diagnostic")
                sink = sink + diagnostic_value_jvp(k)
            end select
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
        if (sink == huge(sink)) print *, sink
    end subroutine time_candidate

    subroutine report(name, values, reps)
        character(*), intent(in) :: name
        real(dp), intent(in) :: values(:)
        integer(int64), intent(in) :: reps
        real(dp) :: ordered(size(values)), deviations(size(values))
        real(dp) :: median, mad

        ordered = values
        call sort_values(ordered)
        median = ordered((size(ordered) + 1)/2)
        deviations = abs(values - median)
        call sort_values(deviations)
        mad = deviations((size(deviations) + 1)/2)
        write (*, "(a,',',i0,',',f12.4,',',f12.4)") &
            name, reps, median, mad
    end subroutine report

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: temporary
        integer :: i, j

        do i = 1, size(values) - 1
            do j = i + 1, size(values)
                if (values(j) < values(i)) then
                    temporary = values(i)
                    values(i) = values(j)
                    values(j) = temporary
                end if
            end do
        end do
    end subroutine sort_values

end program enzyme_ode_forward_sensitivity_jvp
