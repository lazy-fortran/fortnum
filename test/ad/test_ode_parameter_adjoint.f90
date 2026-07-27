module test_ode_parameter_adjoint_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp, ode_integrate_parameter_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    real(dp), parameter :: initial_value = 1.3_dp
    real(dp), parameter :: final_time = 2.0_dp
    real(dp), parameter :: terminal_cotangent = 1.1_dp
    real(dp) :: active_mean = 0.7_dp
    integer :: active_parameter_count = 1

    public :: reverse_parameter_vjp, forward_parameter_vjp
    public :: diagnostic_parameter_vjp, exact_parameter_vjp

contains

    subroutine reverse_parameter_vjp(parameter_count, gradient)
        integer, intent(in) :: parameter_count
        real(dp), allocatable, intent(out) :: gradient(:)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp), allocatable :: initial_vjp(:)

        active_parameter_count = parameter_count
        call make_problem(problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "parameter-adjoint primal failed"
        call ode_integrate_parameter_vjp(problem, adjoint_rhs, parameter_stage_vjp, &
            [terminal_cotangent], parameter_count, solution, initial_vjp, &
            gradient, status)
        if (.not. status_ok(status)) error stop "parameter adjoint failed"
        if (abs(initial_vjp(1) - terminal_cotangent* &
                exp(-active_mean*final_time)) > 2.0e-8_dp) then
            error stop "initial-state adjoint mismatch"
        end if
    end subroutine reverse_parameter_vjp

    subroutine forward_parameter_vjp(parameter_count, gradient)
        integer, intent(in) :: parameter_count
        real(dp), allocatable, intent(out) :: gradient(:)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp), allocatable :: tangent(:)
        integer :: parameter

        active_parameter_count = parameter_count
        call make_problem(problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "parameter-forward primal failed"
        allocate(gradient(parameter_count))
        do parameter = 1, parameter_count
            call ode_integrate_jvp(problem, parameter_tangent_rhs, [0.0_dp], &
                solution, tangent, status)
            if (.not. status_ok(status)) error stop "parameter tangent failed"
            gradient(parameter) = terminal_cotangent*tangent(1)
        end do
    end subroutine forward_parameter_vjp

    subroutine diagnostic_parameter_vjp(parameter_count, gradient)
        integer, intent(in) :: parameter_count
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp), allocatable, intent(out) :: gradient(:)
        real(dp) :: base, plus, minus
        integer :: parameter

        active_parameter_count = parameter_count
        base = primal_objective(active_mean)
        allocate(gradient(parameter_count))
        do parameter = 1, parameter_count
            plus = primal_objective(active_mean + h/real(parameter_count, dp))
            minus = primal_objective(active_mean - h/real(parameter_count, dp))
            gradient(parameter) = (plus - minus)/(2.0_dp*h)
        end do
        if (base == huge(base)) error stop "unreachable parameter objective"
    end subroutine diagnostic_parameter_vjp

    function primal_objective(mean_parameter) result(value)
        real(dp), intent(in) :: mean_parameter
        real(dp) :: value, saved_mean
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        saved_mean = active_mean
        active_mean = mean_parameter
        call make_problem(problem)
        call ode_integrate(problem, workspace, solution, status)
        active_mean = saved_mean
        if (.not. status_ok(status)) error stop "parameter diagnostic failed"
        value = terminal_cotangent*solution%y(1, solution%nsteps + 1)
    end function primal_objective

    subroutine make_problem(problem)
        type(ode_problem_t), intent(out) :: problem

        problem%rhs => primal_rhs
        problem%t0 = 0.0_dp
        problem%t1 = final_time
        problem%y0 = [initial_value]
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
    end subroutine make_problem

    subroutine primal_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt(1) = -active_mean*y(1) + 0.0_dp*t
    end subroutine primal_rhs

    subroutine adjoint_rhs(t, y, cotangent, derivative, ctx)
        real(dp), intent(in) :: t, y(:), cotangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx

        derivative(1) = -active_mean*cotangent(1) + 0.0_dp*(t + y(1))
    end subroutine adjoint_rhs

    subroutine parameter_stage_vjp(t, y, rhs_cotangent, gradient, ctx)
        real(dp), intent(in) :: t, y(:), rhs_cotangent(:)
        real(dp), intent(inout) :: gradient(:)
        class(*), intent(in), optional :: ctx

        gradient = gradient - y(1)*rhs_cotangent(1)/ &
            real(active_parameter_count, dp) + 0.0_dp*t
    end subroutine parameter_stage_vjp

    subroutine parameter_tangent_rhs(t, y, tangent, derivative, ctx)
        real(dp), intent(in) :: t, y(:), tangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx

        derivative(1) = -active_mean*tangent(1) - &
            y(1)/real(active_parameter_count, dp) + 0.0_dp*t
    end subroutine parameter_tangent_rhs

    subroutine exact_parameter_vjp(parameter_count, gradient)
        integer, intent(in) :: parameter_count
        real(dp), allocatable, intent(out) :: gradient(:)

        allocate(gradient(parameter_count))
        gradient = -terminal_cotangent*final_time*initial_value* &
            exp(-active_mean*final_time)/real(parameter_count, dp)
    end subroutine exact_parameter_vjp

end module test_ode_parameter_adjoint_kernel

program test_ode_parameter_adjoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use test_ode_parameter_adjoint_kernel, only: reverse_parameter_vjp, &
        forward_parameter_vjp, diagnostic_parameter_vjp, exact_parameter_vjp
    implicit none

    character(32) :: action, candidate, count_text
    integer :: parameter_count

    call get_environment_variable("FORTNUM_ODE_PARAM_ACTION", action)
    call get_environment_variable("FORTNUM_ODE_PARAM_CANDIDATE", candidate)
    call get_environment_variable("FORTNUM_ODE_PARAM_COUNT", count_text)
    if (trim(action) == "--benchmark" .or. trim(action) == "--peak-rss") then
        read (count_text, *) parameter_count
    else
        parameter_count = 1
    end if

    if (trim(action) == "--benchmark") then
        call run_benchmark(trim(candidate), parameter_count)
    else if (trim(action) == "--peak-rss") then
        call run_peak_rss(trim(candidate), parameter_count)
    else
        call validate_candidates()
    end if

contains

    subroutine validate_candidates()
        integer, parameter :: counts(3) = [1, 4, 16]
        real(dp), allocatable :: reference(:), candidate_value(:)
        real(dp) :: errors(3)
        integer :: count_index

        do count_index = 1, size(counts)
            call exact_parameter_vjp(counts(count_index), reference)
            call reverse_parameter_vjp(counts(count_index), candidate_value)
            errors(1) = maxval(abs(candidate_value - reference))
            call forward_parameter_vjp(counts(count_index), candidate_value)
            errors(2) = maxval(abs(candidate_value - reference))
            call diagnostic_parameter_vjp(counts(count_index), candidate_value)
            errors(3) = maxval(abs(candidate_value - reference))
            if (maxval(errors(1:2)) > 2.0e-8_dp .or. &
                errors(3) > 2.0e-7_dp) then
                print *, "ODE parameter-adjoint mismatch", counts(count_index), &
                    errors
                error stop 1
            end if
        end do
        print *, "PASS ODE parameter adjoint"
    end subroutine validate_candidates

    subroutine run_benchmark(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 200_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_request(name, count)
        do sample = 1, 3
            call time_candidate(name, count, reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, count, reps, elapsed(sample))
        end do
        call report(name, count, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), parameter :: reps = 500_int64
        real(dp) :: elapsed

        call validate_request(name, count)
        call time_candidate(name, count, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_request(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count

        if ((name /= "reverse") .and. (name /= "forward") .and. &
            (name /= "diagnostic")) then
            error stop "parameter candidate must be reverse, forward, or diagnostic"
        end if
        if ((count /= 1) .and. (count /= 4) .and. (count /= 16)) then
            error stop "parameter count must be 1, 4, or 16"
        end if
    end subroutine validate_request

    subroutine time_candidate(name, count, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        real(dp), allocatable :: result(:)
        real(dp) :: sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            select case (name)
            case ("reverse")
                call reverse_parameter_vjp(count, result)
            case ("forward")
                call forward_parameter_vjp(count, result)
            case ("diagnostic")
                call diagnostic_parameter_vjp(count, result)
            end select
            sink = sink + sum(result) + 0.0_dp*real(iteration, dp)
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
        if (sink == huge(sink)) print *, sink
    end subroutine time_candidate

    function peak_rss_bytes() result(bytes)
        integer(int64) :: bytes, kilobytes
        integer :: unit, io_status
        character(256) :: line

        bytes = 0_int64
        open (newunit=unit, file="/proc/self/status", status="old", &
            action="read", iostat=io_status)
        if (io_status /= 0) return
        do
            read (unit, "(a)", iostat=io_status) line
            if (io_status /= 0) exit
            if (index(line, "VmHWM:") == 1) then
                read (line(7:), *, iostat=io_status) kilobytes
                if (io_status == 0) bytes = 1024_int64*kilobytes
                exit
            end if
        end do
        close (unit)
    end function peak_rss_bytes

    subroutine report(name, count, values, reps)
        character(*), intent(in) :: name
        integer, intent(in) :: count
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
        write (*, "(a,',',i0,',',i0,',',f12.4,',',f12.4)") &
            name, count, reps, median, mad
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

end program test_ode_parameter_adjoint
