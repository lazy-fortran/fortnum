program test_ode_event_time_tangent
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ODE_EVENT_RISING
    use fortnum_ode_events, only: ode_event_result_t, ode_event_scan, &
        ode_event_time_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    real(dp) :: current_rate = 1.0_dp
    real(dp) :: current_threshold = 0.8_dp
    character(32) :: action, candidate
    integer :: parameter_count, direction_count

    call get_environment_variable("FORTNUM_EVENT_TIME_ACTION", action)
    call get_environment_variable("FORTNUM_EVENT_TIME_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call read_dimensions(parameter_count, direction_count)
        call run_benchmark(trim(candidate), parameter_count, direction_count)
    else if (trim(action) == "--perf") then
        call read_dimensions(parameter_count, direction_count)
        call run_perf(trim(candidate), parameter_count, direction_count)
    else if (trim(action) == "--peak-rss") then
        call read_dimensions(parameter_count, direction_count)
        call run_peak_rss(trim(candidate), parameter_count, direction_count)
    else
        call validate_product()
    end if

contains

    subroutine validate_product()
        integer, parameter :: nparameter = 2, ndirection = 3
        real(dp) :: parameters(nparameter), parameter_direction(nparameter, ndirection)
        real(dp) :: y0, threshold, y0_direction(ndirection)
        real(dp) :: threshold_direction(ndirection), residual_tangent(ndirection)
        real(dp) :: expected(ndirection), rate_direction, exact_time
        real(dp), allocatable :: tangent(:)
        type(ode_event_result_t) :: event_result
        type(fortnum_status_t) :: status
        integer :: direction

        parameters = [0.8_dp, 1.2_dp]
        y0 = 0.1_dp
        threshold = 0.9_dp
        y0_direction = [0.2_dp, -0.1_dp, 0.3_dp]
        threshold_direction = [-0.2_dp, 0.4_dp, 0.1_dp]
        parameter_direction(:, 1) = [0.3_dp, -0.1_dp]
        parameter_direction(:, 2) = [-0.2_dp, 0.5_dp]
        parameter_direction(:, 3) = [0.1_dp, 0.2_dp]

        call solve_event(y0, parameters, threshold, event_result)
        exact_time = (threshold - y0) / &
            (sum(parameters)/real(nparameter, dp))
        do direction = 1, ndirection
            rate_direction = sum(parameter_direction(:, direction)) / &
                real(nparameter, dp)
            residual_tangent(direction) = y0_direction(direction) + &
                event_result%t_event*rate_direction - &
                threshold_direction(direction)
            expected(direction) = &
                (threshold_direction(direction) - y0_direction(direction)) / &
                current_rate - exact_time*rate_direction/current_rate
        end do

        call ode_event_time_jvp(event_result, residual_tangent, tangent, status)
        if (.not. status_ok(status)) error stop "event-time JVP failed"
        if (maxval(abs(tangent - expected)) > 2.0e-7_dp) then
            print *, "event-time exact-oracle mismatch", tangent, expected
            error stop 1
        end if

        event_result%transversal = .false.
        call ode_event_time_jvp(event_result, residual_tangent, tangent, status)
        if (status_ok(status)) error stop "non-transversal event accepted"
        print *, "PASS analytical event-time tangent"
    end subroutine validate_product

    subroutine rhs_linear(t, y, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        associate (unused_t => t, unused_y => y); end associate
            dydt = current_rate
    end subroutine rhs_linear

    function threshold_event(t, y, ctx) result(value)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: value

        associate (unused_t => t); end associate
            value = y(1) - current_threshold
    end function threshold_event

    subroutine solve_event(y0, parameters, threshold, event_result)
        real(dp), intent(in) :: y0, parameters(:), threshold
        type(ode_event_result_t), intent(out) :: event_result
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        current_rate = sum(parameters)/real(size(parameters), dp)
        current_threshold = threshold
        problem%rhs => rhs_linear
        problem%t0 = 0.0_dp
        problem%t1 = 2.0_dp
        allocate(problem%y0(1))
        problem%y0(1) = y0
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "event primal failed"
        call ode_event_scan(rhs_linear, threshold_event, solution, &
            ODE_EVENT_RISING, 1.0e-12_dp, event_result, status)
        if (.not. status_ok(status)) error stop "event scan failed"
        if (.not. event_result%found) error stop "event not found"
    end subroutine solve_event

    subroutine read_dimensions(nparameter, ndirection)
        integer, intent(out) :: nparameter, ndirection
        character(32) :: text

        call get_environment_variable("FORTNUM_EVENT_TIME_PARAMETERS", text)
        read (text, *) nparameter
        call get_environment_variable("FORTNUM_EVENT_TIME_DIRECTIONS", text)
        read (text, *) ndirection
        if (nparameter < 1) error stop "parameter count must be positive"
        if (ndirection < 1) error stop "direction count must be positive"
    end subroutine read_dimensions

    subroutine run_benchmark(name, nparameter, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, ndirection
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 1000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, nparameter, ndirection, &
                reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, nparameter, ndirection, reps, &
                elapsed(sample))
        end do
        call report(name, nparameter, ndirection, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_perf(name, nparameter, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, ndirection
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, nparameter, ndirection, 100000_int64, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine run_peak_rss(name, nparameter, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, ndirection
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, nparameter, ndirection, 5000_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "analytical" .and. name /= "diagnostic") then
            error stop "candidate must be analytical or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, nparameter, ndirection, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, ndirection
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        real(dp) :: parameters(nparameter), parameter_direction(nparameter, ndirection)
        real(dp) :: y0, threshold, y0_direction(ndirection)
        real(dp) :: threshold_direction(ndirection), residual_tangent(ndirection)
        real(dp), allocatable :: tangent(:)
        type(ode_event_result_t) :: event_result
        type(fortnum_status_t) :: status
        integer(int64) :: iteration, start, finish, rate
        integer :: parameter, direction
        real(dp) :: rate_direction, sink

        do parameter = 1, nparameter
            parameters(parameter) = 0.8_dp + &
                0.4_dp*real(parameter - 1, dp)/real(max(1, nparameter - 1), dp)
            do direction = 1, ndirection
                parameter_direction(parameter, direction) = &
                    0.01_dp*real(mod(parameter + 2*direction, 11) - 5, dp)
            end do
        end do
        do direction = 1, ndirection
            y0_direction(direction) = &
                0.02_dp*real(mod(direction, 7) - 3, dp)
            threshold_direction(direction) = &
                0.015_dp*real(mod(2*direction, 9) - 4, dp)
        end do
        threshold = 0.9_dp
        sink = 0.0_dp

        call system_clock(start, rate)
        do iteration = 1_int64, reps
            y0 = 0.1_dp + 1.0e-8_dp*real(mod(iteration, 101_int64), dp)
            call solve_event(y0, parameters, threshold, event_result)
            if (name == "analytical") then
                do direction = 1, ndirection
                    rate_direction = &
                        sum(parameter_direction(:, direction)) / &
                        real(nparameter, dp)
                    residual_tangent(direction) = y0_direction(direction) + &
                        event_result%t_event*rate_direction - &
                        threshold_direction(direction)
                end do
                call ode_event_time_jvp(event_result, residual_tangent, &
                    tangent, status)
                if (.not. status_ok(status)) error stop "benchmark JVP failed"
            else
                call diagnostic_tangents(y0, parameters, threshold, &
                    y0_direction, parameter_direction, threshold_direction, &
                    tangent)
            end if
            sink = sink + event_result%t_event + sum(tangent)
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = 1.0e9_dp*real(finish - start, dp) / &
            (real(rate, dp)*real(reps, dp))
    end subroutine time_candidate

    subroutine diagnostic_tangents(y0, parameters, threshold, y0_direction, &
            parameter_direction, threshold_direction, tangent)
        real(dp), intent(in) :: y0, parameters(:), threshold
        real(dp), intent(in) :: y0_direction(:), parameter_direction(:,:)
        real(dp), intent(in) :: threshold_direction(:)
        real(dp), allocatable, intent(out) :: tangent(:)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp) :: parameter_work(size(parameters)), plus_time, minus_time
        type(ode_event_result_t) :: event_result
        integer :: direction

        allocate(tangent(size(y0_direction)))
        do direction = 1, size(y0_direction)
            parameter_work = parameters + &
                step*parameter_direction(:, direction)
            call solve_event(y0 + step*y0_direction(direction), parameter_work, &
                threshold + step*threshold_direction(direction), event_result)
            plus_time = event_result%t_event
            parameter_work = parameters - &
                step*parameter_direction(:, direction)
            call solve_event(y0 - step*y0_direction(direction), parameter_work, &
                threshold - step*threshold_direction(direction), event_result)
            minus_time = event_result%t_event
            tangent(direction) = (plus_time - minus_time)/(2.0_dp*step)
        end do
    end subroutine diagnostic_tangents

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

    subroutine report(name, nparameter, ndirection, values, reps)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, ndirection
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
        write (*, "(a,3(',',i0),2(',',f12.4))") name, nparameter, &
            ndirection, reps, median, mad
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

end program test_ode_event_time_tangent
