program test_ode_continuous_sensitivity
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    real(dp), parameter :: initial_value = 1.3_dp
    real(dp), parameter :: final_time = 2.0_dp
    real(dp) :: current_parameter = 0.7_dp
    real(dp) :: current_direction = 1.0_dp
    character(32) :: action, candidate
    integer :: direction_count

    call get_environment_variable("FORTNUM_CONTINUOUS_ACTION", action)
    call get_environment_variable("FORTNUM_CONTINUOUS_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call read_direction_count(direction_count)
        call run_benchmark(trim(candidate), direction_count)
    else if (trim(action) == "--perf") then
        call read_direction_count(direction_count)
        call run_perf(trim(candidate), direction_count)
    else if (trim(action) == "--peak-rss") then
        call read_direction_count(direction_count)
        call run_peak_rss(trim(candidate), direction_count)
    else
        call validate_continuous_contract()
    end if

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        associate (unused_t => t); end associate
            dydt = -current_parameter*y
    end subroutine rhs_decay

    subroutine tangent_decay(t, y, tangent, tangent_rhs, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), tangent(:)
        real(dp), intent(out) :: tangent_rhs(:)
        class(*), intent(in), optional :: ctx

        associate (unused_t => t); end associate
            tangent_rhs = -current_parameter*tangent - current_direction*y
    end subroutine tangent_decay

    subroutine validate_continuous_contract()
        integer, parameter :: refinement_count = 4
        real(dp), parameter :: step_limits(refinement_count) = &
            [0.5_dp, 0.25_dp, 0.125_dp, 0.0625_dp]
        real(dp) :: errors(refinement_count), exact
        real(dp), allocatable :: tangent(:)
        type(ode_solution_t) :: solution
        integer :: refinement

        current_parameter = 0.7_dp
        current_direction = 1.0_dp
        exact = -final_time*initial_value* &
            exp(-current_parameter*final_time)
        do refinement = 1, refinement_count
            call integrate_primal(current_parameter, solution, &
                step_limits(refinement))
            call continuous_jvp(solution, tangent)
            errors(refinement) = abs(tangent(1) - exact)
        end do

        do refinement = 2, refinement_count
            if (errors(refinement) >= errors(refinement - 1)) then
                print *, "continuous sensitivity did not converge", errors
                error stop 1
            end if
        end do
        if (errors(refinement_count) > 1.0e-8_dp) then
            print *, "continuous sensitivity refinement error", errors
            error stop 1
        end if
        print *, "PASS continuous sensitivity contract", errors
    end subroutine validate_continuous_contract

    subroutine integrate_primal(parameter, solution, step_limit)
        real(dp), intent(in) :: parameter
        type(ode_solution_t), intent(out) :: solution
        real(dp), intent(in), optional :: step_limit
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(fortnum_status_t) :: status

        current_parameter = parameter
        problem%rhs => rhs_decay
        problem%t0 = 0.0_dp
        problem%t1 = final_time
        allocate(problem%y0(1))
        problem%y0(1) = initial_value
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
        if (present(step_limit)) then
            problem%h0 = step_limit
            problem%hmax = step_limit
            problem%rtol = 1.0_dp
            problem%atol = 1.0_dp
        end if
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "continuous primal failed"
    end subroutine integrate_primal

    subroutine continuous_jvp(solution, tangent)
        type(ode_solution_t), intent(in) :: solution
        real(dp), allocatable, intent(out) :: tangent(:)
        type(ode_problem_t) :: problem
        type(fortnum_status_t) :: status
        real(dp) :: seed(1)

        problem%rhs => rhs_decay
        seed = 0.0_dp
        call ode_integrate_jvp(problem, tangent_decay, seed, solution, &
            tangent, status)
        if (.not. status_ok(status)) error stop "continuous JVP failed"
    end subroutine continuous_jvp

    subroutine read_direction_count(count)
        integer, intent(out) :: count
        character(32) :: text

        call get_environment_variable("FORTNUM_CONTINUOUS_DIRECTIONS", text)
        read (text, *) count
        if (count < 1) error stop "direction count must be positive"
    end subroutine read_direction_count

    subroutine run_benchmark(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 2000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, count, reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, count, reps, elapsed(sample))
        end do
        call report(name, count, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_perf(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, count, 20000_int64, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine run_peak_rss(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, count, 10000_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "primal" .and. name /= "continuous" .and. &
                name /= "diagnostic") then
            error stop "candidate must be primal, continuous, or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, count, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        type(ode_solution_t) :: solution, perturbed
        real(dp), allocatable :: tangent(:)
        integer(int64) :: iteration, start, finish, rate
        integer :: direction
        real(dp) :: parameter, scale, plus_value, minus_value, sink
        real(dp), parameter :: difference_step = 1.0e-5_dp

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1_int64, reps
            parameter = 0.7_dp + &
                1.0e-8_dp*real(mod(iteration, 101_int64), dp)
            call integrate_primal(parameter, solution)
            if (name == "primal") then
                sink = sink + solution%y(1, solution%nsteps + 1)
            else if (name == "continuous") then
                do direction = 1, count
                    current_direction = direction_scale(direction)
                    call continuous_jvp(solution, tangent)
                    sink = sink + tangent(1)
                end do
            else
                do direction = 1, count
                    scale = direction_scale(direction)
                    call integrate_primal(parameter + &
                        difference_step*scale, perturbed)
                    plus_value = perturbed%y(1, perturbed%nsteps + 1)
                    call integrate_primal(parameter - &
                        difference_step*scale, perturbed)
                    minus_value = perturbed%y(1, perturbed%nsteps + 1)
                    sink = sink + (plus_value - minus_value) / &
                        (2.0_dp*difference_step)
                end do
            end if
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = 1.0e9_dp*real(finish - start, dp) / &
            (real(rate, dp)*real(reps, dp))
    end subroutine time_candidate

    pure real(dp) function direction_scale(direction) result(scale)
        integer, intent(in) :: direction

        scale = 0.1_dp*real(mod(direction, 11) - 5, dp)
    end function direction_scale

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
        write (*, "(a,2(',',i0),2(',',f12.4))") name, count, reps, &
            median, mad
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

end program test_ode_continuous_sensitivity
