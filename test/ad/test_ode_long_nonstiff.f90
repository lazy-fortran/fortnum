program test_ode_long_nonstiff
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_checkpoint_t, ode_recompute_trace_t, ode_integrate, &
        ode_integrate_jvp, ode_integrate_vjp, ode_build_checkpoints, &
        ode_integrate_vjp_checkpointed, ode_build_recompute_trace, &
        ode_integrate_vjp_recomputed
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    real(dp), parameter :: nonstiff_matrix(2, 2) = reshape([ &
        -0.3_dp, 0.7_dp, -1.1_dp, 0.2_dp], [2, 2])
    real(dp), parameter :: stiff_rate = 1000.0_dp
    real(dp) :: system_matrix(2, 2) = nonstiff_matrix
    real(dp) :: final_time = 20.0_dp
    real(dp) :: fixed_step = 0.05_dp
    logical :: stiff_workload = .false.
    character(32) :: action, candidate, workload

    call configure_workload()
    call get_environment_variable("FORTNUM_LONG_ODE_ACTION", action)
    call get_environment_variable("FORTNUM_LONG_ODE_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call run_benchmark(trim(candidate))
    else if (trim(action) == "--perf") then
        call run_perf(trim(candidate))
    else if (trim(action) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        call validate_candidates()
    end if

contains

    subroutine configure_workload()
        character(32) :: requested

        workload = "long_nonstiff"
        call get_environment_variable("FORTNUM_ODE_TRAJECTORY", requested)
        if (trim(requested) == "stiff") then
            workload = "stiff"
            stiff_workload = .true.
            system_matrix = reshape([ &
                -1.0_dp, stiff_rate, 0.0_dp, -stiff_rate], [2, 2])
            final_time = 1.0_dp
            fixed_step = 0.0_dp
        end if
    end subroutine configure_workload

    subroutine primal_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt(1) = system_matrix(1, 1)*y(1) + &
            system_matrix(1, 2)*y(2) + 0.0_dp*t
        dydt(2) = system_matrix(2, 1)*y(1) + &
            system_matrix(2, 2)*y(2)
    end subroutine primal_rhs

    subroutine tangent_rhs(t, y, tangent, derivative, ctx)
        real(dp), intent(in) :: t, y(:), tangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx

        derivative(1) = system_matrix(1, 1)*tangent(1) + &
            system_matrix(1, 2)*tangent(2) + 0.0_dp*(t + sum(y))
        derivative(2) = system_matrix(2, 1)*tangent(1) + &
            system_matrix(2, 2)*tangent(2)
    end subroutine tangent_rhs

    subroutine adjoint_rhs(t, y, cotangent, derivative, ctx)
        real(dp), intent(in) :: t, y(:), cotangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx

        derivative(1) = system_matrix(1, 1)*cotangent(1) + &
            system_matrix(2, 1)*cotangent(2) + 0.0_dp*(t + sum(y))
        derivative(2) = system_matrix(1, 2)*cotangent(1) + &
            system_matrix(2, 2)*cotangent(2)
    end subroutine adjoint_rhs

    subroutine make_problem(y0, problem)
        real(dp), intent(in) :: y0(2)
        type(ode_problem_t), intent(out) :: problem

        problem%rhs => primal_rhs
        problem%t0 = 0.0_dp
        problem%t1 = final_time
        problem%y0 = y0
        if (fixed_step > 0.0_dp) then
            problem%h0 = fixed_step
            problem%hmax = fixed_step
            problem%rtol = 1.0_dp
            problem%atol = 1.0_dp
        else
            problem%rtol = 1.0e-8_dp
            problem%atol = 1.0e-10_dp
        end if
    end subroutine make_problem

    subroutine integrate_primal(y0, problem, solution)
        real(dp), intent(in) :: y0(2)
        type(ode_problem_t), intent(out) :: problem
        type(ode_solution_t), intent(out) :: solution
        type(ode_workspace_t) :: workspace
        type(fortnum_status_t) :: status

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "long ODE primal failed"
    end subroutine integrate_primal

    function full_reverse(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2)
        real(dp), allocatable :: result(:)
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        call integrate_primal(y0, problem, solution)
        call ode_integrate_vjp(problem, adjoint_rhs, cotangent, solution, &
            result, status)
        if (.not. status_ok(status)) error stop "long ODE reverse failed"
        vjp = result
    end function full_reverse

    function forward_reconstruction(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2), seed(2)
        real(dp), allocatable :: tangent(:)
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        integer :: input

        call integrate_primal(y0, problem, solution)
        do input = 1, 2
            seed = 0.0_dp
            seed(input) = 1.0_dp
            call ode_integrate_jvp(problem, tangent_rhs, seed, solution, &
                tangent, status)
            if (.not. status_ok(status)) error stop "long ODE forward failed"
            vjp(input) = dot_product(cotangent, tangent)
        end do
    end function forward_reconstruction

    function checkpointed_reverse(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2)
        real(dp), allocatable :: result(:)
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(ode_checkpoint_t) :: checkpoints
        type(fortnum_status_t) :: status

        call integrate_primal(y0, problem, solution)
        call ode_build_checkpoints(solution, 16, checkpoints, status)
        if (.not. status_ok(status)) error stop "long checkpoint build failed"
        deallocate(solution%y)
        call ode_integrate_vjp_checkpointed(problem, adjoint_rhs, cotangent, &
            checkpoints, result, status)
        if (.not. status_ok(status)) error stop "long checkpoint reverse failed"
        vjp = result
    end function checkpointed_reverse

    function recomputed_reverse(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2)
        real(dp), allocatable :: result(:)
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(ode_recompute_trace_t) :: trace
        type(fortnum_status_t) :: status

        call integrate_primal(y0, problem, solution)
        call ode_build_recompute_trace(solution, trace, status)
        if (.not. status_ok(status)) error stop "long recompute trace failed"
        deallocate(solution%y)
        call ode_integrate_vjp_recomputed(problem, adjoint_rhs, cotangent, &
            trace, result, status)
        if (.not. status_ok(status)) error stop "long recompute reverse failed"
        vjp = result
    end function recomputed_reverse

    function diagnostic_vjp(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp), parameter :: difference_step = 1.0e-5_dp
        real(dp) :: vjp(2), plus(2), minus(2)
        integer :: input

        do input = 1, 2
            plus = y0
            minus = y0
            plus(input) = plus(input) + difference_step
            minus(input) = minus(input) - difference_step
            vjp(input) = (objective(plus, cotangent) - &
                objective(minus, cotangent))/(2.0_dp*difference_step)
        end do
    end function diagnostic_vjp

    function objective(y0, cotangent) result(value)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: value
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution

        call integrate_primal(y0, problem, solution)
        value = dot_product(cotangent, solution%y(:, solution%nsteps + 1))
    end function objective

    pure function exact_vjp(cotangent) result(vjp)
        real(dp), intent(in) :: cotangent(2)
        real(dp) :: vjp(2), trace, determinant, discriminant
        real(dp) :: mean, frequency, scale, c0, c1

        if (stiff_workload) then
            scale = exp(-final_time)
            c1 = exp(-stiff_rate*final_time)
            c0 = stiff_rate*(scale - c1)/(stiff_rate - 1.0_dp)
            vjp(1) = scale*cotangent(1) + c0*cotangent(2)
            vjp(2) = c1*cotangent(2)
            return
        end if

        trace = system_matrix(1, 1) + system_matrix(2, 2)
        determinant = system_matrix(1, 1)*system_matrix(2, 2) - &
            system_matrix(1, 2)*system_matrix(2, 1)
        discriminant = trace*trace - 4.0_dp*determinant
        mean = 0.5_dp*trace
        frequency = 0.5_dp*sqrt(-discriminant)
        scale = exp(mean*final_time)
        c0 = scale*(cos(frequency*final_time) - &
            mean*sin(frequency*final_time)/frequency)
        c1 = scale*sin(frequency*final_time)/frequency
        vjp(1) = c0*cotangent(1) + c1*( &
            system_matrix(1, 1)*cotangent(1) + &
            system_matrix(2, 1)*cotangent(2))
        vjp(2) = c0*cotangent(2) + c1*( &
            system_matrix(1, 2)*cotangent(1) + &
            system_matrix(2, 2)*cotangent(2))
    end function exact_vjp

    subroutine validate_candidates()
        real(dp), parameter :: y0(2) = [0.8_dp, -0.4_dp]
        real(dp), parameter :: cotangent(2) = [0.6_dp, -1.2_dp]
        real(dp) :: reference(2), values(2, 4), diagnostic(2)
        real(dp) :: derivative_error, diagnostic_error
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        integer :: candidate_index

        call integrate_primal(y0, problem, solution)
        if (stiff_workload) then
            if (solution%nsteps < 100) then
                print *, "stiff trajectory has too few steps", solution%nsteps
                error stop 1
            end if
        else
            if (solution%nsteps < 390) then
                print *, "long trajectory has too few steps", solution%nsteps
                error stop 1
            end if
        end if
        reference = exact_vjp(cotangent)
        values(:, 1) = full_reverse(y0, cotangent)
        values(:, 2) = forward_reconstruction(y0, cotangent)
        values(:, 3) = checkpointed_reverse(y0, cotangent)
        values(:, 4) = recomputed_reverse(y0, cotangent)
        diagnostic = diagnostic_vjp(y0, cotangent)
        derivative_error = 0.0_dp
        do candidate_index = 1, 4
            derivative_error = max(derivative_error, &
                maxval(abs(values(:, candidate_index) - reference)))
        end do
        diagnostic_error = maxval(abs(diagnostic - reference))
        if (derivative_error > candidate_tolerance()) then
            print *, "long ODE derivative error", derivative_error
            error stop 1
        end if
        if (diagnostic_error > diagnostic_tolerance()) then
            print *, "long ODE diagnostic error", diagnostic_error
            error stop 1
        end if
        print *, "PASS ODE trajectory", trim(workload), solution%nsteps, &
            derivative_error, diagnostic_error
    end subroutine validate_candidates

    pure real(dp) function candidate_tolerance() result(tolerance)
        tolerance = 2.0e-8_dp
        if (stiff_workload) tolerance = 2.0e-7_dp
    end function candidate_tolerance

    pure real(dp) function diagnostic_tolerance() result(tolerance)
        tolerance = 2.0e-7_dp
        if (stiff_workload) tolerance = 5.0e-7_dp
    end function diagnostic_tolerance

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "full" .and. name /= "forward" .and. &
                name /= "checkpoint" .and. name /= "recompute" .and. &
                name /= "diagnostic") then
            error stop "unknown long ODE candidate"
        end if
    end subroutine validate_name

    subroutine run_benchmark(name)
        character(*), intent(in) :: name
        integer, parameter :: samples = 21
        integer(int64) :: reps
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        reps = candidate_repetitions(name)
        do sample = 1, 3
            call time_candidate(name, max(1_int64, reps/10_int64), sink)
        end do
        do sample = 1, samples
            call time_candidate(name, reps, elapsed(sample))
        end do
        call report(name, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_perf(name)
        character(*), intent(in) :: name
        real(dp) :: elapsed
        integer(int64) :: reps

        call validate_name(name)
        reps = 1000_int64
        if (name == "recompute") reps = 100_int64
        call time_candidate(name, reps, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, 20_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    pure function candidate_repetitions(name) result(reps)
        character(*), intent(in) :: name
        integer(int64) :: reps

        reps = 200_int64
        if (name == "recompute") reps = 10_int64
    end function candidate_repetitions

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        real(dp) :: y0(2), cotangent(2), value(2), sink
        integer(int64) :: iteration, start, finish, rate

        cotangent = [0.6_dp, -1.2_dp]
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1_int64, reps
            y0(1) = 0.8_dp + &
                1.0e-9_dp*real(mod(iteration, 101_int64), dp)
            y0(2) = -0.4_dp
            select case (name)
            case ("full")
                value = full_reverse(y0, cotangent)
            case ("forward")
                value = forward_reconstruction(y0, cotangent)
            case ("checkpoint")
                value = checkpointed_reverse(y0, cotangent)
            case ("recompute")
                value = recomputed_reverse(y0, cotangent)
            case ("diagnostic")
                value = diagnostic_vjp(y0, cotangent)
            end select
            sink = sink + sum(value)
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "long ODE benchmark produced NaN"
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
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
        write (*, "(a,',',i0,2(',',f12.4))") name, reps, median, mad
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

end program test_ode_long_nonstiff
