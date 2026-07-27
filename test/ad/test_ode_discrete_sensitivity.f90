program test_ode_discrete_sensitivity
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp, ode_integrate_vjp
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: state_count = 16
    real(dp), parameter :: final_time = 1.5_dp
    character(32) :: action, candidate
    integer :: product_count

    call get_environment_variable("FORTNUM_DISCRETE_ACTION", action)
    call get_environment_variable("FORTNUM_DISCRETE_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call read_product_count(product_count)
        call run_benchmark(trim(candidate), product_count)
    else if (trim(action) == "--perf") then
        call read_product_count(product_count)
        call run_perf(trim(candidate), product_count)
    else if (trim(action) == "--peak-rss") then
        call read_product_count(product_count)
        call run_peak_rss(trim(candidate), product_count)
    else
        call validate_discrete_contract()
    end if

contains

    subroutine rhs_nonlinear(t, y, dydt, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        integer :: i, previous

        associate (unused_t => t); end associate
        do i = 1, size(y)
            previous = i - 1
            if (previous == 0) previous = size(y)
            dydt(i) = (-0.2_dp - 0.01_dp*real(i, dp))*y(i) + &
                0.03_dp*y(i)*y(i) + 0.02_dp*y(previous)
        end do
    end subroutine rhs_nonlinear

    subroutine tangent_nonlinear(t, y, tangent, derivative, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), tangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx
        integer :: i, previous

        associate (unused_t => t); end associate
        do i = 1, size(y)
            previous = i - 1
            if (previous == 0) previous = size(y)
            derivative(i) = (-0.2_dp - 0.01_dp*real(i, dp) + &
                0.06_dp*y(i))*tangent(i) + 0.02_dp*tangent(previous)
        end do
    end subroutine tangent_nonlinear

    subroutine adjoint_nonlinear(t, y, cotangent, derivative, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:), cotangent(:)
        real(dp), intent(out) :: derivative(:)
        class(*), intent(in), optional :: ctx
        integer :: i, following

        associate (unused_t => t); end associate
        do i = 1, size(y)
            following = i + 1
            if (following > size(y)) following = 1
            derivative(i) = (-0.2_dp - 0.01_dp*real(i, dp) + &
                0.06_dp*y(i))*cotangent(i) + 0.02_dp*cotangent(following)
        end do
    end subroutine adjoint_nonlinear

    subroutine make_problem(y0, problem)
        real(dp), intent(in) :: y0(:)
        type(ode_problem_t), intent(out) :: problem

        problem%rhs => rhs_nonlinear
        problem%t0 = 0.0_dp
        problem%t1 = final_time
        problem%y0 = y0
        problem%rtol = 1.0e-8_dp
        problem%atol = 1.0e-10_dp
    end subroutine make_problem

    subroutine integrate_primal(y0, solution)
        real(dp), intent(in) :: y0(:)
        type(ode_solution_t), intent(out) :: solution
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(fortnum_status_t) :: status

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "discrete primal failed"
    end subroutine integrate_primal

    subroutine frozen_replay(y0, schedule, final_state)
        real(dp), intent(in) :: y0(:)
        type(ode_solution_t), intent(in) :: schedule
        real(dp), intent(out) :: final_state(:)
        real(dp) :: y(size(y0)), y5(size(y0)), yerr(size(y0))
        real(dp) :: ytmp(size(y0)), k1(size(y0)), k2(size(y0))
        real(dp) :: k3(size(y0)), k4(size(y0)), k5(size(y0)), k6(size(y0))
        integer :: step, evaluations

        y = y0
        evaluations = 0
        do step = 1, schedule%nsteps
            call cash_karp_step(rhs_nonlinear, schedule%t(step), y, &
                schedule%h(step), .false., k1, k2, k3, k4, k5, k6, &
                ytmp, y5, yerr, evaluations)
            y = y5
        end do
        final_state = y
    end subroutine frozen_replay

    subroutine apply_jvp(solution, seed, tangent)
        type(ode_solution_t), intent(in) :: solution
        real(dp), intent(in) :: seed(:)
        real(dp), allocatable, intent(out) :: tangent(:)
        type(ode_problem_t) :: problem
        type(fortnum_status_t) :: status

        problem%rhs => rhs_nonlinear
        call ode_integrate_jvp(problem, tangent_nonlinear, seed, solution, &
            tangent, status)
        if (.not. status_ok(status)) error stop "discrete JVP failed"
    end subroutine apply_jvp

    subroutine apply_vjp(solution, cotangent, transpose_product)
        type(ode_solution_t), intent(in) :: solution
        real(dp), intent(in) :: cotangent(:)
        real(dp), allocatable, intent(out) :: transpose_product(:)
        type(ode_problem_t) :: problem
        type(fortnum_status_t) :: status

        problem%rhs => rhs_nonlinear
        call ode_integrate_vjp(problem, adjoint_nonlinear, cotangent, &
            solution, transpose_product, status)
        if (.not. status_ok(status)) error stop "discrete VJP failed"
    end subroutine apply_vjp

    subroutine validate_discrete_contract()
        real(dp), parameter :: difference_step = 1.0e-5_dp
        real(dp) :: y0(state_count), seed(state_count), cotangent(state_count)
        real(dp) :: y_plus(state_count), y_minus(state_count)
        real(dp) :: replay_plus(state_count), replay_minus(state_count)
        real(dp) :: oracle(state_count), jvp_error, vjp_error, dot_error
        real(dp), allocatable :: tangent(:), transpose_product(:)
        type(ode_solution_t) :: solution

        call fill_initial_state(y0)
        call fill_seed(3, seed)
        call fill_cotangent(5, cotangent)
        call integrate_primal(y0, solution)

        y_plus = y0 + difference_step*seed
        y_minus = y0 - difference_step*seed
        call frozen_replay(y_plus, solution, replay_plus)
        call frozen_replay(y_minus, solution, replay_minus)
        oracle = (replay_plus - replay_minus)/(2.0_dp*difference_step)
        call apply_jvp(solution, seed, tangent)
        jvp_error = maxval(abs(tangent - oracle))

        call apply_vjp(solution, cotangent, transpose_product)
        dot_error = abs(dot_product(cotangent, tangent) - &
            dot_product(transpose_product, seed))
        vjp_error = scalar_objective_vjp_error(y0, cotangent, solution, &
            transpose_product, difference_step)

        if (jvp_error > 2.0e-9_dp) then
            print *, "frozen-map JVP error", jvp_error
            error stop 1
        end if
        if (vjp_error > 2.0e-9_dp) then
            print *, "frozen-map VJP error", vjp_error
            error stop 1
        end if
        if (dot_error > 2.0e-12_dp) then
            print *, "discrete adjoint identity error", dot_error
            error stop 1
        end if
        print *, "PASS discrete sensitivity contract", jvp_error, &
            vjp_error, dot_error
    end subroutine validate_discrete_contract

    function scalar_objective_vjp_error(y0, cotangent, solution, &
            transpose_product, difference_step) result(error)
        real(dp), intent(in) :: y0(:), cotangent(:), transpose_product(:)
        type(ode_solution_t), intent(in) :: solution
        real(dp), intent(in) :: difference_step
        real(dp) :: error, perturbed(size(y0)), final_state(size(y0))
        real(dp) :: finite_difference(size(y0)), plus_value, minus_value
        integer :: input

        do input = 1, size(y0)
            perturbed = y0
            perturbed(input) = perturbed(input) + difference_step
            call frozen_replay(perturbed, solution, final_state)
            plus_value = dot_product(cotangent, final_state)
            perturbed(input) = y0(input) - difference_step
            call frozen_replay(perturbed, solution, final_state)
            minus_value = dot_product(cotangent, final_state)
            finite_difference(input) = (plus_value - minus_value)/ &
                (2.0_dp*difference_step)
        end do
        error = maxval(abs(transpose_product - finite_difference))
    end function scalar_objective_vjp_error

    subroutine fill_initial_state(y0)
        real(dp), intent(out) :: y0(:)
        integer :: i

        do i = 1, size(y0)
            y0(i) = 0.3_dp + 0.02_dp*real(i, dp)
        end do
    end subroutine fill_initial_state

    subroutine fill_seed(index, seed)
        integer, intent(in) :: index
        real(dp), intent(out) :: seed(:)
        integer :: i

        do i = 1, size(seed)
            seed(i) = sin(real(i*index, dp))
        end do
    end subroutine fill_seed

    subroutine fill_cotangent(index, cotangent)
        integer, intent(in) :: index
        real(dp), intent(out) :: cotangent(:)
        integer :: i

        do i = 1, size(cotangent)
            cotangent(i) = cos(real(i*index, dp))
        end do
    end subroutine fill_cotangent

    subroutine read_product_count(count)
        integer, intent(out) :: count
        character(32) :: text

        call get_environment_variable("FORTNUM_DISCRETE_PRODUCTS", text)
        read (text, *) count
        if (count < 1 .or. count > state_count) then
            error stop "product count must be in 1:16"
        end if
    end subroutine read_product_count

    subroutine run_benchmark(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 500_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, count, reps/10_int64, sink)
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
        call time_candidate(name, count, 5000_int64, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine run_peak_rss(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, count, 2000_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "primal" .and. name /= "forward" .and. &
                name /= "reverse" .and. name /= "diagnostic") then
            error stop "candidate must be primal, forward, reverse, or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, count, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        type(ode_solution_t) :: solution
        real(dp) :: y0(state_count), seed(state_count), cotangent(state_count)
        real(dp) :: perturbed(state_count), final_state(state_count), sink
        real(dp), allocatable :: product(:)
        integer(int64) :: iteration, start, finish, rate
        integer :: direction
        real(dp), parameter :: difference_step = 1.0e-5_dp

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1_int64, reps
            call fill_initial_state(y0)
            y0(1) = y0(1) + 1.0e-9_dp*real(mod(iteration, 101_int64), dp)
            call integrate_primal(y0, solution)
            if (name == "primal") then
                sink = sink + solution%y(1, solution%nsteps + 1)
            else if (name == "forward") then
                do direction = 1, count
                    call fill_seed(direction, seed)
                    call apply_jvp(solution, seed, product)
                    sink = sink + sum(product)
                end do
            else if (name == "reverse") then
                do direction = 1, count
                    call fill_cotangent(direction, cotangent)
                    call apply_vjp(solution, cotangent, product)
                    sink = sink + sum(product)
                end do
            else
                do direction = 1, count
                    call fill_seed(direction, seed)
                    perturbed = y0 + difference_step*seed
                    call frozen_replay(perturbed, solution, final_state)
                    sink = sink + sum(final_state)
                    perturbed = y0 - difference_step*seed
                    call frozen_replay(perturbed, solution, final_state)
                    sink = sink - sum(final_state)
                end do
            end if
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
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

end program test_ode_discrete_sensitivity
