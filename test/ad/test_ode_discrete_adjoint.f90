module test_ode_discrete_adjoint_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_checkpoint_t, ode_integrate, ode_integrate_jvp, ode_integrate_vjp, &
        ode_build_checkpoints, ode_integrate_vjp_checkpointed
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    real(dp), parameter :: system_matrix(2, 2) = reshape([ &
        -0.3_dp, 0.7_dp, -1.1_dp, 0.2_dp], [2, 2])
    real(dp), parameter :: final_time = 1.3_dp

    public :: reverse_vjp, checkpointed_vjp, forward_reconstruction_vjp
    public :: diagnostic_vjp, exact_vjp

contains

    function reverse_vjp(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp), allocatable :: result(:)

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "adjoint primal failed"
        call ode_integrate_vjp(problem, adjoint_rhs, cotangent, solution, &
            result, status)
        if (.not. status_ok(status)) error stop "discrete adjoint failed"
        vjp = result
    end function reverse_vjp

    function checkpointed_vjp(y0, cotangent, stride) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        integer, intent(in) :: stride
        real(dp) :: vjp(2)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(ode_checkpoint_t) :: checkpoints
        type(fortnum_status_t) :: status
        real(dp), allocatable :: result(:)

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "checkpoint primal failed"
        call ode_build_checkpoints(solution, stride, checkpoints, status)
        if (.not. status_ok(status)) error stop "checkpoint build failed"
        deallocate(solution%y)
        call ode_integrate_vjp_checkpointed(problem, adjoint_rhs, cotangent, &
            checkpoints, result, status)
        if (.not. status_ok(status)) error stop "checkpoint adjoint failed"
        vjp = result
    end function checkpointed_vjp

    function forward_reconstruction_vjp(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: vjp(2), seed(2)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp), allocatable :: tangent(:)
        integer :: input

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "forward primal failed"
        do input = 1, 2
            seed = 0.0_dp
            seed(input) = 1.0_dp
            call ode_integrate_jvp(problem, tangent_rhs, seed, solution, &
                tangent, status)
            if (.not. status_ok(status)) error stop "forward tangent failed"
            vjp(input) = dot_product(cotangent, tangent)
        end do
    end function forward_reconstruction_vjp

    function diagnostic_vjp(y0, cotangent) result(vjp)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: vjp(2), y_plus(2), y_minus(2), base
        integer :: input

        base = objective(y0, cotangent)
        do input = 1, 2
            y_plus = y0
            y_minus = y0
            y_plus(input) = y_plus(input) + h
            y_minus(input) = y_minus(input) - h
            vjp(input) = (objective(y_plus, cotangent) - &
                objective(y_minus, cotangent))/(2.0_dp*h)
        end do
        if (base == huge(base)) error stop "unreachable adjoint objective"
    end function diagnostic_vjp

    function objective(y0, cotangent) result(value)
        real(dp), intent(in) :: y0(2), cotangent(2)
        real(dp) :: value
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        call make_problem(y0, problem)
        call ode_integrate(problem, workspace, solution, status)
        if (.not. status_ok(status)) error stop "diagnostic primal failed"
        value = dot_product(cotangent, &
            solution%y(:, solution%nsteps + 1))
    end function objective

    subroutine make_problem(y0, problem)
        real(dp), intent(in) :: y0(2)
        type(ode_problem_t), intent(out) :: problem

        problem%rhs => primal_rhs
        problem%t0 = 0.0_dp
        problem%t1 = final_time
        problem%y0 = y0
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
    end subroutine make_problem

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

    pure function exact_vjp(cotangent) result(vjp)
        real(dp), intent(in) :: cotangent(2)
        real(dp) :: vjp(2), trace, determinant, discriminant
        real(dp) :: mean, frequency, scale, c0, c1

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

end module test_ode_discrete_adjoint_kernel

program test_ode_discrete_adjoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use test_ode_discrete_adjoint_kernel, only: reverse_vjp, checkpointed_vjp, &
        forward_reconstruction_vjp, diagnostic_vjp, exact_vjp
    implicit none

    character(32) :: action, candidate

    call get_environment_variable("FORTNUM_ODE_ADJOINT_ACTION", action)
    call get_environment_variable("FORTNUM_ODE_ADJOINT_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call run_benchmark(trim(candidate))
    else if (trim(action) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        call validate_candidates()
    end if

contains

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

    subroutine validate_candidates()
        real(dp), parameter :: y0(2) = [1.5_dp, -0.4_dp]
        real(dp), parameter :: cotangent(2) = [0.6_dp, -1.2_dp]
        real(dp) :: reference(2), candidate_value(2), errors(5)

        reference = exact_vjp(cotangent)
        candidate_value = reverse_vjp(y0, cotangent)
        errors(1) = maxval(abs(candidate_value - reference))
        candidate_value = forward_reconstruction_vjp(y0, cotangent)
        errors(2) = maxval(abs(candidate_value - reference))
        candidate_value = diagnostic_vjp(y0, cotangent)
        errors(3) = maxval(abs(candidate_value - reference))
        candidate_value = checkpointed_vjp(y0, cotangent, 4)
        errors(4) = maxval(abs(candidate_value - reference))
        candidate_value = checkpointed_vjp(y0, cotangent, 16)
        errors(5) = maxval(abs(candidate_value - reference))
        if (max(errors(1), errors(2), errors(4), errors(5)) > 2.0e-8_dp .or. &
            errors(3) > 2.0e-7_dp) then
            print *, "ODE discrete-adjoint mismatch", errors
            error stop 1
        end if
        print *, "PASS ODE discrete adjoint", errors
    end subroutine validate_candidates

    subroutine run_benchmark(name)
        character(*), intent(in) :: name
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 1000_int64
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
        integer(int64), parameter :: reps = 2000_int64
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if ((name /= "reverse") .and. (name /= "checkpoint4") .and. &
            (name /= "checkpoint16") .and. (name /= "forward") .and. &
            (name /= "diagnostic")) then
            error stop "unknown ODE adjoint candidate"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        real(dp) :: y0(2), cotangent(2), result(2), sink

        cotangent = [0.6_dp, -1.2_dp]
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            y0(1) = 1.5_dp + &
                1.0e-6_dp*real(mod(iteration, 101_int64), dp)
            y0(2) = -0.4_dp
            select case (name)
            case ("reverse")
                result = reverse_vjp(y0, cotangent)
            case ("checkpoint4")
                result = checkpointed_vjp(y0, cotangent, 4)
            case ("checkpoint16")
                result = checkpointed_vjp(y0, cotangent, 16)
            case ("forward")
                result = forward_reconstruction_vjp(y0, cotangent)
            case ("diagnostic")
                result = diagnostic_vjp(y0, cotangent)
            end select
            sink = sink + sum(result)
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

end program test_ode_discrete_adjoint
