program test_ode_implicit_stage_adjoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: lu_solve
    use fortnum_ode, only: ode_implicit_stage_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    type :: parameter_context_t
        real(dp), allocatable :: derivative(:,:)
    end type parameter_context_t

    character(32) :: action, candidate
    integer :: parameter_count, objective_count

    call get_environment_variable("FORTNUM_IMPLICIT_ADJOINT_ACTION", action)
    call get_environment_variable("FORTNUM_IMPLICIT_ADJOINT_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call read_dimensions(parameter_count, objective_count)
        call run_benchmark(trim(candidate), parameter_count, objective_count)
    else if (trim(action) == "--perf") then
        call read_dimensions(parameter_count, objective_count)
        call run_perf(trim(candidate), parameter_count, objective_count)
    else if (trim(action) == "--peak-rss") then
        call read_dimensions(parameter_count, objective_count)
        call run_peak_rss(trim(candidate), parameter_count, objective_count)
    else
        call validate_product()
    end if

contains

    subroutine validate_product()
        integer, parameter :: n = 3, nparameter = 5, nobjective = 4
        real(dp), parameter :: alpha = 0.2_dp
        real(dp) :: jacobian(n, n), stage(n), cotangent(n, nobjective)
        real(dp) :: expected_base(n, nobjective)
        real(dp) :: expected_parameter(nparameter, nobjective)
        real(dp), allocatable :: base_bar(:,:), parameter_bar(:,:)
        type(parameter_context_t) :: context
        type(fortnum_status_t) :: status
        integer :: i, parameter, objective

        allocate(context%derivative(n, nparameter))
        jacobian = 0.0_dp
        do i = 1, n
            jacobian(i, i) = -real(i, dp)
            stage(i) = 0.3_dp + 0.1_dp*real(i, dp)
            do parameter = 1, nparameter
                context%derivative(i, parameter) = &
                    0.02_dp*real(i + 2*parameter, dp)
            end do
            do objective = 1, nobjective
                cotangent(i, objective) = &
                    0.1_dp*real(2*i - objective, dp)
            end do
        end do

        do objective = 1, nobjective
            do i = 1, n
                expected_base(i, objective) = cotangent(i, objective) / &
                    (1.0_dp + alpha*real(i, dp))
            end do
            do parameter = 1, nparameter
                expected_parameter(parameter, objective) = alpha*sum( &
                    context%derivative(:, parameter) * &
                    expected_base(:, objective))
            end do
        end do

        call ode_implicit_stage_vjp(0.0_dp, stage, alpha, jacobian, &
            parameter_vjp, nparameter, cotangent, base_bar, parameter_bar, &
            status, context)
        if (.not. status_ok(status)) error stop "implicit-stage VJP failed"
        if (maxval(abs(base_bar - expected_base)) > &
                32.0_dp*epsilon(1.0_dp)) then
            error stop "implicit-stage base VJP exact-oracle mismatch"
        end if
        if (maxval(abs(parameter_bar - expected_parameter)) > &
                64.0_dp*epsilon(1.0_dp)) then
            error stop "implicit-stage parameter VJP exact-oracle mismatch"
        end if

        jacobian = 0.0_dp
        jacobian(1, 1) = 1.0_dp/alpha
        call ode_implicit_stage_vjp(0.0_dp, stage, alpha, jacobian, &
            parameter_vjp, nparameter, cotangent, base_bar, parameter_bar, &
            status, context)
        if (status_ok(status)) error stop "singular adjoint stage accepted"
        print *, "PASS analytical implicit-stage adjoint"
    end subroutine validate_product

    subroutine parameter_vjp(t, stage, rhs_cotangent, gradient, ctx)
        real(dp), intent(in) :: t
        real(dp), intent(in) :: stage(:), rhs_cotangent(:)
        real(dp), intent(inout) :: gradient(:)
        class(*), intent(in), optional :: ctx
        integer :: i, parameter

        associate (unused_stage => stage, unused_t => t); end associate
            if (.not. present(ctx)) error stop "parameter context missing"
            select type (context => ctx)
            type is (parameter_context_t)
                do parameter = 1, size(gradient)
                    do i = 1, size(rhs_cotangent)
                        gradient(parameter) = gradient(parameter) + &
                            context%derivative(i, parameter)*rhs_cotangent(i)
                    end do
                end do
            class default
                error stop "wrong parameter context"
            end select
    end subroutine parameter_vjp

    subroutine read_dimensions(nparameter, nobjective)
        integer, intent(out) :: nparameter, nobjective
        character(32) :: text

        call get_environment_variable("FORTNUM_IMPLICIT_ADJOINT_PARAMETERS", text)
        read (text, *) nparameter
        call get_environment_variable("FORTNUM_IMPLICIT_ADJOINT_OBJECTIVES", text)
        read (text, *) nobjective
        if (nparameter < 1) error stop "parameter count must be positive"
        if (nobjective < 1) error stop "objective count must be positive"
    end subroutine read_dimensions

    subroutine run_benchmark(name, nparameter, nobjective)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, nobjective
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 20000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, nparameter, nobjective, &
                reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, nparameter, nobjective, reps, &
                elapsed(sample))
        end do
        call report(name, nparameter, nobjective, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_perf(name, nparameter, nobjective)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, nobjective
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, nparameter, nobjective, 1000000_int64, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine run_peak_rss(name, nparameter, nobjective)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, nobjective
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, nparameter, nobjective, 100000_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "analytical" .and. name /= "diagnostic") then
            error stop "candidate must be analytical or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, nparameter, nobjective, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, nobjective
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer, parameter :: n = 4
        real(dp), parameter :: alpha = 0.03_dp
        real(dp) :: jacobian(n, n), base(n), parameters(nparameter)
        real(dp) :: cotangent(n, nobjective)
        real(dp), allocatable :: stage(:), base_bar(:,:), parameter_bar(:,:)
        type(parameter_context_t) :: context
        type(fortnum_status_t) :: status
        integer(int64) :: iteration, start, finish, rate
        integer :: i, parameter, objective
        real(dp) :: sink

        allocate(context%derivative(n, nparameter))
        jacobian = 0.0_dp
        do i = 1, n
            jacobian(i, i) = -real(i, dp)
            base(i) = 0.2_dp + 0.01_dp*real(i, dp)
            do parameter = 1, nparameter
                context%derivative(i, parameter) = &
                    0.01_dp*real(mod(i + 2*parameter, 11) - 5, dp)
            end do
            do objective = 1, nobjective
                cotangent(i, objective) = &
                    0.02_dp*real(mod(2*i - objective, 9) - 4, dp)
            end do
        end do
        do parameter = 1, nparameter
            parameters(parameter) = 0.05_dp - &
                0.001_dp*real(parameter, dp)
        end do

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1_int64, reps
            base(1) = 0.2_dp + &
                1.0e-8_dp*real(mod(iteration, 101_int64), dp)
            call solve_stage(alpha, jacobian, context%derivative, base, &
                parameters, stage)
            if (name == "analytical") then
                call ode_implicit_stage_vjp(0.0_dp, stage, alpha, jacobian, &
                    parameter_vjp, nparameter, cotangent, base_bar, &
                    parameter_bar, status, context)
                if (.not. status_ok(status)) error stop "benchmark VJP failed"
            else
                call diagnostic_adjoints(alpha, jacobian, context%derivative, &
                    base, parameters, cotangent, base_bar, parameter_bar)
            end if
            sink = sink + stage(1) + sum(base_bar) + sum(parameter_bar)
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = 1.0e9_dp*real(finish - start, dp) / &
            (real(rate, dp)*real(reps, dp))
    end subroutine time_candidate

    subroutine solve_stage(alpha, jacobian, parameter_jacobian, base, &
            parameters, stage)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: jacobian(:,:), parameter_jacobian(:,:)
        real(dp), intent(in) :: base(:), parameters(:)
        real(dp), allocatable, intent(out) :: stage(:)
        real(dp) :: factor(size(base), size(base))
        integer :: i, parameter, info

        factor = -alpha*jacobian
        do i = 1, size(base)
            factor(i, i) = factor(i, i) + 1.0_dp
        end do
        allocate(stage(size(base)))
        stage = base
        do parameter = 1, size(parameters)
            do i = 1, size(base)
                stage(i) = stage(i) + alpha* &
                    parameter_jacobian(i, parameter)*parameters(parameter)
            end do
        end do
        call lu_solve(size(base), factor, stage, info)
        if (info /= 0) error stop "primal stage solve failed"
    end subroutine solve_stage

    subroutine diagnostic_adjoints(alpha, jacobian, parameter_jacobian, base, &
            parameters, cotangent, base_bar, parameter_bar)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: jacobian(:,:), parameter_jacobian(:,:)
        real(dp), intent(in) :: base(:), parameters(:), cotangent(:,:)
        real(dp), allocatable, intent(out) :: base_bar(:,:), parameter_bar(:,:)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp), allocatable :: plus(:), minus(:)
        real(dp) :: base_work(size(base)), parameter_work(size(parameters))
        integer :: input, objective

        allocate(base_bar(size(base), size(cotangent, 2)))
        allocate(parameter_bar(size(parameters), size(cotangent, 2)))
        do input = 1, size(base)
            base_work = base
            base_work(input) = base_work(input) + step
            call solve_stage(alpha, jacobian, parameter_jacobian, base_work, &
                parameters, plus)
            base_work(input) = base_work(input) - 2.0_dp*step
            call solve_stage(alpha, jacobian, parameter_jacobian, base_work, &
                parameters, minus)
            do objective = 1, size(cotangent, 2)
                base_bar(input, objective) = dot_product( &
                    cotangent(:, objective), plus - minus)/(2.0_dp*step)
            end do
        end do
        do input = 1, size(parameters)
            parameter_work = parameters
            parameter_work(input) = parameter_work(input) + step
            call solve_stage(alpha, jacobian, parameter_jacobian, base, &
                parameter_work, plus)
            parameter_work(input) = parameter_work(input) - 2.0_dp*step
            call solve_stage(alpha, jacobian, parameter_jacobian, base, &
                parameter_work, minus)
            do objective = 1, size(cotangent, 2)
                parameter_bar(input, objective) = dot_product( &
                    cotangent(:, objective), plus - minus)/(2.0_dp*step)
            end do
        end do
    end subroutine diagnostic_adjoints

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

    subroutine report(name, nparameter, nobjective, values, reps)
        character(*), intent(in) :: name
        integer, intent(in) :: nparameter, nobjective
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
            nobjective, reps, median, mad
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

end program test_ode_implicit_stage_adjoint
