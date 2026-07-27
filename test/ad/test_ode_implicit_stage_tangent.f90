program test_ode_implicit_stage_tangent
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: lu_solve
    use fortnum_ode, only: ode_implicit_stage_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    character(32) :: action, candidate
    integer :: state_count, direction_count

    call get_environment_variable("FORTNUM_IMPLICIT_STAGE_ACTION", action)
    call get_environment_variable("FORTNUM_IMPLICIT_STAGE_CANDIDATE", candidate)
    if (trim(action) == "--benchmark") then
        call read_dimensions(state_count, direction_count)
        call run_benchmark(trim(candidate), state_count, direction_count)
    else if (trim(action) == "--perf") then
        call read_dimensions(state_count, direction_count)
        call run_perf(trim(candidate), state_count, direction_count)
    else if (trim(action) == "--peak-rss") then
        call read_dimensions(state_count, direction_count)
        call run_peak_rss(trim(candidate), state_count, direction_count)
    else
        call validate_product()
    end if

contains

    subroutine validate_product()
        integer, parameter :: n = 3, ndirection = 4
        real(dp), parameter :: alpha = 0.2_dp
        real(dp) :: jacobian(n, n), base_tangent(n, ndirection)
        real(dp) :: parameter_jvp(n, ndirection), expected(n, ndirection)
        real(dp), allocatable :: tangent(:,:)
        type(fortnum_status_t) :: status
        integer :: i, direction

        jacobian = 0.0_dp
        do i = 1, n
            jacobian(i, i) = -real(i, dp)
        end do
        do direction = 1, ndirection
            do i = 1, n
                base_tangent(i, direction) = &
                    0.1_dp*real(i + 2*direction, dp)
                parameter_jvp(i, direction) = &
                    -0.05_dp*real(2*i - direction, dp)
                expected(i, direction) = &
                    (base_tangent(i, direction) + &
                    alpha*parameter_jvp(i, direction)) / &
                    (1.0_dp + alpha*real(i, dp))
            end do
        end do

        call ode_implicit_stage_jvp(alpha, jacobian, base_tangent, &
            parameter_jvp, tangent, status)
        if (.not. status_ok(status)) error stop "implicit-stage JVP failed"
        if (maxval(abs(tangent - expected)) > 32.0_dp*epsilon(1.0_dp)) then
            print *, "implicit-stage exact-oracle mismatch", &
                maxval(abs(tangent - expected))
            error stop 1
        end if

        jacobian = 0.0_dp
        jacobian(1, 1) = 1.0_dp/alpha
        call ode_implicit_stage_jvp(alpha, jacobian, base_tangent, &
            parameter_jvp, tangent, status)
        if (status_ok(status)) error stop "singular stage accepted"
        print *, "PASS analytical implicit-stage tangent"
    end subroutine validate_product

    subroutine read_dimensions(n, ndirection)
        integer, intent(out) :: n, ndirection
        character(32) :: text

        call get_environment_variable("FORTNUM_IMPLICIT_STAGE_STATES", text)
        read (text, *) n
        call get_environment_variable("FORTNUM_IMPLICIT_STAGE_DIRECTIONS", text)
        read (text, *) ndirection
        if (n < 1 .or. n > 16) error stop "state count must be 1..16"
        if (ndirection < 1) error stop "direction count must be positive"
    end subroutine read_dimensions

    subroutine run_benchmark(name, n, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 20000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_name(name)
        do sample = 1, 3
            call time_candidate(name, n, ndirection, reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, n, ndirection, reps, elapsed(sample))
        end do
        call report(name, n, ndirection, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name, n, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, n, ndirection, 100000_int64, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine run_perf(name, n, ndirection)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
        real(dp) :: elapsed

        call validate_name(name)
        call time_candidate(name, n, ndirection, 1000000_int64, elapsed)
        write (*, "(f0.4)") elapsed
    end subroutine run_perf

    subroutine validate_name(name)
        character(*), intent(in) :: name

        if (name /= "analytical" .and. name /= "diagnostic") then
            error stop "candidate must be analytical or diagnostic"
        end if
    end subroutine validate_name

    subroutine time_candidate(name, n, ndirection, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        real(dp), parameter :: alpha = 0.03_dp
        real(dp) :: jacobian(n, n), base(n), parameter(n)
        real(dp) :: base_tangent(n, ndirection)
        real(dp) :: parameter_jvp(n, ndirection)
        real(dp), allocatable :: stage(:), tangent(:,:)
        type(fortnum_status_t) :: status
        integer(int64) :: iteration, start, finish, rate
        integer :: i, direction
        real(dp) :: sink

        jacobian = 0.0_dp
        do i = 1, n
            jacobian(i, i) = -real(i, dp)
        end do
        do i = 1, n
            base(i) = 0.2_dp + 0.01_dp*real(i, dp)
            parameter(i) = 0.1_dp - 0.002_dp*real(i, dp)
        end do
        do direction = 1, ndirection
            do i = 1, n
                base_tangent(i, direction) = &
                    0.01_dp*real(mod(i + direction, 7) - 3, dp)
                parameter_jvp(i, direction) = &
                    0.02_dp*real(mod(2*i - direction, 9) - 4, dp)
            end do
        end do

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1_int64, reps
            base(1) = 0.2_dp + &
                1.0e-8_dp*real(mod(iteration, 101_int64), dp)
            call solve_stage(alpha, jacobian, base, parameter, stage)
            if (name == "analytical") then
                call ode_implicit_stage_jvp(alpha, jacobian, base_tangent, &
                    parameter_jvp, tangent, status)
                if (.not. status_ok(status)) error stop "benchmark JVP failed"
            else
                call diagnostic_tangents(alpha, jacobian, base, parameter, &
                    base_tangent, parameter_jvp, tangent)
            end if
            sink = sink + stage(1) + sum(tangent)
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = 1.0e9_dp*real(finish - start, dp) / &
            (real(rate, dp)*real(reps, dp))
    end subroutine time_candidate

    subroutine solve_stage(alpha, jacobian, base, parameter, stage)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: jacobian(:,:), base(:), parameter(:)
        real(dp), allocatable, intent(out) :: stage(:)
        real(dp) :: factor(size(base), size(base))
        integer :: i, info

        factor = -alpha*jacobian
        do i = 1, size(base)
            factor(i, i) = factor(i, i) + 1.0_dp
        end do
        stage = base + alpha*parameter
        call lu_solve(size(base), factor, stage, info)
        if (info /= 0) error stop "primal stage solve failed"
    end subroutine solve_stage

    subroutine diagnostic_tangents(alpha, jacobian, base, parameter, &
            base_tangent, parameter_jvp, tangent)
        real(dp), intent(in) :: alpha
        real(dp), intent(in) :: jacobian(:,:), base(:), parameter(:)
        real(dp), intent(in) :: base_tangent(:,:), parameter_jvp(:,:)
        real(dp), allocatable, intent(out) :: tangent(:,:)
        real(dp), parameter :: step = 1.0e-5_dp
        real(dp), allocatable :: plus(:), minus(:)
        real(dp) :: base_work(size(base)), parameter_work(size(parameter))
        integer :: direction

        allocate(tangent(size(base), size(base_tangent, 2)))
        do direction = 1, size(base_tangent, 2)
            base_work = base + step*base_tangent(:, direction)
            parameter_work = parameter + step*parameter_jvp(:, direction)
            call solve_stage(alpha, jacobian, base_work, parameter_work, plus)
            base_work = base - step*base_tangent(:, direction)
            parameter_work = parameter - step*parameter_jvp(:, direction)
            call solve_stage(alpha, jacobian, base_work, parameter_work, minus)
            tangent(:, direction) = (plus - minus)/(2.0_dp*step)
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

    subroutine report(name, n, ndirection, values, reps)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
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
        write (*, "(a,3(',',i0),2(',',f12.4))") name, n, ndirection, reps, &
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

end program test_ode_implicit_stage_tangent
