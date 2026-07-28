module iterative_solver_jvp_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use fortnum_generated_enzyme_iterative_solver_component, only: &
        fortnum_enzyme_iterative_solver_component_jvp
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    real(c_double), parameter :: relaxation = 0.1_c_double
    public :: analytical_jvp, autodiff_jvp, diagnostic_jvp

contains

    pure function iteration_component(a, b, iterations, component) result(value) &
            bind(c, name="fortnum_iterative_solve_component")
        real(c_double), intent(in) :: a(solver_size, solver_size), b(solver_size)
        integer(c_int), value :: iterations, component
        real(c_double) :: value, x(solver_size)

        call iterate(a, b, iterations, x)
        value = x(component)
    end function iteration_component

    pure subroutine iterate(a, b, iterations, x)
        real(c_double), intent(in) :: a(solver_size, solver_size), b(solver_size)
        integer(c_int), intent(in) :: iterations
        real(c_double), intent(out) :: x(solver_size)
        real(c_double) :: next_x(solver_size), residual
        integer :: i, j, iteration

        x = 0.0_c_double
        do iteration = 1, iterations
            do i = 1, solver_size
                residual = b(i)
                do j = 1, solver_size
                    residual = residual - a(i, j)*x(j)
                end do
                next_x(i) = x(i) + relaxation*residual
            end do
            x = next_x
        end do
    end subroutine iterate

    pure subroutine analytical_jvp(a, da, b, db, iterations, x, dx)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        integer(c_int), intent(in) :: iterations
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        real(c_double) :: next_x(solver_size), next_dx(solver_size)
        real(c_double) :: residual, tangent_residual
        integer :: i, j, iteration

        x = 0.0_c_double
        dx = 0.0_c_double
        do iteration = 1, iterations
            do i = 1, solver_size
                residual = b(i)
                tangent_residual = db(i)
                do j = 1, solver_size
                    residual = residual - a(i, j)*x(j)
                    tangent_residual = tangent_residual - da(i, j)*x(j) &
                        - a(i, j)*dx(j)
                end do
                next_x(i) = x(i) + relaxation*residual
                next_dx(i) = dx(i) + relaxation*tangent_residual
            end do
            x = next_x
            dx = next_dx
        end do
    end subroutine analytical_jvp

    subroutine autodiff_jvp(a, da, b, db, iterations, x, dx)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        integer(c_int), intent(in) :: iterations
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        integer :: component

        call iterate(a, b, iterations, x)
        do component = 1, solver_size
            dx(component) = fortnum_enzyme_iterative_solver_component_jvp( &
                a, da, b, db, iterations, int(component, c_int))
        end do
    end subroutine autodiff_jvp

    pure subroutine diagnostic_jvp(a, da, b, db, iterations, x, dx)
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        integer(c_int), intent(in) :: iterations
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        real(c_double) :: plus_a(solver_size, solver_size)
        real(c_double) :: minus_a(solver_size, solver_size)
        real(c_double) :: plus_b(solver_size), minus_b(solver_size)
        real(c_double) :: plus_x(solver_size), minus_x(solver_size)

        plus_a = a + h*da
        minus_a = a - h*da
        plus_b = b + h*db
        minus_b = b - h*db
        call iterate(a, b, iterations, x)
        call iterate(plus_a, plus_b, iterations, plus_x)
        call iterate(minus_a, minus_b, iterations, minus_x)
        dx = (plus_x - minus_x)/(2.0_c_double*h)
    end subroutine diagnostic_jvp

end module iterative_solver_jvp_kernel

program enzyme_iterative_solver_jvp
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use iterative_solver_jvp_kernel, only: analytical_jvp, autodiff_jvp, &
        diagnostic_jvp, solver_size
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none

    real(c_double) :: a(solver_size, solver_size), b(solver_size)
    real(c_double) :: da(solver_size, solver_size, 16), db(solver_size, 16)
    real(c_double) :: x(solver_size), dx(solver_size), reference(solver_size)
    real(c_double) :: samples(fixture_sample_count), sink
    integer :: candidate_kind, direction, direction_count, trace_iterations
    integer :: workloads
    character(32) :: action, candidate
    logical :: valid

    call initialize_inputs()
    call read_fixture_environment("FORTNUM_ITERATIVE_ACTION", "validate", action)
    if (trim(action) == "validate") then
        do trace_iterations = 4, 64, 4
            do direction = 1, 16
                call analytical_jvp(a, da(:, :, direction), b, db(:, direction), &
                    int(trace_iterations, c_int), x, reference)
                call autodiff_jvp(a, da(:, :, direction), b, db(:, direction), &
                    int(trace_iterations, c_int), x, dx)
                if (maxval(abs(dx - reference)) > 3.0e-12_c_double) then
                    error stop "autodiff iterative-solver JVP mismatch"
                end if
                call diagnostic_jvp(a, da(:, :, direction), b, db(:, direction), &
                    int(trace_iterations, c_int), x, dx)
                if (maxval(abs(dx - reference)) > 3.0e-10_c_double) then
                    error stop "diagnostic iterative-solver JVP mismatch"
                end if
            end do
        end do
        write (*, "(a)") "PASS"
        stop
    end if

    call read_fixture_environment("FORTNUM_ITERATIVE_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_integer("FORTNUM_ITERATIVE_DIRECTIONS", 16, &
        direction_count, valid)
    if (.not. valid) error stop "invalid iterative-solver direction count"
    call read_fixture_integer("FORTNUM_ITERATIVE_STEPS", 32, &
        trace_iterations, valid)
    if (.not. valid) error stop "invalid iterative-solver step count"
    call read_fixture_integer("FORTNUM_ITERATIVE_WORKLOADS", 5000, &
        workloads, valid)
    if (.not. valid) error stop "invalid iterative-solver workload count"
    call parse_configuration()

    if (trim(action) == "benchmark" .or. trim(action) == "--benchmark") then
        call benchmark_candidate()
    else if (trim(action) == "peak-rss" .or. trim(action) == "--peak-rss") then
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "iterative-solver benchmark produced NaN"

contains

    subroutine initialize_inputs()
        integer :: i, j, active_direction

        do j = 1, solver_size
            do i = 1, solver_size
                if (i == j) then
                    a(i, j) = 4.0_c_double + 0.5_c_double*real(i, c_double)
                else
                    a(i, j) = 0.1_c_double/real(i + j, c_double)
                end if
                do active_direction = 1, 16
                    da(i, j, active_direction) = 0.01_c_double*real( &
                        mod(3*i + 5*j + active_direction, 11) - 5, c_double)
                end do
            end do
            b(j) = 0.5_c_double*real(j, c_double)
            do active_direction = 1, 16
                db(j, active_direction) = 0.02_c_double*real( &
                    mod(7*j + active_direction, 9) - 4, c_double)
            end do
        end do
    end subroutine initialize_inputs

    subroutine parse_configuration()
        if (direction_count < 1 .or. direction_count > 16) then
            error stop "directions must be 1..16"
        end if
        if (trace_iterations < 1) error stop "steps must be positive"
        if (workloads < 1) error stop "workloads must be positive"
        select case (trim(candidate))
        case ("analytical")
            candidate_kind = 1
        case ("autodiff")
            candidate_kind = 2
        case ("diagnostic")
            candidate_kind = 3
        case default
            error stop "candidate must be analytical, autodiff, or diagnostic"
        end select
    end subroutine parse_configuration

    subroutine benchmark_candidate()
        real(c_double) :: median, mad
        character(96) :: name

        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        write (name, "(a,'_jvp_d',i0,'_s',i0)") trim(candidate), &
            direction_count, trace_iterations
        call write_fixture_result(trim(name), workloads, median, mad)
    end subroutine benchmark_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        integer :: active_direction
        integer :: workload
        real(c_double) :: elapsed_ns, local_sink

        local_sink = 0.0_c_double
        call timer%start()
        do workload = 1, workloads
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(workload, 1024), c_double)
            do active_direction = 1, direction_count
                select case (candidate_kind)
                case (1)
                    call analytical_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(trace_iterations, c_int), x, dx)
                case (2)
                    call autodiff_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(trace_iterations, c_int), x, dx)
                case default
                    call diagnostic_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(trace_iterations, c_int), x, dx)
                end select
                local_sink = local_sink + x(1) + dx(1)
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(workloads, c_double)
        sink = local_sink
    end function measure_candidate

end program enzyme_iterative_solver_jvp
