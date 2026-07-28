module direct_solver_jvp_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use fortnum_generated_enzyme_direct_solver_component, only: &
        fortnum_enzyme_direct_solver_component_jvp
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    public :: analytical_jvp, autodiff_jvp, diagnostic_jvp, solve_direct

contains

    pure function solve_component(a, b, component) result(value) &
            bind(c, name="fortnum_direct_solve_component")
        real(c_double), intent(in) :: a(solver_size, solver_size), b(solver_size)
        integer(c_int), value :: component
        real(c_double) :: value, x(solver_size)

        call solve_direct(a, b, x)
        value = x(component)
    end function solve_component

    pure subroutine solve_direct(a, b, x)
        real(c_double), intent(in) :: a(solver_size, solver_size), b(solver_size)
        real(c_double), intent(out) :: x(solver_size)
        real(c_double) :: factors(solver_size, solver_size), rhs(solver_size)
        real(c_double) :: multiplier, total
        integer :: i, j, k

        factors = a
        rhs = b
        do k = 1, solver_size - 1
            do i = k + 1, solver_size
                multiplier = factors(i, k)/factors(k, k)
                do j = k + 1, solver_size
                    factors(i, j) = factors(i, j) - multiplier*factors(k, j)
                end do
                rhs(i) = rhs(i) - multiplier*rhs(k)
            end do
        end do
        do i = solver_size, 1, -1
            total = rhs(i)
            do j = i + 1, solver_size
                total = total - factors(i, j)*x(j)
            end do
            x(i) = total/factors(i, i)
        end do
    end subroutine solve_direct

    pure subroutine analytical_jvp(a, da, b, db, x, dx)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        real(c_double) :: tangent_rhs(solver_size)
        integer :: i, j

        call solve_direct(a, b, x)
        tangent_rhs = db
        do j = 1, solver_size
            do i = 1, solver_size
                tangent_rhs(i) = tangent_rhs(i) - da(i, j)*x(j)
            end do
        end do
        call solve_direct(a, tangent_rhs, dx)
    end subroutine analytical_jvp

    subroutine autodiff_jvp(a, da, b, db, x, dx)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        integer :: component

        call solve_direct(a, b, x)
        do component = 1, solver_size
            dx(component) = fortnum_enzyme_direct_solver_component_jvp( &
                a, da, b, db, int(component, c_int))
        end do
    end subroutine autodiff_jvp

    pure subroutine diagnostic_jvp(a, da, b, db, x, dx)
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: da(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), db(solver_size)
        real(c_double), intent(out) :: x(solver_size), dx(solver_size)
        real(c_double) :: plus_a(solver_size, solver_size)
        real(c_double) :: minus_a(solver_size, solver_size)
        real(c_double) :: plus_b(solver_size), minus_b(solver_size)
        real(c_double) :: plus_x(solver_size), minus_x(solver_size)

        call solve_direct(a, b, x)
        plus_a = a + h*da
        minus_a = a - h*da
        plus_b = b + h*db
        minus_b = b - h*db
        call solve_direct(plus_a, plus_b, plus_x)
        call solve_direct(minus_a, minus_b, minus_x)
        dx = (plus_x - minus_x)/(2.0_c_double*h)
    end subroutine diagnostic_jvp

end module direct_solver_jvp_kernel

program enzyme_direct_solver_jvp
    use, intrinsic :: iso_c_binding, only: c_double
    use direct_solver_jvp_kernel, only: analytical_jvp, autodiff_jvp, &
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
    integer :: candidate_kind, direction, direction_count, iterations
    character(32) :: action, candidate
    logical :: valid

    call initialize_inputs()
    call read_fixture_environment("FORTNUM_DIRECT_SOLVER_ACTION", &
        "validate", action)
    if (trim(action) == "validate") then
        do direction = 1, 16
            call analytical_jvp(a, da(:, :, direction), b, db(:, direction), &
                x, reference)
            call autodiff_jvp(a, da(:, :, direction), b, db(:, direction), x, dx)
            if (maxval(abs(dx - reference)) > 2.0e-12_c_double) then
                error stop "autodiff direct-solver JVP mismatch"
            end if
            call diagnostic_jvp(a, da(:, :, direction), b, db(:, direction), x, dx)
            if (maxval(abs(dx - reference)) > 2.0e-10_c_double) then
                error stop "diagnostic direct-solver JVP mismatch"
            end if
        end do
        write (*, "(a)") "PASS"
        stop
    end if

    call read_fixture_environment("FORTNUM_DIRECT_SOLVER_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_integer("FORTNUM_DIRECT_SOLVER_DIRECTIONS", 16, &
        direction_count, valid)
    if (.not. valid) error stop "invalid direct-solver direction count"
    call read_fixture_integer("FORTNUM_DIRECT_SOLVER_ITERATIONS", 20000, &
        iterations, valid)
    if (.not. valid) error stop "invalid direct-solver iteration count"
    call parse_configuration()

    if (trim(action) == "benchmark" .or. trim(action) == "--benchmark") then
        call benchmark_candidate()
    else if (trim(action) == "peak-rss" .or. trim(action) == "--peak-rss") then
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "direct-solver JVP benchmark produced NaN"

contains

    subroutine initialize_inputs()
        integer :: i, j, active_direction

        do j = 1, solver_size
            do i = 1, solver_size
                if (i == j) then
                    a(i, j) = 4.0_c_double + real(i, c_double)
                else
                    a(i, j) = 0.2_c_double/real(i + j, c_double)
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
        if (iterations < 1) error stop "iterations must be positive"
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
        write (name, "(a,'_jvp_d',i0)") trim(candidate), direction_count
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine benchmark_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        integer :: active_direction
        integer :: iteration
        real(c_double) :: elapsed_ns, local_sink

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024), c_double)
            do active_direction = 1, direction_count
                select case (candidate_kind)
                case (1)
                    call analytical_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                case (2)
                    call autodiff_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                case default
                    call diagnostic_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                end select
                local_sink = local_sink + x(1) + dx(1)
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

end program enzyme_direct_solver_jvp
