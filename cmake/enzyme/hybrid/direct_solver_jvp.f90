module direct_solver_jvp_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr, c_int
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    public :: analytical_jvp, autodiff_jvp, diagnostic_jvp, solve_direct

    interface
        function enzyme_fwddiff(f, a, da, b, db, component) result(dx) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr, c_int
            type(c_funptr), value :: f
            real(c_double), intent(in) :: a(*), da(*), b(*), db(*)
            integer(c_int), value :: component
            real(c_double) :: dx
        end function enzyme_fwddiff
    end interface

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
            dx(component) = enzyme_fwddiff(c_funloc(solve_component), &
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
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    use direct_solver_jvp_kernel, only: analytical_jvp, autodiff_jvp, &
        diagnostic_jvp, solver_size
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(c_double) :: a(solver_size, solver_size), b(solver_size)
    real(c_double) :: da(solver_size, solver_size, 16), db(solver_size, 16)
    real(c_double) :: x(solver_size), dx(solver_size), reference(solver_size)
    integer :: direction, direction_count
    integer(int64) :: iterations
    character(32) :: action, candidate

    call initialize_inputs()
    call get_environment_variable("FORTNUM_DIRECT_SOLVER_ACTION", action)
    call get_environment_variable("FORTNUM_DIRECT_SOLVER_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_DIRECT_SOLVER_DIRECTIONS", direction_count, 16)
    call read_int64_env("FORTNUM_DIRECT_SOLVER_ITERATIONS", iterations, 20000_int64)
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), direction_count, iterations)
    else if (trim(action) == "--peak-rss") then
        call benchmark_candidate(trim(candidate), direction_count, iterations)
        write (*, "(a,i0)") "peak_rss_bytes=", peak_rss_bytes()
    else
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
    end if

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

    subroutine benchmark_candidate(name, directions, count)
        character(*), intent(in) :: name
        integer, intent(in) :: directions
        integer(int64), intent(in) :: count
        integer :: active_direction
        integer(int64) :: iteration, start, finish, rate
        real(c_double) :: elapsed_ns, sink

        if (directions < 1 .or. directions > 16) error stop "directions must be 1..16"
        if (name /= "analytical" .and. name /= "autodiff" .and. &
            name /= "diagnostic") then
            error stop "candidate must be analytical, autodiff, or diagnostic"
        end if
        sink = 0.0_c_double
        call system_clock(start, rate)
        do iteration = 1, count
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024_int64), c_double)
            do active_direction = 1, directions
                select case (name)
                case ("analytical")
                    call analytical_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                case ("autodiff")
                    call autodiff_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                case default
                    call diagnostic_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), x, dx)
                end select
                sink = sink + x(1) + dx(1)
            end do
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = real(finish - start, c_double)*1.0e9_c_double / &
            (real(rate, c_double)*real(count, c_double))
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " directions=", directions, " iterations=", count, &
            " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

    subroutine read_integer_env(name, value, default_value)
        character(*), intent(in) :: name
        integer, intent(out) :: value
        integer, intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid integer environment value"
    end subroutine read_integer_env

    subroutine read_int64_env(name, value, default_value)
        character(*), intent(in) :: name
        integer(int64), intent(out) :: value
        integer(int64), intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid int64 environment value"
    end subroutine read_int64_env

end program enzyme_direct_solver_jvp
