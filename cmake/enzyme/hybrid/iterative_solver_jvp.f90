module iterative_solver_jvp_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr, c_int
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    real(c_double), parameter :: relaxation = 0.1_c_double
    public :: analytical_jvp, autodiff_jvp, diagnostic_jvp

    interface
        function enzyme_fwddiff(f, a, da, b, db, iterations, component) &
                result(dx) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr, c_int
            type(c_funptr), value :: f
            real(c_double), intent(in) :: a(*), da(*), b(*), db(*)
            integer(c_int), value :: iterations, component
            real(c_double) :: dx
        end function enzyme_fwddiff
    end interface

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
            dx(component) = enzyme_fwddiff(c_funloc(iteration_component), &
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
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    use iterative_solver_jvp_kernel, only: analytical_jvp, autodiff_jvp, &
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
    integer :: direction, direction_count, trace_iterations
    integer(int64) :: workloads
    character(32) :: action, candidate

    call initialize_inputs()
    call get_environment_variable("FORTNUM_ITERATIVE_ACTION", action)
    call get_environment_variable("FORTNUM_ITERATIVE_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_ITERATIVE_DIRECTIONS", direction_count, 16)
    call read_integer_env("FORTNUM_ITERATIVE_STEPS", trace_iterations, 32)
    call read_int64_env("FORTNUM_ITERATIVE_WORKLOADS", workloads, 5000_int64)
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), direction_count, &
            trace_iterations, workloads)
    else if (trim(action) == "--peak-rss") then
        call benchmark_candidate(trim(candidate), direction_count, &
            trace_iterations, workloads)
        write (*, "(a,i0)") "peak_rss_bytes=", peak_rss_bytes()
    else
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
    end if

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

    subroutine benchmark_candidate(name, directions, steps, count)
        character(*), intent(in) :: name
        integer, intent(in) :: directions, steps
        integer(int64), intent(in) :: count
        integer :: active_direction
        integer(int64) :: workload, start, finish, rate
        real(c_double) :: elapsed_ns, sink

        if (directions < 1 .or. directions > 16) error stop "directions must be 1..16"
        if (steps < 1) error stop "steps must be positive"
        if (name /= "analytical" .and. name /= "autodiff" .and. &
            name /= "diagnostic") then
            error stop "candidate must be analytical, autodiff, or diagnostic"
        end if
        sink = 0.0_c_double
        call system_clock(start, rate)
        do workload = 1, count
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(workload, 1024_int64), c_double)
            do active_direction = 1, directions
                select case (name)
                case ("analytical")
                    call analytical_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(steps, c_int), x, dx)
                case ("autodiff")
                    call autodiff_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(steps, c_int), x, dx)
                case default
                    call diagnostic_jvp(a, da(:, :, active_direction), b, &
                        db(:, active_direction), int(steps, c_int), x, dx)
                end select
                sink = sink + x(1) + dx(1)
            end do
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = real(finish - start, c_double)*1.0e9_c_double / &
            (real(rate, c_double)*real(count, c_double))
        write (*, "(a,a,a,i0,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " directions=", directions, " steps=", steps, " workloads=", count, &
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

end program enzyme_iterative_solver_jvp
