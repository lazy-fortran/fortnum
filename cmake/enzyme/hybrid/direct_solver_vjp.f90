module direct_solver_vjp_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    public :: analytical_vjp, autodiff_vjp, diagnostic_vjp

    interface
        function enzyme_autodiff(f, a, abar, b, bbar, u, ubar) result(value) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: a(*), b(*), u(*)
            real(c_double), intent(inout) :: abar(*), bbar(*), ubar(*)
            real(c_double) :: value
        end function enzyme_autodiff
    end interface

contains

    pure function solve_objective(a, b, u) result(value) &
            bind(c, name="fortnum_direct_solve_objective")
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), u(solver_size)
        real(c_double) :: value, x(solver_size)

        call solve_direct(a, b, x)
        value = dot_product(u, x)
    end function solve_objective

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

    pure subroutine analytical_vjp(a, b, u, value, abar, bbar)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), u(solver_size)
        real(c_double), intent(out) :: value
        real(c_double), intent(out) :: abar(solver_size, solver_size)
        real(c_double), intent(out) :: bbar(solver_size)
        real(c_double) :: at(solver_size, solver_size), x(solver_size)
        integer :: i, j

        call solve_direct(a, b, x)
        value = dot_product(u, x)
        at = transpose(a)
        call solve_direct(at, u, bbar)
        do j = 1, solver_size
            do i = 1, solver_size
                abar(i, j) = -bbar(i)*x(j)
            end do
        end do
    end subroutine analytical_vjp

    subroutine autodiff_vjp(a, b, u, value, abar, bbar)
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), u(solver_size)
        real(c_double), intent(out) :: value
        real(c_double), intent(out) :: abar(solver_size, solver_size)
        real(c_double), intent(out) :: bbar(solver_size)
        real(c_double) :: ubar(solver_size)

        abar = 0.0_c_double
        bbar = 0.0_c_double
        ubar = 0.0_c_double
        value = enzyme_autodiff(c_funloc(solve_objective), &
            a, abar, b, bbar, u, ubar)
    end subroutine autodiff_vjp

    pure subroutine diagnostic_vjp(a, b, u, value, abar, bbar)
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double), intent(in) :: a(solver_size, solver_size)
        real(c_double), intent(in) :: b(solver_size), u(solver_size)
        real(c_double), intent(out) :: value
        real(c_double), intent(out) :: abar(solver_size, solver_size)
        real(c_double), intent(out) :: bbar(solver_size)
        real(c_double) :: plus_a(solver_size, solver_size)
        real(c_double) :: minus_a(solver_size, solver_size)
        real(c_double) :: plus_b(solver_size), minus_b(solver_size)
        integer :: i, j

        value = solve_objective(a, b, u)
        do j = 1, solver_size
            do i = 1, solver_size
                plus_a = a
                minus_a = a
                plus_a(i, j) = plus_a(i, j) + h
                minus_a(i, j) = minus_a(i, j) - h
                abar(i, j) = (solve_objective(plus_a, b, u) &
                    - solve_objective(minus_a, b, u))/(2.0_c_double*h)
            end do
        end do
        do i = 1, solver_size
            plus_b = b
            minus_b = b
            plus_b(i) = plus_b(i) + h
            minus_b(i) = minus_b(i) - h
            bbar(i) = (solve_objective(a, plus_b, u) &
                - solve_objective(a, minus_b, u))/(2.0_c_double*h)
        end do
    end subroutine diagnostic_vjp

end module direct_solver_vjp_kernel

program enzyme_direct_solver_vjp
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    use direct_solver_vjp_kernel, only: analytical_vjp, autodiff_vjp, &
        diagnostic_vjp, solver_size
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(c_double) :: a(solver_size, solver_size), b(solver_size)
    real(c_double) :: u(solver_size, 16), abar(solver_size, solver_size)
    real(c_double) :: bbar(solver_size), reference_a(solver_size, solver_size)
    real(c_double) :: reference_b(solver_size), value, reference_value
    integer :: cotangent, cotangent_count
    integer(int64) :: iterations
    character(32) :: action, candidate

    call initialize_inputs()
    call get_environment_variable("FORTNUM_DIRECT_SOLVER_VJP_ACTION", action)
    call get_environment_variable("FORTNUM_DIRECT_SOLVER_VJP_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_DIRECT_SOLVER_VJP_COTANGENTS", &
        cotangent_count, 16)
    call read_int64_env("FORTNUM_DIRECT_SOLVER_VJP_ITERATIONS", iterations, &
        2000_int64)
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), cotangent_count, iterations)
    else if (trim(action) == "--peak-rss") then
        call benchmark_candidate(trim(candidate), cotangent_count, iterations)
        write (*, "(a,i0)") "peak_rss_bytes=", peak_rss_bytes()
    else
        do cotangent = 1, 16
            call analytical_vjp(a, b, u(:, cotangent), reference_value, &
                reference_a, reference_b)
            call autodiff_vjp(a, b, u(:, cotangent), value, abar, bbar)
            if (abs(value - reference_value) > 2.0e-12_c_double) then
                error stop "autodiff direct-solver value mismatch"
            end if
            if (maxval(abs(abar - reference_a)) > 2.0e-12_c_double) then
                error stop "autodiff direct-solver A VJP mismatch"
            end if
            if (maxval(abs(bbar - reference_b)) > 2.0e-12_c_double) then
                error stop "autodiff direct-solver b VJP mismatch"
            end if
            call diagnostic_vjp(a, b, u(:, cotangent), value, abar, bbar)
            if (maxval(abs(abar - reference_a)) > 2.0e-10_c_double) then
                error stop "diagnostic direct-solver A VJP mismatch"
            end if
            if (maxval(abs(bbar - reference_b)) > 2.0e-10_c_double) then
                error stop "diagnostic direct-solver b VJP mismatch"
            end if
        end do
        write (*, "(a)") "PASS"
    end if

contains

    subroutine initialize_inputs()
        integer :: i, j, active_cotangent

        do j = 1, solver_size
            do i = 1, solver_size
                if (i == j) then
                    a(i, j) = 4.0_c_double + real(i, c_double)
                else
                    a(i, j) = 0.2_c_double/real(i + j, c_double)
                end if
            end do
            b(j) = 0.5_c_double*real(j, c_double)
            do active_cotangent = 1, 16
                u(j, active_cotangent) = 0.02_c_double*real( &
                    mod(7*j + active_cotangent, 9) - 4, c_double)
            end do
        end do
    end subroutine initialize_inputs

    subroutine benchmark_candidate(name, cotangents, count)
        character(*), intent(in) :: name
        integer, intent(in) :: cotangents
        integer(int64), intent(in) :: count
        integer :: active_cotangent
        integer(int64) :: iteration, start, finish, rate
        real(c_double) :: elapsed_ns, sink

        if (cotangents < 1 .or. cotangents > 16) then
            error stop "cotangents must be 1..16"
        end if
        if (name /= "analytical" .and. name /= "autodiff" .and. &
            name /= "diagnostic") then
            error stop "candidate must be analytical, autodiff, or diagnostic"
        end if
        sink = 0.0_c_double
        call system_clock(start, rate)
        do iteration = 1, count
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024_int64), c_double)
            do active_cotangent = 1, cotangents
                select case (name)
                case ("analytical")
                    call analytical_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                case ("autodiff")
                    call autodiff_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                case default
                    call diagnostic_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                end select
                sink = sink + value + abar(1, 1) + bbar(1)
            end do
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = real(finish - start, c_double)*1.0e9_c_double / &
            (real(rate, c_double)*real(count, c_double))
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " cotangents=", cotangents, " iterations=", count, &
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

end program enzyme_direct_solver_vjp
