module direct_solver_vjp_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_direct_solver_objective, only: &
        fortnum_enzyme_direct_solver_objective_vjp
    implicit none
    private

    integer, parameter, public :: solver_size = 4
    public :: analytical_vjp, autodiff_vjp, diagnostic_vjp

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

        value = fortnum_enzyme_direct_solver_objective_vjp( &
            a, b, u, 1.0_c_double, abar, bbar, ubar)
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
    use, intrinsic :: iso_c_binding, only: c_double
    use direct_solver_vjp_kernel, only: analytical_vjp, autodiff_vjp, &
        diagnostic_vjp, solver_size
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none

    real(c_double) :: a(solver_size, solver_size), b(solver_size)
    real(c_double) :: u(solver_size, 16), abar(solver_size, solver_size)
    real(c_double) :: bbar(solver_size), reference_a(solver_size, solver_size)
    real(c_double) :: reference_b(solver_size), value, reference_value
    real(c_double) :: samples(fixture_sample_count), sink
    integer :: candidate_kind, cotangent, cotangent_count, iterations
    character(32) :: action, candidate
    logical :: valid

    call initialize_inputs()
    call read_fixture_environment("FORTNUM_DIRECT_SOLVER_VJP_ACTION", &
        "validate", action)
    if (trim(action) == "validate") then
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
        stop
    end if

    call read_fixture_environment("FORTNUM_DIRECT_SOLVER_VJP_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_integer("FORTNUM_DIRECT_SOLVER_VJP_COTANGENTS", 16, &
        cotangent_count, valid)
    if (.not. valid) error stop "invalid direct-solver cotangent count"
    call read_fixture_integer("FORTNUM_DIRECT_SOLVER_VJP_ITERATIONS", 2000, &
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
    if (sink /= sink) error stop "direct-solver VJP benchmark produced NaN"

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

    subroutine parse_configuration()
        if (cotangent_count < 1 .or. cotangent_count > 16) then
            error stop "cotangents must be 1..16"
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
        write (name, "(a,'_vjp_d',i0)") trim(candidate), cotangent_count
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine benchmark_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        integer :: active_cotangent
        integer :: iteration
        real(c_double) :: elapsed_ns, local_sink

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            b(1) = 0.5_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024), c_double)
            do active_cotangent = 1, cotangent_count
                select case (candidate_kind)
                case (1)
                    call analytical_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                case (2)
                    call autodiff_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                case default
                    call diagnostic_vjp(a, b, u(:, active_cotangent), value, &
                        abar, bbar)
                end select
                local_sink = local_sink + value + abar(1, 1) + bbar(1)
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

end program enzyme_direct_solver_vjp
