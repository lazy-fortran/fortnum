program test_linear_solve_jvp_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: lu_factor, lu_solve, linear_solve_jvp_factored, &
        linear_solve_jvp_factored_many, LINALG_MAX_N, LINALG_OK
    implicit none

    integer, parameter :: n = LINALG_MAX_N, max_directions = 16
    real(dp) :: a(n, n), factors(n, n), x(n)
    real(dp) :: da(n, n, max_directions), db(n, max_directions)
    real(dp) :: dx(n, max_directions)
    integer :: pivots(n), info, direction_count
    integer(int64) :: iteration_count
    character(32) :: action, candidate

    call initialize_inputs()
    factors = a
    call lu_factor(n, factors, pivots, info)
    if (info /= LINALG_OK) error stop "factorization failed"
    call get_environment_variable("FORTNUM_JVP_MANY_ACTION", action)
    call get_environment_variable("FORTNUM_JVP_MANY_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_JVP_MANY_DIRECTIONS", direction_count, 16)
    call read_int64_env("FORTNUM_JVP_MANY_ITERATIONS", iteration_count, 20000_int64)
    if (direction_count < 1 .or. direction_count > max_directions) then
        error stop "direction count must be 1..16"
    end if
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), direction_count, iteration_count)
    else
        call validate_products()
    end if

contains

    subroutine initialize_inputs()
        integer :: i, j, direction

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    a(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    a(i, j) = 0.2_dp/real(i + j, dp)
                end if
                do direction = 1, max_directions
                    da(i, j, direction) = 0.001_dp*real( &
                        mod(3*i + 5*j + direction, 11) - 5, dp)
                end do
            end do
            x(j) = real(j, dp)/real(n, dp)
            do direction = 1, max_directions
                db(j, direction) = 0.01_dp*real(mod(7*j + direction, 9) - 4, dp)
            end do
        end do
    end subroutine initialize_inputs

    subroutine validate_products()
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: plus_a(n, n), minus_a(n, n)
        real(dp) :: base_b(n), plus_b(n), minus_b(n), fd(n)
        integer :: direction

        call linear_solve_jvp_factored_many(n, max_directions, factors, pivots, &
            x, da, db, dx, info)
        if (info /= LINALG_OK) error stop "batched JVP failed"
        call matrix_vector(a, x, base_b)
        do direction = 1, max_directions
            plus_a = a + h*da(:, :, direction)
            minus_a = a - h*da(:, :, direction)
            plus_b = base_b + h*db(:, direction)
            minus_b = base_b - h*db(:, direction)
            call lu_solve(n, plus_a, plus_b, info)
            if (info /= LINALG_OK) error stop "plus solve failed"
            call lu_solve(n, minus_a, minus_b, info)
            if (info /= LINALG_OK) error stop "minus solve failed"
            fd = (plus_b - minus_b)/(2.0_dp*h)
            if (maxval(abs(dx(:, direction) - fd)) > 3.0e-9_dp) then
                error stop "batched JVP finite-difference mismatch"
            end if
        end do
    end subroutine validate_products

    subroutine benchmark_candidate(name, directions, iterations)
        character(*), intent(in) :: name
        integer, intent(in) :: directions
        integer(int64), intent(in) :: iterations
        integer :: direction
        integer(int64) :: iteration, start, finish, rate
        real(dp) :: elapsed_ns, sink

        if (name /= "batched" .and. name /= "scalar") then
            error stop "candidate must be batched or scalar"
        end if
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, iterations
            db(1, 1) = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            if (name == "batched") then
                call linear_solve_jvp_factored_many(n, directions, factors, pivots, &
                    x, da(:, :, 1:directions), db(:, 1:directions), &
                    dx(:, 1:directions), info)
            else
                do direction = 1, directions
                    call linear_solve_jvp_factored(n, factors, pivots, x, &
                        da(:, :, direction), db(:, direction), dx(:, direction), info)
                end do
            end if
            sink = sink + dx(1, 1)
        end do
        call system_clock(finish)
        if (info /= LINALG_OK .or. sink /= sink) error stop "benchmark failed"
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(rate, dp)*real(iterations, dp))
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " directions=", directions, " iterations=", iterations, &
            " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

    pure subroutine matrix_vector(matrix, vector, product)
        real(dp), intent(in) :: matrix(n, n), vector(n)
        real(dp), intent(out) :: product(n)
        integer :: i, j

        product = 0.0_dp
        do j = 1, n
            do i = 1, n
                product(i) = product(i) + matrix(i, j)*vector(j)
            end do
        end do
    end subroutine matrix_vector

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

end program test_linear_solve_jvp_many
