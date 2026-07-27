program test_linear_solve_vjp_many
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: lu_factor, lu_solve, linear_solve_vjp_factored, &
        linear_solve_vjp_factored_many, LINALG_MAX_N, LINALG_OK
    implicit none

    integer, parameter :: n = LINALG_MAX_N, max_cotangents = 16
    real(dp) :: a(n, n), transpose_factors(n, n), b(n), x(n)
    real(dp) :: u(n, max_cotangents)
    real(dp) :: abar(n, n, max_cotangents), bbar(n, max_cotangents)
    integer :: pivots(n), info, cotangent_count
    integer(int64) :: iteration_count
    character(32) :: action, candidate

    call initialize_inputs()
    transpose_factors = transpose(a)
    call lu_factor(n, transpose_factors, pivots, info)
    if (info /= LINALG_OK) error stop "factorization failed"
    call get_environment_variable("FORTNUM_VJP_MANY_ACTION", action)
    call get_environment_variable("FORTNUM_VJP_MANY_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_VJP_MANY_COTANGENTS", cotangent_count, 16)
    call read_int64_env("FORTNUM_VJP_MANY_ITERATIONS", iteration_count, 20000_int64)
    if (cotangent_count < 1 .or. cotangent_count > max_cotangents) then
        error stop "cotangent count must be 1..16"
    end if
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), cotangent_count, iteration_count)
    else
        call validate_products()
    end if

contains

    subroutine initialize_inputs()
        integer :: i, j, cotangent

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    a(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    a(i, j) = 0.2_dp/real(i + j, dp)
                end if
            end do
            x(j) = real(j, dp)/real(n, dp)
            do cotangent = 1, max_cotangents
                u(j, cotangent) = 0.01_dp*real( &
                    mod(7*j + cotangent, 9) - 4, dp)
            end do
        end do
        call matrix_vector(a, x, b)
    end subroutine initialize_inputs

    subroutine validate_products()
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: da(n, n), db(n), plus_a(n, n), minus_a(n, n)
        real(dp) :: plus_b(n), minus_b(n), objective_fd, contraction
        integer :: cotangent, i, j

        call linear_solve_vjp_factored_many(n, max_cotangents, transpose_factors, &
            pivots, x, u, abar, bbar, info)
        if (info /= LINALG_OK) error stop "batched VJP failed"
        do cotangent = 1, max_cotangents
            do j = 1, n
                db(j) = 0.01_dp*real(mod(5*j + cotangent, 11) - 5, dp)
                do i = 1, n
                    da(i, j) = 0.001_dp*real( &
                        mod(3*i + 5*j + cotangent, 13) - 6, dp)
                end do
            end do
            plus_a = a + h*da
            minus_a = a - h*da
            plus_b = b + h*db
            minus_b = b - h*db
            call lu_solve(n, plus_a, plus_b, info)
            if (info /= LINALG_OK) error stop "plus solve failed"
            call lu_solve(n, minus_a, minus_b, info)
            if (info /= LINALG_OK) error stop "minus solve failed"
            objective_fd = dot_product(u(:, cotangent), plus_b - minus_b)/(2.0_dp*h)
            contraction = sum(abar(:, :, cotangent)*da) &
                + dot_product(bbar(:, cotangent), db)
            if (abs(contraction - objective_fd) > 3.0e-9_dp) then
                error stop "batched VJP finite-difference mismatch"
            end if
        end do
    end subroutine validate_products

    subroutine benchmark_candidate(name, cotangents, iterations)
        character(*), intent(in) :: name
        integer, intent(in) :: cotangents
        integer(int64), intent(in) :: iterations
        integer :: cotangent
        integer(int64) :: iteration, start, finish, rate
        real(dp) :: elapsed_ns, sink

        if (name /= "batched" .and. name /= "scalar") then
            error stop "candidate must be batched or scalar"
        end if
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, iterations
            u(1, 1) = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            if (name == "batched") then
                call linear_solve_vjp_factored_many(n, cotangents, &
                    transpose_factors, pivots, x, u(:, 1:cotangents), &
                    abar(:, :, 1:cotangents), bbar(:, 1:cotangents), info)
            else
                do cotangent = 1, cotangents
                    call linear_solve_vjp_factored(n, transpose_factors, pivots, x, &
                        u(:, cotangent), abar(:, :, cotangent), &
                        bbar(:, cotangent), info)
                end do
            end if
            sink = sink + abar(1, 1, 1) + bbar(1, 1)
        end do
        call system_clock(finish)
        if (info /= LINALG_OK .or. sink /= sink) error stop "benchmark failed"
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(rate, dp)*real(iterations, dp))
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " cotangents=", cotangents, " iterations=", iterations, &
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

end program test_linear_solve_vjp_many
