program test_lu_factorization_object
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: lu_factorization_t, lu_factor, lu_solve_factored, &
        LINALG_MAX_N, LINALG_OK
    implicit none

    integer, parameter :: n = LINALG_MAX_N
    character(32) :: action, candidate
    integer(int64) :: iteration_count
    real(dp) :: a(n, n), factors(n, n), rhs(n), solution(n)
    integer :: pivots(n), info
    type(lu_factorization_t) :: factorization

    call initialize_matrix()
    factors = a
    call lu_factor(n, factors, pivots, info)
    if (info /= LINALG_OK) error stop "raw factorization failed"
    call factorization%factor(a, info)
    if (info /= LINALG_OK) error stop "object factorization failed"

    call get_environment_variable("FORTNUM_LU_OBJECT_ACTION", action)
    call get_environment_variable("FORTNUM_LU_OBJECT_CANDIDATE", candidate)
    call read_int64_env("FORTNUM_LU_OBJECT_ITERATIONS", iteration_count, 200000_int64)
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), iteration_count)
    else
        call validate_object()
    end if

contains

    subroutine initialize_matrix()
        integer :: i, j

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    a(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    a(i, j) = 0.2_dp/real(i + j, dp)
                end if
            end do
        end do
    end subroutine initialize_matrix

    subroutine validate_object()
        real(dp) :: expected(n), residual(n)
        integer :: i

        do i = 1, n
            expected(i) = real(i, dp)/real(n, dp)
        end do
        rhs = matmul(a, expected)
        call factorization%solve(rhs, info)
        if (info /= LINALG_OK) error stop "object solve failed"
        call matrix_vector(a, rhs, residual)
        call matrix_vector(a, expected, solution)
        residual = residual - solution
        if (maxval(abs(residual)) > 2.0e-14_dp) then
            error stop "object solve residual mismatch"
        end if
    end subroutine validate_object

    subroutine benchmark_candidate(name, iterations)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: iterations
        integer(int64) :: start, finish, rate, iteration
        real(dp) :: elapsed_ns, sink

        if (name /= "object" .and. name /= "raw") then
            error stop "candidate must be object or raw"
        end if
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, iterations
            solution = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            call matrix_vector(a, solution, rhs)
            if (name == "object") then
                call factorization%solve(rhs, info)
            else
                call lu_solve_factored(n, factors, pivots, rhs, info)
            end if
            sink = sink + rhs(1)
        end do
        call system_clock(finish)
        if (info /= LINALG_OK .or. sink /= sink) error stop "benchmark failed"
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(rate, dp)*real(iterations, dp))
        write (*, "(a,a,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " iterations=", iterations, " ns_per_workload=", elapsed_ns, &
            " sink=", sink
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
        if (ios /= 0) error stop "invalid iteration count"
    end subroutine read_int64_env

end program test_lu_factorization_object
