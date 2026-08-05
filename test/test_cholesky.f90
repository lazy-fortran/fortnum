program test_cholesky
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit, int64
    use fortnum_cholesky, only: cholesky_factorization_t
    use fortnum_status, only: fortnum_status_t, status_ok, &
        FORTNUM_CONVERGENCE_ERROR
    implicit none

    character(32) :: action
    integer :: failures

    action = ""
    call get_environment_variable("FORTNUM_CHOLESKY_ACTION", action)
    if (trim(action) == "--benchmark") then
        call benchmark()
        stop
    end if

    failures = 0
    call test_factor_and_solve(failures)
    call test_multiple_rhs(failures)
    call test_non_positive_definite(failures)
    if (failures /= 0) then
        write (error_unit, '(i0,a)') failures, " Cholesky test(s) failed"
        error stop 1
    end if
    write (*, '(a)') "PASS"

contains

    subroutine benchmark()
        integer, parameter :: n = 64, repetitions = 100
        type(cholesky_factorization_t) :: factorization
        type(fortnum_status_t) :: status
        real(dp) :: lower(n, n), matrix(n, n), rhs(n), rhs0(n), exact(n)
        real(dp) :: elapsed_ns, sink
        integer(int64) :: start, finish, rate
        integer :: i, j, iteration

        lower = 0.0_dp
        do j = 1, n
            do i = j, n
                if (i == j) then
                    lower(i, j) = 2.0_dp + real(i, dp)/real(n, dp)
                else
                    lower(i, j) = 0.03_dp/real(i + j, dp)
                end if
            end do
        end do
        matrix = matmul(lower, transpose(lower))
        do i = 1, n
            exact(i) = sin(0.13_dp*real(i, dp))
        end do
        rhs0 = matmul(matrix, exact)
        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, repetitions
            rhs = rhs0
            call factorization%factorize(matrix, status)
            call factorization%solve(rhs, status)
            sink = sink + rhs(1)
        end do
        call system_clock(finish)
        if (.not. status_ok(status) .or. abs(rhs(1) - exact(1)) > 2.0e-12_dp .or. &
            sink /= sink) then
            error stop "Cholesky benchmark failed"
        end if
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(rate, dp)*real(repetitions, dp))
        write (*, '(a,i0,a,i0,a,f0.6)') &
            "matrix_size=", n, " repetitions=", repetitions, &
            " ns_per_factor_and_solve=", elapsed_ns
    end subroutine benchmark

    subroutine test_factor_and_solve(failures)
        integer, intent(inout) :: failures
        type(cholesky_factorization_t) :: factorization
        type(fortnum_status_t) :: status
        real(dp) :: lower(3, 3), matrix(3, 3), rhs(3), exact(3)
        real(dp) :: reconstructed(3, 3), logdet, expected_logdet

        lower = 0.0_dp
        lower(1, 1) = 2.0_dp
        lower(2, 1) = 0.5_dp
        lower(2, 2) = 1.5_dp
        lower(3, 1) = -0.2_dp
        lower(3, 2) = 0.3_dp
        lower(3, 3) = 1.2_dp
        matrix = matmul(lower, transpose(lower))
        exact = [1.2_dp, -0.7_dp, 2.0_dp]
        rhs = matmul(matrix, exact)

        call factorization%factorize(matrix, status)
        call factorization%log_determinant(logdet, status)
        reconstructed = matmul(factorization%lower, transpose(factorization%lower))
        call factorization%solve(rhs, status)
        expected_logdet = 2.0_dp*sum(log([2.0_dp, 1.5_dp, 1.2_dp]))
        if (.not. status_ok(status) .or. maxval(abs(rhs - exact)) > 2.0e-14_dp .or. &
            maxval(abs(reconstructed - matrix)) > 2.0e-14_dp .or. &
            abs(logdet - expected_logdet) > 2.0e-14_dp) then
            write (error_unit, '(a)') "FAIL [factor] reconstruction or solve"
            failures = failures + 1
        end if
    end subroutine test_factor_and_solve

    subroutine test_multiple_rhs(failures)
        integer, intent(inout) :: failures
        type(cholesky_factorization_t) :: factorization
        type(fortnum_status_t) :: status
        real(dp) :: matrix(3, 3), rhs(3, 2), exact(3, 2)

        matrix = reshape([ &
            5.0_dp, 1.0_dp, 0.2_dp, &
            1.0_dp, 4.0_dp, -0.4_dp, &
            0.2_dp, -0.4_dp, 3.0_dp], shape(matrix))
        exact(:, 1) = [1.0_dp, -2.0_dp, 0.5_dp]
        exact(:, 2) = [-0.25_dp, 0.75_dp, 2.0_dp]
        rhs = matmul(matrix, exact)
        call factorization%factorize(matrix, status)
        call factorization%solve(rhs, status)
        if (.not. status_ok(status) .or. maxval(abs(rhs - exact)) > 2.0e-14_dp) then
            write (error_unit, '(a)') "FAIL [multiple_rhs] solve"
            failures = failures + 1
        end if
    end subroutine test_multiple_rhs

    subroutine test_non_positive_definite(failures)
        integer, intent(inout) :: failures
        type(cholesky_factorization_t) :: factorization
        type(fortnum_status_t) :: status
        real(dp) :: matrix(2, 2)

        matrix = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], shape(matrix))
        call factorization%factorize(matrix, status)
        if (status%code /= FORTNUM_CONVERGENCE_ERROR) then
            write (error_unit, '(a)') "FAIL [status] non-positive-definite matrix"
            failures = failures + 1
        end if
    end subroutine test_non_positive_definite

end program test_cholesky
