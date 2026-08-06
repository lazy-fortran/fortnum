program test_fortnum_toeplitz
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_OK
    use fortnum_toeplitz, only: toeplitz_operator_t
    implicit none

    type(toeplitz_operator_t) :: toeplitz
    type(fortnum_status_t) :: status
    real(dp) :: column(5), row_values(5), dense(5, 5)
    real(dp) :: input(5), output(5), expected(5), diagonal(5)
    real(dp) :: inputs(5, 2), outputs(5, 2), expected_matrix(5, 2)
    integer :: i, j, column_index, nfail

    nfail = 0
    column = [2.2_dp, -0.4_dp, 0.15_dp, -0.08_dp, 0.03_dp]
    row_values = [2.2_dp, 0.25_dp, -0.12_dp, 0.06_dp, -0.02_dp]
    call toeplitz%initialize(column, status, row_values)
    call require(status%code == FORTNUM_OK, &
        "nonsymmetric Toeplitz operator initializes", nfail)
    call require(toeplitz%element_count() == 5, &
        "Toeplitz operator reports its size", nfail)

    do i = 1, 5
        do j = 1, 5
            if (i >= j) then
                dense(i, j) = column(i - j + 1)
            else
                dense(i, j) = row_values(j - i + 1)
            end if
        end do
    end do

    input = [0.3_dp, -0.7_dp, 1.1_dp, 0.2_dp, -0.4_dp]
    expected = matmul(dense, input)
    call toeplitz%matvec(input, output)
    call require(maxval(abs(output - expected)) < 2.0e-13_dp, &
        "FFT Toeplitz matrix-vector product matches dense oracle", nfail)

    do column_index = 1, 2
        inputs(:, column_index) = input + &
            0.17_dp*real(column_index, dp)
    end do
    expected_matrix = matmul(dense, inputs)
    call toeplitz%matmat(inputs, outputs)
    call require(maxval(abs(outputs - expected_matrix)) < 2.0e-13_dp, &
        "FFT Toeplitz multi-RHS product matches dense oracle", nfail)

    diagonal = toeplitz%diagonal()
    call require(maxval(abs(diagonal - column(1))) < 2.0e-14_dp, &
        "Toeplitz diagonal uses the first column", nfail)

    call toeplitz%initialize(column, status)
    call require(status%code == FORTNUM_OK, &
        "symmetric Toeplitz default row initializes", nfail)
    do i = 1, 5
        do j = 1, 5
            dense(i, j) = column(abs(i - j) + 1)
        end do
    end do
    expected = matmul(dense, input)
    call toeplitz%matvec(input, output)
    call require(maxval(abs(output - expected)) < 2.0e-13_dp, &
        "symmetric Toeplitz default row matches dense oracle", nfail)

    call toeplitz%initialize(column(:4), status, row_values)
    call require(status%code == FORTNUM_DOMAIN_ERROR, &
        "Toeplitz rejects row and column size mismatch", nfail)

    if (nfail > 0) then
        write (error_unit, '(a,i0)') "FAIL: Toeplitz checks: ", nfail
        error stop 1
    end if
    write (*, '(a)') "PASS: fortnum Toeplitz behavioral tests"

contains

    subroutine require(condition, description, failures)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        integer, intent(inout) :: failures

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            failures = failures + 1
        end if
    end subroutine require

end program test_fortnum_toeplitz
