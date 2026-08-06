program test_fortnum_tensor_product_device
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_tensor_product, only: tensor_factor_t, &
        tensor_product_operator_t
    implicit none

    type(tensor_factor_t) :: factors(2)
    type(tensor_product_operator_t) :: tensor_operator
    type(fortnum_status_t) :: status
    real(dp) :: dense(6, 6), input(6), output(6), expected(6)
    real(dp) :: inputs(6, 2), outputs(6, 2), expected_matrix(6, 2)
    integer :: row, column, rhs, nfail

    nfail = 0
    allocate(factors(1)%values(2, 2), factors(2)%values(3, 3))
    factors(1)%values = reshape([2.0_dp, -0.2_dp, -0.1_dp, 1.7_dp], [2, 2])
    factors(2)%values = reshape([ &
        1.8_dp, -0.3_dp, 0.0_dp, -0.25_dp, 2.1_dp, -0.2_dp, &
        0.0_dp, -0.15_dp, 1.6_dp], [3, 3])
    call tensor_operator%initialize(factors, status)
    call require(status%code == FORTNUM_OK, &
        "device tensor operator initializes", nfail)
    call assemble_dense(factors, dense)
    input = [(0.1_dp + 0.07_dp*real(row, dp), row=1, 6)]
    expected = matmul(dense, input)

    call tensor_operator%enter_data(status, 2)
    call require(status%code == FORTNUM_OK, &
        "tensor factor data enters the accelerator", nfail)
    !$acc data copyin(input) copyout(output)
    call tensor_operator%matvec_device(input, output, status)
    !$acc end data
    call require(status%code == FORTNUM_OK, &
        "accelerator tensor vector product returns success", nfail)
    call require(maxval(abs(output - expected)) < 2.0e-13_dp, &
        "accelerator tensor vector product matches dense oracle", nfail)

    do rhs = 1, 2
        inputs(:, rhs) = input + 0.2_dp*real(rhs, dp)
    end do
    expected_matrix = matmul(dense, inputs)
    !$acc data copyin(inputs) copyout(outputs)
    call tensor_operator%matmat_device(inputs, outputs, status)
    !$acc end data
    call require(status%code == FORTNUM_OK, &
        "accelerator tensor matrix product returns success", nfail)
    call require(maxval(abs(outputs - expected_matrix)) < 2.0e-13_dp, &
        "accelerator tensor matrix product matches dense oracle", nfail)
    call tensor_operator%exit_data(status)
    call require(status%code == FORTNUM_OK, &
        "tensor factor data exits the accelerator", nfail)

    if (nfail > 0) then
        write (error_unit, '(a,i0)') "FAIL: device tensor checks: ", nfail
        error stop 1
    end if
    write (*, '(a)') "PASS: fortnum tensor-product device tests"

contains

    subroutine assemble_dense(input_factors, output_matrix)
        type(tensor_factor_t), intent(in) :: input_factors(:)
        real(dp), intent(out) :: output_matrix(:, :)

        integer :: i1, i2, j1, j2, matrix_row, matrix_column

        do i2 = 1, 3
            do i1 = 1, 2
                matrix_row = i1 + 2*(i2 - 1)
                do j2 = 1, 3
                    do j1 = 1, 2
                        matrix_column = j1 + 2*(j2 - 1)
                        output_matrix(matrix_row, matrix_column) = &
                            input_factors(2)%values(i2, j2)* &
                            input_factors(1)%values(i1, j1)
                    end do
                end do
            end do
        end do
    end subroutine assemble_dense

    subroutine require(condition, description, failures)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        integer, intent(inout) :: failures

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            failures = failures + 1
        end if
    end subroutine require

end program test_fortnum_tensor_product_device
