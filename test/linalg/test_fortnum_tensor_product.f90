program test_fortnum_tensor_product
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_OK
    use fortnum_tensor_product, only: tensor_factor_t, &
        tensor_product_operator_t
    implicit none

    type(tensor_factor_t) :: factors(3), bad_factor(1)
    type(tensor_product_operator_t) :: tensor_operator
    type(fortnum_status_t) :: status
    real(dp) :: input(12), output(12), expected(12), diagonal(12)
    real(dp) :: inputs(12, 3), outputs(12, 3), expected_matrix(12, 3)
    real(dp) :: dense(12, 12)
    integer :: column, i1, i2, i3, j1, j2, j3, row, row_index
    integer :: column_index, nfail

    nfail = 0
    allocate(factors(1)%values(2, 2), factors(2)%values(3, 3))
    allocate(factors(3)%values(2, 2))
    factors(1)%values = reshape([ &
        2.0_dp, -0.5_dp, 0.25_dp, 1.5_dp], [2, 2])
    factors(2)%values = reshape([ &
        1.0_dp, 0.2_dp, -0.1_dp, &
        0.3_dp, 1.7_dp, 0.4_dp, &
        -0.2_dp, 0.5_dp, 1.2_dp], [3, 3])
    factors(3)%values = reshape([ &
        0.8_dp, 0.6_dp, -0.3_dp, 1.1_dp], [2, 2])

    call tensor_operator%initialize(factors, status)
    call require(status%code == FORTNUM_OK, &
        "three-dimensional tensor product initializes", nfail)
    call require(tensor_operator%element_count() == 12, &
        "tensor product reports the flattened grid size", nfail)

    dense = 0.0_dp
    do i3 = 1, 2
        do i2 = 1, 3
            do i1 = 1, 2
                row_index = flatten(i1, i2, i3)
                do j3 = 1, 2
                    do j2 = 1, 3
                        do j1 = 1, 2
                            column_index = flatten(j1, j2, j3)
                            dense(row_index, column_index) = &
                                factors(3)%values(i3, j3)* &
                                factors(2)%values(i2, j2)* &
                                factors(1)%values(i1, j1)
                        end do
                    end do
                end do
            end do
        end do
    end do

    input = [(0.15_dp + 0.07_dp*real(row, dp), row=1, 12)]
    expected = matmul(dense, input)
    call tensor_operator%matvec(input, output, status)
    call require(status%code == FORTNUM_OK, &
        "tensor-product matrix-vector status is successful", nfail)
    call require(maxval(abs(output - expected)) < 2.0e-14_dp, &
        "matrix-free tensor product matches explicit Kronecker oracle", nfail)

    do column = 1, 3
        inputs(:, column) = input + 0.11_dp*real(column, dp)
    end do
    expected_matrix = matmul(dense, inputs)
    call tensor_operator%matmat(inputs, outputs, status)
    call require(status%code == FORTNUM_OK, &
        "tensor-product matrix-matrix status is successful", nfail)
    call require(maxval(abs(outputs - expected_matrix)) < 2.0e-14_dp, &
        "batched tensor product matches explicit Kronecker oracle", nfail)

    diagonal = [(dense(row, row), row=1, 12)]
    call tensor_operator%diagonal(output, status)
    call require(status%code == FORTNUM_OK, &
        "tensor-product diagonal status is successful", nfail)
    call require(maxval(abs(output - diagonal)) < 2.0e-14_dp, &
        "tensor-product diagonal matches explicit matrix diagonal", nfail)

    allocate(bad_factor(1)%values(2, 3))
    call tensor_operator%initialize(bad_factor, status)
    call require(status%code == FORTNUM_DOMAIN_ERROR, &
        "nonsquare tensor factor is rejected", nfail)

    call tensor_operator%matvec(input, output, status)
    call require(status%code == FORTNUM_DOMAIN_ERROR, &
        "uninitialized tensor product rejects a matvec", nfail)

    if (nfail > 0) then
        write (error_unit, '(a,i0)') "FAIL: tensor-product checks: ", nfail
        error stop 1
    end if
    write (*, '(a)') "PASS: fortnum tensor-product behavioral tests"

contains

    integer function flatten(i1, i2, i3) result(index)
        integer, intent(in) :: i1, i2, i3

        index = i1 + 2*(i2 - 1) + 6*(i3 - 1)
    end function flatten

    subroutine require(condition, description, failures)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        integer, intent(inout) :: failures

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            failures = failures + 1
        end if
    end subroutine require

end program test_fortnum_tensor_product
