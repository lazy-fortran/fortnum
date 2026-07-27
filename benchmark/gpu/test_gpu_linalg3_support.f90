module test_gpu_linalg3_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_linalg3_wrapper, only: &
        det3_jvp_batch, det3_vjp_batch, inv3_jvp_batch, inv3_vjp_batch
    implicit none
    private

    public :: run_linalg3_test

contains

    subroutine run_linalg3_test()
        integer, parameter :: batch_size = 4096
        real(dp), parameter :: h = 1.0e-3_dp
        real(dp), parameter :: tolerance = 2.0e-11_dp
        real(dp) :: matrices(batch_size, 9), directions(batch_size, 9)
        real(dp) :: inverses(batch_size, 9), matrix_cotangents(batch_size, 9)
        real(dp) :: det_cotangents(batch_size), det_products(batch_size)
        real(dp) :: det_adjoints(batch_size, 9)
        real(dp) :: inv_products(batch_size, 9), inv_adjoints(batch_size, 9)
        real(dp) :: plus(9), minus(9), expected(9)
        real(dp) :: lhs, rhs, scale, derivative_h, derivative_half_h
        integer :: batch_index

        call initialize_inputs( &
            matrices, directions, inverses, det_cotangents, matrix_cotangents)

        !$acc data copyin(matrices, directions, inverses) &
        !$acc& copyin(det_cotangents, matrix_cotangents) &
        !$acc& create(det_products, det_adjoints, inv_products, inv_adjoints)
        !$omp target data map(to: matrices, directions, inverses) &
        !$omp& map(to: det_cotangents, matrix_cotangents) &
        !$omp& map(alloc: det_products, det_adjoints, inv_products, inv_adjoints)
        call det3_jvp_batch(batch_size, matrices, directions, det_products)
        call det3_vjp_batch( &
            batch_size, matrices, det_cotangents, det_adjoints)
        call inv3_jvp_batch(batch_size, inverses, directions, inv_products)
        call inv3_vjp_batch( &
            batch_size, inverses, matrix_cotangents, inv_adjoints)
        !$omp target update from(det_products, det_adjoints)
        !$omp target update from(inv_products, inv_adjoints)
        !$omp end target data
        !$acc update self(det_products, det_adjoints, inv_products, inv_adjoints)
        !$acc end data

        do batch_index = 1, batch_size
            plus = matrices(batch_index, :) + h*directions(batch_index, :)
            minus = matrices(batch_index, :) - h*directions(batch_index, :)
            derivative_h = &
                (determinant3(plus) - determinant3(minus))/(2.0_dp*h)
            plus = matrices(batch_index, :) + 0.5_dp*h*directions(batch_index, :)
            minus = matrices(batch_index, :) - 0.5_dp*h*directions(batch_index, :)
            derivative_half_h = &
                (determinant3(plus) - determinant3(minus))/h
            lhs = (4.0_dp*derivative_half_h - derivative_h)/3.0_dp
            scale = max(1.0_dp, abs(lhs), abs(det_products(batch_index)))
            if (abs(lhs - det_products(batch_index)) > tolerance*scale) then
                error stop "GPU determinant JVP disagrees with finite difference"
            end if

            call inverse_jvp_oracle( &
                inverses(batch_index, :), directions(batch_index, :), expected)
            scale = max(1.0_dp, maxval(abs(expected)))
            if (maxval(abs(expected - inv_products(batch_index, :))) > &
                tolerance*scale) then
                error stop "GPU inverse JVP disagrees with matrix-product oracle"
            end if
        end do

        lhs = sum(det_cotangents*det_products)
        rhs = sum(directions*det_adjoints)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "GPU determinant products violate adjoint identity"
        end if
        lhs = sum(matrix_cotangents*inv_products)
        rhs = sum(directions*inv_adjoints)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "GPU inverse products violate adjoint identity"
        end if
    end subroutine run_linalg3_test

    subroutine initialize_inputs( &
            matrices, directions, inverses, det_cotangents, matrix_cotangents)
        real(dp), intent(out) :: matrices(:, :), directions(:, :)
        real(dp), intent(out) :: inverses(:, :), det_cotangents(:)
        real(dp), intent(out) :: matrix_cotangents(:, :)
        real(dp) :: shift
        integer :: batch_index, entry

        do batch_index = 1, size(matrices, 1)
            shift = real(mod(batch_index, 29), dp)/113.0_dp
            matrices(batch_index, :) = [ &
                1.8_dp + shift, 0.1_dp, -0.2_dp, &
                0.2_dp, 1.6_dp + shift, 0.15_dp, &
                -0.1_dp, 0.25_dp, 1.9_dp + shift]
            inverses(batch_index, :) = [ &
                0.7_dp + 0.1_dp*shift, -0.04_dp, 0.08_dp, &
                -0.06_dp, 0.8_dp + 0.1_dp*shift, -0.05_dp, &
                0.03_dp, -0.09_dp, 0.65_dp + 0.1_dp*shift]
            det_cotangents(batch_index) = &
                sin(real(batch_index, dp)/17.0_dp)
            do entry = 1, 9
                directions(batch_index, entry) = &
                    cos(real(7*batch_index + 3*entry, dp)/31.0_dp)
                matrix_cotangents(batch_index, entry) = &
                    sin(real(5*batch_index + 2*entry, dp)/23.0_dp)
            end do
        end do
    end subroutine initialize_inputs

    pure function determinant3(matrix) result(value)
        real(dp), intent(in) :: matrix(9)
        real(dp) :: value

        value = matrix(1)*(matrix(5)*matrix(9) - matrix(6)*matrix(8)) - &
            matrix(4)*(matrix(2)*matrix(9) - matrix(3)*matrix(8)) + &
            matrix(7)*(matrix(2)*matrix(6) - matrix(3)*matrix(5))
    end function determinant3

    pure subroutine inverse_jvp_oracle(inverse, direction, product)
        real(dp), intent(in) :: inverse(9), direction(9)
        real(dp), intent(out) :: product(9)
        real(dp) :: first(3, 3), second(3, 3), third(3, 3)
        integer :: row, column

        first = reshape(inverse, [3, 3])
        second = reshape(direction, [3, 3])
        third = -matmul(matmul(first, second), first)
        do column = 1, 3
            do row = 1, 3
                product(row + 3*(column - 1)) = third(row, column)
            end do
        end do
    end subroutine inverse_jvp_oracle

end module test_gpu_linalg3_support
