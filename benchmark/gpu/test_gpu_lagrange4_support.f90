module test_gpu_lagrange4_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_lagrange4_wrapper, only: &
        lagrange4_jvp_batch, lagrange4_vjp_batch
    implicit none
    private

    public :: run_lagrange4_test

contains

    subroutine run_lagrange4_test()
        integer, parameter :: n = 4096
        real(dp), parameter :: tolerance = 8.0e-13_dp
        real(dp), parameter :: nodes(4) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
        real(dp) :: x(n), samples(n, 4), tx(n), sample_tangents(n, 4)
        real(dp) :: cotangents(n), values(n), products(n), adjoint_x(n)
        real(dp) :: adjoint_samples(n, 4)
        real(dp) :: expected_value, expected_product, lhs, rhs, scale
        integer :: i, j

        do i = 1, n
            x(i) = -0.8_dp + 2.6_dp*real(mod(17*i, 4093), dp)/4092.0_dp
            tx(i) = -0.7_dp + 1.4_dp*real(mod(19*i, 4091), dp)/4090.0_dp
            cotangents(i) = -0.9_dp + &
                1.8_dp*real(mod(23*i, 4079), dp)/4078.0_dp
            do j = 1, 4
                samples(i, j) = primal_polynomial(nodes(j), i)
                sample_tangents(i, j) = tangent_polynomial(nodes(j), i)
            end do
        end do

        !$acc data copyin(x, samples, tx, sample_tangents, cotangents) &
        !$acc& create(values, products, adjoint_x, adjoint_samples)
        !$omp target data map(to: x, samples, tx, sample_tangents, cotangents) &
        !$omp& map(alloc: values, products, adjoint_x, adjoint_samples)
        call lagrange4_jvp_batch( &
            n, x, samples, tx, sample_tangents, values, products)
        call lagrange4_vjp_batch( &
            n, x, samples, cotangents, values, adjoint_x, adjoint_samples)
        !$omp target update from(values, products, adjoint_x, adjoint_samples)
        !$omp end target data
        !$acc update self(values, products, adjoint_x, adjoint_samples)
        !$acc end data

        do i = 1, n
            expected_value = primal_polynomial(x(i), i)
            expected_product = primal_derivative(x(i))*tx(i) + &
                tangent_polynomial(x(i), i)
            scale = max(1.0_dp, abs(expected_value), abs(expected_product))
            if (abs(values(i) - expected_value) > tolerance*scale) then
                error stop "GPU Lagrange value disagrees with cubic oracle"
            end if
            if (abs(products(i) - expected_product) > tolerance*scale) then
                error stop "GPU Lagrange JVP disagrees with cubic oracle"
            end if
        end do
        lhs = sum(cotangents*products)
        rhs = sum(tx*adjoint_x) + sum(sample_tangents*adjoint_samples)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "GPU Lagrange products violate adjoint identity"
        end if
    end subroutine run_lagrange4_test

    pure function primal_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 17), dp)/31.0_dp
        value = (0.3_dp + shift) - 0.4_dp*z + 0.2_dp*z*z - 0.05_dp*z*z*z
    end function primal_polynomial

    pure function primal_derivative(z) result(value)
        real(dp), intent(in) :: z
        real(dp) :: value

        value = -0.4_dp + 0.4_dp*z - 0.15_dp*z*z
    end function primal_derivative

    pure function tangent_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 13), dp)/29.0_dp
        value = (-0.2_dp + shift) + 0.3_dp*z - 0.1_dp*z*z + 0.04_dp*z*z*z
    end function tangent_polynomial

end module test_gpu_lagrange4_support
