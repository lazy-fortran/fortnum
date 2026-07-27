module test_gpu_fixed_quadrature_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_gpu_fixed_quadrature_wrapper, only: &
        fixed_quadrature_jvp_batch, fixed_quadrature_vjp_batch
    implicit none
    private

    public :: run_fixed_quadrature_test

contains

    subroutine run_fixed_quadrature_test()
        integer, parameter :: batch_size = 4096, rule_order = 16
        real(dp), parameter :: tolerance = 2.0e-13_dp
        real(dp) :: nodes(rule_order), weights(rule_order)
        real(dp) :: tangents(batch_size, rule_order)
        real(dp) :: cotangents(batch_size), products(batch_size)
        real(dp) :: adjoints(batch_size, rule_order)
        real(dp) :: expected, lhs, rhs, scale, shift
        integer :: batch_index, node_index

        call gauss_legendre(rule_order, nodes, weights)
        do batch_index = 1, batch_size
            shift = real(mod(batch_index, 19), dp)/37.0_dp
            cotangents(batch_index) = -0.9_dp + &
                1.8_dp*real(mod(23*batch_index, 4079), dp)/4078.0_dp
            do node_index = 1, rule_order
                tangents(batch_index, node_index) = &
                    nodes(node_index)**4 + shift
            end do
        end do

        !$acc data copyin(weights, tangents, cotangents) &
        !$acc& create(products, adjoints)
        !$omp target data map(to: weights, tangents, cotangents) &
        !$omp& map(alloc: products, adjoints)
        call fixed_quadrature_jvp_batch( &
            batch_size, rule_order, weights, tangents, products)
        call fixed_quadrature_vjp_batch( &
            batch_size, rule_order, weights, cotangents, adjoints)
        !$omp target update from(products, adjoints)
        !$omp end target data
        !$acc update self(products, adjoints)
        !$acc end data

        do batch_index = 1, batch_size
            shift = real(mod(batch_index, 19), dp)/37.0_dp
            expected = 2.0_dp/5.0_dp + 2.0_dp*shift
            scale = max(1.0_dp, abs(expected))
            if (abs(products(batch_index) - expected) > tolerance*scale) then
                error stop "GPU fixed-quadrature JVP disagrees with exact integral"
            end if
        end do
        lhs = sum(cotangents*products)
        rhs = sum(tangents*adjoints)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "GPU fixed-quadrature products violate adjoint identity"
        end if
    end subroutine run_fixed_quadrature_test

end module test_gpu_fixed_quadrature_support
