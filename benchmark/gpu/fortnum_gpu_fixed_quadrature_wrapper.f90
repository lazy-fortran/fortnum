module fortnum_gpu_fixed_quadrature_wrapper
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: fixed_quadrature_jvp_batch, fixed_quadrature_vjp_batch

contains

    subroutine fixed_quadrature_jvp_batch( &
            batch_size, rule_order, weights, tangents, products)
        integer, intent(in) :: batch_size, rule_order
        real(dp), intent(in) :: weights(rule_order)
        real(dp), intent(in) :: tangents(batch_size, rule_order)
        real(dp), intent(out) :: products(batch_size)
        real(dp) :: product
        integer :: batch_index, node_index

        !$acc parallel loop present(weights, tangents, products) private(product)
        !$omp target teams distribute parallel do &
        !$omp& map(to: weights, tangents) map(from: products) private(product)
        do batch_index = 1, batch_size
            product = 0.0_dp
            do node_index = 1, rule_order
                product = product + &
                    weights(node_index)*tangents(batch_index, node_index)
            end do
            products(batch_index) = product
        end do
    end subroutine fixed_quadrature_jvp_batch

    subroutine fixed_quadrature_vjp_batch( &
            batch_size, rule_order, weights, cotangents, adjoints)
        integer, intent(in) :: batch_size, rule_order
        real(dp), intent(in) :: weights(rule_order)
        real(dp), intent(in) :: cotangents(batch_size)
        real(dp), intent(out) :: adjoints(batch_size, rule_order)
        integer :: batch_index, node_index

        !$acc parallel loop collapse(2) present(weights, cotangents, adjoints)
        !$omp target teams distribute parallel do collapse(2) &
        !$omp& map(to: weights, cotangents) map(from: adjoints)
        do node_index = 1, rule_order
            do batch_index = 1, batch_size
                adjoints(batch_index, node_index) = &
                    cotangents(batch_index)*weights(node_index)
            end do
        end do
    end subroutine fixed_quadrature_vjp_batch

end module fortnum_gpu_fixed_quadrature_wrapper
