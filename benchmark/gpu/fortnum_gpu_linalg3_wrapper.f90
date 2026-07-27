module fortnum_gpu_linalg3_wrapper
    use fortnum_kinds, only: dp
    use fortnum_generated_det3_jvp, only: fortnum_det3_jvp_kernel
    use fortnum_generated_det3_vjp, only: fortnum_det3_vjp_kernel
    use fortnum_generated_inv3_jvp, only: fortnum_inv3_jvp_kernel
    use fortnum_generated_inv3_vjp, only: fortnum_inv3_vjp_kernel
    implicit none
    private

    public :: det3_jvp_batch, det3_vjp_batch
    public :: inv3_jvp_batch, inv3_vjp_batch

contains

    subroutine det3_jvp_batch(batch_size, matrices, directions, products)
        integer, intent(in) :: batch_size
        real(dp), intent(in) :: matrices(batch_size, 9)
        real(dp), intent(in) :: directions(batch_size, 9)
        real(dp), intent(out) :: products(batch_size)
        integer :: i

        !$acc parallel loop present(matrices, directions, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: matrices, directions) map(from: products)
        do i = 1, batch_size
            call fortnum_det3_jvp_kernel( &
                matrices(i, 1), matrices(i, 2), matrices(i, 3), &
                matrices(i, 4), matrices(i, 5), matrices(i, 6), &
                matrices(i, 7), matrices(i, 8), matrices(i, 9), &
                directions(i, 1), directions(i, 2), directions(i, 3), &
                directions(i, 4), directions(i, 5), directions(i, 6), &
                directions(i, 7), directions(i, 8), directions(i, 9), &
                products(i))
        end do
    end subroutine det3_jvp_batch

    subroutine det3_vjp_batch(batch_size, matrices, cotangents, adjoints)
        integer, intent(in) :: batch_size
        real(dp), intent(in) :: matrices(batch_size, 9)
        real(dp), intent(in) :: cotangents(batch_size)
        real(dp), intent(out) :: adjoints(batch_size, 9)
        integer :: i

        !$acc parallel loop present(matrices, cotangents, adjoints)
        !$omp target teams distribute parallel do &
        !$omp& map(to: matrices, cotangents) map(from: adjoints)
        do i = 1, batch_size
            call fortnum_det3_vjp_kernel( &
                matrices(i, 1), matrices(i, 2), matrices(i, 3), &
                matrices(i, 4), matrices(i, 5), matrices(i, 6), &
                matrices(i, 7), matrices(i, 8), matrices(i, 9), &
                cotangents(i), adjoints(i, 1), adjoints(i, 2), &
                adjoints(i, 3), adjoints(i, 4), adjoints(i, 5), &
                adjoints(i, 6), adjoints(i, 7), adjoints(i, 8), &
                adjoints(i, 9))
        end do
    end subroutine det3_vjp_batch

    subroutine inv3_jvp_batch(batch_size, inverses, directions, products)
        integer, intent(in) :: batch_size
        real(dp), intent(in) :: inverses(batch_size, 9)
        real(dp), intent(in) :: directions(batch_size, 9)
        real(dp), intent(out) :: products(batch_size, 9)
        integer :: i

        !$acc parallel loop present(inverses, directions, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: inverses, directions) map(from: products)
        do i = 1, batch_size
            call fortnum_inv3_jvp_kernel( &
                inverses(i, 1), inverses(i, 2), inverses(i, 3), &
                inverses(i, 4), inverses(i, 5), inverses(i, 6), &
                inverses(i, 7), inverses(i, 8), inverses(i, 9), &
                directions(i, 1), directions(i, 2), directions(i, 3), &
                directions(i, 4), directions(i, 5), directions(i, 6), &
                directions(i, 7), directions(i, 8), directions(i, 9), &
                products(i, 1), products(i, 2), products(i, 3), &
                products(i, 4), products(i, 5), products(i, 6), &
                products(i, 7), products(i, 8), products(i, 9))
        end do
    end subroutine inv3_jvp_batch

    subroutine inv3_vjp_batch(batch_size, inverses, cotangents, products)
        integer, intent(in) :: batch_size
        real(dp), intent(in) :: inverses(batch_size, 9)
        real(dp), intent(in) :: cotangents(batch_size, 9)
        real(dp), intent(out) :: products(batch_size, 9)
        integer :: i

        !$acc parallel loop present(inverses, cotangents, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: inverses, cotangents) map(from: products)
        do i = 1, batch_size
            call fortnum_inv3_vjp_kernel( &
                inverses(i, 1), inverses(i, 2), inverses(i, 3), &
                inverses(i, 4), inverses(i, 5), inverses(i, 6), &
                inverses(i, 7), inverses(i, 8), inverses(i, 9), &
                cotangents(i, 1), cotangents(i, 2), cotangents(i, 3), &
                cotangents(i, 4), cotangents(i, 5), cotangents(i, 6), &
                cotangents(i, 7), cotangents(i, 8), cotangents(i, 9), &
                products(i, 1), products(i, 2), products(i, 3), &
                products(i, 4), products(i, 5), products(i, 6), &
                products(i, 7), products(i, 8), products(i, 9))
        end do
    end subroutine inv3_vjp_batch

end module fortnum_gpu_linalg3_wrapper
