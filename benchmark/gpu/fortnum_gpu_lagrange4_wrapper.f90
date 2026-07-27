module fortnum_gpu_lagrange4_wrapper
    use fortnum_kinds, only: dp
    use fortnum_generated_lagrange4_jvp, only: fortnum_lagrange4_jvp_kernel
    use fortnum_generated_lagrange4_vjp, only: fortnum_lagrange4_vjp_kernel
    implicit none
    private

    public :: lagrange4_jvp_batch, lagrange4_vjp_batch

contains

    subroutine lagrange4_jvp_batch( &
            n, x, samples, tx, sample_tangents, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), samples(n, 4), tx(n)
        real(dp), intent(in) :: sample_tangents(n, 4)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, samples, tx, sample_tangents, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, samples, tx, sample_tangents) &
        !$omp& map(from: values, products)
        do i = 1, n
            call fortnum_lagrange4_jvp_kernel( &
                x(i), samples(i, 1), samples(i, 2), samples(i, 3), &
                samples(i, 4), tx(i), sample_tangents(i, 1), &
                sample_tangents(i, 2), sample_tangents(i, 3), &
                sample_tangents(i, 4), values(i), products(i))
        end do
    end subroutine lagrange4_jvp_batch

    subroutine lagrange4_vjp_batch( &
            n, x, samples, cotangents, values, adjoint_x, adjoint_samples)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), samples(n, 4), cotangents(n)
        real(dp), intent(out) :: values(n), adjoint_x(n)
        real(dp), intent(out) :: adjoint_samples(n, 4)
        integer :: i

        !$acc parallel loop present(x, samples, cotangents, values, adjoint_x, adjoint_samples)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, samples, cotangents) &
        !$omp& map(from: values, adjoint_x, adjoint_samples)
        do i = 1, n
            call fortnum_lagrange4_vjp_kernel( &
                x(i), samples(i, 1), samples(i, 2), samples(i, 3), &
                samples(i, 4), cotangents(i), values(i), adjoint_x(i), &
                adjoint_samples(i, 1), adjoint_samples(i, 2), &
                adjoint_samples(i, 3), adjoint_samples(i, 4))
        end do
    end subroutine lagrange4_vjp_batch

end module fortnum_gpu_lagrange4_wrapper
