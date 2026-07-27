module fortnum_gpu_batch_wrappers
    ! Benchmark-only execution wrapper. OpenACC and OpenMP target annotate the
    ! same loop and call the same generated numerical leaf.
    use fortnum_kinds, only: dp
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    implicit none
    private

    public :: dawson_value_jvp_batch

contains

    subroutine dawson_value_jvp_batch(n, x, f, v, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), f(n), v(n)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, f, v, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, f, v) map(from: values, products)
        do i = 1, n
            call fortnum_dawson_outer_kernel( &
                x(i), f(i), v(i), values(i), products(i))
        end do
    end subroutine dawson_value_jvp_batch

end module fortnum_gpu_batch_wrappers
