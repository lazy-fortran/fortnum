module fortnum_gpu_batch_wrappers
    ! Benchmark-only execution wrapper. OpenACC and OpenMP target annotate the
    ! same loop and call the same generated numerical leaf.
    use fortnum_kinds, only: dp
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_generated_dawson_outer_jvp, only: &
        fortnum_dawson_outer_jvp_kernel
    use fortnum_generated_dawson_outer_value, only: &
        fortnum_dawson_outer_value_kernel
    use fortnum_generated_dawson_outer_value_vjp, only: &
        fortnum_dawson_outer_value_vjp_kernel
    use fortnum_generated_dawson_outer_vjp, only: &
        fortnum_dawson_outer_vjp_kernel
    implicit none
    private

    public :: dawson_value_batch, dawson_jvp_batch, dawson_value_jvp_batch
    public :: dawson_vjp_batch, dawson_value_vjp_batch

contains

    subroutine dawson_value_batch(n, f, values)
        integer, intent(in) :: n
        real(dp), intent(in) :: f(n)
        real(dp), intent(out) :: values(n)
        integer :: i

        !$acc parallel loop present(f, values)
        !$omp target teams distribute parallel do &
        !$omp& map(to: f) map(from: values)
        do i = 1, n
            call fortnum_dawson_outer_value_kernel(f(i), values(i))
        end do
    end subroutine dawson_value_batch

    subroutine dawson_jvp_batch(n, x, f, v, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), f(n), v(n)
        real(dp), intent(out) :: products(n)
        integer :: i

        !$acc parallel loop present(x, f, v, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, f, v) map(from: products)
        do i = 1, n
            call fortnum_dawson_outer_jvp_kernel( &
                x(i), f(i), v(i), products(i))
        end do
    end subroutine dawson_jvp_batch

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

    subroutine dawson_vjp_batch(n, x, f, u, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), f(n), u(n)
        real(dp), intent(out) :: products(n)
        integer :: i

        !$acc parallel loop present(x, f, u, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, f, u) map(from: products)
        do i = 1, n
            call fortnum_dawson_outer_vjp_kernel( &
                x(i), f(i), u(i), products(i))
        end do
    end subroutine dawson_vjp_batch

    subroutine dawson_value_vjp_batch(n, x, f, u, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), f(n), u(n)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, f, u, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, f, u) map(from: values, products)
        do i = 1, n
            call fortnum_dawson_outer_value_vjp_kernel( &
                x(i), f(i), u(i), values(i), products(i))
        end do
    end subroutine dawson_value_vjp_batch

end module fortnum_gpu_batch_wrappers
