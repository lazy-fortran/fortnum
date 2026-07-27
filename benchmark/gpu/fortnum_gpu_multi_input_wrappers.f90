module fortnum_gpu_multi_input_wrappers
    use fortnum_kinds, only: dp
    use fortnum_generated_multi_input_p2_jvp, only: &
        fortnum_multi_input_p2_jvp_kernel
    use fortnum_generated_multi_input_p2_vjp, only: &
        fortnum_multi_input_p2_vjp_kernel
    use fortnum_generated_multi_input_p4_jvp, only: &
        fortnum_multi_input_p4_jvp_kernel
    use fortnum_generated_multi_input_p4_vjp, only: &
        fortnum_multi_input_p4_vjp_kernel
    use fortnum_generated_multi_input_p8_jvp, only: &
        fortnum_multi_input_p8_jvp_kernel
    use fortnum_generated_multi_input_p8_vjp, only: &
        fortnum_multi_input_p8_vjp_kernel
    use fortnum_generated_multi_input_p16_jvp, only: &
        fortnum_multi_input_p16_jvp_kernel
    use fortnum_generated_multi_input_p16_vjp, only: &
        fortnum_multi_input_p16_vjp_kernel
    implicit none
    private

    public :: multi_input_p2_jvp_batch, multi_input_p2_vjp_batch
    public :: multi_input_p4_jvp_batch, multi_input_p4_vjp_batch
    public :: multi_input_p8_jvp_batch, multi_input_p8_vjp_batch
    public :: multi_input_p16_jvp_batch, multi_input_p16_vjp_batch

contains

    subroutine multi_input_p2_jvp_batch(n, x, v, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 2), v(n, 2)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, v, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, v) map(from: values, products)
        do i = 1, n
            call fortnum_multi_input_p2_jvp_kernel( &
                x(i, 1), x(i, 2), v(i, 1), v(i, 2), values(i), products(i))
        end do
    end subroutine multi_input_p2_jvp_batch

    subroutine multi_input_p2_vjp_batch(n, x, u, values, adjoints)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 2), u(n)
        real(dp), intent(out) :: values(n), adjoints(n, 2)
        integer :: i

        !$acc parallel loop present(x, u, values, adjoints)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, u) map(from: values, adjoints)
        do i = 1, n
            call fortnum_multi_input_p2_vjp_kernel( &
                x(i, 1), x(i, 2), u(i), values(i), &
                adjoints(i, 1), adjoints(i, 2))
        end do
    end subroutine multi_input_p2_vjp_batch

    subroutine multi_input_p4_jvp_batch(n, x, v, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 4), v(n, 4)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, v, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, v) map(from: values, products)
        do i = 1, n
            call fortnum_multi_input_p4_jvp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                v(i, 1), v(i, 2), v(i, 3), v(i, 4), &
                values(i), products(i))
        end do
    end subroutine multi_input_p4_jvp_batch

    subroutine multi_input_p4_vjp_batch(n, x, u, values, adjoints)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 4), u(n)
        real(dp), intent(out) :: values(n), adjoints(n, 4)
        integer :: i

        !$acc parallel loop present(x, u, values, adjoints)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, u) map(from: values, adjoints)
        do i = 1, n
            call fortnum_multi_input_p4_vjp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), u(i), values(i), &
                adjoints(i, 1), adjoints(i, 2), adjoints(i, 3), &
                adjoints(i, 4))
        end do
    end subroutine multi_input_p4_vjp_batch

    subroutine multi_input_p8_jvp_batch(n, x, v, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 8), v(n, 8)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, v, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, v) map(from: values, products)
        do i = 1, n
            call fortnum_multi_input_p8_jvp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                v(i, 1), v(i, 2), v(i, 3), v(i, 4), &
                v(i, 5), v(i, 6), v(i, 7), v(i, 8), &
                values(i), products(i))
        end do
    end subroutine multi_input_p8_jvp_batch

    subroutine multi_input_p8_vjp_batch(n, x, u, values, adjoints)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 8), u(n)
        real(dp), intent(out) :: values(n), adjoints(n, 8)
        integer :: i

        !$acc parallel loop present(x, u, values, adjoints)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, u) map(from: values, adjoints)
        do i = 1, n
            call fortnum_multi_input_p8_vjp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                x(i, 5), x(i, 6), x(i, 7), x(i, 8), u(i), values(i), &
                adjoints(i, 1), adjoints(i, 2), adjoints(i, 3), &
                adjoints(i, 4), adjoints(i, 5), adjoints(i, 6), &
                adjoints(i, 7), adjoints(i, 8))
        end do
    end subroutine multi_input_p8_vjp_batch

    subroutine multi_input_p16_jvp_batch(n, x, v, values, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 16), v(n, 16)
        real(dp), intent(out) :: values(n), products(n)
        integer :: i

        !$acc parallel loop present(x, v, values, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, v) map(from: values, products)
        do i = 1, n
            call fortnum_multi_input_p16_jvp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                x(i, 9), x(i, 10), x(i, 11), x(i, 12), &
                x(i, 13), x(i, 14), x(i, 15), x(i, 16), &
                v(i, 1), v(i, 2), v(i, 3), v(i, 4), &
                v(i, 5), v(i, 6), v(i, 7), v(i, 8), &
                v(i, 9), v(i, 10), v(i, 11), v(i, 12), &
                v(i, 13), v(i, 14), v(i, 15), v(i, 16), &
                values(i), products(i))
        end do
    end subroutine multi_input_p16_jvp_batch

    subroutine multi_input_p16_vjp_batch(n, x, u, values, adjoints)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n, 16), u(n)
        real(dp), intent(out) :: values(n), adjoints(n, 16)
        integer :: i

        !$acc parallel loop present(x, u, values, adjoints)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, u) map(from: values, adjoints)
        do i = 1, n
            call fortnum_multi_input_p16_vjp_kernel( &
                x(i, 1), x(i, 2), x(i, 3), x(i, 4), &
                x(i, 5), x(i, 6), x(i, 7), x(i, 8), &
                x(i, 9), x(i, 10), x(i, 11), x(i, 12), &
                x(i, 13), x(i, 14), x(i, 15), x(i, 16), u(i), values(i), &
                adjoints(i, 1), adjoints(i, 2), adjoints(i, 3), &
                adjoints(i, 4), adjoints(i, 5), adjoints(i, 6), &
                adjoints(i, 7), adjoints(i, 8), adjoints(i, 9), &
                adjoints(i, 10), adjoints(i, 11), adjoints(i, 12), &
                adjoints(i, 13), adjoints(i, 14), adjoints(i, 15), &
                adjoints(i, 16))
        end do
    end subroutine multi_input_p16_vjp_batch

end module fortnum_gpu_multi_input_wrappers
