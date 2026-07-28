module fortnum_gpu_dawson_variant_wrappers
    use fortnum_kinds, only: dp
    use fortnum_dawson_variant_raw, only: fortnum_dawson_jvp_raw
    use fortnum_dawson_variant_simplified, only: fortnum_dawson_jvp_simplified
    use fortnum_dawson_variant_factored, only: fortnum_dawson_jvp_factored
    implicit none
    private

    public :: dawson_raw_batch, dawson_simplified_batch, dawson_factored_batch

contains

    subroutine dawson_raw_batch(n, x, d, tx, td, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), d(n), tx(n), td(n)
        real(dp), intent(out) :: products(n)
        integer :: i

        !$acc parallel loop present(x, d, tx, td, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, d, tx, td) map(from: products)
        do i = 1, n
            call fortnum_dawson_jvp_raw(x(i), d(i), tx(i), td(i), products(i))
        end do
    end subroutine dawson_raw_batch

    subroutine dawson_simplified_batch(n, x, d, tx, td, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), d(n), tx(n), td(n)
        real(dp), intent(out) :: products(n)
        integer :: i

        !$acc parallel loop present(x, d, tx, td, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, d, tx, td) map(from: products)
        do i = 1, n
            call fortnum_dawson_jvp_simplified( &
                x(i), d(i), tx(i), td(i), products(i))
        end do
    end subroutine dawson_simplified_batch

    subroutine dawson_factored_batch(n, x, d, tx, td, products)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), d(n), tx(n), td(n)
        real(dp), intent(out) :: products(n)
        integer :: i

        !$acc parallel loop present(x, d, tx, td, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: x, d, tx, td) map(from: products)
        do i = 1, n
            call fortnum_dawson_jvp_factored( &
                x(i), d(i), tx(i), td(i), products(i))
        end do
    end subroutine dawson_factored_batch

end module fortnum_gpu_dawson_variant_wrappers
