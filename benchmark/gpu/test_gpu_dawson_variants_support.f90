module test_gpu_dawson_variants_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_dawson_variant_wrappers, only: &
        dawson_raw_batch, dawson_simplified_batch, dawson_factored_batch
    implicit none
    private

    public :: run_dawson_variants_test

contains

    subroutine run_dawson_variants_test()
        integer, parameter :: n = 4096
        real(dp), parameter :: tolerance = 4.0e-14_dp
        real(dp) :: x(n), d(n), tx(n), td(n)
        real(dp) :: raw(n), simplified(n), factored(n)
        real(dp) :: expected, scale
        integer :: i

        do i = 1, n
            x(i) = -2.0_dp + 4.0_dp*real(mod(17*i, 4093), dp)/4092.0_dp
            d(i) = -1.5_dp + 3.0_dp*real(mod(19*i, 4091), dp)/4090.0_dp
            tx(i) = -0.8_dp + 1.6_dp*real(mod(23*i, 4079), dp)/4078.0_dp
            td(i) = -0.7_dp + 1.4_dp*real(mod(29*i, 4057), dp)/4056.0_dp
        end do

        !$acc data copyin(x, d, tx, td) create(raw, simplified, factored)
        !$omp target data map(to: x, d, tx, td) &
        !$omp& map(alloc: raw, simplified, factored)
        call dawson_raw_batch(n, x, d, tx, td, raw)
        call dawson_simplified_batch(n, x, d, tx, td, simplified)
        call dawson_factored_batch(n, x, d, tx, td, factored)
        !$omp target update from(raw, simplified, factored)
        !$omp end target data
        !$acc update self(raw, simplified, factored)
        !$acc end data

        do i = 1, n
            expected = d(i)*tx(i) + &
                td(i)*(x(i) - 2.0_dp*d(i)*exp(-d(i)*d(i)))
            scale = max(1.0_dp, abs(expected))
            if (abs(raw(i) - expected) > tolerance*scale) then
                error stop "raw Dawson variant disagrees with oracle"
            end if
            if (abs(simplified(i) - expected) > tolerance*scale) then
                error stop "simplified Dawson variant disagrees with oracle"
            end if
            if (abs(factored(i) - expected) > tolerance*scale) then
                error stop "factored Dawson variant disagrees with oracle"
            end if
        end do
    end subroutine run_dawson_variants_test

end module test_gpu_dawson_variants_support
