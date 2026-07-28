program bench_dawson_variants_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_dawson_variant_wrappers, only: &
        dawson_raw_batch, dawson_simplified_batch, dawson_factored_batch
    implicit none

    integer, parameter :: sample_count = 31
    real(dp), allocatable :: x(:), d(:), tx(:), td(:), products(:)
    character(16) :: variant, argument, residency, option
    integer :: i, n, sample_index, repetition, repetitions
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, variant)
    call get_command_argument(2, argument)
    read (argument, *) n
    call get_command_argument(3, residency)
    call get_command_argument(4, option)
    if (trim(variant) /= "raw" .and. trim(variant) /= "simplified" .and. &
        trim(variant) /= "factored") then
        error stop "variant must be raw, simplified, or factored"
    end if
    if (n < 1) error stop "batch size must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"
    repetitions = max(1, (65536 + n - 1)/n)

    allocate (x(n), d(n), tx(n), td(n), products(n))
    do i = 1, n
        x(i) = -2.0_dp + 4.0_dp*real(mod(17*i, 4093), dp)/4092.0_dp
        d(i) = -1.5_dp + 3.0_dp*real(mod(19*i, 4091), dp)/4090.0_dp
        tx(i) = -0.8_dp + 1.6_dp*real(mod(23*i, 4079), dp)/4078.0_dp
        td(i) = -0.7_dp + 1.4_dp*real(mod(29*i, 4057), dp)/4056.0_dp
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    sink = products(n)
    if (sink /= sink) error stop "Dawson-variant benchmark produced NaN"

contains

    subroutine run_resident()
        !$acc data copyin(x, d, tx, td) create(products)
        !$omp target data map(to: x, d, tx, td) map(alloc: products)
        call run_benchmark()
        !$omp target update from(products)
        !$omp end target data
        !$acc update self(products)
        !$acc end data
    end subroutine run_resident

    subroutine run_benchmark()
        sink = run_sample()
        if (peak_only) then
            write (*, "(i0)") peak_rss_bytes()
            return
        end if
        do sample_index = 1, 3
            sink = run_sample()
        end do
        do sample_index = 1, sample_count
            write (*, "(f0.6)") run_sample()
        end do
    end subroutine run_benchmark

    function run_sample() result(milliseconds)
        real(dp) :: milliseconds
        integer(int64) :: tick0, tick1, rate

        call system_clock(tick0, rate)
        do repetition = 1, repetitions
            if (trim(residency) == "transfer") then
                call execute_with_transfer()
            else
                call execute_variant()
            end if
        end do
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/ &
            (real(rate, dp)*real(repetitions, dp))
    end function run_sample

    subroutine execute_with_transfer()
        !$acc data copyin(x, d, tx, td) copyout(products)
        !$omp target data map(to: x, d, tx, td) map(from: products)
        call execute_variant()
        !$omp end target data
        !$acc end data
    end subroutine execute_with_transfer

    subroutine execute_variant()
        select case (trim(variant))
        case ("raw")
            call dawson_raw_batch(n, x, d, tx, td, products)
        case ("simplified")
            call dawson_simplified_batch(n, x, d, tx, td, products)
        case default
            call dawson_factored_batch(n, x, d, tx, td, products)
        end select
    end subroutine execute_variant

end program bench_dawson_variants_gpu
