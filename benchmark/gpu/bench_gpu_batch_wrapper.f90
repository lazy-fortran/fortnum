program bench_gpu_batch_wrapper
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_gpu_batch_wrappers, only: dawson_value_jvp_batch
    implicit none

    integer, parameter :: n = 65536
    integer, parameter :: samples = 31
    integer(int64), parameter :: repetitions = 100_int64
    real(dp), allocatable :: x(:), f(:), v(:), values(:), products(:)
    character(16) :: candidate, mode
    real(dp) :: warmup
    integer :: i, sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "wrapper" .and. trim(candidate) /= "direct") then
        error stop "candidate must be wrapper or direct"
    end if

    allocate (x(n), f(n), v(n), values(n), products(n))
    do i = 1, n
        x(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
        f(i) = 0.2_dp + 0.1_dp*real(mod(i, 31), dp)/30.0_dp
        v(i) = -0.5_dp + real(mod(i, 17), dp)/16.0_dp
    end do

    warmup = run_sample(5_int64)
    if (trim(mode) == "--peak-rss") then
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if

    do sample = 1, 3
        warmup = run_sample(10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(repetitions)
    end do

contains

    function run_sample(count) result(ns_per_element)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_element
        integer(int64) :: iteration, tick0, tick1, rate
        integer :: j
        real(dp) :: sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            if (trim(candidate) == "wrapper") then
                call dawson_value_jvp_batch(n, x, f, v, values, products)
            else
                do j = 1, n
                    call fortnum_dawson_outer_kernel( &
                        x(j), f(j), v(j), values(j), products(j))
                end do
            end if
            sink = sink + values(1) + products(n)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark produced NaN"
        ns_per_element = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp)*real(n, dp))
    end function run_sample

end program bench_gpu_batch_wrapper
