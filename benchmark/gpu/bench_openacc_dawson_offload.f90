program bench_openacc_dawson_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_gpu_batch_wrappers, only: dawson_value_jvp_batch
    implicit none

    integer, parameter :: n = 1048576
    integer, parameter :: samples = 31
    real(dp), allocatable :: x(:), f(:), v(:), values(:), products(:)
    character(16) :: candidate, mode
    integer :: i, sample
    real(dp) :: sink

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "openacc" .and. trim(candidate) /= "cpu") then
        error stop "candidate must be openacc or cpu"
    end if
    allocate (x(n), f(n), v(n), values(n), products(n))
    do i = 1, n
        x(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
        f(i) = 0.2_dp + 0.1_dp*real(mod(i, 31), dp)/30.0_dp
        v(i) = -0.5_dp + real(mod(i, 17), dp)/16.0_dp
    end do

    sink = run_sample()
    if (trim(mode) == "--peak-rss") then
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if
    do sample = 1, 3
        sink = run_sample()
    end do
    do sample = 1, samples
        write (*, "(f0.6)") run_sample()
    end do

contains

    function run_sample() result(milliseconds)
        real(dp) :: milliseconds
        integer(int64) :: tick0, tick1, rate

        call system_clock(tick0, rate)
        if (trim(candidate) == "openacc") then
            !$acc data copyin(x, f, v) copyout(values, products)
            call dawson_value_jvp_batch(n, x, f, v, values, products)
            !$acc end data
        else
            do i = 1, n
                call fortnum_dawson_outer_kernel( &
                    x(i), f(i), v(i), values(i), products(i))
            end do
        end if
        call system_clock(tick1)
        if (values(1) + products(n) /= values(1) + products(n)) then
            error stop "benchmark produced NaN"
        end if
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

end program bench_openacc_dawson_offload
