program bench_stable_sqrt_difference
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_generated_sqrt1pm1_raw, only: fortnum_sqrt1pm1_raw_kernel
    use fortnum_generated_sqrt1pm1_stable, only: fortnum_sqrt1pm1_stable_kernel
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: batch_size = 1048576
    integer, parameter :: repetitions = 32
    real(dp), allocatable :: x(:), y(:)
    real(dp) :: elapsed_ns, sink
    integer(int64) :: tick0, tick1, rate
    character(32) :: candidate, action
    integer :: i, repetition

    call get_command_argument(1, candidate)
    call get_command_argument(2, action)
    if (trim(candidate) /= "raw" .and. trim(candidate) /= "stable") then
        error stop "usage: bench_stable_sqrt_difference raw|stable [--peak-rss]"
    end if

    allocate (x(batch_size), y(batch_size))
    do i = 1, batch_size
        x(i) = epsilon(1.0_dp)*real(1 + mod(i, 1024), dp)
    end do
    sink = 0.0_dp
    call system_clock(tick0, rate)
    do repetition = 1, repetitions
        if (trim(candidate) == "raw") then
            call fortnum_sqrt1pm1_raw_kernel(x, y)
        else
            call fortnum_sqrt1pm1_stable_kernel(x, y)
        end if
        sink = sink + sum(y)
    end do
    call system_clock(tick1)
    if (sink /= sink) error stop "square-root benchmark produced NaN"

    if (trim(action) == "--peak-rss") then
        write (*, "(i0)") peak_rss_bytes()
    else
        elapsed_ns = real(tick1 - tick0, dp)*1.0e9_dp/real(rate, dp)
        write (*, "(f0.6)") elapsed_ns/real(batch_size*repetitions, dp)
    end if

end program bench_stable_sqrt_difference
