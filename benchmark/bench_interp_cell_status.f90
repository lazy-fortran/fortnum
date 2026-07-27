program bench_interp_cell_status
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_interp, only: grid_search, grid_search_derivative_status
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: max_n = 4096, samples = 15
    integer(int64), parameter :: reps = 2000000_int64
    real(dp) :: grid(max_n), warmup
    integer :: n, sample
    character(16) :: candidate, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, size_arg)
    call get_command_argument(3, mode)
    read (size_arg, *) n
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "plain") .and. &
        (trim(candidate) /= "status")) error stop "invalid candidate"
    if ((n /= 16) .and. (n /= 256) .and. (n /= 4096)) error stop "invalid size"
    call initialize_grid()

    if (memory_only) then
        warmup = run_sample(reps/10_int64)
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if
    do sample = 1, 3
        warmup = run_sample(reps/10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(reps)
    end do

contains

    subroutine initialize_grid()
        integer :: i

        do i = 1, n
            grid(i) = real(i - 1, dp)/real(n - 1, dp)
        end do
    end subroutine initialize_grid

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, x, sink
        integer(int64) :: iteration, tick0, tick1, rate
        integer :: cell
        type(fortnum_status_t) :: status

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 0.25_dp + 0.5_dp*real(mod(iteration, 997_int64), dp)/997.0_dp
            if (trim(candidate) == "plain") then
                call grid_search(grid(:n), 1, n, x, cell)
                sink = sink + real(cell, dp)
            else
                call grid_search_derivative_status(grid(:n), 1, n, x, &
                    1.0_dp, 1.0e-8_dp, status)
                sink = sink + real(status%code, dp)
            end if
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark failed"
        if (trim(candidate) == "status") then
            if (status%code /= FORTNUM_OK) error stop "unexpected crossing"
        end if
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_interp_cell_status
