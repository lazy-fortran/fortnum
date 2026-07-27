program bench_bspline_span_status
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_span_index, bspline_span_derivative_status
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: max_breaks = 4096, samples = 15
    integer(int64), parameter :: reps = 2000000_int64
    type(bspline_workspace_t) :: ws
    real(dp) :: breakpoints(max_breaks), warmup
    integer :: nbreak, sample
    character(16) :: candidate, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, size_arg)
    call get_command_argument(3, mode)
    read (size_arg, *) nbreak
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "plain") .and. &
        (trim(candidate) /= "status")) error stop "invalid candidate"
    if ((nbreak /= 16) .and. (nbreak /= 256) .and. &
        (nbreak /= 4096)) error stop "invalid size"
    call initialize_workspace()

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

    subroutine initialize_workspace()
        type(fortnum_status_t) :: status
        integer :: i

        do i = 1, nbreak
            breakpoints(i) = real(i - 1, dp)/real(nbreak - 1, dp)
        end do
        call bspline_init(ws, 4, nbreak, status)
        call bspline_set_knots(ws, breakpoints(:nbreak), status)
        if (status%code /= FORTNUM_OK) error stop "workspace setup failed"
    end subroutine initialize_workspace

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, x, sink
        integer(int64) :: iteration, tick0, tick1, rate
        integer :: span
        type(fortnum_status_t) :: status

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 0.25_dp + 0.5_dp*real(mod(iteration, 997_int64), dp)/997.0_dp
            if (trim(candidate) == "plain") then
                span = bspline_span_index(ws, x)
                sink = sink + real(span, dp)
            else
                call bspline_span_derivative_status(ws, x, 1.0_dp, &
                    1.0e-8_dp, status)
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

end program bench_bspline_span_status
