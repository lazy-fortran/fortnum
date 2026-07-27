program bench_bspline_combined
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_jvp, bspline_eval_vjp, &
        bspline_eval_coef_jvp, bspline_eval_coef_vjp, &
        bspline_eval_combined_jvp, bspline_eval_combined_vjp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: max_n = 16, samples = 15
    integer(int64), parameter :: reps = 500000_int64
    type(bspline_workspace_t) :: ws
    real(dp) :: coef(max_n), vcoef(max_n), coefbar(max_n), breakpoints(max_n)
    real(dp) :: x, vx, u, warmup
    integer :: n, sample
    character(16) :: candidate, product, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, size_arg)
    call get_command_argument(4, mode)
    read (size_arg, *) n
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "separate") .and. &
        (trim(candidate) /= "fused")) error stop "invalid candidate"
    if ((trim(product) /= "jvp") .and. &
        (trim(product) /= "vjp")) error stop "invalid product"
    if ((n /= 4) .and. (n /= 8) .and. (n /= 16)) error stop "invalid size"
    call initialize_inputs()

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

    subroutine initialize_inputs()
        type(fortnum_status_t) :: status
        integer :: i, nbreak

        nbreak = n - 2
        do i = 1, nbreak
            breakpoints(i) = real(i - 1, dp)/real(nbreak - 1, dp)
        end do
        call bspline_init(ws, 4, nbreak, status)
        call bspline_set_knots(ws, breakpoints(:nbreak), status)
        if (status%code /= FORTNUM_OK) error stop "workspace setup failed"
        x = 0.37_dp
        vx = -0.4_dp
        u = 1.7_dp
        do i = 1, n
            coef(i) = sin(0.6_dp*real(i, dp))
            vcoef(i) = cos(0.4_dp*real(i, dp))
        end do
    end subroutine initialize_inputs

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, sink, jv, xbar, xpart, coefpart
        integer(int64) :: iteration, tick0, tick1, rate
        type(fortnum_status_t) :: status

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            vcoef(1) = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            if (trim(product) == "jvp") then
                if (trim(candidate) == "fused") then
                    call bspline_eval_combined_jvp(ws, x, coef(:n), vx, &
                        vcoef(:n), jv, status)
                else
                    call bspline_eval_jvp(ws, x, coef(:n), vx, xpart, status)
                    call bspline_eval_coef_jvp(ws, x, vcoef(:n), coefpart, status)
                    jv = xpart + coefpart
                end if
                sink = sink + jv
            else
                if (trim(candidate) == "fused") then
                    call bspline_eval_combined_vjp(ws, x, coef(:n), u, xbar, &
                        coefbar(:n), status)
                else
                    call bspline_eval_vjp(ws, x, coef(:n), u, xbar, status)
                    call bspline_eval_coef_vjp(ws, x, u, coefbar(:n), status)
                end if
                sink = sink + xbar + coefbar(1)
            end if
        end do
        call system_clock(tick1)
        if (status%code /= FORTNUM_OK .or. sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_bspline_combined
