program bench_bspline_knots
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_basis, bspline_eval_knots_jvp, &
        bspline_eval_knots_vjp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: max_breaks = 18
    integer, parameter :: samples = 15
    integer(int64), parameter :: analytical_reps = 500000_int64
    integer(int64), parameter :: diagnostic_jvp_reps = 250000_int64
    integer(int64), parameter :: diagnostic_vjp_reps = 30000_int64
    real(dp), parameter :: h = 1.0e-5_dp
    type(bspline_workspace_t) :: ws, ws_plus, ws_minus
    real(dp), allocatable :: coef(:)
    real(dp) :: breaks(max_breaks), vbreak(max_breaks), breakbar(max_breaks)
    real(dp) :: warmup
    integer :: nbreak, sample
    integer(int64) :: reps
    character(16) :: candidate, product, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, size_arg)
    call get_command_argument(4, mode)
    read (size_arg, *) nbreak
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
        (trim(candidate) /= "diagnostic")) error stop "invalid candidate"
    if ((trim(product) /= "jvp") .and. &
        (trim(product) /= "vjp")) error stop "invalid product"
    if ((nbreak < 4) .or. (nbreak > max_breaks)) error stop "invalid nbreak"

    call initialize_inputs()
    if (trim(candidate) == "analytical") then
        reps = analytical_reps
    else if (trim(product) == "jvp") then
        reps = diagnostic_jvp_reps
    else
        reps = diagnostic_vjp_reps
    end if

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
        integer :: i

        do i = 1, nbreak
            breaks(i) = real(i - 1, dp)/real(nbreak - 1, dp)
            vbreak(i) = 0.1_dp*cos(real(i, dp))
        end do
        call bspline_init(ws, 4, nbreak, status)
        call bspline_set_knots(ws, breaks(:nbreak), status)
        call bspline_init(ws_plus, 4, nbreak, status)
        call bspline_init(ws_minus, 4, nbreak, status)
        if (status%code /= FORTNUM_OK) error stop "benchmark setup failed"
        allocate (coef(ws%ncoef))
        do i = 1, ws%ncoef
            coef(i) = sin(0.3_dp*real(i, dp))
        end do
    end subroutine initialize_inputs

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, value, sink
        integer(int64) :: i, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do i = 1_int64, count
            vbreak(1) = 0.01_dp*real(mod(i, 17_int64) - 8_int64, dp)
            if (trim(product) == "jvp") then
                call run_jvp(value)
                sink = sink + value
            else
                call run_vjp(breakbar)
                sink = sink + breakbar(1)
            end if
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine run_jvp(value)
        real(dp), intent(out) :: value
        type(fortnum_status_t) :: status
        real(dp) :: bp(max_breaks), bm(max_breaks), vp, vm

        if (trim(candidate) == "analytical") then
            call bspline_eval_knots_jvp(ws, 0.53_dp, coef, vbreak(:nbreak), &
                value, status)
        else
            bp(:nbreak) = breaks(:nbreak) + h*vbreak(:nbreak)
            bm(:nbreak) = breaks(:nbreak) - h*vbreak(:nbreak)
            call spline_value(bp(:nbreak), ws_plus, vp)
            call spline_value(bm(:nbreak), ws_minus, vm)
            value = (vp - vm)/(2.0_dp*h)
        end if
    end subroutine run_jvp

    subroutine run_vjp(result)
        real(dp), intent(out) :: result(max_breaks)
        type(fortnum_status_t) :: status
        real(dp) :: work(max_breaks), vp, vm
        integer :: i

        result = 0.0_dp
        if (trim(candidate) == "analytical") then
            call bspline_eval_knots_vjp(ws, 0.53_dp, coef, 1.3_dp, &
                result(:nbreak), status)
            return
        end if
        do i = 1, nbreak
            work(:nbreak) = breaks(:nbreak)
            work(i) = work(i) + h
            call spline_value(work(:nbreak), ws_plus, vp)
            work(i) = work(i) - 2.0_dp*h
            call spline_value(work(:nbreak), ws_minus, vm)
            result(i) = 1.3_dp*(vp - vm)/(2.0_dp*h)
        end do
    end subroutine run_vjp

    subroutine spline_value(breakpoints, local_ws, value)
        real(dp), intent(in) :: breakpoints(:)
        type(bspline_workspace_t), intent(inout) :: local_ws
        real(dp), intent(out) :: value
        type(fortnum_status_t) :: status
        real(dp) :: basis(size(coef))

        call bspline_set_knots(local_ws, breakpoints, status)
        call bspline_eval_basis(local_ws, 0.53_dp, basis, status)
        if (status%code /= FORTNUM_OK) error stop "benchmark evaluation failed"
        value = dot_product(coef, basis)
    end subroutine spline_value

end program bench_bspline_knots
