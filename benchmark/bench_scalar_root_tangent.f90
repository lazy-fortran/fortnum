program bench_scalar_root_tangent
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_roots, only: root_implicit_jvp, root_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 2000000_int64
    real(dp), parameter :: xstar = 1.213721857123380_dp
    real(dp) :: p(2), tp(2), warmup
    type(fortnum_status_t) :: status
    integer :: sample
    character(16) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "assembled") .and. &
            (trim(candidate) /= "boundary")) then
        error stop "usage: bench_scalar_root_tangent assembled|boundary"
    end if

    p = [0.7_dp, xstar**3 + 0.7_dp*xstar]
    tp = [0.4_dp, -0.6_dp]

    if (memory_only) then
        warmup = run_sample(candidate, reps)
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if

    do sample = 1, 3
        warmup = run_sample(candidate, reps/10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(candidate, reps)
    end do

contains

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        real(dp) :: dx, f_p(2), sink

        sink = 0.0_dp
        f_p = [xstar, -1.0_dp]
        call system_clock(tick0, rate)
        do k = 1_int64, count
            tp(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "boundary") then
                call root_implicit_jvp(residual_jvp, xstar, p, tp, dx, status)
            else
                call root_jvp(3.0_dp*xstar*xstar + p(1), f_p, tp, &
                    dx, status)
            end if
            sink = sink + dx
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine residual_jvp(x, parameters, direction, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, parameters(:), direction(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        f_x = 3.0_dp*x*x + parameters(1)
        f_p_tp = x*direction(1) - direction(2)
    end subroutine residual_jvp

end program bench_scalar_root_tangent
