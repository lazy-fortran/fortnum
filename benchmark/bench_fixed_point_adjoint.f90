program bench_fixed_point_adjoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_fixed_point, only: fixed_point_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: samples = 15
    integer(int64), parameter :: analytical_reps = 2000000_int64
    integer(int64), parameter :: reference_reps = 25000_int64
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), xstar(2), map_x(2, 2), map_p(2, 2), u(2), warmup
    type(fortnum_status_t) :: status
    integer(int64) :: reps
    integer :: sample
    character(32) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
        (trim(candidate) /= "reference")) then
        error stop "usage: bench_fixed_point_adjoint analytical|reference"
    end if

    p = [0.1_dp, -0.2_dp]
    u = [1.3_dp, -0.4_dp]
    call solve_fixed_point(p, xstar)
    call map_derivatives(xstar, p, map_x, map_p)
    if (trim(mode) == "--validation-error") then
        call report_validation_error()
        stop
    end if
    if (trim(candidate) == "analytical") then
        reps = analytical_reps
    else
        reps = reference_reps
    end if

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

    subroutine report_validation_error()
        real(dp) :: analytical(2), reference(2), pp(2), pm(2), xp(2), xm(2)
        real(dp) :: error
        integer :: i

        call fixed_point_vjp(map_x, map_p, u, analytical, status)
        if (.not. status_ok(status)) error stop "fixed-point VJP failed"
        do i = 1, 2
            pp = p
            pm = p
            pp(i) = pp(i) + h
            pm(i) = pm(i) - h
            call solve_fixed_point(pp, xp)
            call solve_fixed_point(pm, xm)
            reference(i) = dot_product(u, xp - xm)/(2.0_dp*h)
        end do
        error = maxval(abs(analytical - reference))
        if (error > 1.0e-9_dp) error stop "fixed-point VJP validation failed"
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        integer :: i
        real(dp) :: jtu(2), pp(2), pm(2), xp(2), xm(2), sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            u(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "analytical") then
                call fixed_point_vjp(map_x, map_p, u, jtu, status)
            else
                do i = 1, 2
                    pp = p
                    pm = p
                    pp(i) = pp(i) + h
                    pm(i) = pm(i) - h
                    call solve_fixed_point(pp, xp)
                    call solve_fixed_point(pm, xm)
                    jtu(i) = dot_product(u, xp - xm)/(2.0_dp*h)
                end do
            end if
            sink = sink + jtu(1)
        end do
        call system_clock(tick1)
        if (trim(name) == "analytical") then
            if (.not. status_ok(status)) error stop "benchmark failed"
        end if
        if (sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine solve_fixed_point(parameters, x)
        real(dp), intent(in) :: parameters(2)
        real(dp), intent(out) :: x(2)
        real(dp) :: next(2)
        integer :: iteration

        x = 0.0_dp
        do iteration = 1, 1000
            call fixed_point_map(x, parameters, next)
            if (maxval(abs(next - x)) <= 1.0e-14_dp) then
                x = next
                return
            end if
            x = next
        end do
        error stop "fixed-point benchmark did not converge"
    end subroutine solve_fixed_point

    pure subroutine fixed_point_map(x, parameters, next)
        real(dp), intent(in) :: x(2), parameters(2)
        real(dp), intent(out) :: next(2)

        next(1) = tanh(0.2_dp*x(1) + 0.1_dp*x(2) + parameters(1))
        next(2) = tanh(0.05_dp*x(1) + 0.25_dp*x(2) + parameters(2))
    end subroutine fixed_point_map

    pure subroutine map_derivatives(x, parameters, jac_x, jac_p)
        real(dp), intent(in) :: x(2), parameters(2)
        real(dp), intent(out) :: jac_x(2, 2), jac_p(2, 2)
        real(dp) :: value(2), scale(2)

        call fixed_point_map(x, parameters, value)
        scale = 1.0_dp - value*value
        jac_x(1, :) = scale(1)*[0.2_dp, 0.1_dp]
        jac_x(2, :) = scale(2)*[0.05_dp, 0.25_dp]
        jac_p = 0.0_dp
        jac_p(1, 1) = scale(1)
        jac_p(2, 2) = scale(2)
    end subroutine map_derivatives

end program bench_fixed_point_adjoint
