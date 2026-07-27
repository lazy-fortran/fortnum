program bench_integrate_fixed_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_integrate, only: integrate, integrate_fixed_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    type :: parameter_t
        real(dp) :: value
    end type parameter_t

    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 10000_int64
    type(fortnum_status_t) :: status
    character(16) :: candidate, mode
    logical :: memory_only
    real(dp) :: warmup
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
            (trim(candidate) /= "diagnostic")) then
        error stop "usage: bench_integrate_fixed_jvp analytical|diagnostic"
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

    function run_sample(name, count) result(ns_per_jvp)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_jvp
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: derivative, minus_value, plus_value, sink
        type(parameter_t) :: parameter

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            parameter%value = 2.0_dp &
                + 0.001_dp*real(mod(iteration, 101_int64), dp)
            if (name == "analytical") then
                call integrate_fixed_jvp(integrand_tangent, 0.0_dp, 1.0_dp, &
                    derivative, status, epsrel=1.0e-10_dp, ctx=parameter)
            else
                parameter%value = parameter%value + h
                call integrate(primal_integrand, 0.0_dp, 1.0_dp, plus_value, &
                    status, epsrel=1.0e-10_dp, ctx=parameter)
                parameter%value = parameter%value - 2.0_dp*h
                call integrate(primal_integrand, 0.0_dp, 1.0_dp, minus_value, &
                    status, epsrel=1.0e-10_dp, ctx=parameter)
                derivative = (plus_value - minus_value)/(2.0_dp*h)
            end if
            sink = sink + derivative
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_jvp = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    function primal_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
        type is (parameter_t)
            value = exp(context%value*x)
        class default
            error stop "missing benchmark parameter"
        end select
    end function primal_integrand

    function integrand_tangent(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
        type is (parameter_t)
            value = x*exp(context%value*x)
        class default
            error stop "missing benchmark parameter"
        end select
    end function integrand_tangent

end program bench_integrate_fixed_jvp
