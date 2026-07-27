program bench_integrate_frozen_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_integrate, only: integrate_workspace_t, integrate_result_t, &
        integrate_qag, integrate_qag_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    type :: parameter_t
        real(dp) :: value
    end type parameter_t

    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 10000_int64
    character(16) :: candidate, mode
    logical :: memory_only
    real(dp) :: warmup
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
        (trim(candidate) /= "diagnostic")) then
        error stop "usage: bench_integrate_frozen_jvp analytical|diagnostic"
    end if

    call validate_candidate()
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

    subroutine validate_candidate()
        real(dp), parameter :: parameter_value = 12.0_dp
        real(dp) :: derivative, exact
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = parameter_value
        call integrate_qag(primal_integrand, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0e-12_dp, workspace, result, status, key=15, ctx=parameter)
        call integrate_qag_jvp(integrand_tangent, result, derivative, status, &
            ctx=parameter)
        exact = (exp(parameter_value)*(parameter_value - 1.0_dp) + 1.0_dp) &
            /(parameter_value*parameter_value)
        if ((.not. status_ok(status)) .or. result%nsub <= 1 .or. &
            abs(derivative - exact) > 1.0e-7_dp*abs(exact)) then
            error stop "frozen-trace benchmark validation failed"
        end if
    end subroutine validate_candidate

    function run_sample(name, count) result(ns_per_workload)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: derivative, minus_value, plus_value, sink
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            parameter%value = 12.0_dp &
                + 0.001_dp*real(mod(iteration, 101_int64), dp)
            call integrate_qag(primal_integrand, 0.0_dp, 1.0_dp, 0.0_dp, &
                1.0e-12_dp, workspace, result, status, key=15, ctx=parameter)
            if (name == "analytical") then
                call integrate_qag_jvp(integrand_tangent, result, derivative, &
                    status, ctx=parameter)
            else
                parameter%value = parameter%value + h
                call integrate_qag_jvp(primal_integrand, result, plus_value, &
                    status, ctx=parameter)
                parameter%value = parameter%value - 2.0_dp*h
                call integrate_qag_jvp(primal_integrand, result, minus_value, &
                    status, ctx=parameter)
                derivative = (plus_value - minus_value)/(2.0_dp*h)
            end if
            sink = sink + result%value + derivative
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp &
            /(real(rate, dp)*real(count, dp))
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

end program bench_integrate_frozen_jvp
