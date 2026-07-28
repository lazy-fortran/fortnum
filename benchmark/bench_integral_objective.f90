program bench_integral_objective
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_integrate, only: integrate, integrate_moving_upper_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    type :: parameter_t
        real(dp) :: value
    end type parameter_t

    integer, parameter :: samples = 15
    integer(int64), parameter :: analytical_reps = 20000_int64
    integer(int64), parameter :: diagnostic_reps = 12000_int64
    real(dp), parameter :: benchmark_parameter = 2.0_dp
    real(dp), parameter :: bound_offset = 0.5_dp
    real(dp), parameter :: bound_rate = 0.1_dp
    real(dp), parameter :: difference_step = 1.0e-5_dp
    character(32) :: candidate, mode
    integer(int64) :: reps
    integer :: sample
    real(dp) :: warmup

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "analytical" .and. &
        trim(candidate) /= "diagnostic") then
        error stop "usage: bench_integral_objective analytical|diagnostic"
    end if
    if (trim(mode) == "--validation-error") then
        call report_validation_error()
        stop
    end if

    if (trim(candidate) == "analytical") then
        reps = analytical_reps
    else
        reps = diagnostic_reps
    end if
    if (trim(mode) == "--peak-rss") then
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
        real(dp) :: value, gradient, integral, derivative, bound, exact_integral
        real(dp) :: exact_derivative, exact_value, exact_gradient, error

        call analytical_objective(benchmark_parameter, value, gradient)
        bound = bound_offset + bound_rate*benchmark_parameter
        exact_integral = (exp(benchmark_parameter*bound) - 1.0_dp)/ &
            benchmark_parameter
        exact_derivative = exp(benchmark_parameter*bound)*bound_rate + &
            ((benchmark_parameter*bound - 1.0_dp)* &
            exp(benchmark_parameter*bound) + 1.0_dp)/benchmark_parameter**2
        exact_value = 0.5_dp*exact_integral**2
        exact_gradient = exact_integral*exact_derivative
        error = max(abs(value - exact_value), abs(gradient - exact_gradient))
        if (error > 2.0e-10_dp) error stop "integral objective validation failed"

        call integral_and_tangent(benchmark_parameter, integral, derivative)
        if (abs(integral - exact_integral) > 2.0e-11_dp) then
            error stop "integral value validation failed"
        end if
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp) :: parameter, value, gradient, sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            parameter = benchmark_parameter + &
                1.0e-3_dp*real(mod(iteration, 101_int64), dp)
            if (trim(name) == "analytical") then
                call analytical_objective(parameter, value, gradient)
            else
                call diagnostic_objective(parameter, value, gradient)
            end if
            sink = sink + value + gradient
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "integral objective benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine analytical_objective(parameter, value, gradient)
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: value, gradient
        real(dp) :: integral, derivative

        call integral_and_tangent(parameter, integral, derivative)
        value = 0.5_dp*integral**2
        gradient = integral*derivative
    end subroutine analytical_objective

    subroutine diagnostic_objective(parameter, value, gradient)
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: value, gradient
        real(dp) :: integral, plus_integral, minus_integral

        integral = primal_integral(parameter)
        plus_integral = primal_integral(parameter + difference_step)
        minus_integral = primal_integral(parameter - difference_step)
        value = 0.5_dp*integral**2
        gradient = (0.5_dp*plus_integral**2 - &
            0.5_dp*minus_integral**2)/(2.0_dp*difference_step)
    end subroutine diagnostic_objective

    subroutine integral_and_tangent(parameter, integral, derivative)
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: integral, derivative
        type(parameter_t) :: context
        type(fortnum_status_t) :: status
        real(dp) :: bound

        context%value = parameter
        bound = bound_offset + bound_rate*parameter
        call integrate(primal_integrand, 0.0_dp, bound, integral, status, &
            epsrel=1.0e-10_dp, ctx=context)
        if (.not. status_ok(status)) error stop "primal integration failed"
        call integrate_moving_upper_jvp(primal_integrand, integrand_tangent, &
            0.0_dp, bound, bound_rate, derivative, status, &
            epsrel=1.0e-10_dp, ctx=context)
        if (.not. status_ok(status)) error stop "tangent integration failed"
    end subroutine integral_and_tangent

    function primal_integral(parameter) result(integral)
        real(dp), intent(in) :: parameter
        real(dp) :: integral, bound
        type(parameter_t) :: context
        type(fortnum_status_t) :: status

        context%value = parameter
        bound = bound_offset + bound_rate*parameter
        call integrate(primal_integrand, 0.0_dp, bound, integral, status, &
            epsrel=1.0e-10_dp, ctx=context)
        if (.not. status_ok(status)) error stop "primal integration failed"
    end function primal_integral

    function primal_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = exp(context%value*x)
        class default
            error stop "missing integral parameter"
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
            error stop "missing integral parameter"
        end select
    end function integrand_tangent

end program bench_integral_objective
