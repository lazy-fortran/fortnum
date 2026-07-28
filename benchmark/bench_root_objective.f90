program bench_root_objective
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_roots, only: root_newton
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: samples = 15
    integer(int64), parameter :: analytical_reps = 1000000_int64
    integer(int64), parameter :: diagnostic_reps = 200000_int64
    real(dp), parameter :: difference_step = 1.0e-5_dp
    real(dp), parameter :: benchmark_parameters(2) = [0.7_dp, 1.2_dp]
    real(dp) :: callback_parameters(2), warmup
    integer(int64) :: reps
    integer :: sample
    character(32) :: candidate, mode

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "analytical" .and. &
        trim(candidate) /= "diagnostic") then
        error stop "usage: bench_root_objective analytical|diagnostic"
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
        real(dp) :: value, gradient(2), reference(2), pp(2), pm(2)
        real(dp) :: error
        integer :: direction

        call analytical_objective(benchmark_parameters, value, gradient)
        do direction = 1, 2
            pp = benchmark_parameters
            pm = benchmark_parameters
            pp(direction) = pp(direction) + difference_step
            pm(direction) = pm(direction) - difference_step
            reference(direction) = (reference_objective(pp) - &
                reference_objective(pm))/(2.0_dp*difference_step)
        end do
        error = maxval(abs(gradient - reference))
        if (error > 2.0e-10_dp) error stop "root objective validation failed"
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        real(dp) :: parameters(2), value, gradient(2), sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            parameters = benchmark_parameters
            if (trim(name) == "analytical") then
                call analytical_objective(parameters, value, gradient)
            else
                call diagnostic_objective(parameters, value, gradient)
            end if
            sink = sink + value + gradient(1)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "root objective benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine analytical_objective(parameters, value, gradient)
        real(dp), intent(in) :: parameters(2)
        real(dp), intent(out) :: value, gradient(2)
        real(dp) :: root, residual_derivative

        call solve_primal(parameters, root)
        value = 0.5_dp*root*root + 0.1_dp*parameters(1)**2
        residual_derivative = 3.0_dp*root*root + parameters(1)
        gradient(1) = -root*root/residual_derivative + &
            0.2_dp*parameters(1)
        gradient(2) = root/residual_derivative
    end subroutine analytical_objective

    subroutine diagnostic_objective(parameters, value, gradient)
        real(dp), intent(in) :: parameters(2)
        real(dp), intent(out) :: value, gradient(2)
        real(dp) :: pp(2), pm(2)
        integer :: direction

        value = primal_objective(parameters)
        do direction = 1, 2
            pp = parameters
            pm = parameters
            pp(direction) = pp(direction) + difference_step
            pm(direction) = pm(direction) - difference_step
            gradient(direction) = (primal_objective(pp) - &
                primal_objective(pm))/(2.0_dp*difference_step)
        end do
    end subroutine diagnostic_objective

    function primal_objective(parameters) result(value)
        real(dp), intent(in) :: parameters(2)
        real(dp) :: value, root

        call solve_primal(parameters, root)
        value = 0.5_dp*root*root + 0.1_dp*parameters(1)**2
    end function primal_objective

    subroutine solve_primal(parameters, root)
        real(dp), intent(in) :: parameters(2)
        real(dp), intent(out) :: root
        type(fortnum_status_t) :: status

        callback_parameters = parameters
        call root_newton(residual_and_derivative, 0.0_dp, 2.0_dp, 1.0_dp, &
            root, status, xtol=1.0e-12_dp, ftol=1.0e-12_dp)
        if (.not. status_ok(status)) error stop "root solve failed"
    end subroutine solve_primal

    pure subroutine residual_and_derivative(x, value, derivative)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative

        value = x**3 + callback_parameters(1)*x - callback_parameters(2)
        derivative = 3.0_dp*x*x + callback_parameters(1)
    end subroutine residual_and_derivative

    function reference_objective(parameters) result(value)
        real(dp), intent(in) :: parameters(2)
        real(dp) :: value, lo, hi, midpoint, residual
        integer :: iteration

        lo = 0.0_dp
        hi = 2.0_dp
        do iteration = 1, 100
            midpoint = lo + 0.5_dp*(hi - lo)
            residual = midpoint**3 + parameters(1)*midpoint - parameters(2)
            if (residual > 0.0_dp) then
                hi = midpoint
            else
                lo = midpoint
            end if
        end do
        midpoint = lo + 0.5_dp*(hi - lo)
        value = 0.5_dp*midpoint*midpoint + 0.1_dp*parameters(1)**2
    end function reference_objective

end program bench_root_objective
