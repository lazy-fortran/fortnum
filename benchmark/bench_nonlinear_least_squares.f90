program bench_nonlinear_least_squares
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: input_count = 8, residual_count = 24
    integer, parameter :: samples = 15
    integer(int64), parameter :: fast_reps = 200000_int64
    integer(int64), parameter :: diagnostic_reps = 20000_int64
    real(dp), parameter :: difference_step = 1.0e-5_dp
    real(dp) :: design(residual_count, input_count), target(residual_count)
    real(dp) :: parameters(input_count), warmup
    character(32) :: candidate, mode
    integer(int64) :: reps
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "analytical_reverse" .and. &
        trim(candidate) /= "analytical_forward" .and. &
        trim(candidate) /= "diagnostic") then
        error stop "usage: bench_nonlinear_least_squares " // &
            "analytical_reverse|analytical_forward|diagnostic"
    end if
    call initialize_problem()
    if (trim(mode) == "--validation-error") then
        call report_validation_error()
        stop
    end if
    reps = merge(diagnostic_reps, fast_reps, trim(candidate) == "diagnostic")
    if (trim(mode) == "--peak-rss") then
        warmup = run_sample(reps)
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

    subroutine initialize_problem()
        integer :: i, j

        do j = 1, input_count
            parameters(j) = 0.03_dp*real(j - 4, dp)
            do i = 1, residual_count
                design(i, j) = 0.08_dp*sin(0.17_dp*real(3*i + 5*j, dp))
            end do
        end do
        do i = 1, residual_count
            target(i) = 0.2_dp*cos(0.11_dp*real(2*i - 1, dp))
        end do
    end subroutine initialize_problem

    subroutine report_validation_error()
        real(dp) :: value, gradient(input_count), reference(input_count), error
        integer :: j

        call analytical_reverse(parameters, value, gradient)
        do j = 1, input_count
            reference(j) = complex_step_gradient(parameters, j)
        end do
        error = maxval(abs(gradient - reference))
        if (error > 2.0e-12_dp) error stop "least-squares validation failed"
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, x(input_count), value, gradient(input_count)
        real(dp) :: sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = parameters
            x(1) = x(1) + 1.0e-4_dp* &
                real(mod(iteration, 17_int64) - 8_int64, dp)
            select case (trim(candidate))
            case ("analytical_reverse")
                call analytical_reverse(x, value, gradient)
            case ("analytical_forward")
                call analytical_forward(x, value, gradient)
            case default
                call diagnostic_gradient(x, value, gradient)
            end select
            sink = sink + value + gradient(1)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "least-squares benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine analytical_reverse(x, value, gradient)
        real(dp), intent(in) :: x(input_count)
        real(dp), intent(out) :: value, gradient(input_count)
        real(dp) :: phase(residual_count), residual(residual_count), weight
        integer :: i, j

        phase = matmul(design, x)
        residual = sin(phase) - target
        value = 0.5_dp*dot_product(residual, residual)
        gradient = 0.0_dp
        do i = 1, residual_count
            weight = residual(i)*cos(phase(i))
            do j = 1, input_count
                gradient(j) = gradient(j) + design(i, j)*weight
            end do
        end do
    end subroutine analytical_reverse

    subroutine analytical_forward(x, value, gradient)
        real(dp), intent(in) :: x(input_count)
        real(dp), intent(out) :: value, gradient(input_count)
        real(dp) :: phase(residual_count), residual(residual_count)
        integer :: j

        phase = matmul(design, x)
        residual = sin(phase) - target
        value = 0.5_dp*dot_product(residual, residual)
        do j = 1, input_count
            gradient(j) = dot_product(residual, &
                cos(phase)*design(:, j))
        end do
    end subroutine analytical_forward

    subroutine diagnostic_gradient(x, value, gradient)
        real(dp), intent(in) :: x(input_count)
        real(dp), intent(out) :: value, gradient(input_count)
        real(dp) :: plus(input_count), minus(input_count)
        integer :: j

        value = objective(x)
        do j = 1, input_count
            plus = x
            minus = x
            plus(j) = plus(j) + difference_step
            minus(j) = minus(j) - difference_step
            gradient(j) = (objective(plus) - objective(minus))/ &
                (2.0_dp*difference_step)
        end do
    end subroutine diagnostic_gradient

    function objective(x) result(value)
        real(dp), intent(in) :: x(input_count)
        real(dp) :: value, residual(residual_count)

        residual = sin(matmul(design, x)) - target
        value = 0.5_dp*dot_product(residual, residual)
    end function objective

    function complex_step_gradient(x, direction) result(gradient)
        real(dp), intent(in) :: x(input_count)
        integer, intent(in) :: direction
        real(dp) :: gradient
        complex(dp) :: active(input_count), residual(residual_count), value
        real(dp), parameter :: step = 1.0e-30_dp

        active = cmplx(x, 0.0_dp, dp)
        active(direction) = active(direction) + cmplx(0.0_dp, step, dp)
        residual = sin(matmul(cmplx(design, 0.0_dp, dp), active)) - target
        value = 0.5_dp*sum(residual*residual)
        gradient = aimag(value)/step
    end function complex_step_gradient

end program bench_nonlinear_least_squares
