program bench_least_squares_hvp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: input_count = 8, residual_count = 24, samples = 15
    integer(int64), parameter :: analytical_reps = 300000_int64
    integer(int64), parameter :: diagnostic_reps = 50000_int64
    real(dp), parameter :: difference_step = 1.0e-5_dp
    real(dp) :: design(residual_count, input_count), target(residual_count)
    real(dp) :: parameters(input_count), direction(input_count), warmup
    character(32) :: candidate, mode
    integer(int64) :: reps
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "analytical" .and. &
        trim(candidate) /= "diagnostic") then
        error stop "usage: bench_least_squares_hvp analytical|diagnostic"
    end if
    call initialize_problem()
    if (trim(mode) == "--validation-error") then
        call report_validation_error()
        stop
    end if
    reps = merge(analytical_reps, diagnostic_reps, &
        trim(candidate) == "analytical")
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
            direction(j) = 0.04_dp*real(mod(5*j, 9) - 4, dp)
            do i = 1, residual_count
                design(i, j) = 0.08_dp*sin(0.17_dp*real(3*i + 5*j, dp))
            end do
        end do
        do i = 1, residual_count
            target(i) = 0.2_dp*cos(0.11_dp*real(2*i - 1, dp))
        end do
    end subroutine initialize_problem

    subroutine report_validation_error()
        real(dp) :: product(input_count), reference(input_count), error

        call analytical_hvp(parameters, direction, product)
        call complex_step_hvp(parameters, direction, reference)
        error = maxval(abs(product - reference))
        if (error > 2.0e-13_dp) error stop "least-squares HVP validation failed"
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, x(input_count), product(input_count), sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = parameters
            x(1) = x(1) + 1.0e-4_dp* &
                real(mod(iteration, 17_int64) - 8_int64, dp)
            if (trim(candidate) == "analytical") then
                call analytical_hvp(x, direction, product)
            else
                call diagnostic_hvp(x, direction, product)
            end if
            sink = sink + product(1)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "least-squares HVP benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine analytical_hvp(x, v, product)
        real(dp), intent(in) :: x(input_count), v(input_count)
        real(dp), intent(out) :: product(input_count)
        real(dp) :: phase(residual_count), residual(residual_count)
        real(dp) :: tangent(residual_count), weight(residual_count)

        phase = matmul(design, x)
        residual = sin(phase) - target
        tangent = matmul(design, v)
        weight = (cos(phase)**2 - residual*sin(phase))*tangent
        product = matmul(transpose(design), weight)
    end subroutine analytical_hvp

    subroutine diagnostic_hvp(x, v, product)
        real(dp), intent(in) :: x(input_count), v(input_count)
        real(dp), intent(out) :: product(input_count)
        real(dp) :: plus(input_count), minus(input_count)
        real(dp) :: plus_gradient(input_count), minus_gradient(input_count)

        plus = x + difference_step*v
        minus = x - difference_step*v
        call analytical_gradient(plus, plus_gradient)
        call analytical_gradient(minus, minus_gradient)
        product = (plus_gradient - minus_gradient)/(2.0_dp*difference_step)
    end subroutine diagnostic_hvp

    subroutine analytical_gradient(x, gradient)
        real(dp), intent(in) :: x(input_count)
        real(dp), intent(out) :: gradient(input_count)
        real(dp) :: phase(residual_count), residual(residual_count)

        phase = matmul(design, x)
        residual = sin(phase) - target
        gradient = matmul(transpose(design), residual*cos(phase))
    end subroutine analytical_gradient

    subroutine complex_step_hvp(x, v, product)
        real(dp), intent(in) :: x(input_count), v(input_count)
        real(dp), intent(out) :: product(input_count)
        complex(dp) :: active(input_count), phase(residual_count)
        complex(dp) :: residual(residual_count), gradient(input_count)
        real(dp), parameter :: step = 1.0e-30_dp

        active = cmplx(x, step*v, dp)
        phase = matmul(cmplx(design, 0.0_dp, dp), active)
        residual = sin(phase) - target
        gradient = matmul(transpose(cmplx(design, 0.0_dp, dp)), &
            residual*cos(phase))
        product = aimag(gradient)/step
    end subroutine complex_step_hvp

end program bench_least_squares_hvp
