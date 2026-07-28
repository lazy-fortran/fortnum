program bench_nonlinear_poisson
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: state_count = 64, samples = 15
    integer(int64), parameter :: reverse_reps = 200000_int64
    integer(int64), parameter :: reconstruction_reps = 5000_int64
    integer(int64), parameter :: diagnostic_reps = 2000_int64
    real(dp), parameter :: difference_step = 1.0e-5_dp
    real(dp) :: state(state_count), forcing(state_count), warmup
    character(32) :: candidate, mode
    integer(int64) :: reps
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    if (trim(candidate) /= "analytical_reverse" .and. &
        trim(candidate) /= "analytical_forward" .and. &
        trim(candidate) /= "diagnostic") then
        error stop "invalid nonlinear Poisson candidate"
    end if
    call initialize_problem()
    if (trim(mode) == "--validation-error") then
        call report_validation_error()
        stop
    end if
    select case (trim(candidate))
    case ("analytical_reverse")
        reps = reverse_reps
    case ("analytical_forward")
        reps = reconstruction_reps
    case default
        reps = diagnostic_reps
    end select
    if (trim(mode) == "--peak-rss") then
        warmup = run_sample(reps)
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if
    do sample = 1, 3
        warmup = run_sample(max(1_int64, reps/10_int64))
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(reps)
    end do

contains

    subroutine initialize_problem()
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: i

        do i = 1, state_count
            state(i) = 0.1_dp*sin(pi*real(i, dp)/real(state_count + 1, dp))
            forcing(i) = sin(2.0_dp*pi*real(i, dp)/ &
                real(state_count + 1, dp))
        end do
    end subroutine initialize_problem

    subroutine report_validation_error()
        real(dp) :: value, gradient(state_count), reference(state_count), error
        integer :: i

        call reverse_gradient(state, value, gradient)
        do i = 1, state_count
            reference(i) = complex_step_gradient(state, i)
        end do
        error = maxval(abs(gradient - reference))
        if (error > 2.0e-10_dp) error stop "Poisson gradient validation failed"
        write (*, "(es24.16)") error
    end subroutine report_validation_error

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, x(state_count), value
        real(dp) :: gradient(state_count), sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = state
            x(1) = x(1) + 1.0e-6_dp* &
                real(mod(iteration, 17_int64) - 8_int64, dp)
            select case (trim(candidate))
            case ("analytical_reverse")
                call reverse_gradient(x, value, gradient)
            case ("analytical_forward")
                call forward_reconstruction(x, value, gradient)
            case default
                call diagnostic_gradient(x, value, gradient)
            end select
            sink = sink + value + gradient(1)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "Poisson benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine reverse_gradient(x, value, gradient)
        real(dp), intent(in) :: x(state_count)
        real(dp), intent(out) :: value, gradient(state_count)
        real(dp) :: residual(state_count), scale
        integer :: i

        call residual_value(x, residual)
        value = 0.5_dp*dot_product(residual, residual)
        scale = real(state_count + 1, dp)**2
        do i = 1, state_count
            gradient(i) = (2.0_dp*scale + cos(x(i)))*residual(i)
            if (i > 1) gradient(i) = gradient(i) - scale*residual(i - 1)
            if (i < state_count) then
                gradient(i) = gradient(i) - scale*residual(i + 1)
            end if
        end do
    end subroutine reverse_gradient

    subroutine forward_reconstruction(x, value, gradient)
        real(dp), intent(in) :: x(state_count)
        real(dp), intent(out) :: value, gradient(state_count)
        real(dp) :: residual(state_count), direction(state_count)
        real(dp) :: product(state_count)
        integer :: i

        call residual_value(x, residual)
        value = 0.5_dp*dot_product(residual, residual)
        do i = 1, state_count
            direction = 0.0_dp
            direction(i) = 1.0_dp
            call residual_jvp(x, direction, product)
            gradient(i) = dot_product(residual, product)
        end do
    end subroutine forward_reconstruction

    subroutine diagnostic_gradient(x, value, gradient)
        real(dp), intent(in) :: x(state_count)
        real(dp), intent(out) :: value, gradient(state_count)
        real(dp) :: plus(state_count), minus(state_count)
        integer :: i

        value = objective(x)
        do i = 1, state_count
            plus = x
            minus = x
            plus(i) = plus(i) + difference_step
            minus(i) = minus(i) - difference_step
            gradient(i) = (objective(plus) - objective(minus))/ &
                (2.0_dp*difference_step)
        end do
    end subroutine diagnostic_gradient

    subroutine residual_value(x, residual)
        real(dp), intent(in) :: x(state_count)
        real(dp), intent(out) :: residual(state_count)
        real(dp) :: left, right, scale
        integer :: i

        scale = real(state_count + 1, dp)**2
        do i = 1, state_count
            left = 0.0_dp
            right = 0.0_dp
            if (i > 1) left = x(i - 1)
            if (i < state_count) right = x(i + 1)
            residual(i) = scale*(2.0_dp*x(i) - left - right) + &
                sin(x(i)) - forcing(i)
        end do
    end subroutine residual_value

    subroutine residual_jvp(x, direction, product)
        real(dp), intent(in) :: x(state_count), direction(state_count)
        real(dp), intent(out) :: product(state_count)
        real(dp) :: left, right, scale
        integer :: i

        scale = real(state_count + 1, dp)**2
        do i = 1, state_count
            left = 0.0_dp
            right = 0.0_dp
            if (i > 1) left = direction(i - 1)
            if (i < state_count) right = direction(i + 1)
            product(i) = scale*(2.0_dp*direction(i) - left - right) + &
                cos(x(i))*direction(i)
        end do
    end subroutine residual_jvp

    function objective(x) result(value)
        real(dp), intent(in) :: x(state_count)
        real(dp) :: value, residual(state_count)

        call residual_value(x, residual)
        value = 0.5_dp*dot_product(residual, residual)
    end function objective

    function complex_step_gradient(x, direction) result(gradient)
        real(dp), intent(in) :: x(state_count)
        integer, intent(in) :: direction
        real(dp) :: gradient
        complex(dp) :: active(state_count), residual(state_count)
        complex(dp) :: left, right, value
        real(dp), parameter :: step = 1.0e-30_dp
        real(dp) :: scale
        integer :: i

        active = cmplx(x, 0.0_dp, dp)
        active(direction) = active(direction) + cmplx(0.0_dp, step, dp)
        scale = real(state_count + 1, dp)**2
        do i = 1, state_count
            left = cmplx(0.0_dp, 0.0_dp, dp)
            right = cmplx(0.0_dp, 0.0_dp, dp)
            if (i > 1) left = active(i - 1)
            if (i < state_count) right = active(i + 1)
            residual(i) = scale*(2.0_dp*active(i) - left - right) + &
                sin(active(i)) - forcing(i)
        end do
        value = 0.5_dp*sum(residual*residual)
        gradient = aimag(value)/step
    end function complex_step_gradient

end program bench_nonlinear_poisson
