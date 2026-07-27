program bench_multiroot_jacobian_reuse
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 8
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 5000_int64
    real(dp) :: parameters(n)
    character(16) :: candidate, mode
    logical :: memory_only
    integer :: sample
    real(dp) :: warmup

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "reevaluate") .and. &
            (trim(candidate) /= "reuse")) then
        error stop "usage: bench_multiroot_jacobian_reuse reevaluate|reuse"
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

    function run_sample(name, count) result(ns_per_workload)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp) :: x(n), x0(n), target(n), f(n), jacobian(n, n), sink
        type(fortnum_status_t) :: status
        integer :: i

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            do i = 1, n
                target(i) = 0.5_dp + 0.04_dp*real(i, dp) &
                    + 1.0e-4_dp*real(mod(iteration + i, 17_int64) - 8_int64, dp)
            end do
            call set_parameters(target)
            x0 = target + 0.05_dp
            if (name == "reuse") then
                call multiroot_hybrid(residual, n, x0, x, &
                    status, ftol=1.0e-14_dp, jacobian=jacobian)
            else
                call multiroot_hybrid(residual, n, x0, x, &
                    status, ftol=1.0e-14_dp)
                call residual(x, f, jacobian)
            end if
            sink = sink + x(1) + jacobian(1, 1)
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine set_parameters(x)
        real(dp), intent(in) :: x(n)
        integer :: i, j

        do i = 1, n
            parameters(i) = sin(x(i))
            do j = 1, n
                parameters(i) = parameters(i) &
                    + 0.05_dp*x(j)*x(j)/real(i + j, dp)
            end do
        end do
    end subroutine set_parameters

    subroutine residual(x, f, jacobian, context)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:), jacobian(:, :)
        class(*), intent(in), optional :: context
        integer :: i, j

        do i = 1, n
            f(i) = sin(x(i)) - parameters(i)
            do j = 1, n
                f(i) = f(i) + 0.05_dp*x(j)*x(j)/real(i + j, dp)
                jacobian(i, j) = 0.1_dp*x(j)/real(i + j, dp)
            end do
            jacobian(i, i) = jacobian(i, i) + cos(x(i))
        end do
    end subroutine residual

end program bench_multiroot_jacobian_reuse
