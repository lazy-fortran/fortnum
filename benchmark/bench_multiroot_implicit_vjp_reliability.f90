program bench_multiroot_implicit_vjp_reliability
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_multiroot, only: multiroot_implicit_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 10000_int64
    real(dp) :: state(n), parameters(n), cotangent(n), parameter_adjoint(n)
    real(dp) :: warmup
    type(fortnum_status_t) :: status
    integer :: sample
    character(16) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "plain") .and. &
            (trim(candidate) /= "reliability")) then
        error stop "usage: bench_multiroot_implicit_vjp_reliability " // &
            "plain|reliability"
    end if
    call initialize_inputs()

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

    subroutine initialize_inputs()
        integer :: i

        do i = 1, n
            state(i) = 4.0_dp + real(i, dp)/real(n, dp)
            parameters(i) = 0.0_dp
            cotangent(i) = 0.02_dp*real(mod(7*i, 9) - 4, dp)
        end do
    end subroutine initialize_inputs

    function run_sample(name, count) result(ns_per_vjp)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_vjp
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp) :: reciprocal_condition, sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            cotangent(1) = 0.02_dp &
                *real(mod(iteration, 17_int64) - 8_int64, dp)
            if (name == "reliability") then
                call multiroot_implicit_vjp(state_jacobian, parameter_vjp, &
                    state, parameters, cotangent, parameter_adjoint, status, &
                    reciprocal_condition=reciprocal_condition, &
                    minimum_reciprocal_condition=1.0e-12_dp)
                sink = sink + reciprocal_condition
            else
                call multiroot_implicit_vjp(state_jacobian, parameter_vjp, &
                    state, parameters, cotangent, parameter_adjoint, status)
            end if
            sink = sink + parameter_adjoint(1)
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_vjp = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine state_jacobian(x, p, jac_x, context)
        real(dp), intent(in) :: x(:), p(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        class(*), intent(inout), optional :: context
        integer :: i, j

        jac_x = 0.0_dp
        do j = 1, n
            do i = 1, n
                if (i == j) then
                    jac_x(i, j) = x(i) + 0.0_dp*p(i)
                else
                    jac_x(i, j) = 0.2_dp/real(i + j, dp)
                end if
            end do
        end do
    end subroutine state_jacobian

    subroutine parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        f_p_t_u = -u + 0.0_dp*x + 0.0_dp*p
    end subroutine parameter_vjp

end program bench_multiroot_implicit_vjp_reliability
