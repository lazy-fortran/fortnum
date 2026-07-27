program bench_vector_root_adjoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_multiroot, only: multiroot_implicit_vjp, multiroot_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 2000000_int64
    real(dp) :: xstar(2), p(2), jac_x(2, 2), f_p(2, 2), u(2), warmup
    type(fortnum_status_t) :: status
    integer :: sample
    character(16) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "assembled") .and. &
            (trim(candidate) /= "boundary")) then
        error stop "usage: bench_vector_root_adjoint assembled|boundary"
    end if

    xstar = [1.2_dp, 0.8_dp]
    p = [0.7_dp, 0.3_dp]
    u = [0.4_dp, -0.6_dp]
    call assemble_products(xstar, p, jac_x, f_p)

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

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        real(dp) :: jtu(2), sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            u(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "boundary") then
                call multiroot_implicit_vjp(state_jacobian, parameter_vjp, &
                    xstar, p, u, jtu, status)
            else
                call multiroot_vjp(jac_x, f_p, u, jtu, status)
            end if
            sink = sink + jtu(1)
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine state_jacobian(x, parameters, jacobian, context)
        real(dp), intent(in) :: x(:), parameters(:)
        real(dp), intent(out) :: jacobian(size(x), size(x))
        class(*), intent(inout), optional :: context

        jacobian(1, 1) = 2.0_dp*x(1) + parameters(1)
        jacobian(1, 2) = 1.0_dp
        jacobian(2, 1) = 1.0_dp
        jacobian(2, 2) = 2.0_dp*x(2) + parameters(2)
    end subroutine state_jacobian

    subroutine parameter_vjp(x, parameters, residual_u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), parameters(:), residual_u(:)
        real(dp), intent(out) :: f_p_t_u(size(parameters))
        class(*), intent(inout), optional :: context

        f_p_t_u(1) = (x(1) - 1.0_dp)*residual_u(1)
        f_p_t_u(2) = (x(2) - 1.0_dp)*residual_u(2)
    end subroutine parameter_vjp

    subroutine assemble_products(x, parameters, jacobian, parameter_jacobian)
        real(dp), intent(in) :: x(2), parameters(2)
        real(dp), intent(out) :: jacobian(2, 2), parameter_jacobian(2, 2)

        jacobian(1, 1) = 2.0_dp*x(1) + parameters(1)
        jacobian(1, 2) = 1.0_dp
        jacobian(2, 1) = 1.0_dp
        jacobian(2, 2) = 2.0_dp*x(2) + parameters(2)
        parameter_jacobian = 0.0_dp
        parameter_jacobian(1, 1) = x(1) - 1.0_dp
        parameter_jacobian(2, 2) = x(2) - 1.0_dp
    end subroutine assemble_products

end program bench_vector_root_adjoint
