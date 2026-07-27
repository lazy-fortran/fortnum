program bench_multiroot_condition
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_multiroot, only: multiroot_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 10000_int64
    real(dp) :: jacobian(n, n), parameter_jacobian(n, n)
    real(dp) :: direction(n), tangent(n), warmup
    type(fortnum_status_t) :: status
    integer :: sample
    character(16) :: candidate

    call get_command_argument(1, candidate)
    if ((trim(candidate) /= "plain") .and. (trim(candidate) /= "diagnostic")) then
        error stop "usage: bench_multiroot_condition plain|diagnostic"
    end if
    call initialize_inputs()

    do sample = 1, 3
        warmup = run_sample(candidate, reps/10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(candidate, reps)
    end do

contains

    subroutine initialize_inputs()
        integer :: i, j

        parameter_jacobian = 0.0_dp
        do j = 1, n
            do i = 1, n
                if (i == j) then
                    jacobian(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    jacobian(i, j) = 0.2_dp/real(i + j, dp)
                end if
            end do
            parameter_jacobian(j, j) = -1.0_dp
            direction(j) = 0.02_dp*real(mod(7*j, 9) - 4, dp)
        end do
    end subroutine initialize_inputs

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        real(dp) :: reciprocal_condition, sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            direction(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "diagnostic") then
                call multiroot_jvp(jacobian, parameter_jacobian, direction, &
                    tangent, status, reciprocal_condition=reciprocal_condition)
                sink = sink + reciprocal_condition
            else
                call multiroot_jvp(jacobian, parameter_jacobian, direction, &
                    tangent, status)
            end if
            sink = sink + tangent(1)
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_multiroot_condition
