program bench_multiroot_preconditioner
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: LINALG_OK, lu_solve
    use fortnum_multiroot, only: multiroot_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 200000_int64

    type :: diagonal_preconditioner_t
        real(dp) :: inverse_diagonal(n)
    end type diagonal_preconditioner_t

    real(dp) :: jacobian(n, n), parameter_jacobian(n, n)
    real(dp) :: direction(n), tangent(n), warmup
    type(diagonal_preconditioner_t) :: preconditioner
    type(fortnum_status_t) :: status
    integer :: sample
    character(16) :: candidate

    call get_command_argument(1, candidate)
    if ((trim(candidate) /= "default") .and. &
            (trim(candidate) /= "preconditioned")) then
        error stop "usage: bench_multiroot_preconditioner default|preconditioned"
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
            preconditioner%inverse_diagonal(j) = 1.0_dp/jacobian(j, j)
        end do
    end subroutine initialize_inputs

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        real(dp) :: sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            direction(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "preconditioned") then
                call multiroot_jvp(jacobian, parameter_jacobian, direction, &
                    tangent, status, diagonal_solve, preconditioner)
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

    subroutine diagonal_solve(a, b, x, info, context)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        class(*), intent(inout), optional :: context
        real(dp) :: work_a(n, n), work_b(n)
        integer :: i

        info = 1
        if (.not. present(context)) return
        select type (preconditioner_context => context)
        type is (diagonal_preconditioner_t)
            do i = 1, n
                work_a(i, :) = preconditioner_context%inverse_diagonal(i)*a(i, :)
                work_b(i) = preconditioner_context%inverse_diagonal(i)*b(i)
            end do
        class default
            return
        end select
        call lu_solve(n, work_a, work_b, info)
        if (info == LINALG_OK) x = work_b
    end subroutine diagonal_solve

end program bench_multiroot_preconditioner
