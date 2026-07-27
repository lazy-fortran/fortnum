program bench_linear_solve_vjp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: LINALG_MAX_N, LINALG_OK, linear_solve_vjp, &
        linear_solve_vjp_factored, lu_factor
    implicit none

    integer, parameter :: n = LINALG_MAX_N
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 200000_int64
    real(dp) :: a(n, n), transpose_factors(n, n), x(n), u(n)
    real(dp) :: abar(n, n), bbar(n), warmup
    integer :: ipiv(n), info, sample
    character(16) :: candidate

    call get_command_argument(1, candidate)
    if ((trim(candidate) /= "refactor") .and. (trim(candidate) /= "reuse")) then
        error stop "usage: bench_linear_solve_vjp refactor|reuse"
    end if

    call initialize_inputs()
    transpose_factors = transpose(a)
    call lu_factor(n, transpose_factors, ipiv, info)
    if (info /= LINALG_OK) error stop "benchmark factorization failed"

    do sample = 1, 3
        warmup = run_sample(candidate, reps/10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(candidate, reps)
    end do

contains

    subroutine initialize_inputs()
        integer :: i, j

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    a(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    a(i, j) = 0.2_dp/real(i + j, dp)
                end if
            end do
            x(j) = real(j, dp)/real(n, dp)
            u(j) = 0.02_dp*real(mod(7*j, 9) - 4, dp)
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
            u(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "reuse") then
                call linear_solve_vjp_factored(n, transpose_factors, ipiv, x, &
                    u, abar, bbar, info)
            else
                call linear_solve_vjp(n, a, x, u, abar, bbar, info)
            end if
            sink = sink + abar(1, 1) + bbar(1)
        end do
        call system_clock(tick1)
        if ((info /= LINALG_OK) .or. (sink /= sink)) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_linear_solve_vjp
