program bench_linear_solve_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: LINALG_MAX_N, LINALG_OK, linear_solve_jvp, &
        linear_solve_jvp_factored, lu_factor
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: n = LINALG_MAX_N
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 200000_int64
    real(dp) :: a(n, n), factors(n, n), x(n), da(n, n), db(n), dx(n)
    real(dp) :: warmup
    integer :: ipiv(n), info, sample
    character(16) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "refactor") .and. (trim(candidate) /= "reuse")) then
        error stop "usage: bench_linear_solve_jvp refactor|reuse"
    end if

    call initialize_inputs()
    factors = a
    call lu_factor(n, factors, ipiv, info)
    if (info /= LINALG_OK) error stop "benchmark factorization failed"

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
        integer :: i, j

        do j = 1, n
            do i = 1, n
                if (i == j) then
                    a(i, j) = 4.0_dp + real(i, dp)/real(n, dp)
                else
                    a(i, j) = 0.2_dp/real(i + j, dp)
                end if
                da(i, j) = 0.01_dp*real(mod(3*i + 5*j, 11) - 5, dp)
            end do
            x(j) = real(j, dp)/real(n, dp)
            db(j) = 0.02_dp*real(mod(7*j, 9) - 4, dp)
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
            db(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(name) == "reuse") then
                call linear_solve_jvp_factored(n, factors, ipiv, x, da, db, dx, info)
            else
                call linear_solve_jvp(n, a, x, da, db, dx, info)
            end if
            sink = sink + dx(1)
        end do
        call system_clock(tick1)
        if ((info /= LINALG_OK) .or. (sink /= sink)) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_linear_solve_jvp
