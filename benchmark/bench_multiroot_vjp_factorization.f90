program bench_multiroot_vjp_factorization
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_linalg, only: LINALG_OK, lu_factor
    use fortnum_multiroot, only: multiroot_vjp, multiroot_vjp_factored
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 16
    integer, parameter :: samples = 15
    integer(int64), parameter :: reps = 200000_int64
    real(dp) :: jacobian(n, n), transpose_factors(n, n)
    real(dp) :: parameter_jacobian(n, n), cotangent(n), parameter_vjp(n)
    real(dp) :: warmup
    integer :: ipiv(n), info, sample
    type(fortnum_status_t) :: status
    character(16) :: candidate, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "refactor") .and. &
            (trim(candidate) /= "reuse")) then
        error stop "usage: bench_multiroot_vjp_factorization refactor|reuse"
    end if
    call initialize_inputs()
    transpose_factors = transpose(jacobian)
    call lu_factor(n, transpose_factors, ipiv, info)
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
            cotangent(j) = 0.02_dp*real(mod(7*j, 9) - 4, dp)
        end do
    end subroutine initialize_inputs

    function run_sample(name, count) result(ns_per_vjp)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_vjp
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp) :: sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            cotangent(1) = 0.02_dp &
                *real(mod(iteration, 17_int64) - 8_int64, dp)
            if (name == "reuse") then
                call multiroot_vjp_factored(transpose_factors, ipiv, &
                    parameter_jacobian, cotangent, parameter_vjp, status)
            else
                call multiroot_vjp(jacobian, parameter_jacobian, cotangent, &
                    parameter_vjp, status)
            end if
            sink = sink + parameter_vjp(1)
        end do
        call system_clock(tick1)
        if ((.not. status_ok(status)) .or. (sink /= sink)) then
            error stop "benchmark failed"
        end if
        ns_per_vjp = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_multiroot_vjp_factorization
