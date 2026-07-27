program bench_lagrange_combined
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_polynomial, only: lagrange_weights_jvp, lagrange_weights_vjp, &
        lagrange_fval_jvp, lagrange_fval_vjp, lagrange_combined_jvp, &
        lagrange_combined_vjp
    implicit none

    integer, parameter :: max_n = 16, samples = 15
    integer(int64), parameter :: reps = 100000_int64
    real(dp) :: xp(max_n), f(max_n), vf(max_n), fbar(max_n), warmup
    real(dp) :: x, vx, u
    integer :: n, sample
    character(16) :: candidate, product, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, size_arg)
    call get_command_argument(4, mode)
    read (size_arg, *) n
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "separate") .and. &
        (trim(candidate) /= "fused")) error stop "invalid candidate"
    if ((trim(product) /= "jvp") .and. &
        (trim(product) /= "vjp")) error stop "invalid product"
    if ((n /= 4) .and. (n /= 8) .and. (n /= 16)) error stop "invalid size"
    call initialize_inputs()

    if (memory_only) then
        warmup = run_sample(reps/10_int64)
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

    subroutine initialize_inputs()
        integer :: i

        x = 0.37_dp
        vx = -0.4_dp
        u = 1.7_dp
        do i = 1, n
            xp(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
            f(i) = sin(0.6_dp*real(i, dp))
            vf(i) = cos(0.4_dp*real(i, dp))
        end do
    end subroutine initialize_inputs

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, sink, jv, xbar
        real(dp) :: xpart, fpart
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            vf(1) = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            if (trim(product) == "jvp") then
                if (trim(candidate) == "fused") then
                    call lagrange_combined_jvp(n, x, xp(:n), f(:n), vx, &
                        vf(:n), jv)
                else
                    call lagrange_weights_jvp(n, x, xp(:n), f(:n), vx, xpart)
                    call lagrange_fval_jvp(n, x, xp(:n), vf(:n), fpart)
                    jv = xpart + fpart
                end if
                sink = sink + jv
            else
                if (trim(candidate) == "fused") then
                    call lagrange_combined_vjp(n, x, xp(:n), f(:n), u, xbar, &
                        fbar(:n))
                else
                    call lagrange_weights_vjp(n, x, xp(:n), f(:n), u, xbar)
                    call lagrange_fval_vjp(n, x, xp(:n), u, fbar(:n))
                end if
                sink = sink + xbar + fbar(1)
            end if
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_lagrange_combined
