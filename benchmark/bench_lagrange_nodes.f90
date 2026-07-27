program bench_lagrange_nodes
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_polynomial, only: lagrange_nodes_jvp, lagrange_nodes_vjp, &
        lagrange_weights
    implicit none

    integer, parameter :: max_n = 16
    integer, parameter :: samples = 15
    integer(int64), parameter :: analytical_reps = 1000000_int64
    integer(int64), parameter :: diagnostic_jvp_reps = 500000_int64
    integer(int64), parameter :: diagnostic_vjp_reps = 50000_int64
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: xp(max_n), f(max_n), vxp(max_n), xpbar(max_n), warmup
    integer :: n, sample
    integer(int64) :: reps
    character(16) :: candidate, product, size_arg, mode
    logical :: memory_only

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, size_arg)
    call get_command_argument(4, mode)
    read (size_arg, *) n
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
        (trim(candidate) /= "diagnostic")) then
        error stop "candidate must be analytical or diagnostic"
    end if
    if ((trim(product) /= "jvp") .and. (trim(product) /= "vjp")) then
        error stop "product must be jvp or vjp"
    end if
    if ((n < 2) .or. (n > max_n)) error stop "n must be in [2, 16]"

    call initialize_inputs()
    if (trim(candidate) == "analytical") then
        reps = analytical_reps
    else if (trim(product) == "jvp") then
        reps = diagnostic_jvp_reps
    else
        reps = diagnostic_vjp_reps
    end if
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

        do i = 1, n
            xp(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
            f(i) = sin(0.7_dp*xp(i)) + 0.1_dp*xp(i)**2
            vxp(i) = 0.2_dp*cos(real(i, dp))
        end do
    end subroutine initialize_inputs

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        real(dp) :: jv, sink
        integer(int64) :: k, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            vxp(1) = 0.02_dp*real(mod(k, 17_int64) - 8_int64, dp)
            if (trim(product) == "jvp") then
                call run_jvp(jv)
                sink = sink + jv
            else
                call run_vjp(xpbar)
                sink = sink + xpbar(1)
            end if
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine run_jvp(jv)
        real(dp), intent(out) :: jv
        real(dp) :: cp(max_n), cm(max_n), xp_plus(max_n), xp_minus(max_n)

        if (trim(candidate) == "analytical") then
            call lagrange_nodes_jvp(n, 0.17_dp, xp(:n), f(:n), vxp(:n), jv)
        else
            xp_plus(:n) = xp(:n) + h*vxp(:n)
            xp_minus(:n) = xp(:n) - h*vxp(:n)
            call lagrange_weights(n, 0.17_dp, xp_plus(:n), cp(:n))
            call lagrange_weights(n, 0.17_dp, xp_minus(:n), cm(:n))
            jv = dot_product(f(:n), cp(:n) - cm(:n))/(2.0_dp*h)
        end if
    end subroutine run_jvp

    subroutine run_vjp(result)
        real(dp), intent(out) :: result(max_n)
        real(dp) :: cp(max_n), cm(max_n), xp_work(max_n)
        integer :: i

        result = 0.0_dp
        if (trim(candidate) == "analytical") then
            call lagrange_nodes_vjp(n, 0.17_dp, xp(:n), f(:n), 1.3_dp, result(:n))
            return
        end if
        do i = 1, n
            xp_work(:n) = xp(:n)
            xp_work(i) = xp_work(i) + h
            call lagrange_weights(n, 0.17_dp, xp_work(:n), cp(:n))
            xp_work(i) = xp_work(i) - 2.0_dp*h
            call lagrange_weights(n, 0.17_dp, xp_work(:n), cm(:n))
            result(i) = 1.3_dp*dot_product(f(:n), cp(:n) - cm(:n))/(2.0_dp*h)
        end do
    end subroutine run_vjp

end program bench_lagrange_nodes
