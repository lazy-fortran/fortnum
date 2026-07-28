program bench_dawson_generated_family
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_generated_dawson_outer_jvp, only: &
        fortnum_dawson_outer_jvp_kernel
    use fortnum_generated_dawson_outer_value, only: &
        fortnum_dawson_outer_value_kernel
    use fortnum_generated_dawson_outer_value_vjp, only: &
        fortnum_dawson_outer_value_vjp_kernel
    use fortnum_generated_dawson_outer_vjp, only: &
        fortnum_dawson_outer_vjp_kernel
    use fortnum_special_dawson, only: dawson
    implicit none

    integer, parameter :: samples = 31
    integer(int64), parameter :: reps = 10000000_int64
    character(16) :: candidate, product, mode
    real(dp) :: warmup
    integer :: sample

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, mode)
    if (trim(candidate) /= "fused" .and. trim(candidate) /= "separate") then
        error stop "candidate must be fused or separate"
    end if
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if

    call validate_candidates()
    if (trim(mode) == "--validate") stop

    if (trim(mode) == "--peak-rss") then
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

    subroutine validate_candidates()
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: x = 0.9_dp
        real(dp), parameter :: direction = -0.4_dp
        real(dp), parameter :: cotangent = 1.7_dp
        real(dp) :: f, f_plus, f_minus, value, separate_value
        real(dp) :: jvp, separate_jvp, vjp, separate_vjp, reference

        f = dawson(x)
        f_plus = dawson(x + h*direction)
        f_minus = dawson(x - h*direction)
        reference = ((sin(f_plus) + f_plus**2) - &
            (sin(f_minus) + f_minus**2))/(2.0_dp*h)
        call fortnum_dawson_outer_kernel(x, f, direction, value, jvp)
        call fortnum_dawson_outer_value_kernel(f, separate_value)
        call fortnum_dawson_outer_jvp_kernel(x, f, direction, separate_jvp)
        call fortnum_dawson_outer_value_vjp_kernel( &
            x, f, cotangent, value, vjp)
        call fortnum_dawson_outer_vjp_kernel( &
            x, f, cotangent, separate_vjp)

        if (abs(value - separate_value) > 1.0e-13_dp) &
            error stop "fused and separate values disagree"
        if (abs(jvp - separate_jvp) > 1.0e-13_dp) &
            error stop "fused and separate JVPs disagree"
        if (abs(jvp - reference) > 1.0e-7_dp) &
            error stop "JVP disagrees with central difference"
        if (abs(vjp - separate_vjp) > 1.0e-13_dp) &
            error stop "fused and separate VJPs disagree"
        if (abs(cotangent*jvp - direction*vjp) > 1.0e-13_dp) &
            error stop "JVP and VJP fail the adjoint identity"
    end subroutine validate_candidates

    function run_sample(count) result(ns_per_workload)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload

        if (trim(product) == "jvp") then
            if (trim(candidate) == "fused") then
                ns_per_workload = time_fused_jvp(count)
            else
                ns_per_workload = time_separate_jvp(count)
            end if
        else
            if (trim(candidate) == "fused") then
                ns_per_workload = time_fused_vjp(count)
            else
                ns_per_workload = time_separate_vjp(count)
            end if
        end if
    end function run_sample

    function time_fused_jvp(count) result(ns_per_workload)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        real(dp) :: x, value, product_value, sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 1.0_dp + 1.0e-12_dp*real(mod(iteration, 1024_int64), dp)
            call fortnum_dawson_outer_kernel( &
                x, 0.5380795069127684_dp, -0.4_dp, value, product_value)
            sink = sink + value + product_value
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark produced NaN"
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function time_fused_jvp

    function time_separate_jvp(count) result(ns_per_workload)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        real(dp) :: x, value, product_value, sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 1.0_dp + 1.0e-12_dp*real(mod(iteration, 1024_int64), dp)
            call fortnum_dawson_outer_value_kernel( &
                0.5380795069127684_dp, value)
            call fortnum_dawson_outer_jvp_kernel( &
                x, 0.5380795069127684_dp, -0.4_dp, product_value)
            sink = sink + value + product_value
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark produced NaN"
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function time_separate_jvp

    function time_fused_vjp(count) result(ns_per_workload)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        real(dp) :: x, value, product_value, sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 1.0_dp + 1.0e-12_dp*real(mod(iteration, 1024_int64), dp)
            call fortnum_dawson_outer_value_vjp_kernel( &
                x, 0.5380795069127684_dp, 1.7_dp, value, product_value)
            sink = sink + value + product_value
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark produced NaN"
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function time_fused_vjp

    function time_separate_vjp(count) result(ns_per_workload)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_workload
        real(dp) :: x, value, product_value, sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            x = 1.0_dp + 1.0e-12_dp*real(mod(iteration, 1024_int64), dp)
            call fortnum_dawson_outer_value_kernel( &
                0.5380795069127684_dp, value)
            call fortnum_dawson_outer_vjp_kernel( &
                x, 0.5380795069127684_dp, 1.7_dp, product_value)
            sink = sink + value + product_value
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark produced NaN"
        ns_per_workload = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function time_separate_vjp

end program bench_dawson_generated_family
