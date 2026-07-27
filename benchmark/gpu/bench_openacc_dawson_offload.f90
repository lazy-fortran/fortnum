program bench_openacc_dawson_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_generated_dawson_outer_vjp, only: &
        fortnum_dawson_outer_vjp_kernel
    use fortnum_gpu_batch_wrappers, only: dawson_value_jvp_batch, &
        dawson_vjp_batch
    implicit none

    integer, parameter :: n = 1048576
    integer, parameter :: samples = 31
    real(dp), allocatable :: x(:), f(:), v(:), u(:)
    real(dp), allocatable :: values(:), products(:), adjoints(:)
    character(16) :: candidate, mode
    integer :: i, sample
    logical :: peak_only, resident, vjp_workload
    real(dp) :: sink

    call get_command_argument(1, candidate)
    call get_command_argument(2, mode)
    select case (trim(candidate))
    case ("openacc", "cpu")
        resident = .false.
        vjp_workload = .false.
    case ("openacc_resident")
        resident = .true.
        vjp_workload = .false.
    case ("vjp", "cpu_vjp")
        resident = .false.
        vjp_workload = .true.
    case ("vjp_resident")
        resident = .true.
        vjp_workload = .true.
    case default
        error stop "unknown OpenACC benchmark candidate"
    end select
    peak_only = trim(mode) == "--peak-rss"
    allocate (x(n), f(n), v(n), u(n), values(n), products(n), adjoints(n))
    do i = 1, n
        x(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(n - 1, dp)
        f(i) = 0.2_dp + 0.1_dp*real(mod(i, 31), dp)/30.0_dp
        v(i) = -0.5_dp + real(mod(i, 17), dp)/16.0_dp
        u(i) = -0.8_dp + 1.6_dp*real(mod(i, 19), dp)/18.0_dp
    end do

    if (resident) then
        if (vjp_workload) then
            !$acc data copyin(x, f, u) create(adjoints)
            call run_benchmark()
            !$acc update self(adjoints)
            !$acc end data
        else
            !$acc data copyin(x, f, v) create(values, products)
            call run_benchmark()
            !$acc update self(values, products)
            !$acc end data
        end if
    else
        call run_benchmark()
    end if
    if (vjp_workload) then
        if (adjoints(1) + adjoints(n) /= adjoints(1) + adjoints(n)) then
            error stop "benchmark produced NaN"
        end if
    else
        if (values(1) + products(n) /= values(1) + products(n)) then
            error stop "benchmark produced NaN"
        end if
    end if

contains

    subroutine run_benchmark()
        sink = run_sample()
        if (peak_only) then
            write (*, "(i0)") peak_rss_bytes()
            return
        end if
        do sample = 1, 3
            sink = run_sample()
        end do
        do sample = 1, samples
            write (*, "(f0.6)") run_sample()
        end do
    end subroutine run_benchmark

    function run_sample() result(milliseconds)
        real(dp) :: milliseconds
        integer(int64) :: tick0, tick1, rate

        call system_clock(tick0, rate)
        if (trim(candidate) == "openacc") then
            !$acc data copyin(x, f, v) copyout(values, products)
            call dawson_value_jvp_batch(n, x, f, v, values, products)
            !$acc end data
        else if (trim(candidate) == "openacc_resident") then
            call dawson_value_jvp_batch(n, x, f, v, values, products)
        else if (trim(candidate) == "cpu") then
            do i = 1, n
                call fortnum_dawson_outer_kernel( &
                    x(i), f(i), v(i), values(i), products(i))
            end do
        else if (trim(candidate) == "vjp") then
            !$acc data copyin(x, f, u) copyout(adjoints)
            call dawson_vjp_batch(n, x, f, u, adjoints)
            !$acc end data
        else if (trim(candidate) == "vjp_resident") then
            call dawson_vjp_batch(n, x, f, u, adjoints)
        else
            do i = 1, n
                call fortnum_dawson_outer_vjp_kernel( &
                    x(i), f(i), u(i), adjoints(i))
            end do
        end if
        call system_clock(tick1)
        if (.not. resident) then
            if (vjp_workload) then
                if (adjoints(1) + adjoints(n) /= adjoints(1) + adjoints(n)) then
                    error stop "benchmark produced NaN"
                end if
            else
                if (values(1) + products(n) /= values(1) + products(n)) then
                    error stop "benchmark produced NaN"
                end if
            end if
        end if
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

end program bench_openacc_dawson_offload
