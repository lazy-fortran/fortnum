program bench_dawson_gpu_fusion
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_batch_wrappers, only: dawson_value_batch, &
        dawson_jvp_batch, dawson_value_jvp_batch, dawson_vjp_batch, &
        dawson_value_vjp_batch
    implicit none

    integer, parameter :: samples = 31
    real(dp), allocatable :: x(:), f(:), v(:), u(:), values(:), products(:)
    character(16) :: candidate, residency, size_argument, option
    integer :: i, n, sample
    logical :: peak_only, vjp_workload
    real(dp) :: sink

    call get_command_argument(1, candidate)
    call get_command_argument(2, residency)
    call get_command_argument(3, size_argument)
    call get_command_argument(4, option)
    select case (trim(candidate))
    case ("fused_jvp", "separate_jvp")
        vjp_workload = .false.
    case ("fused_vjp", "separate_vjp")
        vjp_workload = .true.
    case default
        error stop "unknown fusion candidate"
    end select
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    read (size_argument, *) n
    if (n < 1) error stop "batch size must be positive"
    peak_only = trim(option) == "--peak-rss"

    allocate (x(n), f(n), v(n), u(n), values(n), products(n))
    do i = 1, n
        x(i) = -1.0_dp + 2.0_dp*real(i - 1, dp)/real(max(1, n - 1), dp)
        f(i) = 0.2_dp + 0.1_dp*real(mod(i, 31), dp)/30.0_dp
        v(i) = -0.5_dp + real(mod(i, 17), dp)/16.0_dp
        u(i) = -0.8_dp + 1.6_dp*real(mod(i, 19), dp)/18.0_dp
    end do

    if (trim(residency) == "resident") then
        if (vjp_workload) then
            !$acc data copyin(x, f, u) create(values, products)
            !$omp target data map(to: x, f, u) map(alloc: values, products)
            call run_benchmark()
            !$omp target update from(values, products)
            !$omp end target data
            !$acc update self(values, products)
            !$acc end data
        else
            !$acc data copyin(x, f, v) create(values, products)
            !$omp target data map(to: x, f, v) map(alloc: values, products)
            call run_benchmark()
            !$omp target update from(values, products)
            !$omp end target data
            !$acc update self(values, products)
            !$acc end data
        end if
    else
        call run_benchmark()
    end if
    if (values(1) + products(n) /= values(1) + products(n)) then
        error stop "benchmark produced NaN"
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
        if (trim(residency) == "transfer") then
            if (vjp_workload) then
                !$acc data copyin(x, f, u) copyout(values, products)
                !$omp target data map(to: x, f, u) &
                !$omp& map(from: values, products)
                call execute_candidate()
                !$omp end target data
                !$acc end data
            else
                !$acc data copyin(x, f, v) copyout(values, products)
                !$omp target data map(to: x, f, v) &
                !$omp& map(from: values, products)
                call execute_candidate()
                !$omp end target data
                !$acc end data
            end if
        else
            call execute_candidate()
        end if
        call system_clock(tick1)
        if (trim(residency) == "transfer") then
            if (values(1) + products(n) /= values(1) + products(n)) then
                error stop "benchmark produced NaN"
            end if
        end if
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

    subroutine execute_candidate()
        select case (trim(candidate))
        case ("fused_jvp")
            call dawson_value_jvp_batch(n, x, f, v, values, products)
        case ("separate_jvp")
            call dawson_value_batch(n, f, values)
            call dawson_jvp_batch(n, x, f, v, products)
        case ("fused_vjp")
            call dawson_value_vjp_batch(n, x, f, u, values, products)
        case ("separate_vjp")
            call dawson_value_batch(n, f, values)
            call dawson_vjp_batch(n, x, f, u, products)
        end select
    end subroutine execute_candidate

end program bench_dawson_gpu_fusion
