program bench_multi_input_scaling
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_multi_input_wrappers, only: &
        multi_input_p2_jvp_batch, multi_input_p2_vjp_batch, &
        multi_input_p4_jvp_batch, multi_input_p4_vjp_batch, &
        multi_input_p8_jvp_batch, multi_input_p8_vjp_batch, &
        multi_input_p16_jvp_batch, multi_input_p16_vjp_batch
    implicit none

    integer, parameter :: samples = 31
    real(dp), allocatable :: x(:, :), directions(:, :, :)
    real(dp), allocatable :: cotangents(:, :), values(:, :), products(:, :)
    real(dp), allocatable :: adjoints(:, :, :)
    character(16) :: product, residency, argument, option
    integer :: i, j, k, n, nactive, nproducts, sample
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, product)
    call get_command_argument(2, argument)
    read (argument, *) nactive
    call get_command_argument(3, argument)
    read (argument, *) n
    call get_command_argument(4, argument)
    read (argument, *) nproducts
    call get_command_argument(5, residency)
    call get_command_argument(6, option)
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if
    select case (nactive)
    case (2, 4, 8, 16)
    case default
        error stop "active inputs must be 2, 4, 8, or 16"
    end select
    if (n < 1) error stop "batch size must be positive"
    if (nproducts < 1) error stop "product count must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"

    allocate (x(n, nactive), values(n, nproducts))
    if (trim(product) == "jvp") then
        allocate (directions(n, nactive, nproducts), products(n, nproducts))
        allocate (cotangents(0, 0), adjoints(0, 0, 0))
    else
        allocate (directions(0, 0, 0), products(0, 0))
        allocate (cotangents(n, nproducts), adjoints(n, nactive, nproducts))
    end if
    do j = 1, nactive
        do i = 1, n
            x(i, j) = -0.4_dp + 0.8_dp* &
                real(mod(17*i + 11*j, 257), dp)/256.0_dp
            if (trim(product) == "jvp") then
                do k = 1, nproducts
                    directions(i, j, k) = -0.7_dp + 1.4_dp* &
                        real(mod(13*i + 19*j + 23*k, 251), dp)/250.0_dp
                end do
            end if
        end do
    end do
    if (trim(product) == "vjp") then
        do k = 1, nproducts
            do i = 1, n
                cotangents(i, k) = -0.8_dp + 1.6_dp* &
                    real(mod(23*i + 29*k, 263), dp)/262.0_dp
            end do
        end do
    end if

    if (trim(residency) == "resident") then
        if (trim(product) == "jvp") then
            !$acc data copyin(x, directions) create(values, products)
            !$omp target data map(to: x, directions) &
            !$omp& map(alloc: values, products)
            call run_benchmark()
            !$omp target update from(values, products)
            !$omp end target data
            !$acc update self(values, products)
            !$acc end data
        else
            !$acc data copyin(x, cotangents) create(values, adjoints)
            !$omp target data map(to: x, cotangents) &
            !$omp& map(alloc: values, adjoints)
            call run_benchmark()
            !$omp target update from(values, adjoints)
            !$omp end target data
            !$acc update self(values, adjoints)
            !$acc end data
        end if
    else
        call run_benchmark()
    end if
    if (trim(product) == "jvp") then
        sink = values(1, 1) + products(n, nproducts)
    else
        sink = values(1, 1) + adjoints(n, nactive, nproducts)
    end if
    if (sink /= sink) error stop "benchmark produced NaN"

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
            if (trim(product) == "jvp") then
                !$acc data copyin(x, directions) copyout(values, products)
                !$omp target data map(to: x, directions) &
                !$omp& map(from: values, products)
                call execute_products()
                !$omp end target data
                !$acc end data
            else
                !$acc data copyin(x, cotangents) copyout(values, adjoints)
                !$omp target data map(to: x, cotangents) &
                !$omp& map(from: values, adjoints)
                call execute_products()
                !$omp end target data
                !$acc end data
            end if
        else
            call execute_products()
        end if
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

    subroutine execute_products()
        do k = 1, nproducts
            select case (nactive)
            case (2)
                if (trim(product) == "jvp") then
                    call multi_input_p2_jvp_batch( &
                        n, x, directions(:, :, k), values(:, k), products(:, k))
                else
                    call multi_input_p2_vjp_batch( &
                        n, x, cotangents(:, k), values(:, k), adjoints(:, :, k))
                end if
            case (4)
                if (trim(product) == "jvp") then
                    call multi_input_p4_jvp_batch( &
                        n, x, directions(:, :, k), values(:, k), products(:, k))
                else
                    call multi_input_p4_vjp_batch( &
                        n, x, cotangents(:, k), values(:, k), adjoints(:, :, k))
                end if
            case (8)
                if (trim(product) == "jvp") then
                    call multi_input_p8_jvp_batch( &
                        n, x, directions(:, :, k), values(:, k), products(:, k))
                else
                    call multi_input_p8_vjp_batch( &
                        n, x, cotangents(:, k), values(:, k), adjoints(:, :, k))
                end if
            case (16)
                if (trim(product) == "jvp") then
                    call multi_input_p16_jvp_batch( &
                        n, x, directions(:, :, k), values(:, k), products(:, k))
                else
                    call multi_input_p16_vjp_batch( &
                        n, x, cotangents(:, k), values(:, k), adjoints(:, :, k))
                end if
            end select
        end do
    end subroutine execute_products

end program bench_multi_input_scaling
