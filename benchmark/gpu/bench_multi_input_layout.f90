program bench_multi_input_layout
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_gpu_multi_input_wrappers, only: &
        multi_input_p8_jvp_batch, multi_input_p8_vjp_batch, &
        multi_input_p8_jvp_active_first_batch, &
        multi_input_p8_vjp_active_first_batch
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: p = 8, samples = 31
    real(dp), allocatable :: x_batch(:, :), direction_batch(:, :)
    real(dp), allocatable :: x_active(:, :), direction_active(:, :)
    real(dp), allocatable :: cotangent(:), values(:), products(:)
    real(dp), allocatable :: adjoints_batch(:, :), adjoints_active(:, :)
    character(16) :: product, residency, layout, argument, option
    integer :: i, j, n, sample
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, product)
    call get_command_argument(2, argument)
    read (argument, *) n
    call get_command_argument(3, residency)
    call get_command_argument(4, layout)
    call get_command_argument(5, option)
    peak_only = trim(option) == "--peak-rss"
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if
    if (n < 1) error stop "batch size must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    if (trim(layout) /= "batch_first" .and. &
        trim(layout) /= "active_first") then
        error stop "layout must be batch_first or active_first"
    end if

    allocate (values(n))
    if (trim(product) == "jvp") then
        allocate (cotangent(0))
    else
        allocate (cotangent(n))
    end if
    if (trim(layout) == "batch_first") then
        allocate (x_batch(n, p), x_active(0, 0))
        if (trim(product) == "jvp") then
            allocate (direction_batch(n, p), products(n))
            allocate (direction_active(0, 0), adjoints_batch(0, 0))
        else
            allocate (direction_batch(0, 0), products(0))
            allocate (direction_active(0, 0), adjoints_batch(n, p))
        end if
        allocate (adjoints_active(0, 0))
    else
        allocate (x_batch(0, 0), x_active(p, n))
        if (trim(product) == "jvp") then
            allocate (direction_active(p, n), products(n))
            allocate (direction_batch(0, 0), adjoints_active(0, 0))
        else
            allocate (direction_active(0, 0), products(0))
            allocate (direction_batch(0, 0), adjoints_active(p, n))
        end if
        allocate (adjoints_batch(0, 0))
    end if

    do i = 1, n
        if (trim(product) == "vjp") then
            cotangent(i) = -0.8_dp + 1.6_dp* &
                real(mod(23*i, 263), dp)/262.0_dp
        end if
        do j = 1, p
            if (trim(layout) == "batch_first") then
                x_batch(i, j) = input_value(i, j)
                if (trim(product) == "jvp") then
                    direction_batch(i, j) = direction_value(i, j)
                end if
            else
                x_active(j, i) = input_value(i, j)
                if (trim(product) == "jvp") then
                    direction_active(j, i) = direction_value(i, j)
                end if
            end if
        end do
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    if (trim(product) == "jvp") then
        sink = values(1) + products(n)
    else if (trim(layout) == "batch_first") then
        sink = values(1) + adjoints_batch(n, p)
    else
        sink = values(1) + adjoints_active(p, n)
    end if
    if (sink /= sink) error stop "benchmark produced NaN"

contains

    pure function input_value(batch_index, active_index) result(value)
        integer, intent(in) :: batch_index, active_index
        real(dp) :: value

        value = -0.4_dp + 0.8_dp* &
            real(mod(17*batch_index + 11*active_index, 257), dp)/256.0_dp
    end function input_value

    pure function direction_value(batch_index, active_index) result(value)
        integer, intent(in) :: batch_index, active_index
        real(dp) :: value

        value = -0.7_dp + 1.4_dp* &
            real(mod(13*batch_index + 19*active_index, 251), dp)/250.0_dp
    end function direction_value

    subroutine run_resident()
        if (trim(layout) == "batch_first") then
            if (trim(product) == "jvp") then
                !$acc data copyin(x_batch, direction_batch) &
                !$acc& create(values, products)
                !$omp target data map(to: x_batch, direction_batch) &
                !$omp& map(alloc: values, products)
                call run_benchmark()
                !$omp target update from(values, products)
                !$omp end target data
                !$acc update self(values, products)
                !$acc end data
            else
                !$acc data copyin(x_batch, cotangent) &
                !$acc& create(values, adjoints_batch)
                !$omp target data map(to: x_batch, cotangent) &
                !$omp& map(alloc: values, adjoints_batch)
                call run_benchmark()
                !$omp target update from(values, adjoints_batch)
                !$omp end target data
                !$acc update self(values, adjoints_batch)
                !$acc end data
            end if
        else
            if (trim(product) == "jvp") then
                !$acc data copyin(x_active, direction_active) &
                !$acc& create(values, products)
                !$omp target data map(to: x_active, direction_active) &
                !$omp& map(alloc: values, products)
                call run_benchmark()
                !$omp target update from(values, products)
                !$omp end target data
                !$acc update self(values, products)
                !$acc end data
            else
                !$acc data copyin(x_active, cotangent) &
                !$acc& create(values, adjoints_active)
                !$omp target data map(to: x_active, cotangent) &
                !$omp& map(alloc: values, adjoints_active)
                call run_benchmark()
                !$omp target update from(values, adjoints_active)
                !$omp end target data
                !$acc update self(values, adjoints_active)
                !$acc end data
            end if
        end if
    end subroutine run_resident

    subroutine run_benchmark()
        if (peak_only) then
            sink = run_sample()
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
            call execute_with_transfer()
        else
            call execute_product()
        end if
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

    subroutine execute_with_transfer()
        if (trim(layout) == "batch_first") then
            if (trim(product) == "jvp") then
                !$acc data copyin(x_batch, direction_batch) &
                !$acc& copyout(values, products)
                !$omp target data map(to: x_batch, direction_batch) &
                !$omp& map(from: values, products)
                call execute_product()
                !$omp end target data
                !$acc end data
            else
                !$acc data copyin(x_batch, cotangent) &
                !$acc& copyout(values, adjoints_batch)
                !$omp target data map(to: x_batch, cotangent) &
                !$omp& map(from: values, adjoints_batch)
                call execute_product()
                !$omp end target data
                !$acc end data
            end if
        else
            if (trim(product) == "jvp") then
                !$acc data copyin(x_active, direction_active) &
                !$acc& copyout(values, products)
                !$omp target data map(to: x_active, direction_active) &
                !$omp& map(from: values, products)
                call execute_product()
                !$omp end target data
                !$acc end data
            else
                !$acc data copyin(x_active, cotangent) &
                !$acc& copyout(values, adjoints_active)
                !$omp target data map(to: x_active, cotangent) &
                !$omp& map(from: values, adjoints_active)
                call execute_product()
                !$omp end target data
                !$acc end data
            end if
        end if
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(layout) == "batch_first") then
            if (trim(product) == "jvp") then
                call multi_input_p8_jvp_batch( &
                    n, x_batch, direction_batch, values, products)
            else
                call multi_input_p8_vjp_batch( &
                    n, x_batch, cotangent, values, adjoints_batch)
            end if
        else
            if (trim(product) == "jvp") then
                call multi_input_p8_jvp_active_first_batch( &
                    n, x_active, direction_active, values, products)
            else
                call multi_input_p8_vjp_active_first_batch( &
                    n, x_active, cotangent, values, adjoints_active)
            end if
        end if
    end subroutine execute_product

end program bench_multi_input_layout
