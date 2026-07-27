program bench_fixed_quadrature_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_gpu_fixed_quadrature_wrapper, only: &
        fixed_quadrature_jvp_batch, fixed_quadrature_vjp_batch
    implicit none

    integer, parameter :: rule_order = 16, sample_count = 31
    real(dp) :: nodes(rule_order), weights(rule_order)
    real(dp), allocatable :: tangents(:, :), cotangents(:), products(:)
    real(dp), allocatable :: adjoints(:, :)
    character(16) :: product, argument, residency, option
    integer :: batch_index, node_index, batch_size, sample_index
    integer :: repetition, repetitions
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, product)
    call get_command_argument(2, argument)
    read (argument, *) batch_size
    call get_command_argument(3, residency)
    call get_command_argument(4, option)
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if
    if (batch_size < 1) error stop "batch size must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"
    repetitions = max(1, (65536 + batch_size - 1)/batch_size)

    call gauss_legendre(rule_order, nodes, weights)
    if (trim(product) == "jvp") then
        allocate (tangents(batch_size, rule_order), products(batch_size))
        allocate (cotangents(0), adjoints(0, 0))
        do node_index = 1, rule_order
            do batch_index = 1, batch_size
                tangents(batch_index, node_index) = &
                    sin(nodes(node_index) + real(mod(batch_index, 31), dp))
            end do
        end do
    else
        allocate (tangents(0, 0), products(0))
        allocate (cotangents(batch_size), adjoints(batch_size, rule_order))
        do batch_index = 1, batch_size
            cotangents(batch_index) = &
                cos(real(mod(batch_index, 37), dp)/13.0_dp)
        end do
    end if

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    if (trim(product) == "jvp") then
        sink = products(batch_size)
    else
        sink = adjoints(batch_size, rule_order)
    end if
    if (sink /= sink) error stop "fixed-quadrature benchmark produced NaN"

contains

    subroutine run_resident()
        if (trim(product) == "jvp") then
            !$acc data copyin(weights, tangents) create(products)
            !$omp target data map(to: weights, tangents) map(alloc: products)
            call run_benchmark()
            !$omp target update from(products)
            !$omp end target data
            !$acc update self(products)
            !$acc end data
        else
            !$acc data copyin(weights, cotangents) create(adjoints)
            !$omp target data map(to: weights, cotangents) map(alloc: adjoints)
            call run_benchmark()
            !$omp target update from(adjoints)
            !$omp end target data
            !$acc update self(adjoints)
            !$acc end data
        end if
    end subroutine run_resident

    subroutine run_benchmark()
        sink = run_sample()
        if (peak_only) then
            write (*, "(i0)") peak_rss_bytes()
            return
        end if
        do sample_index = 1, 3
            sink = run_sample()
        end do
        do sample_index = 1, sample_count
            write (*, "(f0.6)") run_sample()
        end do
    end subroutine run_benchmark

    function run_sample() result(milliseconds)
        real(dp) :: milliseconds
        integer(int64) :: tick0, tick1, rate

        call system_clock(tick0, rate)
        do repetition = 1, repetitions
            if (trim(residency) == "transfer") then
                call execute_with_transfer()
            else
                call execute_product()
            end if
        end do
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/ &
            (real(rate, dp)*real(repetitions, dp))
    end function run_sample

    subroutine execute_with_transfer()
        if (trim(product) == "jvp") then
            !$acc data copyin(weights, tangents) copyout(products)
            !$omp target data map(to: weights, tangents) map(from: products)
            call execute_product()
            !$omp end target data
            !$acc end data
        else
            !$acc data copyin(weights, cotangents) copyout(adjoints)
            !$omp target data map(to: weights, cotangents) map(from: adjoints)
            call execute_product()
            !$omp end target data
            !$acc end data
        end if
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(product) == "jvp") then
            call fixed_quadrature_jvp_batch( &
                batch_size, rule_order, weights, tangents, products)
        else
            call fixed_quadrature_vjp_batch( &
                batch_size, rule_order, weights, cotangents, adjoints)
        end if
    end subroutine execute_product

end program bench_fixed_quadrature_gpu
