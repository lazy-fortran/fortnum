program bench_lagrange4_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_lagrange4_wrapper, only: &
        lagrange4_jvp_batch, lagrange4_vjp_batch
    implicit none

    integer, parameter :: samples_count = 31
    real(dp), parameter :: nodes(4) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    real(dp), allocatable :: x(:), samples(:, :), tx(:)
    real(dp), allocatable :: sample_tangents(:, :), cotangents(:)
    real(dp), allocatable :: values(:), products(:), adjoint_x(:)
    real(dp), allocatable :: adjoint_samples(:, :)
    character(16) :: product, argument, residency, option
    integer :: i, j, n, sample_index
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, product)
    call get_command_argument(2, argument)
    read (argument, *) n
    call get_command_argument(3, residency)
    call get_command_argument(4, option)
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if
    if (n < 1) error stop "batch size must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"

    allocate (x(n), samples(n, 4), values(n))
    if (trim(product) == "jvp") then
        allocate (tx(n), sample_tangents(n, 4), products(n))
        allocate (cotangents(0), adjoint_x(0), adjoint_samples(0, 0))
    else
        allocate (tx(0), sample_tangents(0, 0), products(0))
        allocate (cotangents(n), adjoint_x(n), adjoint_samples(n, 4))
    end if
    do i = 1, n
        x(i) = -0.8_dp + 2.6_dp*real(mod(17*i, 4093), dp)/4092.0_dp
        if (trim(product) == "jvp") then
            tx(i) = -0.7_dp + 1.4_dp*real(mod(19*i, 4091), dp)/4090.0_dp
        else
            cotangents(i) = -0.9_dp + &
                1.8_dp*real(mod(23*i, 4079), dp)/4078.0_dp
        end if
        do j = 1, 4
            samples(i, j) = primal_polynomial(nodes(j), i)
            if (trim(product) == "jvp") then
                sample_tangents(i, j) = tangent_polynomial(nodes(j), i)
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
    else
        sink = values(1) + adjoint_x(n) + adjoint_samples(n, 4)
    end if
    if (sink /= sink) error stop "Lagrange benchmark produced NaN"

contains

    pure function primal_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 17), dp)/31.0_dp
        value = (0.3_dp + shift) - 0.4_dp*z + 0.2_dp*z*z - 0.05_dp*z*z*z
    end function primal_polynomial

    pure function tangent_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 13), dp)/29.0_dp
        value = (-0.2_dp + shift) + 0.3_dp*z - 0.1_dp*z*z + 0.04_dp*z*z*z
    end function tangent_polynomial

    subroutine run_resident()
        if (trim(product) == "jvp") then
            !$acc data copyin(x, samples, tx, sample_tangents) &
            !$acc& create(values, products)
            !$omp target data map(to: x, samples, tx, sample_tangents) &
            !$omp& map(alloc: values, products)
            call run_benchmark()
            !$omp target update from(values, products)
            !$omp end target data
            !$acc update self(values, products)
            !$acc end data
        else
            !$acc data copyin(x, samples, cotangents) &
            !$acc& create(values, adjoint_x, adjoint_samples)
            !$omp target data map(to: x, samples, cotangents) &
            !$omp& map(alloc: values, adjoint_x, adjoint_samples)
            call run_benchmark()
            !$omp target update from(values, adjoint_x, adjoint_samples)
            !$omp end target data
            !$acc update self(values, adjoint_x, adjoint_samples)
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
        do sample_index = 1, samples_count
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
        if (trim(product) == "jvp") then
            !$acc data copyin(x, samples, tx, sample_tangents) &
            !$acc& copyout(values, products)
            !$omp target data map(to: x, samples, tx, sample_tangents) &
            !$omp& map(from: values, products)
            call execute_product()
            !$omp end target data
            !$acc end data
        else
            !$acc data copyin(x, samples, cotangents) &
            !$acc& copyout(values, adjoint_x, adjoint_samples)
            !$omp target data map(to: x, samples, cotangents) &
            !$omp& map(from: values, adjoint_x, adjoint_samples)
            call execute_product()
            !$omp end target data
            !$acc end data
        end if
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(product) == "jvp") then
            call lagrange4_jvp_batch( &
                n, x, samples, tx, sample_tangents, values, products)
        else
            call lagrange4_vjp_batch( &
                n, x, samples, cotangents, values, adjoint_x, adjoint_samples)
        end if
    end subroutine execute_product

end program bench_lagrange4_gpu
