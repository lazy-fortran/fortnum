program bench_linalg3_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_linalg3_wrapper, only: &
        det3_jvp_batch, det3_vjp_batch, inv3_jvp_batch, inv3_vjp_batch
    implicit none

    integer, parameter :: sample_count = 31
    real(dp), allocatable :: primary(:, :), secondary(:, :), outputs(:, :)
    real(dp), allocatable :: scalar_inputs(:), scalar_outputs(:)
    character(16) :: operation, product, argument, residency, option
    integer :: batch_index, entry, batch_size, sample_index
    integer :: repetition, repetitions
    logical :: peak_only
    real(dp) :: shift, sink

    call get_command_argument(1, operation)
    call get_command_argument(2, product)
    call get_command_argument(3, argument)
    read (argument, *) batch_size
    call get_command_argument(4, residency)
    call get_command_argument(5, option)
    if (trim(operation) /= "det" .and. trim(operation) /= "inv") then
        error stop "operation must be det or inv"
    end if
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

    allocate (primary(batch_size, 9))
    if (trim(operation) == "det" .and. trim(product) == "jvp") then
        allocate (secondary(batch_size, 9), scalar_outputs(batch_size))
        allocate (outputs(0, 0), scalar_inputs(0))
    else if (trim(operation) == "det") then
        allocate (scalar_inputs(batch_size), outputs(batch_size, 9))
        allocate (secondary(0, 0), scalar_outputs(0))
    else
        allocate (secondary(batch_size, 9), outputs(batch_size, 9))
        allocate (scalar_inputs(0), scalar_outputs(0))
    end if
    do batch_index = 1, batch_size
        shift = real(mod(batch_index, 29), dp)/113.0_dp
        primary(batch_index, :) = [ &
            1.8_dp + shift, 0.1_dp, -0.2_dp, &
            0.2_dp, 1.6_dp + shift, 0.15_dp, &
            -0.1_dp, 0.25_dp, 1.9_dp + shift]
        if (trim(operation) == "inv") then
            primary(batch_index, :) = [ &
                0.7_dp + 0.1_dp*shift, -0.04_dp, 0.08_dp, &
                -0.06_dp, 0.8_dp + 0.1_dp*shift, -0.05_dp, &
                0.03_dp, -0.09_dp, 0.65_dp + 0.1_dp*shift]
        end if
        if (trim(operation) == "det" .and. trim(product) == "vjp") then
            scalar_inputs(batch_index) = sin(real(batch_index, dp)/17.0_dp)
        else
            do entry = 1, 9
                secondary(batch_index, entry) = &
                    cos(real(7*batch_index + 3*entry, dp)/31.0_dp)
            end do
        end if
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    if (trim(operation) == "det" .and. trim(product) == "jvp") then
        sink = scalar_outputs(batch_size)
    else
        sink = outputs(batch_size, 9)
    end if
    if (sink /= sink) error stop "linalg3 benchmark produced NaN"

contains

    subroutine run_resident()
        if (trim(operation) == "det" .and. trim(product) == "jvp") then
            !$acc data copyin(primary, secondary) create(scalar_outputs)
            !$omp target data map(to: primary, secondary) map(alloc: scalar_outputs)
            call run_benchmark()
            !$omp target update from(scalar_outputs)
            !$omp end target data
            !$acc update self(scalar_outputs)
            !$acc end data
        else if (trim(operation) == "det") then
            !$acc data copyin(primary, scalar_inputs) create(outputs)
            !$omp target data map(to: primary, scalar_inputs) map(alloc: outputs)
            call run_benchmark()
            !$omp target update from(outputs)
            !$omp end target data
            !$acc update self(outputs)
            !$acc end data
        else
            !$acc data copyin(primary, secondary) create(outputs)
            !$omp target data map(to: primary, secondary) map(alloc: outputs)
            call run_benchmark()
            !$omp target update from(outputs)
            !$omp end target data
            !$acc update self(outputs)
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
        if (trim(operation) == "det" .and. trim(product) == "jvp") then
            !$acc data copyin(primary, secondary) copyout(scalar_outputs)
            !$omp target data map(to: primary, secondary) map(from: scalar_outputs)
            call execute_product()
            !$omp end target data
            !$acc end data
        else if (trim(operation) == "det") then
            !$acc data copyin(primary, scalar_inputs) copyout(outputs)
            !$omp target data map(to: primary, scalar_inputs) map(from: outputs)
            call execute_product()
            !$omp end target data
            !$acc end data
        else
            !$acc data copyin(primary, secondary) copyout(outputs)
            !$omp target data map(to: primary, secondary) map(from: outputs)
            call execute_product()
            !$omp end target data
            !$acc end data
        end if
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(operation) == "det" .and. trim(product) == "jvp") then
            call det3_jvp_batch( &
                batch_size, primary, secondary, scalar_outputs)
        else if (trim(operation) == "det") then
            call det3_vjp_batch( &
                batch_size, primary, scalar_inputs, outputs)
        else if (trim(product) == "jvp") then
            call inv3_jvp_batch(batch_size, primary, secondary, outputs)
        else
            call inv3_vjp_batch(batch_size, primary, secondary, outputs)
        end if
    end subroutine execute_product

end program bench_linalg3_gpu
