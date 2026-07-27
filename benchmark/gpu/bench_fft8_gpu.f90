program bench_fft8_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_fft8_wrapper, only: fft8_jvp_batch, fft8_vjp_batch
    implicit none

    integer, parameter :: sample_count = 31
    complex(dp), allocatable :: inputs(:, :), outputs(:, :)
    character(16) :: product, argument, residency, option
    integer :: batch_index, entry, batch_size, sample_index
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

    allocate (inputs(batch_size, 8), outputs(batch_size, 8))
    do entry = 1, 8
        do batch_index = 1, batch_size
            inputs(batch_index, entry) = cmplx( &
                sin(real(3*batch_index + entry, dp)/19.0_dp), &
                cos(real(batch_index + 5*entry, dp)/23.0_dp), dp)
        end do
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    sink = real(outputs(batch_size, 8), dp)
    if (sink /= sink) error stop "FFT8 benchmark produced NaN"

contains

    subroutine run_resident()
        !$acc data copyin(inputs) create(outputs)
        !$omp target data map(to: inputs) map(alloc: outputs)
        call run_benchmark()
        !$omp target update from(outputs)
        !$omp end target data
        !$acc update self(outputs)
        !$acc end data
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
        !$acc data copyin(inputs) copyout(outputs)
        !$omp target data map(to: inputs) map(from: outputs)
        call execute_product()
        !$omp end target data
        !$acc end data
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(product) == "jvp") then
            call fft8_jvp_batch(batch_size, inputs, outputs)
        else
            call fft8_vjp_batch(batch_size, inputs, outputs)
        end if
    end subroutine execute_product

end program bench_fft8_gpu
