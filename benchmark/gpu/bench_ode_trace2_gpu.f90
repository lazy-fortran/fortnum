program bench_ode_trace2_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_ode_trace2_wrapper, only: &
        ode_trace2_jvp_batch, ode_trace2_vjp_batch
    implicit none

    integer, parameter :: sample_count = 31
    real(dp), allocatable :: transitions(:, :), inputs(:, :), outputs(:, :)
    character(16) :: product, argument, residency, option
    integer :: batch_index, batch_size, step_count, step, sample_index
    integer :: repetition, repetitions
    logical :: peak_only
    real(dp) :: angle, scale, sink

    call get_command_argument(1, product)
    call get_command_argument(2, argument)
    read (argument, *) batch_size
    call get_command_argument(3, argument)
    read (argument, *) step_count
    call get_command_argument(4, residency)
    call get_command_argument(5, option)
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if
    if (batch_size < 1 .or. step_count < 1) then
        error stop "batch size and step count must be positive"
    end if
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"
    repetitions = max(1, (65536 + batch_size - 1)/batch_size)

    allocate (transitions(step_count, 4))
    allocate (inputs(batch_size, 2), outputs(batch_size, 2))
    do step = 1, step_count
        angle = 0.0005_dp*real(1 + mod(7*step, 17), dp)
        scale = 1.0_dp + 0.0001_dp*real(mod(5*step, 11) - 5, dp)
        transitions(step, :) = [ &
            scale*cos(angle), scale*sin(angle), &
            -scale*sin(angle), scale*cos(angle)]
    end do
    do batch_index = 1, batch_size
        inputs(batch_index, :) = [ &
            sin(real(batch_index, dp)/19.0_dp), &
            cos(real(batch_index, dp)/23.0_dp)]
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    sink = outputs(batch_size, 2)
    if (sink /= sink) error stop "ODE-trace benchmark produced NaN"

contains

    subroutine run_resident()
        !$acc data copyin(transitions, inputs) create(outputs)
        !$omp target data map(to: transitions, inputs) map(alloc: outputs)
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
        !$acc data copyin(transitions, inputs) copyout(outputs)
        !$omp target data map(to: transitions, inputs) map(from: outputs)
        call execute_product()
        !$omp end target data
        !$acc end data
    end subroutine execute_with_transfer

    subroutine execute_product()
        if (trim(product) == "jvp") then
            call ode_trace2_jvp_batch( &
                batch_size, step_count, transitions, inputs, outputs)
        else
            call ode_trace2_vjp_batch( &
                batch_size, step_count, transitions, inputs, outputs)
        end if
    end subroutine execute_product

end program bench_ode_trace2_gpu
