program bench_ode_objective_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_ode_objective_wrapper, only: &
        ode_trace2_terminal_objective_batch
    implicit none

    integer, parameter :: sample_count = 31
    real(dp), allocatable :: transitions(:, :)
    real(dp), allocatable :: initial_states(:, :), targets(:, :)
    real(dp), allocatable :: terminal_states(:, :), cotangents(:, :)
    real(dp), allocatable :: losses(:), gradients(:, :)
    character(16) :: argument, residency, option
    integer :: batch_index, batch_size, step_count, step, sample_index
    integer :: repetition, repetitions
    logical :: peak_only
    real(dp) :: angle, scale, sink

    call get_command_argument(1, argument)
    read (argument, *) batch_size
    call get_command_argument(2, argument)
    read (argument, *) step_count
    call get_command_argument(3, residency)
    call get_command_argument(4, option)
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
    allocate (initial_states(batch_size, 2), targets(batch_size, 2))
    allocate (terminal_states(batch_size, 2), cotangents(batch_size, 2))
    allocate (losses(batch_size), gradients(batch_size, 2))
    do step = 1, step_count
        angle = 0.0005_dp*real(1 + mod(7*step, 17), dp)
        scale = 1.0_dp + 0.0001_dp*real(mod(5*step, 11) - 5, dp)
        transitions(step, :) = [ &
            scale*cos(angle), scale*sin(angle), &
            -scale*sin(angle), scale*cos(angle)]
    end do
    do batch_index = 1, batch_size
        initial_states(batch_index, :) = [ &
            sin(real(batch_index, dp)/19.0_dp), &
            cos(real(batch_index, dp)/23.0_dp)]
        targets(batch_index, :) = [ &
            0.2_dp*cos(real(batch_index, dp)/31.0_dp), &
            -0.3_dp*sin(real(batch_index, dp)/37.0_dp)]
    end do

    if (trim(residency) == "resident") then
        call run_resident()
    else
        call run_benchmark()
    end if
    sink = losses(batch_size) + gradients(batch_size, 2)
    if (sink /= sink) error stop "ODE-objective benchmark produced NaN"

contains

    subroutine run_resident()
        !$acc data copyin(transitions, initial_states, targets) &
        !$acc& create(terminal_states, cotangents, losses, gradients)
        !$omp target data map(to: transitions, initial_states, targets) &
        !$omp& map(alloc: terminal_states, cotangents, losses, gradients)
        call run_benchmark()
        !$omp target update from(terminal_states, losses, gradients)
        !$omp end target data
        !$acc update self(terminal_states, losses, gradients)
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
                call execute_application()
            end if
        end do
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/ &
            (real(rate, dp)*real(repetitions, dp))
    end function run_sample

    subroutine execute_with_transfer()
        !$acc data copyin(transitions, initial_states, targets) &
        !$acc& create(cotangents) copyout(terminal_states, losses, gradients)
        !$omp target data map(to: transitions, initial_states, targets) &
        !$omp& map(alloc: cotangents) map(from: terminal_states, losses, gradients)
        call execute_application()
        !$omp end target data
        !$acc end data
    end subroutine execute_with_transfer

    subroutine execute_application()
        call ode_trace2_terminal_objective_batch( &
            batch_size, step_count, transitions, initial_states, targets, &
            terminal_states, cotangents, losses, gradients)
    end subroutine execute_application

end program bench_ode_objective_gpu
