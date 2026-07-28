module test_gpu_ode_objective_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_ode_objective_wrapper, only: &
        ode_trace2_terminal_objective_batch
    implicit none
    private

    public :: run_ode_objective_test

contains

    subroutine run_ode_objective_test()
        integer, parameter :: batch_size = 4096, step_count = 64
        real(dp), parameter :: tolerance = 4.0e-13_dp
        real(dp) :: transitions(step_count, 4)
        real(dp) :: initial_states(batch_size, 2), targets(batch_size, 2)
        real(dp) :: terminal_states(batch_size, 2)
        real(dp) :: cotangents(batch_size, 2), gradients(batch_size, 2)
        real(dp) :: losses(batch_size)
        real(dp) :: total_angle, total_scale, angle, scale
        real(dp) :: terminal1, terminal2, delta1, delta2
        real(dp) :: expected_loss, gradient1, gradient2, magnitude
        integer :: batch_index, step

        total_angle = 0.0_dp
        total_scale = 1.0_dp
        do step = 1, step_count
            angle = 0.0005_dp*real(1 + mod(7*step, 17), dp)
            scale = 1.0_dp + 0.0001_dp*real(mod(5*step, 11) - 5, dp)
            transitions(step, :) = [ &
                scale*cos(angle), scale*sin(angle), &
                -scale*sin(angle), scale*cos(angle)]
            total_angle = total_angle + angle
            total_scale = total_scale*scale
        end do
        do batch_index = 1, batch_size
            initial_states(batch_index, :) = [ &
                sin(real(batch_index, dp)/19.0_dp), &
                cos(real(batch_index, dp)/23.0_dp)]
            targets(batch_index, :) = [ &
                0.2_dp*cos(real(batch_index, dp)/31.0_dp), &
                -0.3_dp*sin(real(batch_index, dp)/37.0_dp)]
        end do

        !$acc data copyin(transitions, initial_states, targets) &
        !$acc& create(terminal_states, cotangents, losses, gradients)
        !$omp target data map(to: transitions, initial_states, targets) &
        !$omp& map(alloc: terminal_states, cotangents, losses, gradients)
        call ode_trace2_terminal_objective_batch( &
            batch_size, step_count, transitions, initial_states, targets, &
            terminal_states, cotangents, losses, gradients)
        !$omp target update from(terminal_states, losses, gradients)
        !$omp end target data
        !$acc update self(terminal_states, losses, gradients)
        !$acc end data

        do batch_index = 1, batch_size
            call rotated( &
                initial_states(batch_index, 1), initial_states(batch_index, 2), &
                total_angle, total_scale, terminal1, terminal2)
            delta1 = terminal1 - targets(batch_index, 1)
            delta2 = terminal2 - targets(batch_index, 2)
            expected_loss = 0.5_dp*(delta1*delta1 + delta2*delta2)
            call rotated( &
                delta1, delta2, -total_angle, total_scale, gradient1, gradient2)
            magnitude = max(1.0_dp, abs(expected_loss), abs(gradient1), abs(gradient2))
            if (abs(losses(batch_index) - expected_loss) > tolerance*magnitude) then
                error stop "GPU ODE objective disagrees with closed-form loss"
            end if
            if (max( &
                abs(gradients(batch_index, 1) - gradient1), &
                abs(gradients(batch_index, 2) - gradient2)) > &
                tolerance*magnitude) then
                error stop "GPU ODE objective disagrees with closed-form gradient"
            end if
        end do
    end subroutine run_ode_objective_test

    pure subroutine rotated(x, y, angle, scale, output1, output2)
        real(dp), intent(in) :: x, y, angle, scale
        real(dp), intent(out) :: output1, output2

        output1 = scale*(cos(angle)*x - sin(angle)*y)
        output2 = scale*(sin(angle)*x + cos(angle)*y)
    end subroutine rotated

end module test_gpu_ode_objective_support
