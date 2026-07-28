module fortnum_gpu_ode_objective_wrapper
    use fortnum_kinds, only: dp
    use fortnum_gpu_ode_trace2_wrapper, only: &
        ode_trace2_jvp_batch, ode_trace2_vjp_batch
    implicit none
    private

    public :: ode_trace2_terminal_objective_batch

contains

    subroutine ode_trace2_terminal_objective_batch( &
            batch_size, step_count, transitions, initial_states, targets, &
            terminal_states, cotangents, losses, gradients)
        integer, intent(in) :: batch_size, step_count
        real(dp), intent(in) :: transitions(step_count, 4)
        real(dp), intent(in) :: initial_states(batch_size, 2)
        real(dp), intent(in) :: targets(batch_size, 2)
        real(dp), intent(out) :: terminal_states(batch_size, 2)
        real(dp), intent(out) :: cotangents(batch_size, 2)
        real(dp), intent(out) :: losses(batch_size)
        real(dp), intent(out) :: gradients(batch_size, 2)

        call ode_trace2_jvp_batch( &
            batch_size, step_count, transitions, initial_states, &
            terminal_states)
        call terminal_objective_seed_batch( &
            batch_size, terminal_states, targets, losses, cotangents)
        call ode_trace2_vjp_batch( &
            batch_size, step_count, transitions, cotangents, gradients)
    end subroutine ode_trace2_terminal_objective_batch

    subroutine terminal_objective_seed_batch( &
            batch_size, terminal_states, targets, losses, cotangents)
        integer, intent(in) :: batch_size
        real(dp), intent(in) :: terminal_states(batch_size, 2)
        real(dp), intent(in) :: targets(batch_size, 2)
        real(dp), intent(out) :: losses(batch_size)
        real(dp), intent(out) :: cotangents(batch_size, 2)
        real(dp) :: delta1, delta2
        integer :: batch_index

        !$acc parallel loop present(terminal_states, targets, losses, cotangents) &
        !$acc& private(delta1, delta2)
        !$omp target teams distribute parallel do &
        !$omp& map(to: terminal_states, targets) map(from: losses, cotangents) &
        !$omp& private(delta1, delta2)
        do batch_index = 1, batch_size
            delta1 = terminal_states(batch_index, 1) - targets(batch_index, 1)
            delta2 = terminal_states(batch_index, 2) - targets(batch_index, 2)
            cotangents(batch_index, 1) = delta1
            cotangents(batch_index, 2) = delta2
            losses(batch_index) = 0.5_dp*(delta1*delta1 + delta2*delta2)
        end do
    end subroutine terminal_objective_seed_batch

end module fortnum_gpu_ode_objective_wrapper
