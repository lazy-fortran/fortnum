module fortnum_gpu_ode_trace2_wrapper
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: ode_trace2_jvp_batch, ode_trace2_vjp_batch

contains

    subroutine ode_trace2_jvp_batch( &
            batch_size, step_count, transitions, directions, products)
        integer, intent(in) :: batch_size, step_count
        real(dp), intent(in) :: transitions(step_count, 4)
        real(dp), intent(in) :: directions(batch_size, 2)
        real(dp), intent(out) :: products(batch_size, 2)
        integer :: batch_index

        !$acc parallel loop present(transitions, directions, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: transitions, directions) map(from: products)
        do batch_index = 1, batch_size
            call ode_trace2_product( &
                step_count, transitions, directions(batch_index, 1), &
                directions(batch_index, 2), .false., &
                products(batch_index, 1), products(batch_index, 2))
        end do
    end subroutine ode_trace2_jvp_batch

    subroutine ode_trace2_vjp_batch( &
            batch_size, step_count, transitions, cotangents, products)
        integer, intent(in) :: batch_size, step_count
        real(dp), intent(in) :: transitions(step_count, 4)
        real(dp), intent(in) :: cotangents(batch_size, 2)
        real(dp), intent(out) :: products(batch_size, 2)
        integer :: batch_index

        !$acc parallel loop present(transitions, cotangents, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: transitions, cotangents) map(from: products)
        do batch_index = 1, batch_size
            call ode_trace2_product( &
                step_count, transitions, cotangents(batch_index, 1), &
                cotangents(batch_index, 2), .true., &
                products(batch_index, 1), products(batch_index, 2))
        end do
    end subroutine ode_trace2_vjp_batch

    pure subroutine ode_trace2_product( &
            step_count, transitions, input1, input2, transpose_product, &
            output1, output2)
        !$omp declare target
        !$acc routine seq
        integer, intent(in) :: step_count
        real(dp), intent(in) :: transitions(step_count, 4)
        real(dp), intent(in) :: input1, input2
        logical, intent(in) :: transpose_product
        real(dp), intent(out) :: output1, output2
        real(dp) :: value1, value2, next1, next2
        integer :: step

        value1 = input1
        value2 = input2
        if (transpose_product) then
            do step = step_count, 1, -1
                next1 = transitions(step, 1)*value1 + &
                    transitions(step, 2)*value2
                next2 = transitions(step, 3)*value1 + &
                    transitions(step, 4)*value2
                value1 = next1
                value2 = next2
            end do
        else
            do step = 1, step_count
                next1 = transitions(step, 1)*value1 + &
                    transitions(step, 3)*value2
                next2 = transitions(step, 2)*value1 + &
                    transitions(step, 4)*value2
                value1 = next1
                value2 = next2
            end do
        end if
        output1 = value1
        output2 = value2
    end subroutine ode_trace2_product

end module fortnum_gpu_ode_trace2_wrapper
