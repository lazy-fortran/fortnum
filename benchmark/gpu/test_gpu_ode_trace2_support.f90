module test_gpu_ode_trace2_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_ode_trace2_wrapper, only: &
        ode_trace2_jvp_batch, ode_trace2_vjp_batch
    implicit none
    private

    public :: run_ode_trace2_test

contains

    subroutine run_ode_trace2_test()
        integer, parameter :: batch_size = 4096, step_count = 64
        real(dp), parameter :: tolerance = 3.0e-13_dp
        real(dp) :: transitions(step_count, 4)
        real(dp) :: directions(batch_size, 2), cotangents(batch_size, 2)
        real(dp) :: jvps(batch_size, 2), vjps(batch_size, 2)
        real(dp) :: total_angle, total_scale, angle, scale
        real(dp) :: expected1, expected2, lhs, rhs, magnitude
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
            directions(batch_index, :) = [ &
                sin(real(batch_index, dp)/19.0_dp), &
                cos(real(batch_index, dp)/23.0_dp)]
            cotangents(batch_index, :) = [ &
                cos(real(batch_index, dp)/17.0_dp), &
                sin(real(batch_index, dp)/29.0_dp)]
        end do

        !$acc data copyin(transitions, directions, cotangents) create(jvps, vjps)
        !$omp target data map(to: transitions, directions, cotangents) &
        !$omp& map(alloc: jvps, vjps)
        call ode_trace2_jvp_batch( &
            batch_size, step_count, transitions, directions, jvps)
        call ode_trace2_vjp_batch( &
            batch_size, step_count, transitions, cotangents, vjps)
        !$omp target update from(jvps, vjps)
        !$omp end target data
        !$acc update self(jvps, vjps)
        !$acc end data

        do batch_index = 1, batch_size
            call rotated( &
                directions(batch_index, 1), directions(batch_index, 2), &
                total_angle, total_scale, expected1, expected2)
            magnitude = max(1.0_dp, abs(expected1), abs(expected2))
            if (max( &
                abs(jvps(batch_index, 1) - expected1), &
                abs(jvps(batch_index, 2) - expected2)) > &
                tolerance*magnitude) then
                error stop "GPU ODE-trace JVP disagrees with closed form"
            end if
            call rotated( &
                cotangents(batch_index, 1), cotangents(batch_index, 2), &
                -total_angle, total_scale, expected1, expected2)
            magnitude = max(1.0_dp, abs(expected1), abs(expected2))
            if (max( &
                abs(vjps(batch_index, 1) - expected1), &
                abs(vjps(batch_index, 2) - expected2)) > &
                tolerance*magnitude) then
                error stop "GPU ODE-trace VJP disagrees with closed form"
            end if
        end do

        lhs = sum(cotangents*jvps)
        rhs = sum(vjps*directions)
        magnitude = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*magnitude) then
            error stop "GPU ODE-trace products violate adjoint identity"
        end if
    end subroutine run_ode_trace2_test

    pure subroutine rotated(x, y, angle, scale, output1, output2)
        real(dp), intent(in) :: x, y, angle, scale
        real(dp), intent(out) :: output1, output2

        output1 = scale*(cos(angle)*x - sin(angle)*y)
        output2 = scale*(sin(angle)*x + cos(angle)*y)
    end subroutine rotated

end module test_gpu_ode_trace2_support
