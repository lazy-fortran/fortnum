module fortnum_gpu_rk54_wrapper
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode_rk54_device, only: rk54_controls4_t, rk54_state4_t, &
        rk54_initialize4, rk54_request4, rk54_supply4, &
        RK54_DORMAND_PRINCE, RK54_NEED_RHS
    implicit none
    private

    public :: rk54_dopri_linear_batch

contains

    subroutine rk54_dopri_linear_batch(n, y0, h, yout, result)
        integer, intent(in) :: n
        real(dp), intent(in) :: y0(n, 4), h
        real(dp), intent(out) :: yout(n, 4)
        integer, intent(out) :: result(n)
        real(dp), parameter :: lambda(4) = [-1.0_dp, 0.5_dp, 1.0_dp, -2.0_dp]
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4), initial(4)
        integer :: i, request

        !$acc parallel loop present(y0, yout, result) &
        !$acc& private(controls, state, t_eval, y_eval, derivative, initial, request)
        !$omp target teams distribute parallel do map(to: y0, h) &
        !$omp& map(from: yout, result) &
        !$omp& private(controls, state, t_eval, y_eval, derivative, initial, request)
        do i = 1, n
            controls%method = RK54_DORMAND_PRINCE
            controls%rtol = 1.0e-6_dp
            controls%atol = 1.0e-6_dp
            controls%hmin = 0.0_dp
            controls%hmax = h
            initial = [y0(i, 1), y0(i, 2), y0(i, 3), y0(i, 4)]
            call rk54_initialize4(state, 0.0_dp, initial, h)
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                derivative = lambda*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, y_eval, &
                    request)
            end do
            yout(i, 1) = state%y(1)
            yout(i, 2) = state%y(2)
            yout(i, 3) = state%y(3)
            yout(i, 4) = state%y(4)
            result(i) = request
        end do
    end subroutine rk54_dopri_linear_batch

end module fortnum_gpu_rk54_wrapper
