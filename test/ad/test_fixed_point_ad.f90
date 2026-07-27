program test_fixed_point_ad
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_fixed_point, only: fixed_point_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    real(dp) :: p(2), pp(2), pm(2), tp(2)
    real(dp) :: xstar(2), xp(2), xm(2), dx(2), dx_fd(2)
    real(dp) :: map_x(2, 2), map_p(2, 2), h
    type(fortnum_status_t) :: status

    p = [0.1_dp, -0.2_dp]
    tp = [0.4_dp, -0.6_dp]
    h = 1.0e-5_dp
    call solve_fixed_point(p, xstar)
    call map_derivatives(xstar, p, map_x, map_p)
    call fixed_point_jvp(map_x, map_p, tp, dx, status)

    pp = p + h*tp
    pm = p - h*tp
    call solve_fixed_point(pp, xp)
    call solve_fixed_point(pm, xm)
    dx_fd = (xp - xm)/(2.0_dp*h)

    if (.not. status_ok(status)) then
        write (error_unit, '(a)') "FAIL fixed-point tangent status"
        stop 1
    end if
    if (maxval(abs(dx - dx_fd)) > 1.0e-9_dp) then
        write (error_unit, '(a,2es24.16)') "FAIL fixed-point implicit=", dx
        write (error_unit, '(a,2es24.16)') "FAIL fixed-point finite difference=", dx_fd
        stop 1
    end if

    write (*, '(a)') "PASS"

contains

    subroutine solve_fixed_point(parameters, x)
        real(dp), intent(in) :: parameters(2)
        real(dp), intent(out) :: x(2)
        real(dp) :: next(2)
        integer :: iteration

        x = 0.0_dp
        do iteration = 1, 1000
            call fixed_point_map(x, parameters, next)
            if (maxval(abs(next - x)) <= 1.0e-14_dp) then
                x = next
                return
            end if
            x = next
        end do
        error stop "fixed-point oracle did not converge"
    end subroutine solve_fixed_point

    pure subroutine fixed_point_map(x, parameters, next)
        real(dp), intent(in) :: x(2), parameters(2)
        real(dp), intent(out) :: next(2)

        next(1) = tanh(0.2_dp*x(1) + 0.1_dp*x(2) + parameters(1))
        next(2) = tanh(0.05_dp*x(1) + 0.25_dp*x(2) + parameters(2))
    end subroutine fixed_point_map

    pure subroutine map_derivatives(x, parameters, map_x, map_p)
        real(dp), intent(in) :: x(2), parameters(2)
        real(dp), intent(out) :: map_x(2, 2), map_p(2, 2)
        real(dp) :: value(2), scale(2)

        call fixed_point_map(x, parameters, value)
        scale = 1.0_dp - value*value
        map_x(1, :) = scale(1)*[0.2_dp, 0.1_dp]
        map_x(2, :) = scale(2)*[0.05_dp, 0.25_dp]
        map_p = 0.0_dp
        map_p(1, 1) = scale(1)
        map_p(2, 2) = scale(2)
    end subroutine map_derivatives

end program test_fixed_point_ad
