module fortnum_periodize
    ! Periodic continuation of functions on a window [x_mid - L/2, x_mid + L/2),
    ! L = 2 (a + w): equal to the function on the core |x - x_mid| <= a and
    ! blended smoothly into the neighbouring periods across transition bands
    ! of width w, so that the continuation is L-periodic and C-infinity when
    ! the function is. With x0 the point of the window congruent to x mod L,
    !   g(x) = w1 (1 - w2) f(x0) + (1 - w1)(1 - w2) f(x0 - L) + w2 f(x0 + L),
    ! w1 = step((x0 - x_mid - a)/(2 w)), w2 = step((x0 - x_mid + a + 2 w)/(2 w)),
    ! with the smooth step step(t) = exp(-2 pi/(1 - t) exp(-sqrt(2)/t)) on
    ! (0, 1), 1 below and 0 above, whose derivatives of every order vanish at
    ! both ends. periodize_points returns the three points and weights; the
    ! caller evaluates its own functions there, so any number of functions is
    ! continued with one call. The derivative of g follows from the derivatives
    ! of f at the same points and the weight derivatives dws.
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: periodize_points, smooth_step

    real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp
    real(dp), parameter :: sqrt2 = 1.414213562373095048801688724210_dp

contains

    ! The smooth step and its derivative with respect to t.
    elemental subroutine smooth_step(t, s, ds)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: s, ds
        real(dp) :: e, a
        if (t <= 0) then
            s = 1
            ds = 0
        else if (t >= 1) then
            s = 0
            ds = 0
        else
            e = exp(-sqrt2/t)
            a = two_pi/(1 - t)*e
            s = exp(-a)
            ! da/dt = 2 pi e (1/(1 - t)^2 + sqrt(2)/(t^2 (1 - t)))
            ds = -s*two_pi*e*(1/(1 - t)**2 + sqrt2/(t*t*(1 - t)))
        end if
    end subroutine smooth_step

    ! Points xs and weights ws (sum ws f(xs) is the continuation at x) and
    ! the derivatives of the weights with respect to x.
    pure subroutine periodize_points(x, x_mid, a, w, xs, ws, dws)
        real(dp), intent(in) :: x, x_mid, a, w
        real(dp), intent(out) :: xs(3), ws(3), dws(3)
        real(dp) :: period, x0, s1, d1, s2, d2, w1, w2, dw1, dw2
        period = 2*(a + w)
        x0 = x_mid - period/2 + modulo(x - (x_mid - period/2), period)
        xs = [x0, x0 - period, x0 + period]
        call smooth_step((x0 - x_mid - a)/(2*w), s1, d1)
        call smooth_step((x0 - x_mid + a + 2*w)/(2*w), s2, d2)
        w1 = s1
        dw1 = d1/(2*w)
        w2 = s2
        dw2 = d2/(2*w)
        ws = [w1*(1 - w2), (1 - w1)*(1 - w2), w2]
        dws = [dw1*(1 - w2) - w1*dw2, -dw1*(1 - w2) - (1 - w1)*dw2, dw2]
    end subroutine periodize_points

end module fortnum_periodize
