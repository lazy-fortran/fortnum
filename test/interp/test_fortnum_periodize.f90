program test_fortnum_periodize
    ! Periodic continuation (fortnum_periodize): unchanged on the core, value
    ! and first derivative continuous across the window edge for a smooth
    ! aperiodic function, weights summing to one, weight derivatives against
    ! central differences, and the smooth step's limits and derivative.
    use fortnum_kinds, only: dp
    use fortnum_periodize, only: periodize_points, smooth_step
    implicit none
    real(dp), parameter :: xm = 2.0_dp, a = 1.5_dp, w = 0.25_dp
    real(dp) :: period, x, g, gp, gl, glp, gr, grp, xs(3), ws(3), dws(3), h, s, ds, sp, sm, dum
    real(dp) :: wsp(3), wsm(3)
    integer :: nfail, i
    nfail = 0
    period = 2*(a + w)
    ! Core: the function itself.
    do i = 0, 10
        x = xm - a + 2*a*i/10.0_dp
        call cont(x, g, gp)
        call check('core value', g, f(x), 1e-15_dp)
        call check('core derivative', gp, fp(x), 1e-14_dp)
    end do
    ! Across the window edge x_mid + L/2 the continuation is continuous.
    h = 1e-9_dp
    call cont(xm + period/2 - h, gl, glp)
    call cont(xm + period/2 + h, gr, grp)
    call check('edge value', gl, gr, 1e-7_dp)
    call check('edge derivative', glp, grp, 1e-6_dp)
    ! Weights sum to one; their derivatives match central differences.
    do i = 1, 9
        x = xm - period/2 + period*i/10.0_dp
        call periodize_points(x, xm, a, w, xs, ws, dws)
        call check('weights sum', sum(ws), 1.0_dp, 1e-15_dp)
        h = 1e-6_dp
        call periodize_points(x + h, xm, a, w, xs, wsp, dws)
        call periodize_points(x - h, xm, a, w, xs, wsm, dws)
        call periodize_points(x, xm, a, w, xs, ws, dws)
        call check_abs('weight derivative', dws(1), (wsp(1) - wsm(1))/(2*h), 1e-6_dp)
        call check_abs('weight derivative', dws(3), (wsp(3) - wsm(3))/(2*h), 1e-6_dp)
    end do
    call smooth_step(-0.1_dp, s, ds)
    call check('step below', s, 1.0_dp, 0.0_dp)
    call smooth_step(1.1_dp, s, ds)
    call check_abs('step above', s, 0.0_dp, 0.0_dp)
    call smooth_step(0.4_dp + 1e-7_dp, sp, dum)
    call smooth_step(0.4_dp - 1e-7_dp, sm, dum)
    call smooth_step(0.4_dp, s, ds)
    call check('step derivative', ds, (sp - sm)/2e-7_dp, 1e-6_dp)
    if (nfail > 0) error stop 'FAIL'
    print '(a)', 'PASS'
contains
    pure real(dp) function f(x)
        real(dp), intent(in) :: x
        f = exp(0.3_dp*x)*sin(1.7_dp*x)
    end function f
    pure real(dp) function fp(x)
        real(dp), intent(in) :: x
        fp = exp(0.3_dp*x)*(0.3_dp*sin(1.7_dp*x) + 1.7_dp*cos(1.7_dp*x))
    end function fp
    subroutine cont(x, g, gp)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: g, gp
        real(dp) :: xs(3), ws(3), dws(3)
        integer :: k
        call periodize_points(x, xm, a, w, xs, ws, dws)
        g = 0
        gp = 0
        do k = 1, 3
            g = g + ws(k)*f(xs(k))
            gp = gp + ws(k)*fp(xs(k)) + dws(k)*f(xs(k))
        end do
    end subroutine cont
    subroutine check(name, got, want, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        if (abs(got - want) <= tol*max(abs(want), 1e-300_dp)) return
        nfail = nfail + 1
        print '(a,a,2es24.15)', 'FAIL ', name, got, want
    end subroutine check
    subroutine check_abs(name, got, want, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, want, tol
        if (abs(got - want) <= tol) return
        nfail = nfail + 1
        print '(a,a,2es24.15)', 'FAIL ', name, got, want
    end subroutine check_abs
end program test_fortnum_periodize
