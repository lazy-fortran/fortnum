!> Interval forward AD: on random boxes, the value and gradient enclosures of
!> f(x, y) = sin(x) exp(y)/(1 + x^2) + sqrt(x + y) log(2 + y) + (x y)^3
!>           - cos(x)/y
!> must contain the hand-derived real128 value and gradient at sampled
!> points of the box (point boxes make the enclosures tight, so a wrong
!> derivative rule fails).
program test_fortnum_idual
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, interval
    use fortnum_idual, only: idual_t, idual_var, idual_const, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sqrt, exp, log, &
        sin, cos, sqr
    implicit none

    integer :: i, j, nfail, nbad, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: s(4), x0, y0, hx, hy
    real(qp) :: xq, yq, f, fx, fy
    type(idual_t) :: x, y, g, c

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 99
    call random_seed(put=seed)

    nbad = 0
    do i = 1, 3000
        call random_number(s)
        x0 = 0.2_dp + 2.0_dp*s(1)
        y0 = 0.5_dp + 2.0_dp*s(2)
        hx = 1.0e-3_dp*s(3)
        hy = 1.0e-3_dp*s(4)
        if (modulo(i, 3) == 0) then
            hx = 0.0_dp
            hy = 0.0_dp
        end if
        x = idual_var(interval(x0, x0 + hx), 1, 2)
        y = idual_var(interval(y0, y0 + hy), 2, 2)
        g = sin(x)*exp(y)/(1.0_dp + sqr(x)) + sqrt(x + y)*log(2.0_dp + y) &
            + (x*y)**3 - cos(x)/y
        do j = 0, 4
            xq = real(x0, qp) + real(hx, qp)*real(j, qp)/4.0_qp
            yq = real(y0, qp) + real(hy, qp)*real(4 - j, qp)/4.0_qp
            f = sin(xq)*exp(yq)/(1 + xq**2) + sqrt(xq + yq)*log(2 + yq) &
                + (xq*yq)**3 - cos(xq)/yq
            fx = (cos(xq)*(1 + xq**2) - 2*xq*sin(xq))*exp(yq)/(1 + xq**2)**2 &
                + log(2 + yq)/(2*sqrt(xq + yq)) + 3*xq**2*yq**3 + sin(xq)/yq
            fy = sin(xq)*exp(yq)/(1 + xq**2) + log(2 + yq)/(2*sqrt(xq + yq)) &
                + sqrt(xq + yq)/(2 + yq) + 3*xq**3*yq**2 + cos(xq)/yq**2
            if (.not. encl(g%v, f)) nbad = nbad + 1
            if (.not. encl(g%d(1), fx)) nbad = nbad + 1
            if (.not. encl(g%d(2), fy)) nbad = nbad + 1
        end do
    end do
    call require(nbad == 0, "value and gradient enclosures contain real128 values", &
        nfail)

    c = idual_const(interval(3.0_dp), 2)
    g = c*x - 2.0_dp
    call require(encl(g%d(1), 3.0_qp) .and. encl(g%d(2), 0.0_qp), &
        "constants have zero gradient", nfail)
    g = x**0
    call require(encl(g%v, 1.0_qp) .and. encl(g%d(1), 0.0_qp), "x**0 = 1", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " idual test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_idual: all tests passed"

contains

    logical function encl(v, e)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e

        encl = real(v%lo, qp) <= e .and. e <= real(v%hi, qp)
    end function encl

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_idual
