!> Interval arithmetic against real128 oracles: containment of the exact
!> value for point and interval arguments, tightness (a few ulps for point
!> arguments), and containment of exact identities.
program test_fortnum_interval
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, cinterval_t, interval, cinterval, &
        operator(+), operator(-), operator(*), operator(/), operator(**), &
        sqrt, exp, log, sin, cos, sinh, cosh, abs, sqr, interval_pi, hull, &
        intersect, mid, width, mag, contains, contains_zero, disjoint, &
        is_empty, cabs_up, cabs_down
    implicit none

    integer, parameter :: npt = 20000
    integer :: i, j, nfail, nbad, nwide, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: x, a, b, lo, hi, r(4)
    real(qp) :: xq, pq
    type(interval_t) :: ia, ib, ic, one
    type(cinterval_t) :: za, zb, zc

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 17
    call random_seed(put=seed)

    pq = 4.0_qp*atan(1.0_qp)
    ia = interval_pi()
    call require(encl(ia, pq) .and. width(ia) <= 1.0e-15_dp, &
        "interval_pi encloses pi to one ulp", nfail)

    ! Point arguments: exp, log, sin, cos, sinh, cosh, sqrt.
    nbad = 0
    nwide = 0
    do i = 1, npt
        call random_number(r)
        x = (r(1) - 0.5_dp)*2.0_dp*10.0_dp**(4.0_dp*r(2) - 2.0_dp)
        if (r(3) < 0.02_dp) x = (r(1) - 0.5_dp)*1400.0_dp
        xq = real(x, qp)
        call check1(exp(interval(x)), exp(xq), 16)
        call check1(sinh(interval(x)), sinh(xq), 32)
        call check1(cosh(interval(x)), cosh(xq), 32)
        if (abs(x) < 1.0e6_dp) then
            call check1(sin(interval(x)), sin(xq), 32)
            call check1(cos(interval(x)), cos(xq), 32)
        end if
        call check1(log(interval(abs(x))), log(abs(xq)), 32)
        call check1(sqrt(interval(abs(x))), sqrt(abs(xq)), 2)
    end do
    call require(nbad == 0, "point transcendentals enclose real128 values", nfail)
    call require(nwide == 0, "point transcendentals are within 32 ulps", nfail)

    ! Large and tiny arguments.
    nbad = 0
    nwide = 0
    call check1(exp(interval(-700.0_dp)), exp(-700.0_qp), 16)
    call check1(exp(interval(709.0_dp)), exp(709.0_qp), 16)
    call check1(sin(interval(1.0e5_dp)), sin(1.0e5_qp), 32)
    call check1(sin(interval(1.0e-300_dp)), sin(1.0e-300_qp), 32)
    call check1(sinh(interval(1.0e-300_dp)), sinh(1.0e-300_qp), 32)
    call check1(log(interval(1.0_dp)), 0.0_qp, 4)
    call check1(log(interval(huge(1.0_dp))), log(real(huge(1.0_dp), qp)), 32)
    call check1(log(interval(tiny(1.0_dp))), log(real(tiny(1.0_dp), qp)), 32)
    ia = sin(interval(1.0e12_dp))
    call require(encl(ia, sin(1.0e12_qp)), "sin(1e12) enclosed", nfail)
    call require(nbad == 0 .and. nwide == 0, "extreme arguments", nfail)
    ia = exp(interval(800.0_dp))
    call require(ia%lo >= huge(1.0_dp) .and. ia%hi > huge(1.0_dp), &
        "exp overflow encloses", nfail)
    ia = exp(interval(-800.0_dp))
    call require(ia%lo >= 0.0_dp .and. ia%hi > 0.0_dp, "exp underflow encloses", &
        nfail)

    ! Interval arguments: sampled points must lie in the enclosure.
    nbad = 0
    do i = 1, 4000
        call random_number(r)
        a = (r(1) - 0.5_dp)*20.0_dp
        b = a + r(2)*r(2)*8.0_dp
        ia = interval(a, b)
        ib = interval(r(3) - 2.0_dp, r(4) + 0.5_dp)
        do j = 0, 8
            xq = real(a, qp) + (real(b, qp) - real(a, qp))*real(j, qp)/8.0_qp
            call check_in(sin(ia), sin(xq))
            call check_in(cos(ia), cos(xq))
            call check_in(exp(ia), exp(xq))
            call check_in(sinh(ia), sinh(xq))
            call check_in(cosh(ia), cosh(xq))
            call check_in(sqr(ia), xq*xq)
            call check_in(ia**3, xq**3)
            call check_in(ia**(-2), 1.0_qp/xq**2)
            call check_in(abs(ia), abs(xq))
            call check_in(ia*ib, xq*real(ib%lo, qp))
            call check_in(ia/ib, xq/real(ib%hi, qp))
            call check_in(ia - ib, xq - real(ib%lo, qp))
            if (xq > 0.0_qp) then
                call check_in(log(ia), log(xq))
                call check_in(sqrt(ia), sqrt(xq))
            end if
        end do
    end do
    call require(nbad == 0, "interval functions contain sampled real128 values", &
        nfail)
    ia = sin(interval(1.0_dp, 2.0_dp))
    call require(ia%hi >= 1.0_dp, "sin range includes the maximum at pi/2", nfail)
    ia = cos(interval(3.0_dp, 3.2_dp))
    call require(ia%lo <= -1.0_dp, "cos range includes the minimum at pi", nfail)

    ! Division by an interval containing zero, emptiness, set operations.
    ia = interval(1.0_dp) / interval(-1.0_dp, 1.0_dp)
    call require(ia%lo < -huge(1.0_dp) .and. ia%hi > huge(1.0_dp), &
        "division by zero-containing interval gives the entire line", nfail)
    ia = intersect(interval(0.0_dp, 1.0_dp), interval(2.0_dp, 3.0_dp))
    call require(is_empty(ia), "disjoint intersection is empty", nfail)
    ia = hull(interval(0.0_dp, 1.0_dp), interval(2.0_dp, 3.0_dp))
    call require(contains(ia, 1.5_dp) .and. contains_zero(ia), "hull", nfail)
    call require(disjoint(interval(0.0_dp, 1.0_dp), interval(2.0_dp, 3.0_dp)), &
        "disjoint", nfail)

    ! Exact identities.
    nbad = 0
    one = interval(1.0_dp)
    do i = 1, 2000
        call random_number(r)
        x = (r(1) - 0.5_dp)*40.0_dp
        ia = interval(x, x + r(2)*1.0e-6_dp)
        if (.not. contains(sqr(sin(ia)) + sqr(cos(ia)), 1.0_dp)) nbad = nbad + 1
        ib = interval(abs(x) + 0.1_dp)
        if (.not. contains(exp(log(ib)), ib%lo)) nbad = nbad + 1
        if (.not. contains(sqr(cosh(ia)) - sqr(sinh(ia)), 1.0_dp)) nbad = nbad + 1
        ic = intersect(exp(ia + ib), exp(ia)*exp(ib))
        if (is_empty(ic)) nbad = nbad + 1
    end do
    ia = sin(interval_pi())
    if (.not. contains_zero(ia)) nbad = nbad + 1
    ia = cos(interval_pi()/2)
    if (.not. contains_zero(ia)) nbad = nbad + 1
    call require(nbad == 0, "identities sin^2+cos^2, exp(log), cosh^2-sinh^2, &
        &exp(a+b) hold by containment", nfail)

    ! Complex boxes against real128 complex arithmetic at sampled points.
    nbad = 0
    do i = 1, 2000
        call random_number(r)
        za = cinterval(interval(r(1), r(1) + 1.0e-3_dp), &
            interval(r(2) - 0.5_dp, r(2) - 0.5_dp + 1.0e-3_dp))
        zb = cinterval(interval(r(3) + 0.1_dp), interval(r(4) - 0.5_dp))
        call check_c(za*zb, cmplx(r(1), r(2) - 0.5_dp, qp)*cmplx(r(3) + 0.1_dp, &
            r(4) - 0.5_dp, qp))
        call check_c(za/zb, cmplx(r(1), r(2) - 0.5_dp, qp)/cmplx(r(3) + 0.1_dp, &
            r(4) - 0.5_dp, qp))
        call check_c(sin(za), sin(cmplx(r(1), r(2) - 0.5_dp, qp)))
        call check_c(cos(za), cos(cmplx(r(1), r(2) - 0.5_dp, qp)))
        call check_c(exp(za), exp(cmplx(r(1), r(2) - 0.5_dp, qp)))
        call check_c(za**3, cmplx(r(1), r(2) - 0.5_dp, qp)**3)
        call check_c(sqrt(zb), sqrt(cmplx(r(3) + 0.1_dp, r(4) - 0.5_dp, qp)))
        zc = zb
        lo = cabs_down(zc)
        hi = cabs_up(zc)
        xq = abs(cmplx(r(3) + 0.1_dp, r(4) - 0.5_dp, qp))
        if (.not. (real(lo, qp) <= xq .and. xq <= real(hi, qp))) nbad = nbad + 1
    end do
    call require(nbad == 0, "complex boxes contain real128 values", nfail)
    ia = interval(mid(interval(1.0_dp, 3.0_dp)) + mag(interval(-4.0_dp, 1.0_dp)))
    call require(contains(ia, 6.0_dp), "mid and mag", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " interval test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_interval: all tests passed"

contains

    logical function encl(v, e)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e

        encl = real(v%lo, qp) <= e .and. e <= real(v%hi, qp)
    end function encl

    !> Containment plus tightness: width at most nulp ulps of the value.
    subroutine check1(v, e, nulp)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e
        integer, intent(in) :: nulp

        if (.not. encl(v, e)) then
            nbad = nbad + 1
            write (error_unit, *) "not enclosed:", v%lo, v%hi, real(e, dp)
        end if
        if (real(v%hi, qp) - real(v%lo, qp) > &
            real(nulp, qp)*spacing(max(abs(real(e, dp)), tiny(1.0_dp)))) then
            nwide = nwide + 1
            write (error_unit, *) "too wide:", v%lo, v%hi, real(e, dp)
        end if
    end subroutine check1

    subroutine check_in(v, e)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e

        if (.not. encl(v, e)) nbad = nbad + 1
    end subroutine check_in

    subroutine check_c(v, e)
        type(cinterval_t), intent(in) :: v
        complex(qp), intent(in) :: e

        if (.not. (encl(v%re, real(e, qp)) .and. encl(v%im, aimag(e)))) &
            nbad = nbad + 1
    end subroutine check_c

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_interval
