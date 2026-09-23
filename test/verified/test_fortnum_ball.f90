!> Complex ball arithmetic against real128 complex arithmetic: for random
!> balls and random points inside them, the real128 result of each operation
!> must lie in the result ball; moduli bounds over the whole exponent range;
!> and containment of exact identities.
program test_fortnum_ball
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_ball, only: ball_t, ball, badd, bsub, bmul, bdiv, bneg, binv, &
        bsqrt, bpowi, bscale, bpoint, bcpoint, benclose, babs_hi, babs_lo, &
        bre_lo, bre_hi, bim_lo, bim_hi, bcontains, bcontains_zero, cabs_hi, &
        cabs_lo, ball_from_cinterval, cinterval_from_ball, operator(+), &
        operator(-), operator(*), operator(/), operator(**)
    use fortnum_interval, only: cinterval_t, interval, cinterval
    implicit none

    integer :: i, j, nfail, nbad, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: r(8), x, y
    real(qp) :: aq
    complex(qp) :: za, zb
    type(ball_t) :: a, b, c
    type(cinterval_t) :: box

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 4242
    call random_seed(put=seed)

    ! Modulus bounds across the exponent range, no libm hypot.
    nbad = 0
    do i = 1, 200000
        call random_number(r)
        x = (r(1) - 0.5_dp)*2.0_dp**nint(2000.0_dp*r(2) - 1000.0_dp)
        y = (r(3) - 0.5_dp)*2.0_dp**nint(2000.0_dp*r(4) - 1000.0_dp)
        if (r(5) < 0.1_dp) y = x*(1.0_dp + r(6))
        aq = abs(cmplx(real(x, qp), real(y, qp), qp))
        if (.not. (real(cabs_lo(cmplx(x, y, dp)), qp) <= aq .and. &
            aq <= real(cabs_hi(cmplx(x, y, dp)), qp))) nbad = nbad + 1
    end do
    call require(nbad == 0, "cabs_lo <= abs(z) <= cabs_hi over all exponents", &
        nfail)

    ! Operations on random balls, checked at random interior points.
    nbad = 0
    do i = 1, 20000
        a = random_ball()
        b = random_ball()
        do j = 1, 4
            za = interior(a)
            zb = interior(b)
            call check(badd(a, b), za + zb)
            call check(bsub(a, b), za - zb)
            call check(bmul(a, b), za*zb)
            call check(bscale(a, 0.37_dp), za*0.37_qp)
            call check(bneg(a), -za)
            call check(bpowi(a, 5), za**5)
            call check(a**(-2), za**(-2))
            if (.not. bcontains_zero(b)) then
                call check(bdiv(a, b), za/zb)
                call check(binv(b), 1.0_qp/zb)
            end if
            if (real(a%c, dp) > 0.0_dp .or. abs(aimag(a%c)) > a%r) then
                call check_sqrt(bsqrt(a), za)
            end if
            call check(a + 2.0_dp*b - (1.0_dp, 3.0_dp), za + 2*zb - (1, 3))
        end do
    end do
    call require(nbad == 0, "ball operations contain real128 results", nfail)

    ! Real and imaginary bounds, modulus bounds of balls, conversions.
    nbad = 0
    do i = 1, 5000
        a = random_ball()
        za = interior(a)
        if (.not. (real(bre_lo(a), qp) <= real(za, qp) .and. &
            real(za, qp) <= real(bre_hi(a), qp))) nbad = nbad + 1
        if (.not. (real(bim_lo(a), qp) <= aimag(za) .and. &
            aimag(za) <= real(bim_hi(a), qp))) nbad = nbad + 1
        if (.not. (real(babs_lo(a), qp) <= abs(za) .and. &
            abs(za) <= real(babs_hi(a), qp))) nbad = nbad + 1
        box = cinterval_from_ball(a)
        c = ball_from_cinterval(box)
        call check(c, za)
    end do
    call require(nbad == 0, "component, modulus bounds and conversions", nfail)

    ! Exact identities by containment.
    nbad = 0
    do i = 1, 5000
        a = random_ball()
        if (bcontains_zero(a)) cycle
        if (.not. bcontains(binv(binv(a)), a%c)) nbad = nbad + 1
        if (.not. bcontains(bmul(a, binv(a)), (1.0_dp, 0.0_dp))) nbad = nbad + 1
        if (real(a%c, dp) > 0.0_dp) then
            if (.not. bcontains(bsqrt(a)**2, a%c)) nbad = nbad + 1
        end if
    end do
    c = bsqrt(bpoint(-4.0_dp))
    call require(c%r > 1.0e300_dp, "sqrt centre on the branch cut is uninformative", &
        nfail)
    c = bsqrt(bcpoint(-4.0_dp, 1.0e-3_dp))
    if (.not. bcontains(c, cmplx(sqrt(cmplx(-4.0_qp, 1.0e-3_qp, qp)), kind=dp))) &
        nbad = nbad + 1
    c = binv(benclose(0.1_dp, 0.2_dp))
    call require(c%r > 1.0e300_dp, "inverse of a ball containing 0", nfail)
    call require(nbad == 0, "identities 1/(1/a) = a, a/a = 1, sqrt(a)^2 = a", &
        nfail)
    box = cinterval(interval(1.0_dp, 2.0_dp), interval(-1.0_dp, 0.5_dp))
    c = ball_from_cinterval(box)
    call require(bcontains(c, (1.0_dp, -1.0_dp)) .and. &
        bcontains(c, (2.0_dp, 0.5_dp)), "ball of a box contains its corners", nfail)
    c = ball((1.0_dp, 0.0_dp), 0.0_dp)
    call require(bcontains(c, (1.0_dp, 0.0_dp)), "point ball contains its point", &
        nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " ball test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_ball: all tests passed"

contains

    type(ball_t) function random_ball() result(bb)
        real(dp) :: s(5)

        call random_number(s)
        bb%c = cmplx(s(1) - 0.5_dp, s(2) - 0.5_dp, dp)*10.0_dp**(4.0_dp*s(3) - 2.0_dp)
        bb%r = s(4)*abs(bb%c)*10.0_dp**(-8.0_dp*s(5))
    end function random_ball

    !> A random point of the ball (in real128, inside by a relative margin).
    complex(qp) function interior(bb) result(z)
        type(ball_t), intent(in) :: bb
        real(dp) :: s(2)
        real(qp) :: rho, th

        call random_number(s)
        rho = real(bb%r, qp)*real(s(1), qp)*(1.0_qp - 1.0e-20_qp)
        th = 2.0_qp*acos(-1.0_qp)*real(s(2), qp)
        z = cmplx(real(bb%c, qp), aimag(cmplx(bb%c, kind=qp)), qp) + &
            rho*cmplx(cos(th), sin(th), qp)
    end function interior

    subroutine check(bb, z)
        type(ball_t), intent(in) :: bb
        complex(qp), intent(in) :: z

        if (.not. (abs(z - cmplx(bb%c, kind=qp)) <= real(bb%r, qp))) &
            nbad = nbad + 1
    end subroutine check

    !> The principal square root: points just across the branch cut map to
    !> the other half-plane, so compare with the root of either sign
    !> convention the real128 library picks near the cut.
    subroutine check_sqrt(bb, z)
        type(ball_t), intent(in) :: bb
        complex(qp), intent(in) :: z
        complex(qp) :: w

        w = sqrt(z)
        if (.not. (abs(w - cmplx(bb%c, kind=qp)) <= real(bb%r, qp))) &
            nbad = nbad + 1
    end subroutine check_sqrt

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_ball
