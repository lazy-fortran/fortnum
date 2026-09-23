!> Exhaustive random test of the successor-formula rounding primitives against
!> the exact neighbours from the nearest intrinsic: round_up(x) must be
!> succ(x) or succ(succ(x)), round_down(x) pred(x) or pred(pred(x)), for
!> finite x over all binades including subnormals.
program test_fortnum_rounding
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        int64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_positive_inf
    use fortnum_rounding, only: round_up, round_down, sum_up, sum_down, &
        gamma_up, mul_up, mul_down, div_up, div_down, sqrt_up, sqrt_down, &
        unit_roundoff
    implicit none

    integer, parameter :: nrand = 2000000
    real(dp) :: x, s1, s2, inf, a, b, y(1000)
    real(qp) :: exact
    integer :: i, nfail, nbad, seed_size
    integer, allocatable :: seed(:)

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 20260923
    call random_seed(put=seed)

    nbad = 0
    do i = 1, nrand
        x = random_finite()
        if (.not. neighbour_ok(x)) nbad = nbad + 1
    end do
    call require(nbad == 0, "random bit patterns: round_up/down within two ulps", &
        nfail)

    nbad = 0
    do i = -1074, 1023
        x = 2.0_dp**i
        if (.not. neighbour_ok(x)) nbad = nbad + 1
        if (.not. neighbour_ok(-x)) nbad = nbad + 1
        if (.not. neighbour_ok(nearest(x, 1.0_dp))) nbad = nbad + 1
        if (.not. neighbour_ok(nearest(x, -1.0_dp))) nbad = nbad + 1
    end do
    if (.not. neighbour_ok(0.0_dp)) nbad = nbad + 1
    if (.not. neighbour_ok(huge(1.0_dp))) nbad = nbad + 1
    if (.not. neighbour_ok(-huge(1.0_dp))) nbad = nbad + 1
    if (.not. neighbour_ok(tiny(1.0_dp))) nbad = nbad + 1
    call require(nbad == 0, "powers of two, binade edges, subnormals, huge", nfail)

    inf = ieee_value(1.0_dp, ieee_positive_inf)
    call require(round_up(inf) > huge(1.0_dp) .and. &
        round_down(-inf) < -huge(1.0_dp), &
        "infinities stay infinite in the outward direction", nfail)
    call require(same(round_down(inf), huge(1.0_dp)) .and. &
        same(round_up(-inf), -huge(1.0_dp)), &
        "an overflowed endpoint becomes the largest finite inward bound", nfail)
    call require(round_up(huge(1.0_dp)) > huge(1.0_dp), "round_up(huge) overflows", &
        nfail)

    ! Directed operations against real128 (products and sqrt of binary64
    ! values are exact or far more accurate in real128).
    nbad = 0
    do i = 1, 200000
        a = random_moderate()
        b = random_moderate()
        exact = real(a, qp)*real(b, qp)
        if (.not. (real(mul_down(a, b), qp) <= exact .and. &
            exact <= real(mul_up(a, b), qp))) nbad = nbad + 1
        exact = real(a, qp)/real(b, qp)
        if (.not. (real(div_down(a, b), qp) <= exact .and. &
            exact <= real(div_up(a, b), qp))) nbad = nbad + 1
        exact = sqrt(real(abs(a), qp))
        if (.not. (real(sqrt_down(abs(a)), qp) <= exact .and. &
            exact <= real(sqrt_up(abs(a)), qp))) nbad = nbad + 1
    end do
    call require(nbad == 0, "directed mul, div, sqrt enclose real128 results", &
        nfail)

    ! Sums with cancellation: the real128 sum of 1000 binary64 values of
    ! moderate range is exact.
    nbad = 0
    do i = 1, 200
        call random_number(y)
        y = (y - 0.5_dp)*10.0_dp**nint(6.0_dp*y(1))
        exact = sum(real(y, qp))
        s1 = sum_down(y)
        s2 = sum_up(y)
        if (.not. (real(s1, qp) <= exact .and. exact <= real(s2, qp))) &
            nbad = nbad + 1
    end do
    call require(nbad == 0, "sum_down <= exact sum <= sum_up", nfail)

    exact = 10.0_qp*real(unit_roundoff, qp)/(1.0_qp - 10.0_qp*real(unit_roundoff, qp))
    call require(real(gamma_up(10), qp) >= exact .and. &
        real(gamma_up(10), qp) <= exact*(1.0_qp + 1.0e-14_qp), &
        "gamma_up(10) is a tight upper bound", nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " rounding test(s) failed"
        error stop 1
    end if
    deallocate (seed)
    print '(a)', "fortnum_rounding: all tests passed"

contains

    logical function same(p, q)
        real(dp), intent(in) :: p, q

        same = p <= q .and. q <= p
    end function same

    logical function neighbour_ok(v) result(ok)
        real(dp), intent(in) :: v
        real(dp) :: u1, u2, d1, d2, ru, rd

        ru = round_up(v)
        rd = round_down(v)
        u1 = nearest(v, 1.0_dp)
        d1 = nearest(v, -1.0_dp)
        ok = ru >= u1 .and. rd <= d1
        if (ieee_is_finite(u1)) then
            u2 = nearest(u1, 1.0_dp)
            ok = ok .and. ru <= u2
        end if
        if (ieee_is_finite(d1)) then
            d2 = nearest(d1, -1.0_dp)
            ok = ok .and. rd >= d2
        end if
    end function neighbour_ok

    !> A uniformly random finite binary64 bit pattern.
    real(dp) function random_finite() result(v)
        real(dp) :: r(2)
        integer(int64) :: hi, lo, bits

        do
            call random_number(r)
            hi = int(r(1)*4294967296.0_dp, int64)
            lo = int(r(2)*4294967296.0_dp, int64)
            bits = ior(ishft(hi, 32), lo)
            v = transfer(bits, 1.0_dp)
            if (ieee_is_finite(v)) exit
        end do
    end function random_finite

    !> A random value with exponent in [-60, 60] and random sign.
    real(dp) function random_moderate() result(v)
        real(dp) :: r(3)

        call random_number(r)
        v = (0.5_dp + r(1))*2.0_dp**nint(120.0_dp*r(2) - 60.0_dp)
        if (r(3) < 0.5_dp) v = -v
    end function random_moderate

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_rounding
