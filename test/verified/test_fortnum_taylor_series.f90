!> Order-by-order Taylor arithmetic: coefficients of exp(t), 1/(1-t),
!> sin/cos(t), and sqrt(t) expanded about a random center t0, checked against
!> the closed-form derivative formulas exp(t0)/k!, 1/(1-t0)^(k+1),
!> sin/cos(t0 + k*pi/2)/k!, and sqrt(t0)*binomial(1/2,k)/t0^k. The ODE
!> y' = y^2, y(0) = y0 is integrated coefficient by coefficient with ts_mul
!> and checked against the closed-form solution y(t) = y0/(1 - y0 t), whose
!> Taylor coefficients are y0^(k+1); for idual and cdual coefficients the
!> seeded derivative/holomorphic-derivative of y_k with respect to y0,
!> (k+1) y0^k, is checked too.
program test_fortnum_taylor_series
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, interval, cinterval, real_part, &
        imag_part, operator(+), operator(-), operator(*), operator(/)
    use fortnum_idual, only: idual_t, idual_var, idual_const, operator(+), &
        operator(-), operator(*), operator(/)
    use fortnum_cdual, only: cdual_t, cdual_var, cdual_const, operator(+), &
        operator(-), operator(*), operator(/)
    use fortnum_taylor_series, only: ts_mul, ts_div, ts_sincos, ts_exp, &
        ts_sqrt, ts_lin, ts_zero
    implicit none

    integer, parameter :: q = 6
    integer :: i, k, nfail, nbad, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: s(2), t0, x0, y0
    complex(dp) :: zy0

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 41
    call random_seed(put=seed)

    nbad = 0
    do i = 1, 500
        call random_number(s)
        t0 = 0.3_dp + 0.5_dp*s(1)
        call check_known_functions(t0, nbad)
    end do
    call require(nbad == 0, "exp/1-1-t/sincos/sqrt series match closed-form "// &
        "derivative formulas", nfail)

    nbad = 0
    do i = 1, 500
        call random_number(s)
        y0 = 0.1_dp + 0.3_dp*s(1)
        call check_ode_interval(y0, nbad)
        call check_ode_idual(y0, nbad)
        zy0 = cmplx(y0, 0.2_dp + 0.3_dp*s(2), dp)
        call check_ode_cdual(zy0, nbad)
    end do
    call require(nbad == 0, "y'=y^2 Taylor coefficients and their y0-"// &
        "derivatives match the closed-form flow y0/(1-y0 t)", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " taylor_series test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_taylor_series: all tests passed"

contains

    subroutine check_known_functions(t0, nbad)
        real(dp), intent(in) :: t0
        integer, intent(inout) :: nbad
        type(interval_t) :: a(0:q), e(0:q), one(0:q), oneminus(0:q), c(0:q)
        type(interval_t) :: sn(0:q), cs(0:q), sq(0:q)
        real(qp) :: t0q, fact, binom, ref
        integer :: k

        t0q = real(t0, qp)
        a(0) = interval(t0)
        a(1) = interval(1.0_dp)
        do k = 2, q
            call ts_zero(k, a)
        end do
        do k = 0, q
            call ts_exp(k, a, e)
        end do
        fact = 1.0_qp
        do k = 0, q
            if (k > 0) fact = fact*real(k, qp)
            ref = exp(t0q)/fact
            if (.not. encl(e(k), ref)) nbad = nbad + 1
        end do

        one(0) = interval(1.0_dp)
        do k = 1, q
            call ts_zero(k, one)
        end do
        oneminus(0) = interval(1.0_dp) - a(0)
        oneminus(1) = interval(0.0_dp) - a(1)
        do k = 2, q
            call ts_zero(k, oneminus)
        end do
        do k = 0, q
            call ts_div(k, one, oneminus, c)
            ref = 1.0_qp/(1.0_qp - t0q)**(k + 1)
            if (.not. encl(c(k), ref)) nbad = nbad + 1
        end do

        do k = 0, q
            call ts_sincos(k, a, sn, cs)
        end do
        fact = 1.0_qp
        do k = 0, q
            if (k > 0) fact = fact*real(k, qp)
            ref = sin(t0q + real(k, qp)*qacos(-1.0_qp)/2.0_qp)/fact
            if (.not. encl(sn(k), ref)) nbad = nbad + 1
            ref = cos(t0q + real(k, qp)*qacos(-1.0_qp)/2.0_qp)/fact
            if (.not. encl(cs(k), ref)) nbad = nbad + 1
        end do

        do k = 0, q
            call ts_sqrt(k, a, sq)
        end do
        binom = 1.0_qp
        do k = 0, q
            if (k > 0) binom = binom*(0.5_qp - real(k - 1, qp))/real(k, qp)
            ref = sqrt(t0q)*binom/t0q**k
            if (.not. encl(sq(k), ref)) nbad = nbad + 1
        end do

        ! Linear combination: 3 exp(t) - 2/(1-t), coefficient by coefficient.
        block
            type(interval_t) :: lc(0:q)
            fact = 1.0_qp
            do k = 0, q
                call ts_lin(k, 3.0_dp, e, -2.0_dp, c, lc)
                if (k > 0) fact = fact*real(k, qp)
                ref = 3.0_qp*exp(t0q)/fact - 2.0_qp/(1.0_qp - t0q)**(k + 1)
                if (.not. encl(lc(k), ref)) nbad = nbad + 1
            end do
        end block
    end subroutine check_known_functions

    pure function qacos(x) result(r)
        real(qp), intent(in) :: x
        real(qp) :: r
        r = acos(x)
    end function qacos

    subroutine check_ode_interval(y0, nbad)
        real(dp), intent(in) :: y0
        integer, intent(inout) :: nbad
        type(interval_t) :: y(0:q), g(0:q)
        real(qp) :: y0q, ref
        integer :: k

        y(0) = interval(y0)
        do k = 1, q
            call ts_mul(k - 1, y, y, g)
            y(k) = g(k - 1)/interval(real(k, dp))
        end do
        y0q = real(y0, qp)
        do k = 0, q
            ref = y0q**(k + 1)
            if (.not. encl(y(k), ref)) nbad = nbad + 1
        end do
    end subroutine check_ode_interval

    subroutine check_ode_idual(y0, nbad)
        real(dp), intent(in) :: y0
        integer, intent(inout) :: nbad
        type(idual_t) :: y(0:q), g(0:q)
        real(qp) :: y0q, ref, refd
        integer :: k

        y(0) = idual_var(interval(y0), 1, 1)
        do k = 1, q
            call ts_mul(k - 1, y, y, g)
            y(k) = g(k - 1)/idual_const(interval(real(k, dp)), 1)
        end do
        y0q = real(y0, qp)
        do k = 0, q
            ref = y0q**(k + 1)
            refd = real(k + 1, qp)*y0q**k
            if (.not. encl(y(k)%v, ref)) nbad = nbad + 1
            if (.not. encl(y(k)%d(1), refd)) nbad = nbad + 1
        end do
    end subroutine check_ode_idual

    subroutine check_ode_cdual(zy0, nbad)
        complex(dp), intent(in) :: zy0
        integer, intent(inout) :: nbad
        type(cdual_t) :: y(0:q), g(0:q)
        complex(qp) :: zy0q, ref, refd
        integer :: k

        y(0) = cdual_var(cinterval(real(zy0, dp), aimag(zy0)), 1, 1)
        do k = 1, q
            call ts_mul(k - 1, y, y, g)
            y(k) = g(k - 1)/cdual_const(cinterval(real(k, dp), 0.0_dp), 1)
        end do
        zy0q = cmplx(real(real(zy0, dp), qp), real(aimag(zy0), qp), qp)
        do k = 0, q
            ref = zy0q**(k + 1)
            refd = real(k + 1, qp)*zy0q**k
            if (.not. encl_c(y(k)%v, ref)) nbad = nbad + 1
            if (.not. encl_c(y(k)%d(1), refd)) nbad = nbad + 1
        end do
    end subroutine check_ode_cdual

    logical function encl(v, e)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e
        encl = real(v%lo, qp) <= e .and. e <= real(v%hi, qp)
    end function encl

    logical function encl_c(v, e)
        use fortnum_interval, only: cinterval_t
        type(cinterval_t), intent(in) :: v
        complex(qp), intent(in) :: e
        type(interval_t) :: re, im
        re = real_part(v)
        im = imag_part(v)
        encl_c = real(re%lo, qp) <= real(e, qp) .and. &
            real(e, qp) <= real(re%hi, qp) .and. &
            real(im%lo, qp) <= aimag(e) .and. aimag(e) <= real(im%hi, qp)
    end function encl_c

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail
        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_taylor_series
