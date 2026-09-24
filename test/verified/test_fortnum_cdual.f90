!> Complex-dual forward AD: on random points z1, z2, the value and derivative
!> enclosures of
!>   g(z1, z2) = sin(z1) exp(z2)/(1 + z1) + sqrt(z1 + z2 + 4) + z1/z2 - cos(z2)
!> must contain complex128 central-difference estimates of dg/dz1 and dg/dz2
!> taken along both the real and the imaginary axis (Cauchy-Riemann
!> consistency: an analytic function has the same complex derivative from
!> either direction, so a directionally wrong chain rule fails this check
!> even though the plain value enclosure would still pass).
program test_fortnum_cdual
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, interval, cinterval_t, cinterval, &
        real_part, imag_part
    use fortnum_cdual, only: cdual_t, cdual_var, cdual_const, operator(+), &
        operator(-), operator(*), operator(/), sqrt, exp, sin, cos
    implicit none

    integer :: i, nfail, nbad, seed_size
    integer, allocatable :: seed(:)
    real(dp) :: s(4), x1, y1, x2, y2
    complex(qp) :: z1, z2, h, fp, fm, d1, d2
    type(cdual_t) :: a, b, g, scaled

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 17
    call random_seed(put=seed)
    h = cmplx(1.0e-6_qp, 0.0_qp, qp)

    nbad = 0
    do i = 1, 2000
        call random_number(s)
        x1 = 0.3_dp + 1.5_dp*s(1)
        y1 = 0.2_dp + 1.5_dp*s(2)
        x2 = 0.6_dp + 1.5_dp*s(3)
        y2 = 0.3_dp + 1.5_dp*s(4)
        z1 = cmplx(real(x1, qp), real(y1, qp), qp)
        z2 = cmplx(real(x2, qp), real(y2, qp), qp)

        a = cdual_var(cinterval(x1, y1), 1, 2)
        b = cdual_var(cinterval(x2, y2), 2, 2)
        g = eval(a, b)

        ! Real-axis central difference for d/dz1, d/dz2.
        fp = feval(z1 + h, z2)
        fm = feval(z1 - h, z2)
        d1 = (fp - fm)/(2.0_qp*h)
        fp = feval(z1, z2 + h)
        fm = feval(z1, z2 - h)
        d2 = (fp - fm)/(2.0_qp*h)
        if (.not. encl_c(g%d(1), d1)) nbad = nbad + 1
        if (.not. encl_c(g%d(2), d2)) nbad = nbad + 1

        ! Imaginary-axis central difference: dg/dz = (f(z+ih)-f(z-ih))/(2ih).
        fp = feval(z1 + cmplx(0.0_qp, 1.0_qp, qp)*h, z2)
        fm = feval(z1 - cmplx(0.0_qp, 1.0_qp, qp)*h, z2)
        d1 = (fp - fm)/(2.0_qp*cmplx(0.0_qp, 1.0_qp, qp)*h)
        fp = feval(z1, z2 + cmplx(0.0_qp, 1.0_qp, qp)*h)
        fm = feval(z1, z2 - cmplx(0.0_qp, 1.0_qp, qp)*h)
        d2 = (fp - fm)/(2.0_qp*cmplx(0.0_qp, 1.0_qp, qp)*h)
        if (.not. encl_c(g%d(1), d1)) nbad = nbad + 1
        if (.not. encl_c(g%d(2), d2)) nbad = nbad + 1

        if (.not. encl_c(g%v, feval(z1, z2))) nbad = nbad + 1
    end do
    call require(nbad == 0, "value and derivative enclosures contain "// &
        "real- and imaginary-axis central differences (Cauchy-Riemann)", nfail)

    ! Constant has zero derivative.
    block
        type(interval_t) :: dre, dim
        a = cdual_const(cinterval(2.0_dp, -1.0_dp), 2)
        dre = real_part(a%d(1))
        dim = imag_part(a%d(1))
        call require(dre%lo == 0.0_dp .and. dre%hi == 0.0_dp .and. &
            dim%lo == 0.0_dp .and. dim%hi == 0.0_dp, &
            "constant has zero derivative", nfail)
    end block

    ! Multiplication by a real interval must scale the value and derivative
    ! enclosures in either operand order.
    block
        type(interval_t) :: scale
        scale = interval(0.25_dp, 0.75_dp)
        a = cdual_var(cinterval(1.0_dp, 2.0_dp), 1, 1)
        scaled = a*scale
        call require(encl_c(scaled%v, cmplx(0.25_qp, 0.5_qp, qp)) .and. &
            encl_c(scaled%v, cmplx(0.75_qp, 1.5_qp, qp)), &
            "complex-dual times real interval encloses endpoint values", nfail)
        call require(scaled%d(1)%re%lo <= 0.25_dp .and. &
            scaled%d(1)%re%hi >= 0.75_dp .and. &
            scaled%d(1)%im%lo <= 0.0_dp .and. scaled%d(1)%im%hi >= 0.0_dp, &
            "complex-dual times real interval encloses derivative", nfail)

        scaled = scale*a
        call require(encl_c(scaled%v, cmplx(0.25_qp, 0.5_qp, qp)) .and. &
            encl_c(scaled%v, cmplx(0.75_qp, 1.5_qp, qp)), &
            "real interval times complex-dual encloses endpoint values", nfail)
        call require(scaled%d(1)%re%lo <= 0.25_dp .and. &
            scaled%d(1)%re%hi >= 0.75_dp .and. &
            scaled%d(1)%im%lo <= 0.0_dp .and. scaled%d(1)%im%hi >= 0.0_dp, &
            "real interval times complex-dual encloses derivative", nfail)
    end block

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " cdual test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_cdual: all tests passed"

contains

    !> The Fortran-dual and complex128-reference evaluations of the same
    !> analytic expression, kept textually parallel on purpose.
    pure function eval(z1, z2) result(g)
        type(cdual_t), intent(in) :: z1, z2
        type(cdual_t) :: g

        g = sin(z1)*exp(z2)/(1.0_dp + z1) + sqrt(z1 + z2 + 4.0_dp) &
            + z1/z2 - cos(z2)
    end function eval

    pure function feval(z1, z2) result(g)
        complex(qp), intent(in) :: z1, z2
        complex(qp) :: g

        g = sin(z1)*exp(z2)/(1.0_qp + z1) + sqrt(z1 + z2 + 4.0_qp) &
            + z1/z2 - cos(z2)
    end function feval

    logical function encl_c(v, e)
        type(cinterval_t), intent(in) :: v
        complex(qp), intent(in) :: e
        real(qp) :: pad

        ! Central-difference truncation error is O(h^2) ~ 1e-12 at h=1e-6;
        ! the enclosure itself is a point (width at machine epsilon), so the
        ! comparison pad only needs to absorb the finite-difference scheme.
        pad = 2.0e-9_qp
        encl_c = real(v%re%lo, qp) - pad <= real(e, qp) .and. &
            real(e, qp) <= real(v%re%hi, qp) + pad .and. &
            real(v%im%lo, qp) - pad <= aimag(e) .and. &
            aimag(e) <= real(v%im%hi, qp) + pad
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
end program test_fortnum_cdual
