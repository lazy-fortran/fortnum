!> binary128 intervals: containment of exact rationals, square roots and pi
!> (residual of the nearest binary128 value from an independent 400-bit
!> computation), and the binary64 conversion.
program test_fortnum_interval_qp
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_interval, only: interval_t
    use fortnum_interval_qp, only: qp, qinterval_t, qinterval, qrat, &
        qinterval_pi, qsqrt, qsqr, to_interval, qmid, qrad, operator(+), &
        operator(-), operator(*), operator(/)
    implicit none

    ! pi - (binary128 value nearest pi), from mpmath at 400 bits.
    real(qp), parameter :: pi_residual = 8.6718101301237810248e-35_qp
    type(qinterval_t) :: a, b
    type(interval_t) :: d
    integer :: nfail, k

    nfail = 0
    a = qinterval_pi()
    call require(a%lo < 3.14159265358979323846264338327950288_qp .and. &
        a%hi - 3.14159265358979323846264338327950288_qp >= pi_residual, &
        "qinterval_pi encloses pi", nfail)

    do k = 1, 50
        a = qrat(k, 7)*qinterval(7) - qinterval(k)
        if (.not. (a%lo <= 0.0_qp .and. a%hi >= 0.0_qp)) nfail = nfail + 1
        a = qsqr(qsqrt(qinterval(real(k, qp))))
        if (.not. (a%lo <= real(k, qp) .and. a%hi >= real(k, qp))) nfail = nfail + 1
        if (qrad(qsqrt(qinterval(real(k, qp)))) > 1.0e-32_qp) nfail = nfail + 1
    end do
    call require(nfail == 0, "rationals and square roots", nfail)

    b = qsqrt(qinterval(2)) - qinterval(1.4142135623730950488016887242096980786_qp)
    call require(abs(qmid(b)) < 1.0e-33_qp, "sqrt 2 to binary128 accuracy", nfail)

    d = to_interval(qrat(1, 3))
    call require(real(d%lo, qp) <= 1.0_qp/3.0_qp .and. &
        real(d%hi, qp) >= 1.0_qp/3.0_qp .and. d%hi - d%lo < 1.0e-15_dp, &
        "binary64 enclosure of 1/3", nfail)
    a = qinterval(1)/qinterval(-1.0_qp, 1.0_qp)
    call require(a%hi > huge(1.0_qp), "division by zero-containing interval", nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " interval_qp test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_interval_qp: all tests passed"

contains

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_interval_qp
