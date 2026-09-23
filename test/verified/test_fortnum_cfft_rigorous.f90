!> Rigorous FFT against real128 oracles: certified twiddles against
!> real128 cos/sin, the Higham bound against a direct real128 DFT, and the
!> convolution bounds against exact real128 convolutions, in 1D and 2D.
program test_fortnum_cfft_rigorous
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_rounding, only: unit_roundoff
    use fortnum_cfft_rigorous, only: rigorous_fft_plan_t, rigorous_fft_plan_init, &
        rigorous_fft_apply, rigorous_fft2_apply, conv_error_bound, &
        conv_error_bound_2d, &
        conv_nonneg, conv_nonneg_2d, next_pow2, fft_kappa
    implicit none

    integer :: nfail, p, n, k, j, na, nb, seed_size
    integer, allocatable :: seed(:)
    type(rigorous_fft_plan_t) :: plan, plan_m, plan_n
    type(fortnum_status_t) :: status
    real(qp) :: pi_q, errmax, ex(2)
    real(dp) :: e
    logical :: ok

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 777
    call random_seed(put=seed)
    pi_q = 4.0_qp*atan(1.0_qp)

    ! Certified twiddles.
    ok = .true.
    do p = 1, 14
        n = 2**p
        call rigorous_fft_plan_init(plan, n)
        errmax = 0.0_qp
        do k = 0, n/2 - 1
            ex(1) = cos(2.0_qp*pi_q*real(k, qp)/real(n, qp))
            ex(2) = -sin(2.0_qp*pi_q*real(k, qp)/real(n, qp))
            errmax = max(errmax, abs(cmplx(plan%w(k), kind=qp) - &
                cmplx(ex(1), ex(2), qp)))
        end do
        if (.not. (errmax <= real(plan%mu, qp))) ok = .false.
        if (plan%mu > unit_roundoff) ok = .false.
    end do
    call require(ok, "twiddle errors are bounded by the certified mu <= u", nfail)
    call rigorous_fft_plan_init(plan, 12, status)
    call require(status%code /= FORTNUM_OK, "non power of two is rejected", nfail)

    ! Forward and inverse transforms against a direct real128 DFT.
    ok = .true.
    do p = 1, 8
        n = 2**p
        call rigorous_fft_plan_init(plan, n)
        if (.not. dft_ok(plan, .false.)) ok = .false.
        if (.not. dft_ok(plan, .true.)) ok = .false.
    end do
    call require(ok, "FFT error within the Higham bound kappa norm2(F x)", nfail)
    call require(fft_kappa(2**20, 0.0_dp) < 1.0e-13_dp, "kappa is small", nfail)

    ! 1D convolution bound.
    ok = .true.
    do j = 1, 30
        na = 1 + mod(7*j, 40)
        nb = 1 + mod(11*j, 30)
        call rigorous_fft_plan_init(plan, next_pow2(2*(na + nb) + 2))
        if (.not. conv_ok(plan, na, nb)) ok = .false.
    end do
    call require(ok, "FFT convolution error within conv_error_bound", nfail)

    ok = .true.
    do j = 1, 30
        na = 1 + mod(5*j, 25)
        nb = 1 + mod(3*j, 20)
        call rigorous_fft_plan_init(plan, next_pow2(2*(na + nb) + 2))
        if (.not. nonneg_ok(plan, na, nb)) ok = .false.
    end do
    call require(ok, "conv_nonneg is an upper bound within twice the FFT margin", nfail)

    ! 2D.
    call rigorous_fft_plan_init(plan_m, 32)
    call rigorous_fft_plan_init(plan_n, 16)
    call require(conv2_ok(plan_m, plan_n, 6, 3, 7, 4), &
        "2D FFT convolution error within conv_error_bound_2d and conv_nonneg_2d &
        &bounds", nfail)
    e = conv_error_bound_2d(plan_m, plan_n, 1.0_dp, 1.0_dp)
    call require(e > conv_error_bound(plan_m, 1.0_dp, 1.0_dp), &
        "2D bound exceeds the 1D bound of one factor", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " cfft test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_cfft_rigorous: all tests passed"

contains

    logical function dft_ok(pl, inverse) result(good)
        type(rigorous_fft_plan_t), intent(in) :: pl
        logical, intent(in) :: inverse
        complex(dp) :: x(0:pl%n - 1)
        complex(qp) :: y(0:pl%n - 1)
        real(dp) :: re(pl%n), im(pl%n)
        real(qp) :: sgn, err, nrm
        integer :: jj, kk, m

        m = pl%n
        call random_number(re)
        call random_number(im)
        x = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
        sgn = merge(1.0_qp, -1.0_qp, inverse)
        do kk = 0, m - 1
            y(kk) = (0.0_qp, 0.0_qp)
            do jj = 0, m - 1
                y(kk) = y(kk) + cmplx(x(jj), kind=qp)*exp(cmplx(0.0_qp, &
                    sgn*2.0_qp*pi_q*real(mod(jj*kk, m), qp)/real(m, qp), qp))
            end do
        end do
        call rigorous_fft_apply(pl, x, inverse)
        if (inverse) then
            x = x*real(m, dp)
        end if
        err = sqrt(sum(abs(cmplx(x, kind=qp) - y)**2))
        nrm = sqrt(sum(abs(y)**2))
        good = err <= real(pl%kappa, qp)*nrm*(1.0_qp + 1.0e-20_qp)
        ! The rescaled inverse check is exact only up to the 1/n scaling,
        ! which is exact for a power of two.
    end function dft_ok

    logical function conv_ok(pl, ka, kb) result(good)
        type(rigorous_fft_plan_t), intent(in) :: pl
        integer, intent(in) :: ka, kb
        complex(dp) :: a(-ka:ka), b(-kb:kb), xa(0:pl%n - 1), xb(0:pl%n - 1)
        complex(qp) :: c
        real(dp) :: r1(2*ka + 1), r2(2*ka + 1), s1(2*kb + 1), s2(2*kb + 1)
        real(dp) :: l1a, l1b, bound
        integer :: i, m, kk

        m = pl%n
        call random_number(r1)
        call random_number(r2)
        call random_number(s1)
        call random_number(s2)
        a = cmplx(r1 - 0.5_dp, r2 - 0.5_dp, dp)
        b = cmplx(s1 - 0.5_dp, s2 - 0.5_dp, dp)*1.0e3_dp
        xa = (0.0_dp, 0.0_dp)
        xb = (0.0_dp, 0.0_dp)
        do i = -ka, ka
            xa(modulo(i, m)) = a(i)
        end do
        do i = -kb, kb
            xb(modulo(i, m)) = b(i)
        end do
        call rigorous_fft_apply(pl, xa, .false.)
        call rigorous_fft_apply(pl, xb, .false.)
        xa = xa*xb
        call rigorous_fft_apply(pl, xa, .true.)
        l1a = sum(abs(a))*(1.0_dp + 1.0e-10_dp)
        l1b = sum(abs(b))*(1.0_dp + 1.0e-10_dp)
        bound = conv_error_bound(pl, l1a, l1b)
        good = .true.
        do kk = -(ka + kb), ka + kb
            c = (0.0_qp, 0.0_qp)
            do i = max(-ka, kk - kb), min(ka, kk + kb)
                c = c + cmplx(a(i), kind=qp)*cmplx(b(kk - i), kind=qp)
            end do
            if (abs(cmplx(xa(modulo(kk, m)), kind=qp) - c) > real(bound, qp)) &
                good = .false.
        end do
    end function conv_ok

    logical function nonneg_ok(pl, ka, kb) result(good)
        type(rigorous_fft_plan_t), intent(in) :: pl
        integer, intent(in) :: ka, kb
        real(dp) :: a(-ka:ka), b(-kb:kb), cc(0:pl%n - 1)
        real(qp) :: c, margin
        integer :: i, kk

        call random_number(a)
        call random_number(b)
        call conv_nonneg(pl, a, ka, b, kb, cc)
        margin = real(conv_error_bound(pl, sum(a)*1.000001_dp, sum(b)*1.000001_dp), qp)
        good = .true.
        do kk = -(ka + kb), ka + kb
            c = 0.0_qp
            do i = max(-ka, kk - kb), min(ka, kk + kb)
                c = c + real(a(i), qp)*real(b(kk - i), qp)
            end do
            if (real(cc(modulo(kk, pl%n)), qp) < c) good = .false.
            if (real(cc(modulo(kk, pl%n)), qp) > c + 2*margin + 1.0e-15_qp*c) &
                good = .false.
        end do
    end function nonneg_ok

    logical function conv2_ok(pm, pn, ma, na2, mb, nb2) result(good)
        type(rigorous_fft_plan_t), intent(in) :: pm, pn
        integer, intent(in) :: ma, na2, mb, nb2
        real(dp) :: a(-ma:ma, -na2:na2), b(-mb:mb, -nb2:nb2)
        real(dp) :: cc(0:pm%n - 1, 0:pn%n - 1)
        complex(dp) :: xa(0:pm%n - 1, 0:pn%n - 1), xb(0:pm%n - 1, 0:pn%n - 1)
        real(qp) :: c
        real(dp) :: bound
        integer :: i1, i2, k1, k2

        call random_number(a)
        call random_number(b)
        a = a - 0.3_dp
        xa = (0.0_dp, 0.0_dp)
        xb = (0.0_dp, 0.0_dp)
        do i2 = -na2, na2
            do i1 = -ma, ma
                xa(modulo(i1, pm%n), modulo(i2, pn%n)) = a(i1, i2)
            end do
        end do
        do i2 = -nb2, nb2
            do i1 = -mb, mb
                xb(modulo(i1, pm%n), modulo(i2, pn%n)) = b(i1, i2)
            end do
        end do
        call rigorous_fft2_apply(pm, pn, xa, .false.)
        call rigorous_fft2_apply(pm, pn, xb, .false.)
        xa = xa*xb
        call rigorous_fft2_apply(pm, pn, xa, .true.)
        bound = conv_error_bound_2d(pm, pn, sum(abs(a))*1.000001_dp, &
            sum(abs(b))*1.000001_dp)
        call conv_nonneg_2d(pm, pn, abs(a), ma, na2, b, mb, nb2, cc)
        good = .true.
        do k2 = -(na2 + nb2), na2 + nb2
            do k1 = -(ma + mb), ma + mb
                c = 0.0_qp
                do i2 = max(-na2, k2 - nb2), min(na2, k2 + nb2)
                    do i1 = max(-ma, k1 - mb), min(ma, k1 + mb)
                        c = c + real(a(i1, i2), qp)*real(b(k1 - i1, k2 - i2), qp)
                    end do
                end do
                if (abs(cmplx(xa(modulo(k1, pm%n), modulo(k2, pn%n)), kind=qp) &
                    - c) > real(bound, qp)) good = .false.
                c = 0.0_qp
                do i2 = max(-na2, k2 - nb2), min(na2, k2 + nb2)
                    do i1 = max(-ma, k1 - mb), min(ma, k1 + mb)
                        c = c + abs(real(a(i1, i2), qp))*real(b(k1 - i1, k2 - i2), qp)
                    end do
                end do
                if (real(cc(modulo(k1, pm%n), modulo(k2, pn%n)), qp) < c) &
                    good = .false.
            end do
        end do
    end function conv2_ok

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_cfft_rigorous
