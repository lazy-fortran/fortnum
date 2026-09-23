!> Verified linear algebra against exact and real128 oracles: eigenvalue
!> and norm bounds of the tridiagonal Toeplitz matrix tridiag(-1, 2, -1)
!> (closed-form eigenvalues 2 - 2 cos(k pi/(n + 1))) and of an exactly
!> unitarily similar Hermitian matrix, product error bounds against real128
!> products, and verified inverses against a real128 Gauss-Jordan inverse.
program test_fortnum_verified_linalg
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, interval
    use fortnum_verified_linalg, only: mball_t, fro_up, norm2_up, norm2_tight_up, &
        matmul_err, verified_inverse, sym_lower_bound, eig_lower_bound, &
        eig_upper_bound, interval_matvec
    implicit none

    integer, parameter :: nt = 40
    real(dp) :: hh(32, 32), dd(32, 32), ad(32, 32)
    real(dp) :: t(nt, nt), a(7, 5), b(5, 6), p(7, 6), err, lo, hi, s, rn
    complex(dp) :: h(nt, nt), z(8, 8), zp(8, 8), zc(6, 4), zd(4, 3), pc(6, 3)
    complex(qp) :: zq(8, 8), zi(8, 8)
    real(qp) :: lmin, lmax, pi_q, fro_err
    real(dp) :: re(8, 8), im(8, 8), x(8), y(8)
    type(mball_t) :: g
    type(interval_t) :: ai(3, 3), xi(3), yi(3)
    logical :: ok, allok
    integer :: nfail, i, j, k, seed_size
    integer, allocatable :: seed(:)

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 2718
    call random_seed(put=seed)
    pi_q = 4.0_qp*atan(1.0_qp)

    t = 0.0_dp
    do i = 1, nt
        t(i, i) = 2.0_dp
        if (i > 1) t(i, i - 1) = -1.0_dp
        if (i < nt) t(i, i + 1) = -1.0_dp
    end do
    lmin = 2.0_qp - 2.0_qp*cos(pi_q/real(nt + 1, qp))
    lmax = 2.0_qp + 2.0_qp*cos(pi_q/real(nt + 1, qp))

    call eig_lower_bound(t, lo, ok)
    call require(ok .and. real(lo, qp) <= lmin .and. &
        real(lo, qp) >= lmin - 1.0e-12_qp, &
        "lambda_min lower bound of tridiag(-1,2,-1) is rigorous and tight", nfail)
    call eig_upper_bound(t, hi, ok)
    call require(ok .and. real(hi, qp) >= lmax .and. &
        real(hi, qp) <= lmax + 1.0e-12_qp, &
        "lambda_max upper bound", nfail)
    call eig_lower_bound(t, lo, ok, rho=1.0e-3_dp)
    call require(real(lo, qp) <= lmin - 1.0e-3_qp, "matrix ball lowers the bound", &
        nfail)
    s = norm2_tight_up(t)
    call require(real(s, qp) >= lmax .and. real(s, qp) <= lmax*(1 + 1.0e-10_qp), &
        "tight spectral norm bound", nfail)
    call require(real(norm2_up(t), qp) >= lmax .and. real(fro_up(t), qp) >= lmax, &
        "cheap spectral norm bounds", nfail)

    ! A shift above lambda_min must never produce a wrong certificate.
    allok = .true.
    do k = -20, 20
        call sym_lower_bound(t, real(lmin, dp) + real(k, dp)*1.0e-4_dp, lo, ok)
        if (ok .and. real(lo, qp) > lmin) allok = .false.
    end do
    ! Shifts within the rounding noise of the factorization, where a
    ! completed Cholesky alone proves nothing and the gamma_(n+1) term must
    ! carry the certificate.
    do k = -200, 200
        call sym_lower_bound(t, real(lmin, dp) + real(k, dp)*1.0e-16_dp, lo, ok)
        if (ok .and. real(lo, qp) > lmin) allok = .false.
    end do
    call require(allok, "Cholesky certificates never exceed lambda_min", nfail)

    ! Dense A = H D H with the exact Householder reflection H = I - (1/16) 1 1^T
    ! (n = 32): every entry of A is exact in binary64 and the eigenvalues are
    ! exactly diag(D) = 1, 2, ..., 32 scaled by 30, with lambda_min = 1.
    hh = -1.0_dp/16.0_dp
    dd = 0.0_dp
    do i = 1, 32
        hh(i, i) = hh(i, i) + 1.0_dp
        dd(i, i) = real(30*i, dp)
    end do
    dd(1, 1) = 1.0_dp
    ad = matmul(hh, matmul(dd, hh))
    call eig_lower_bound(ad, lo, ok)
    call require(ok .and. lo <= 1.0_dp .and. lo >= 1.0_dp - 1.0e-9_dp, &
        "dense exact-spectrum lambda_min lower bound", nfail)
    allok = .true.
    do k = -300, 300
        call sym_lower_bound(ad, 1.0_dp + real(k, dp)*1.0e-13_dp, lo, ok)
        if (ok .and. lo > 1.0_dp) allok = .false.
    end do
    call require(allok, "dense certificates never exceed lambda_min", nfail)

    ! Hermitian H = D T D^H, D = diag(i^k): exact entries, same eigenvalues.
    do j = 1, nt
        do i = 1, nt
            h(i, j) = t(i, j)*(0.0_dp, 1.0_dp)**modulo(i - j, 4)
        end do
    end do
    call eig_lower_bound(h, lo, ok)
    call require(ok .and. real(lo, qp) <= lmin .and. &
        real(lo, qp) >= lmin - 1.0e-11_qp, &
        "Hermitian lambda_min lower bound", nfail)
    call eig_upper_bound(h, hi, ok)
    call require(ok .and. real(hi, qp) >= lmax .and. &
        real(hi, qp) <= lmax + 1.0e-11_qp, &
        "Hermitian lambda_max upper bound", nfail)
    s = norm2_tight_up(h)
    call require(real(s, qp) >= lmax .and. real(s, qp) <= lmax*(1 + 1.0e-10_qp), &
        "complex tight spectral norm bound", nfail)

    ! Product error bounds (Frobenius of the real128 error bounds its norm2).
    call random_number(a)
    call random_number(b)
    a = a - 0.5_dp
    call matmul_err(a, b, p, err)
    fro_err = sqrt(sum((real(p, qp) - matmul(real(a, qp), real(b, qp)))**2))
    call require(fro_err <= real(err, qp), "real product error bound", nfail)
    call random_number(re(1:6, 1:4))
    call random_number(im(1:6, 1:4))
    zc = cmplx(re(1:6, 1:4), im(1:6, 1:4) - 0.5_dp, dp)
    call random_number(re(1:4, 1:3))
    call random_number(im(1:4, 1:3))
    zd = cmplx(re(1:4, 1:3) - 0.5_dp, im(1:4, 1:3), dp)
    call matmul_err(zc, zd, pc, err)
    fro_err = sqrt(sum(abs(cmplx(pc, kind=qp) - matmul(cmplx(zc, kind=qp), &
        cmplx(zd, kind=qp)))**2))
    call require(fro_err <= real(err, qp), "complex product error bound", nfail)

    ! Verified inverse of a well-conditioned complex matrix, with and
    ! without a perturbation of norm eta.
    call random_number(re)
    call random_number(im)
    z = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
    do i = 1, 8
        z(i, i) = z(i, i) + (4.0_dp, 1.0_dp)
    end do
    call verified_inverse(z, 0.0_dp, 1.0e-3_dp, g, ok)
    call require(ok, "verified inverse certifies", nfail)
    zq = cmplx(z, kind=qp)
    call qinverse(zq, zi)
    fro_err = sqrt(sum(abs(zi - cmplx(g%c, kind=qp))**2))
    call require(fro_err <= real(g%rho, qp), "exact inverse lies in the ball", nfail)
    call random_number(x)
    call random_number(y)
    x = x/norm2(x)
    y = y/norm2(y)
    do j = 1, 8
        do i = 1, 8
            zp(i, j) = z(i, j) + 0.999e-3_dp*x(i)*y(j)
        end do
    end do
    call qinverse(cmplx(zp, kind=qp), zi)
    fro_err = sqrt(sum(abs(zi - cmplx(g%c, kind=qp))**2))
    call require(fro_err <= real(g%rho, qp), "perturbed inverse lies in the ball", &
        nfail)
    rn = g%rho
    call verified_inverse(z, 0.0_dp, 10.0_dp, g, ok)
    call require(.not. ok, "an overlarge perturbation is refused", nfail)
    call require(rn < 1.0e-2_dp, "inverse ball radius is small", nfail)

    ! Interval matrix-vector product.
    do j = 1, 3
        do i = 1, 3
            ai(i, j) = interval(real(i - j, dp), real(i - j, dp) + 0.1_dp)
        end do
        xi(j) = interval(1.0_dp, 1.5_dp)
    end do
    yi = interval_matvec(ai, xi)
    call require(yi(3)%lo <= 3.0_dp .and. yi(3)%hi >= 3.0_dp*1.5_dp + 0.45_dp, &
        "interval matrix-vector product encloses its corners", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " verified linalg test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_verified_linalg: all tests passed"

contains

    !> Gauss-Jordan inverse in real128 (the oracle).
    subroutine qinverse(a, ainv)
        complex(qp), intent(in) :: a(:, :)
        complex(qp), intent(out) :: ainv(size(a, 1), size(a, 1))
        complex(qp) :: w(size(a, 1), 2*size(a, 1)), f
        integer :: n, kk, ii

        n = size(a, 1)
        w = (0.0_qp, 0.0_qp)
        w(:, 1:n) = a
        do ii = 1, n
            w(ii, n + ii) = (1.0_qp, 0.0_qp)
        end do
        do kk = 1, n
            w(kk, :) = w(kk, :)/w(kk, kk)
            do ii = 1, n
                if (ii == kk) cycle
                f = w(ii, kk)
                w(ii, :) = w(ii, :) - f*w(kk, :)
            end do
        end do
        ainv = w(:, n + 1:2*n)
    end subroutine qinverse

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_verified_linalg
