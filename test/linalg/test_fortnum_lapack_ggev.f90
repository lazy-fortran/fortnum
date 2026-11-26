program test_fortnum_lapack_ggev
    ! The complex generalized eigenproblem by its defining equations in quad
    ! precision: beta A v = alpha B v and u^H (beta A - alpha B) = 0 with
    ! LAPACK's normalization, spectra known by construction including
    ! infinite eigenvalues of a singular B, triangular pencils, extreme
    ! magnitudes and argument errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_ggev, only: fortnum_zggev
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: eps = epsilon(1.0_dp)

    call general(1, 0.3_dp, 1.0_dp, 1.0_dp)
    call general(2, 0.9_dp, 1.0_dp, 1.0_dp)
    call general(7, 1.7_dp, 1.0_dp, 1.0_dp)
    call general(40, 0.29_dp, 1.0_dp, 1.0_dp)
    call general(80, 0.61_dp, 1.0_dp, 1.0_dp)
    call general(20, 0.53_dp, 1.0e-300_dp, 1.0e300_dp)
    call general(20, 0.53_dp, 1.0e300_dp, 1.0e-300_dp)
    call known_spectrum(16, .false.)
    call known_spectrum(16, .true.)
    call triangular_pencil(9, .false.)
    call triangular_pencil(9, .true.)
    call zero_column(12)
    call skew_definite(2, 0.3_dp)
    call skew_definite(6, 0.9_dp)
    call skew_definite(10, 1.7_dp)
    call argument_errors()
    print '(a)', 'generalized eigenproblem matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function cgeneral(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: a(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                a(i, j) = cmplx(sin(seed*i + 1.91_dp*j + 0.37_dp*i*j), &
                    cos(0.7_dp*seed*i - 0.53_dp*j + 0.11_dp*i*j), dp)
            end do
        end do
    end function cgeneral

    ! unitary I - 2 v v^H / v^H v
    function householder(n, seed) result(u)
        integer, intent(in) :: n
        real(qp), intent(in) :: seed
        complex(qp) :: u(n, n)
        complex(qp) :: v(n)
        integer :: i, j
        do i = 1, n
            v(i) = cmplx(cos(seed*i), sin(1.3_qp*i + seed), qp)
        end do
        do j = 1, n
            do i = 1, n
                u(i, j) = -2*v(i)*conjg(v(j))/sum(abs(v)**2)
            end do
            u(j, j) = u(j, j) + 1
        end do
    end function householder

    ! residuals of all pairs and LAPACK's normalization
    subroutine check_pairs(n, a, b, alpha, beta, vl, vr, tol, what)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a(n, n), b(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        real(dp), intent(in) :: tol
        character(*), intent(in) :: what
        complex(qp) :: aq(n, n), bq(n, n), x(n), al, be
        real(qp) :: anorm, bnorm, res
        integer :: j
        aq = a
        bq = b
        anorm = sqrt(sum(abs(aq)**2))
        bnorm = sqrt(sum(abs(bq)**2))
        do j = 1, n
            al = alpha(j)
            be = beta(j)
            call require(aimag(beta(j)) == 0 .and. real(beta(j), dp) >= 0, what//': beta real')
            x = vr(:, j)
            res = sqrt(sum(abs(be*matmul(aq, x) - al*matmul(bq, x))**2))
            call require(res <= tol*(abs(be)*anorm + abs(al)*bnorm)*sqrt(sum(abs(x)**2)), &
                what//': right residual')
            call require(abs(maxval(abs(real(x, qp)) + abs(aimag(x))) - 1) <= 4*eps, &
                what//': right normalized')
            x = vl(:, j)
            res = sqrt(sum(abs(conjg(be)*matmul(conjg(transpose(aq)), x) &
                - conjg(al)*matmul(conjg(transpose(bq)), x))**2))
            call require(res <= tol*(abs(be)*anorm + abs(al)*bnorm)*sqrt(sum(abs(x)**2)), &
                what//': left residual')
            call require(abs(maxval(abs(real(x, qp)) + abs(aimag(x))) - 1) <= 4*eps, &
                what//': left normalized')
        end do
    end subroutine check_pairs

    subroutine general(n, seed, sa, sb)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed, sa, sb
        complex(dp) :: a(n, n), b(n, n), aw(n, n), bw(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        complex(dp) :: alpha2(n), beta2(n), work(2*n), dum(1, 1)
        real(dp) :: rwork(8*n)
        integer :: info
        a = cgeneral(n, seed)*sa
        b = cgeneral(n, seed + 0.77_dp)*sb
        aw = a
        bw = b
        call fortnum_zggev('V', 'V', n, aw, n, bw, n, alpha, beta, vl, n, vr, n, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev general: info')
        call check_pairs(n, a, b, alpha, beta, vl, vr, 50*n*eps, 'zggev general')
        ! eigenvalues only: the same pairs
        aw = a
        bw = b
        call fortnum_zggev('N', 'N', n, aw, n, bw, n, alpha2, beta2, dum, 1, dum, 1, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev eigenvalues only: info')
        call require(all(abs(alpha2*beta - alpha*beta2) <= 1.0e-10_dp*abs(alpha*beta2)), &
            'zggev eigenvalues only: same eigenvalues')
    end subroutine general

    ! A = U diag(d) V, B = U diag(e) V with U, V unitary: lambda = d/e,
    ! infinite where e = 0 (singular B)
    subroutine known_spectrum(n, singular)
        integer, intent(in) :: n
        logical, intent(in) :: singular
        complex(dp) :: a(n, n), b(n, n), aw(n, n), bw(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        complex(dp) :: work(2*n)
        complex(qp) :: u(n, n), v(n, n), d(n, n), e(n, n), lam, exact(n)
        real(dp) :: rwork(8*n)
        logical :: finite(n)
        integer :: info, i, ninf
        u = householder(n, 0.7_qp)
        v = householder(n, 1.9_qp)
        d = 0
        e = 0
        do i = 1, n
            d(i, i) = cmplx(0.5_qp*i - 3, 0.2_qp*(-1)**i*i, qp)
            e(i, i) = cmplx(1 + 0.1_qp*i, 0.3_qp*sin(1.0_qp*i), qp)
            ! a strictly upper part keeps the problem nonnormal
            if (i > 1) d(i - 1, i) = 0.3_qp
        end do
        if (singular) then
            e(3, 3) = 0
            e(7, 7) = 0
        end if
        a = cmplx(matmul(matmul(u, d), v), kind=dp)
        b = cmplx(matmul(matmul(u, e), v), kind=dp)
        do i = 1, n
            finite(i) = abs(e(i, i)) > 0
            exact(i) = d(i, i)
            if (finite(i)) exact(i) = d(i, i)/e(i, i)
        end do
        aw = a
        bw = b
        call fortnum_zggev('V', 'V', n, aw, n, bw, n, alpha, beta, vl, n, vr, n, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev known spectrum: info')
        call check_pairs(n, a, b, alpha, beta, vl, vr, 50*n*eps, 'zggev known spectrum')
        ninf = 0
        do i = 1, n
            if (abs(beta(i)) <= 1.0e-12_dp*abs(alpha(i))) then
                ninf = ninf + 1
                cycle
            end if
            lam = cmplx(alpha(i)/beta(i), kind=qp)
            call require(minval(abs(exact - lam), mask=finite) <= 1.0e-10_qp*abs(lam), &
                'zggev known spectrum: eigenvalue')
        end do
        if (singular) then
            call require(ninf == 2, 'zggev singular B: two infinite eigenvalues')
        else
            call require(ninf == 0, 'zggev regular B: no infinite eigenvalue')
        end if
    end subroutine known_spectrum

    ! (A, B) already upper triangular: lambda_i = a(i,i)/b(i,i). Variant
    ! special: b(1,1) = b(5,5) = 0 give infinite eigenvalues whose zero is
    ! moved down from the top of a block, and a repeated eigenvalue
    ! a(3,3)/b(3,3) = a(7,7)/b(7,7) needs the perturbed divisor.
    subroutine triangular_pencil(n, special)
        integer, intent(in) :: n
        logical, intent(in) :: special
        complex(dp) :: a(n, n), b(n, n), aw(n, n), bw(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        complex(dp) :: work(2*n), lam
        real(dp) :: rwork(8*n)
        logical :: finite(n), infinite(n)
        integer :: info, i, j
        a = 0
        b = 0
        do j = 1, n
            do i = 1, j
                a(i, j) = cmplx(cos(1.0_dp*i + 0.5_dp*j), sin(0.3_dp*i*j), dp)
                b(i, j) = cmplx(0.2_dp*sin(1.0_dp*i - j), 0.1_dp*cos(2.0_dp*j), dp)
            end do
            b(j, j) = cmplx(1 + 0.1_dp*j, -0.2_dp*j, dp)
        end do
        finite = .true.
        if (special) then
            b(1, 1) = 0
            b(5, 5) = 0
            finite([1, 5]) = .false.
            a(7, 7) = a(3, 3)
            b(7, 7) = b(3, 3)
        end if
        aw = a
        bw = b
        call fortnum_zggev('V', 'V', n, aw, n, bw, n, alpha, beta, vl, n, vr, n, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev triangular: info')
        call check_pairs(n, a, b, alpha, beta, vl, vr, 50*n*eps, 'zggev triangular')
        infinite = abs(beta) <= 1.0e-13_dp*abs(alpha)
        call require(count(infinite) == count(.not. finite), 'zggev triangular: infinite eigenvalues')
        do i = 1, n
            if (.not. finite(i)) cycle
            lam = a(i, i)/b(i, i)
            call require(minval(abs(alpha/beta - lam), mask=.not. infinite) <= 1.0e-7_dp*abs(lam), &
                'zggev triangular: eigenvalue')
        end do
    end subroutine triangular_pencil

    ! B e1 = 0: one infinite eigenvalue, whose zero on the diagonal of the
    ! triangular factor sits at the top of the Hessenberg block
    subroutine zero_column(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), b(n, n), aw(n, n), bw(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        complex(dp) :: work(2*n)
        real(dp) :: rwork(8*n)
        integer :: info
        a = cgeneral(n, 0.47_dp)
        b = cgeneral(n, 1.13_dp)
        b(:, 1) = 0
        aw = a
        bw = b
        call fortnum_zggev('V', 'V', n, aw, n, bw, n, alpha, beta, vl, n, vr, n, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev zero column: info')
        call check_pairs(n, a, b, alpha, beta, vl, vr, 50*n*eps, 'zggev zero column')
        call require(count(abs(beta) <= 1.0e-13_dp*abs(alpha)) == 1, &
            'zggev zero column: one infinite eigenvalue')
    end subroutine zero_column

    ! real skew-symmetric A and symmetric positive definite B: B^-1 A has a
    ! purely imaginary spectrum, which a real Wilkinson shift cannot reach
    subroutine skew_definite(n, seed)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: a(n, n), b(n, n), aw(n, n), bw(n, n), alpha(n), beta(n), vl(n, n), vr(n, n)
        complex(dp) :: work(2*n)
        real(dp) :: g(n, n), rwork(8*n)
        integer :: info, i
        g = real(cgeneral(n, seed), dp)
        a = g - transpose(g)
        g = aimag(cgeneral(n, seed + 0.77_dp))
        b = matmul(g, transpose(g))
        do i = 1, n
            b(i, i) = b(i, i) + n
        end do
        aw = a
        bw = b
        call fortnum_zggev('V', 'V', n, aw, n, bw, n, alpha, beta, vl, n, vr, n, work, &
            size(work), rwork, info)
        call require(info == 0, 'zggev skew-definite: info')
        call check_pairs(n, a, b, alpha, beta, vl, vr, 50*n*eps, 'zggev skew-definite')
        call require(all(abs(real(alpha, dp)) <= 1.0e-10_dp*abs(alpha) + 50*n*eps*abs(beta)), &
            'zggev skew-definite: imaginary eigenvalues')
    end subroutine skew_definite

    subroutine argument_errors()
        complex(dp) :: a(2, 2), b(2, 2), alpha(2), beta(2), v(2, 2), work(4)
        real(dp) :: rwork(16)
        integer :: info
        a = 0
        b = 0
        call fortnum_zggev('X', 'N', 2, a, 2, b, 2, alpha, beta, v, 2, v, 2, work, 4, rwork, info)
        call require(info == -1, 'zggev: jobvl')
        call fortnum_zggev('N', 'X', 2, a, 2, b, 2, alpha, beta, v, 2, v, 2, work, 4, rwork, info)
        call require(info == -2, 'zggev: jobvr')
        call fortnum_zggev('N', 'N', -1, a, 2, b, 2, alpha, beta, v, 2, v, 2, work, 4, rwork, info)
        call require(info == -3, 'zggev: n')
        call fortnum_zggev('N', 'N', 2, a, 1, b, 2, alpha, beta, v, 2, v, 2, work, 4, rwork, info)
        call require(info == -5, 'zggev: lda')
        call fortnum_zggev('N', 'N', 2, a, 2, b, 1, alpha, beta, v, 2, v, 2, work, 4, rwork, info)
        call require(info == -7, 'zggev: ldb')
        call fortnum_zggev('V', 'N', 2, a, 2, b, 2, alpha, beta, v, 1, v, 2, work, 4, rwork, info)
        call require(info == -11, 'zggev: ldvl')
        call fortnum_zggev('N', 'V', 2, a, 2, b, 2, alpha, beta, v, 2, v, 1, work, 4, rwork, info)
        call require(info == -13, 'zggev: ldvr')
        call fortnum_zggev('N', 'N', 2, a, 2, b, 2, alpha, beta, v, 2, v, 2, work, 3, rwork, info)
        call require(info == -15, 'zggev: lwork')
        call fortnum_zggev('V', 'V', 2, a, 2, b, 2, alpha, beta, v, 2, v, 2, work, -1, rwork, info)
        call require(info == 0 .and. int(real(work(1), dp)) >= 4, 'zggev: workspace query')
    end subroutine argument_errors

end program test_fortnum_lapack_ggev
