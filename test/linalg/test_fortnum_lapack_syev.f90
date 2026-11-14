program test_fortnum_lapack_syev
    ! Symmetric and Hermitian eigensolvers by their defining equations in quad
    ! precision: A V = V diag(w) with orthonormal V and ascending w, the
    ! generalized problems of every type, tridiagonal and band drivers.
    use fortnum_kinds, only: dp
    use fortnum_lapack_syev
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: eps = real(epsilon(1.0_dp), qp)

    call real_case(1, 'L', 1.0_dp)
    call real_case(5, 'U', 1.0_dp)
    call real_case(80, 'L', 1.0_dp)
    call real_case(80, 'U', 1.0e-290_dp)
    call real_case(40, 'L', 1.0e290_dp)
    call complex_case(60, 'L')
    call complex_case(60, 'U')
    call generalized_real(40)
    call generalized_complex(30)
    call tridiagonal(50)
    call band(60, 3)
    print '(a)', 'symmetric and Hermitian eigensolvers match quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function sym(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        real(dp) :: a(n, n)
        integer :: i, j
        do j = 1, n
            do i = j, n
                a(i, j) = sin(seed*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j))
                a(j, i) = a(i, j)
            end do
        end do
    end function sym

    function identity(n) result(e)
        integer, intent(in) :: n
        real(qp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function identity

    subroutine check_eigen(a, v, w, what)
        real(dp), intent(in) :: a(:, :), v(:, :), w(:)
        character(*), intent(in) :: what
        integer :: n, j
        real(qp) :: anorm
        n = size(a, 1)
        anorm = maxval(sum(abs(real(a, qp)), 1))
        call require(maxval(abs(matmul(transpose(real(v, qp)), real(v, qp)) - identity(n))) <= 32*n*eps, &
            what//': V orthonormal')
        do j = 1, n
            call require(maxval(abs(matmul(real(a, qp), real(v(:, j), qp)) - w(j)*real(v(:, j), qp))) <= &
                32*n*eps*anorm, what//': A v = w v')
            if (j > 1) call require(w(j) >= w(j - 1), what//': ascending')
        end do
    end subroutine check_eigen

    subroutine real_case(n, uplo, scale)
        integer, intent(in) :: n
        character, intent(in) :: uplo
        real(dp), intent(in) :: scale
        real(dp) :: a(n, n), v(n + 1, n), w(n), w2(n), work(max(1, 3*n)), wq(1)
        integer :: info, iwork(3 + 5*n)
        a = scale*sym(n, 0.37_dp)
        v = 9
        v(1:n, :) = a
        call fortnum_dsyev('V', uplo, n, v, n + 1, w, wq, -1, info)
        call require(info == 0 .and. wq(1) >= 3*n - 1, 'dsyev workspace query')
        call fortnum_dsyev('V', uplo, n, v, n + 1, w, work, size(work), info)
        call require(info == 0, 'dsyev info')
        call check_eigen(a/scale, v(1:n, :), w/scale, 'dsyev')
        v(1:n, :) = a
        call fortnum_dsyev('N', uplo, n, v, n + 1, w2, work, size(work), info)
        call require(maxval(abs(w2 - w)) <= 1e-12_dp*maxval(abs(w)), 'values only agree')
        v(1:n, :) = a
        block
            real(dp) :: dwork(1 + 6*n + 2*n**2)
            call fortnum_dsyevd('V', uplo, n, v, n + 1, w2, dwork, size(dwork), iwork, size(iwork), info)
        end block
        call require(info == 0 .and. maxval(abs(w2 - w)) <= 1e-12_dp*maxval(abs(w)), 'dsyevd')
    end subroutine real_case

    subroutine complex_case(n, uplo)
        integer, intent(in) :: n
        character, intent(in) :: uplo
        complex(dp) :: a(n, n), v(n, n), work(2*n)
        real(dp) :: w(n), rwork(3*n), re(n, n)
        complex(qp) :: av(n), vq(n, n)
        integer :: info, i, j
        re = sym(n, 0.51_dp)
        do j = 1, n
            do i = 1, n
                a(i, j) = cmplx(re(i, j), 0.3_dp*sin(real(i - j, dp)*0.7_dp), dp)
            end do
        end do
        v = a
        call fortnum_zheev('V', uplo, n, v, n, w, work, size(work), rwork, info)
        call require(info == 0, 'zheev info')
        vq = cmplx(v, kind=qp)
        call require(maxval(abs(matmul(conjg(transpose(vq)), vq) - identity(n))) <= 64*n*eps, 'zheev V unitary')
        do j = 1, n
            av = matmul(cmplx(a, kind=qp), vq(:, j))
            call require(maxval(abs(av - w(j)*vq(:, j))) <= 64*n*eps*maxval(sum(abs(cmplx(a, kind=qp)), 1)), &
                'zheev A v = w v')
        end do
    end subroutine complex_case

    ! A x = w B x (1), A B x = w x (2), B A x = w x (3)
    subroutine generalized_real(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), b(n, n), v(n, n), bf(n, n), w(n), work(3*n)
        real(qp) :: lhs(n), rhs(n), scale
        integer :: info, itype, j, k
        character :: uplo
        a = sym(n, 0.29_dp)
        b = matmul(sym(n, 0.83_dp), sym(n, 0.83_dp))
        do j = 1, n
            b(j, j) = b(j, j) + n
        end do
        do k = 1, 6
            itype = 1 + mod(k - 1, 3)
            uplo = merge('U', 'L', k > 3)
            v = a
            bf = b
            call fortnum_dsygv(itype, 'V', uplo, n, v, n, bf, n, w, work, size(work), info)
            call require(info == 0, 'dsygv info')
            scale = maxval(abs(real(a, qp)))*maxval(abs(real(b, qp)))*maxval(abs(real(v, qp)))*n
            do j = 1, n
                select case (itype)
                case (1)
                    lhs = matmul(real(a, qp), real(v(:, j), qp))
                    rhs = w(j)*matmul(real(b, qp), real(v(:, j), qp))
                case (2)
                    lhs = matmul(real(a, qp), matmul(real(b, qp), real(v(:, j), qp)))
                    rhs = w(j)*real(v(:, j), qp)
                case (3)
                    lhs = matmul(real(b, qp), matmul(real(a, qp), real(v(:, j), qp)))
                    rhs = w(j)*real(v(:, j), qp)
                end select
                call require(maxval(abs(lhs - rhs)) <= 1e-11_qp*scale*max(1.0_qp, abs(real(w(j), qp))), &
                    'dsygv residual')
            end do
        end do
        bf = -b
        v = a
        call fortnum_dsygv(1, 'N', 'L', n, v, n, bf, n, w, work, size(work), info)
        call require(info == n + 1, 'dsygv reports B not positive definite')
    end subroutine generalized_real

    subroutine generalized_complex(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), b(n, n), v(n, n), bf(n, n), work(2*n), m(n, n)
        real(dp) :: w(n), rwork(3*n), re(n, n), im(n, n)
        complex(qp) :: lhs(n), rhs(n)
        integer :: info, j
        re = sym(n, 0.41_dp)
        im = sym(n, 0.77_dp)
        a = cmplx(re, 0.0_dp, dp)
        do j = 1, n
            a(j + 1:n, j) = a(j + 1:n, j) + cmplx(0.0_dp, 0.2_dp, dp)
            a(j, j + 1:n) = conjg(a(j + 1:n, j))
        end do
        m = cmplx(re, im, dp)
        b = matmul(m, conjg(transpose(m)))
        do j = 1, n
            b(j, j) = b(j, j) + n
        end do
        v = a
        bf = b
        call fortnum_zhegv(1, 'V', 'L', n, v, n, bf, n, w, work, size(work), rwork, info)
        call require(info == 0, 'zhegv info')
        do j = 1, n
            lhs = matmul(cmplx(a, kind=qp), cmplx(v(:, j), kind=qp))
            rhs = w(j)*matmul(cmplx(b, kind=qp), cmplx(v(:, j), kind=qp))
            call require(maxval(abs(lhs - rhs)) <= 1e-10_qp*maxval(abs(rhs)) + 1e-10_qp, 'zhegv residual')
        end do
    end subroutine generalized_complex

    subroutine tridiagonal(n)
        integer, intent(in) :: n
        real(dp) :: d(n), e(n), d0(n), e0(n), z(n, n), t(n, n), work(1)
        integer :: info, i
        do i = 1, n
            d0(i) = cos(real(i, dp))
            e0(i) = 0.5_dp*sin(1.7_dp*i)
        end do
        t = 0
        do i = 1, n
            t(i, i) = d0(i)
            if (i < n) then
                t(i + 1, i) = e0(i)
                t(i, i + 1) = e0(i)
            end if
        end do
        d = d0
        e = e0
        call fortnum_dstev('V', n, d, e, z, n, work, info)
        call require(info == 0, 'dstev info')
        call check_eigen(t, z, d, 'dstev')
        d = d0
        e = e0
        call fortnum_dsterf(n, d, e, info)
        call require(info == 0, 'dsterf info')
        block
            real(dp) :: zz(n, n), wv(n)
            zz = 0
            wv = 0
            d0 = d
            d = [(cos(real(i, dp)), i = 1, n)]
            e(1:n - 1) = [(0.5_dp*sin(1.7_dp*i), i = 1, n - 1)]
            call fortnum_dstev('N', n, d, e, zz, n, work, info)
            call require(maxval(abs(d - d0)) <= 1e-13_dp, 'dsterf and dstev values agree')
        end block
    end subroutine tridiagonal

    subroutine band(n, kd)
        integer, intent(in) :: n, kd
        real(dp) :: a(n, n), ab(kd + 1, n), w(n), z(n, n), q(n, n), wx(n), zx(n, n)
        real(dp) :: work(1 + 5*n + 2*n**2)
        integer :: iwork(3 + 5*n), ifail(n), info, i, j, m
        a = sym(n, 0.63_dp)
        do j = 1, n
            do i = 1, n
                if (abs(i - j) > kd) a(i, j) = 0
            end do
        end do
        do j = 1, n
            do i = j, min(n, j + kd)
                ab(1 + i - j, j) = a(i, j)
            end do
        end do
        call fortnum_dsbevd('V', 'L', n, kd, ab, kd + 1, w, z, n, work, size(work), iwork, size(iwork), info)
        call require(info == 0, 'dsbevd info')
        call check_eigen(a, z, w, 'dsbevd')
        call fortnum_dsbevx('V', 'I', 'L', n, kd, ab, kd + 1, q, n, 0.0_dp, 0.0_dp, 5, 9, 0.0_dp, m, wx, &
            zx, n, work, iwork, ifail, info)
        call require(info == 0 .and. m == 5 .and. maxval(abs(wx(1:5) - w(5:9))) <= 1e-12_dp, 'dsbevx indices')
        call fortnum_dsbevx('N', 'V', 'L', n, kd, ab, kd + 1, q, n, w(10) - 1e-9_dp, w(20), 0, 0, 0.0_dp, m, wx, &
            zx, n, work, iwork, ifail, info)
        call require(info == 0 .and. m == 11, 'dsbevx value range (vl, vu]')
    end subroutine band

end program test_fortnum_lapack_syev
