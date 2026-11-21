program test_fortnum_lapack_lu
    ! LU, solves, inverse and band LU against quad-precision residuals and
    ! reconstructions; singular and illegal inputs report LAPACK info codes.
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: eps = real(epsilon(1.0_dp), qp)

    call dense_real(5)
    call dense_real(150)
    call dense_complex(7)
    call dense_complex(130)
    call rectangular(90, 40)
    call rectangular(40, 90)
    call singular()
    call band_real(40, 3, 2)
    call band_real(9, 0, 4)
    call band_complex(37, 2, 3)
    call interchanges(50, 70, 1)
    call interchanges(50, 70, -1)
    call interchanges(1500, 1400, 1)
    print '(a)', 'LU family matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    ! DLASWP against the interchanges applied one by one to row labels
    subroutine interchanges(m, n, incx)
        integer, intent(in) :: m, n, incx
        real(dp), allocatable :: a(:, :), a0(:, :)
        integer, allocatable :: ipiv(:), rows(:)
        integer :: k, k1, k2, t
        allocate (a(m + 2, n), ipiv(m), rows(m))
        call fill(a, 0.53_dp)
        a0 = a
        k1 = 3
        k2 = m - 4
        do k = 1, m
            ipiv(k) = k + mod(7*k*k + 3, m - k + 1)
        end do
        do k = 1, m
            rows(k) = k
        end do
        if (incx > 0) then
            do k = k1, k2
                t = rows(k); rows(k) = rows(ipiv(k)); rows(ipiv(k)) = t
            end do
        else
            ! reverse order: ipiv(k) for k = k2 down to k1
            do k = k2, k1, -1
                t = rows(k); rows(k) = rows(ipiv(k)); rows(ipiv(k)) = t
            end do
        end if
        call fortnum_dlaswp(n, a, m + 2, k1, k2, ipiv, incx)
        call require(all(a(1:m, :) == a0(rows, :)) .and. all(a(m + 1:, :) == a0(m + 1:, :)), 'dlaswp')
    end subroutine interchanges

    subroutine fill(x, seed)
        real(dp), intent(out) :: x(:, :)
        real(dp), intent(in) :: seed
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(seed*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j + seed))
            end do
        end do
    end subroutine fill

    ! normwise backward error of op(A) X = B in quad
    real(qp) function backward(a, x, b, trans)
        real(dp), intent(in) :: a(:, :), x(:, :), b(:, :)
        character, intent(in) :: trans
        real(qp), allocatable :: oa(:, :)
        oa = real(a, qp)
        if (trans /= 'N') oa = transpose(oa)
        backward = maxval(abs(matmul(oa, real(x, qp)) - real(b, qp)))/ &
            (maxval(sum(abs(oa), 2))*maxval(abs(real(x, qp))) + maxval(abs(real(b, qp))))
    end function backward

    real(qp) function zbackward(a, x, b, trans)
        complex(dp), intent(in) :: a(:, :), x(:, :), b(:, :)
        character, intent(in) :: trans
        complex(qp), allocatable :: oa(:, :)
        oa = cmplx(a, kind=qp)
        if (trans == 'T') oa = transpose(oa)
        if (trans == 'C') oa = conjg(transpose(oa))
        zbackward = maxval(abs(matmul(oa, cmplx(x, kind=qp)) - cmplx(b, kind=qp)))/ &
            (maxval(sum(abs(oa), 2))*maxval(abs(cmplx(x, kind=qp))) + maxval(abs(cmplx(b, kind=qp))))
    end function zbackward

    subroutine dense_real(n)
        integer, intent(in) :: n
        real(dp) :: a(n + 3, n), lu(n + 3, n), b(n + 1, 4), x(n + 1, 4), work(n)
        integer :: ipiv(n), info, k
        character :: tr
        call fill(a, 0.37_dp)
        call fill(b, 1.1_dp)
        lu = a
        x = b
        call fortnum_dgesv(n, 4, lu, n + 3, ipiv, x, n + 1, info)
        call require(info == 0, 'dgesv info')
        call require(backward(a(1:n, :), x(1:n, :), b(1:n, :), 'N') <= 4*n*eps, 'dgesv backward error')
        do k = 1, 2
            tr = 'NT'(k:k)
            x = b
            call fortnum_dgetrs(tr, n, 4, lu, n + 3, ipiv, x, n + 1, info)
            call require(info == 0 .and. backward(a(1:n, :), x(1:n, :), b(1:n, :), tr) <= 4*n*eps, 'dgetrs '//tr)
        end do
        call fortnum_dgetri(n, lu, n + 3, ipiv, work, n, info)
        call require(info == 0, 'dgetri info')
        call require(maxval(abs(matmul(real(a(1:n, :), qp), real(lu(1:n, :), qp)) - identity(n))) <= &
            8*n*eps*maxval(sum(abs(real(a(1:n, :), qp)), 2))*maxval(abs(real(lu(1:n, :), qp))), 'dgetri inverse')
        call fortnum_dgetri(n, lu, n + 3, ipiv, work, -1, info)
        call require(info == 0 .and. nint(work(1)) == n, 'dgetri workspace query')
        call fortnum_dgesv(n, 4, lu, n - 1, ipiv, x, n + 1, info)
        call require(info == -4, 'dgesv reports an illegal lda')
    end subroutine dense_real

    function identity(n) result(e)
        integer, intent(in) :: n
        real(qp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function identity

    subroutine dense_complex(n)
        integer, intent(in) :: n
        real(dp) :: re(n, n), im(n, n), rb(n, 3), ib(n, 3)
        complex(dp) :: a(n, n), lu(n, n), b(n, 3), x(n, 3), work(n)
        integer :: ipiv(n), info, k
        character :: tr
        call fill(re, 0.37_dp)
        call fill(im, 0.71_dp)
        call fill(rb, 1.1_dp)
        call fill(ib, 1.7_dp)
        a = cmplx(re, im, dp)
        b = cmplx(rb, ib, dp)
        lu = a
        x = b
        call fortnum_zgesv(n, 3, lu, n, ipiv, x, n, info)
        call require(info == 0 .and. zbackward(a, x, b, 'N') <= 8*n*eps, 'zgesv backward error')
        do k = 1, 3
            tr = 'NTC'(k:k)
            x = b
            call fortnum_zgetrs(tr, n, 3, lu, n, ipiv, x, n, info)
            call require(info == 0 .and. zbackward(a, x, b, tr) <= 8*n*eps, 'zgetrs '//tr)
        end do
        call fortnum_zgetri(n, lu, n, ipiv, work, n, info)
        call require(info == 0 .and. maxval(abs(matmul(cmplx(a, kind=qp), cmplx(lu, kind=qp)) - identity(n))) <= &
            16*n*eps*maxval(sum(abs(cmplx(a, kind=qp)), 2))*maxval(abs(cmplx(lu, kind=qp))), 'zgetri inverse')
    end subroutine dense_complex

    ! P L U reproduces A for tall and wide matrices
    subroutine rectangular(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), f(m, n)
        real(qp) :: l(m, min(m, n)), u(min(m, n), n), plu(m, n), t
        integer :: ipiv(min(m, n)), info, i, j
        call fill(a, 0.53_dp)
        f = a
        call fortnum_dgetrf(m, n, f, m, ipiv, info)
        call require(info == 0, 'dgetrf rectangular info')
        l = 0
        u = 0
        do j = 1, min(m, n)
            l(j, j) = 1
            l(j + 1:m, j) = f(j + 1:m, j)
        end do
        do j = 1, n
            do i = 1, min(j, min(m, n))
                u(i, j) = f(i, j)
            end do
        end do
        plu = matmul(l, u)
        ! undo the interchanges in reverse order
        do j = min(m, n), 1, -1
            do i = 1, n
                t = plu(j, i)
                plu(j, i) = plu(ipiv(j), i)
                plu(ipiv(j), i) = t
            end do
        end do
        call require(maxval(abs(plu - real(a, qp))) <= 4*max(m, n)*eps*maxval(abs(real(a, qp)))*min(m, n), &
            'P L U = A')
    end subroutine rectangular

    subroutine singular()
        real(dp) :: a(4, 4)
        integer :: ipiv(4), info
        a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 2.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [4, 4])
        call fortnum_dgetrf(4, 4, a, 4, ipiv, info)
        call require(info == 2, 'dgetrf reports the zero pivot of the dependent column')
    end subroutine singular

    ! band LU and solves against the dense solution of the same matrix
    subroutine band_real(n, kl, ku)
        integer, intent(in) :: n, kl, ku
        real(dp) :: a(n, n), dense(n, n), ab(2*kl + ku + 2, n), b(n, 2), x(n, 2), y(n, 2)
        integer :: ipiv(n), info, i, j, k
        character :: tr
        call fill(a, 0.29_dp)
        do j = 1, n
            do i = 1, n
                if (i - j > kl .or. j - i > ku) a(i, j) = 0
            end do
        end do
        call fill(b, 0.9_dp)
        do k = 1, 2
            tr = 'NT'(k:k)
            ab = 777
            do j = 1, n
                do i = max(1, j - ku), min(n, j + kl)
                    ab(kl + ku + 1 + i - j, j) = a(i, j)
                end do
            end do
            call fortnum_dgbtrf(n, n, kl, ku, ab, size(ab, 1), ipiv, info)
            call require(info == 0, 'dgbtrf info')
            x = b
            call fortnum_dgbtrs(tr, n, kl, ku, 2, ab, size(ab, 1), ipiv, x, n, info)
            call require(info == 0 .and. backward(a, x, b, tr) <= 8*n*eps, 'dgbtrs '//tr)
            dense = a
            y = b
            if (tr == 'T') dense = transpose(a)
            call fortnum_dgesv(n, 2, dense, n, ipiv, y, n, info)
            call require(maxval(abs(x - y)) <= 1e-10_dp*maxval(abs(y)), 'band and dense solutions agree')
        end do
        ab = 0
        do j = 1, n
            do i = max(1, j - ku), min(n, j + kl)
                ab(kl + ku + 1 + i - j, j) = a(i, j)
            end do
        end do
        x = b
        call fortnum_dgbsv(n, kl, ku, 2, ab, size(ab, 1), ipiv, x, n, info)
        call require(info == 0 .and. backward(a, x, b, 'N') <= 8*n*eps, 'dgbsv')
        call fortnum_dgbsv(n, kl, ku, 2, ab, 2*kl + ku, ipiv, x, n, info)
        call require(info == -6, 'dgbsv reports a short ldab')
    end subroutine band_real

    subroutine band_complex(n, kl, ku)
        integer, intent(in) :: n, kl, ku
        real(dp) :: re(n, n), im(n, n), rb(n, 1)
        complex(dp) :: a(n, n), ab(2*kl + ku + 1, n), b(n, 1), x(n, 1)
        integer :: ipiv(n), info, i, j, k
        character :: tr
        call fill(re, 0.31_dp)
        call fill(im, 0.83_dp)
        call fill(rb, 1.3_dp)
        a = cmplx(re, im, dp)
        do j = 1, n
            do i = 1, n
                if (i - j > kl .or. j - i > ku) a(i, j) = 0
            end do
        end do
        b = cmplx(rb, -rb, dp)
        do k = 1, 3
            tr = 'NTC'(k:k)
            ab = 0
            do j = 1, n
                do i = max(1, j - ku), min(n, j + kl)
                    ab(kl + ku + 1 + i - j, j) = a(i, j)
                end do
            end do
            call fortnum_zgbtrf(n, n, kl, ku, ab, size(ab, 1), ipiv, info)
            x = b
            call fortnum_zgbtrs(tr, n, kl, ku, 1, ab, size(ab, 1), ipiv, x, n, info)
            call require(info == 0 .and. zbackward(a, x, b, tr) <= 16*n*eps, 'zgbtrs '//tr)
        end do
    end subroutine band_complex

end program test_fortnum_lapack_lu
