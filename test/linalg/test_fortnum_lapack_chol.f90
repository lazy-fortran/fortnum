program test_fortnum_lapack_chol
    ! Cholesky, triangular inverses, condition estimate and pivoted Cholesky
    ! against quad-precision reconstructions and residuals.
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: eps = real(epsilon(1.0_dp), qp)

    call spd_real(5, 'L')
    call spd_real(150, 'L')
    call spd_real(150, 'U')
    call hpd_complex(70, 'U')
    call hpd_complex(70, 'L')
    call not_definite()
    call triangular(90)
    call pivoted(40, 23)
    print '(a)', 'Cholesky family matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

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

    function spd(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        real(dp) :: a(n, n), m(n, n)
        integer :: i
        call fill(m, seed)
        a = matmul(m, transpose(m))
        do i = 1, n
            a(i, i) = a(i, i) + 0.1_dp*n
        end do
    end function spd

    subroutine spd_real(n, uplo)
        integer, intent(in) :: n
        character, intent(in) :: uplo
        real(dp) :: a(n, n), f(n + 2, n), b(n, 3), x(n, 3), inv(n, n), work(3*n), rcond, anorm
        real(qp) :: l(n, n), exact
        integer :: info, iwork(n), i, j
        a = spd(n, 0.43_dp)
        f = 7
        f(1:n, :) = a
        call fortnum_dpotrf(uplo, n, f, n + 2, info)
        call require(info == 0, 'dpotrf info')
        l = 0
        do j = 1, n
            do i = 1, n
                if (uplo == 'L' .and. i >= j) l(i, j) = f(i, j)
                if (uplo == 'U' .and. i <= j) l(j, i) = f(i, j)
            end do
        end do
        call require(maxval(abs(matmul(l, transpose(l)) - real(a, qp))) <= 4*n*eps*maxval(abs(real(a, qp))), &
            'dpotrf reconstructs A')
        call fill(b, 1.3_dp)
        x = b
        call fortnum_dpotrs(uplo, n, 3, f, n + 2, x, n, info)
        call require(info == 0 .and. maxval(abs(matmul(real(a, qp), real(x, qp)) - real(b, qp))) <= &
            8*n*eps*(maxval(sum(abs(real(a, qp)), 2))*maxval(abs(real(x, qp)))), 'dpotrs residual')
        ! condition estimate bracketed by the exact value from the inverse
        anorm = maxval(sum(abs(a), 1))
        call fortnum_dpocon(uplo, n, f, n + 2, anorm, rcond, work, iwork, info)
        inv = f(1:n, :)
        call fortnum_dpotri(uplo, n, inv, n, info)
        call require(info == 0, 'dpotri info')
        do j = 1, n
            do i = 1, n
                if ((uplo == 'L' .and. i < j) .or. (uplo == 'U' .and. i > j)) inv(i, j) = inv(j, i)
            end do
        end do
        call require(maxval(abs(matmul(real(a, qp), real(inv, qp)) - identity(n))) <= &
            16*n*eps*maxval(sum(abs(real(a, qp)), 2))*maxval(abs(real(inv, qp))), 'dpotri inverse')
        exact = 1/(real(anorm, qp)*maxval(sum(abs(real(inv, qp)), 1)))
        call require(rcond >= exact*(1 - 1e-8_qp) .and. rcond <= 10*exact, 'dpocon within a factor 10')
        x = b
        f(1:n, :) = a
        call fortnum_dposv(uplo, n, 3, f, n + 2, x, n, info)
        call require(info == 0 .and. maxval(abs(matmul(real(a, qp), real(x, qp)) - real(b, qp))) <= &
            8*n*eps*(maxval(sum(abs(real(a, qp)), 2))*maxval(abs(real(x, qp)))), 'dposv')
    end subroutine spd_real

    function identity(n) result(e)
        integer, intent(in) :: n
        real(qp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function identity

    subroutine hpd_complex(n, uplo)
        integer, intent(in) :: n
        character, intent(in) :: uplo
        real(dp) :: re(n, n), im(n, n)
        complex(dp) :: m(n, n), a(n, n), f(n, n), b(n, 2), x(n, 2), inv(n, n)
        complex(qp) :: l(n, n)
        integer :: info, i, j
        call fill(re, 0.29_dp)
        call fill(im, 0.61_dp)
        m = cmplx(re, im, dp)
        a = matmul(m, conjg(transpose(m)))
        do i = 1, n
            a(i, i) = a(i, i) + 0.1_dp*n
        end do
        f = a
        call fortnum_zpotrf(uplo, n, f, n, info)
        call require(info == 0, 'zpotrf info')
        l = 0
        do j = 1, n
            do i = 1, n
                if (uplo == 'L' .and. i >= j) l(i, j) = f(i, j)
                if (uplo == 'U' .and. i <= j) l(j, i) = conjg(cmplx(f(i, j), kind=qp))
            end do
        end do
        call require(maxval(abs(matmul(l, conjg(transpose(l))) - cmplx(a, kind=qp))) <= &
            8*n*eps*maxval(abs(cmplx(a, kind=qp))), 'zpotrf reconstructs A')
        b(:, 1) = cmplx(re(:, 1), -im(:, 2), dp)
        b(:, 2) = cmplx(im(:, 3), re(:, 4), dp)
        x = b
        f = a
        call fortnum_zposv(uplo, n, 2, f, n, x, n, info)
        call require(info == 0 .and. maxval(abs(matmul(cmplx(a, kind=qp), cmplx(x, kind=qp)) - b)) <= &
            16*n*eps*maxval(sum(abs(cmplx(a, kind=qp)), 2))*maxval(abs(cmplx(x, kind=qp))), 'zposv')
        inv = f
        call fortnum_zpotri(uplo, n, inv, n, info)
        do j = 1, n
            do i = 1, n
                if ((uplo == 'L' .and. i < j) .or. (uplo == 'U' .and. i > j)) inv(i, j) = conjg(inv(j, i))
            end do
        end do
        call require(info == 0 .and. maxval(abs(matmul(cmplx(a, kind=qp), cmplx(inv, kind=qp)) - identity(n))) <= &
            32*n*eps*maxval(sum(abs(cmplx(a, kind=qp)), 2))*maxval(abs(cmplx(inv, kind=qp))), 'zpotri inverse')
    end subroutine hpd_complex

    subroutine not_definite()
        real(dp) :: a(3, 3)
        integer :: info
        a = reshape([4.0_dp, 2.0_dp, 0.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp], [3, 3])
        call fortnum_dpotrf('L', 3, a, 3, info)
        call require(info == 2, 'dpotrf reports the first non-positive pivot')
    end subroutine not_definite

    subroutine triangular(n)
        integer, intent(in) :: n
        real(dp) :: t(n, n), inv(n, n), b(n, 2), x(n, 2)
        complex(dp) :: zt(n, n), zinv(n, n)
        real(qp) :: tq(n, n)
        integer :: info, i, j, k
        call fill(t, 0.77_dp)
        do i = 1, n
            t(i, i) = 3 + cos(real(i, dp))
        end do
        t = t/sqrt(real(n, dp))
        do k = 1, 2
            tq = 0
            do j = 1, n
                do i = 1, n
                    if ((k == 1 .and. i <= j) .or. (k == 2 .and. i >= j)) tq(i, j) = t(i, j)
                end do
            end do
            inv = t
            call fortnum_dtrtri('UL'(k:k), 'N', n, inv, n, info)
            call require(info == 0, 'dtrtri info')
            do j = 1, n
                do i = 1, n
                    if ((k == 1 .and. i > j) .or. (k == 2 .and. i < j)) inv(i, j) = 0
                end do
            end do
            call require(maxval(abs(matmul(tq, real(inv, qp)) - identity(n))) <= &
                16*n*eps*maxval(sum(abs(tq), 2))*maxval(abs(real(inv, qp))), 'dtrtri inverse')
            call fill(b, 0.5_dp)
            x = b
            call fortnum_dtrtrs('UL'(k:k), 'T', 'N', n, 2, t, n, x, n, info)
            call require(info == 0 .and. maxval(abs(matmul(transpose(tq), real(x, qp)) - real(b, qp))) <= &
                16*n*eps*maxval(sum(abs(tq), 1))*maxval(abs(real(x, qp))), 'dtrtrs T')
        end do
        zt = cmplx(t, transpose(t), dp)
        zinv = zt
        call fortnum_ztrtri('U', 'U', n, zinv, n, info)
        call require(info == 0 .and. all(zinv == zinv), 'ztrtri unit')
        do i = 1, n
            call require(zinv(i, i) == zt(i, i), 'ztrtri keeps a unit diagonal untouched')
        end do
        t(5, 5) = 0
        call fortnum_dtrtri('L', 'N', n, t, n, info)
        call require(info == 5, 'dtrtri reports a zero diagonal element')
    end subroutine triangular

    ! rank-r positive semidefinite matrix: the pivoted factor finds r and
    ! P^T A P = L L^T on its first r columns
    subroutine pivoted(n, r)
        integer, intent(in) :: n, r
        real(dp) :: m(n, r), a(n, n), f(n, n), work(2*n)
        real(qp) :: l(n, r), pap(n, n)
        integer :: piv(n), rank, info, i, j
        call fill(m, 0.39_dp)
        a = matmul(m, transpose(m))
        f = a
        call fortnum_dpstrf('L', n, f, n, piv, rank, -1.0_dp, work, info)
        call require(info == 1 .and. rank == r, 'dpstrf finds the rank')
        l = 0
        do j = 1, r
            l(j:n, j) = f(j:n, j)
        end do
        do j = 1, n
            do i = 1, n
                pap(i, j) = a(piv(i), piv(j))
            end do
        end do
        call require(maxval(abs(matmul(l, transpose(l)) - pap)) <= 1e-10_qp*maxval(abs(real(a, qp))), &
            'dpstrf reconstructs P^T A P')
        do j = 2, r
            call require(f(j, j) <= f(j - 1, j - 1)*(1 + 1e-12_dp), 'pivots decrease')
        end do
    end subroutine pivoted

end program test_fortnum_lapack_chol
