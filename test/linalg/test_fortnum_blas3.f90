program test_fortnum_blas3
    ! Level-3 BLAS against quad precision: TRSM by its backward error, TRMM,
    ! SYRK, SYR2K and HERK by the exact product; every side, triangle,
    ! transpose and unit diagonal, sizes across the 64 x 64 diagonal blocks.
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm, &
        fortnum_dsyrk, fortnum_dsyr2k, fortnum_zherk
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    character, parameter :: sides(2) = ['L', 'R'], uplos(2) = ['U', 'L'], diags(2) = ['N', 'U']
    character, parameter :: rtrans(2) = ['N', 'T'], ztrans(3) = ['N', 'T', 'C']
    integer, parameter :: shapes(2, 3) = reshape([5, 3, 70, 33, 9, 130], [2, 3])
    integer :: is, iu, it, id, sh

    do sh = 1, size(shapes, 2)
        do is = 1, 2
            do iu = 1, 2
                do id = 1, 2
                    do it = 1, 2
                        call check_real(sides(is), uplos(iu), rtrans(it), diags(id), shapes(1, sh), shapes(2, sh))
                    end do
                    do it = 1, 3
                        call check_complex(sides(is), uplos(iu), ztrans(it), diags(id), shapes(1, sh), &
                            shapes(2, sh))
                    end do
                end do
            end do
        end do
    end do
    ! real triangles of order 300: two diagonal blocks of the packed solve,
    ! the update of the rows below them and partial register tiles
    do iu = 1, 2
        do id = 1, 2
            do it = 1, 2
                call check_real('L', uplos(iu), rtrans(it), diags(id), 300, 7)
                call check_real('R', uplos(iu), rtrans(it), diags(id), 7, 300)
            end do
        end do
    end do
    do iu = 1, 2
        do it = 1, 2
            call check_rank_k(uplos(iu), rtrans(it), 150, 37)
            call check_rank_k(uplos(iu), rtrans(it), 7, 0)
        end do
        call check_herk(uplos(iu), 'N', 140, 21)
        call check_herk(uplos(iu), 'C', 9, 30)
    end do
    print '(a)', 'level-3 BLAS match quad precision'

contains

    ! Well-conditioned triangle: diagonal 2 + something, off-diagonal small.
    subroutine fill_tri(a)
        real(dp), intent(out) :: a(:, :)
        integer :: i, j
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                a(i, j) = 0.3_dp*sin(0.7_dp*i + 1.3_dp*j)/sqrt(real(size(a, 1), dp))
            end do
            a(j, j) = 2.0_dp + cos(real(j, dp))
        end do
    end subroutine fill_tri

    subroutine fill(x)
        real(dp), intent(out) :: x(:, :)
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(0.37_dp*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j))
            end do
        end do
    end subroutine fill

    ! op(T) as a dense quad matrix, T the uplo triangle of a (unit diagonal if asked)
    function optri_q(a, uplo, trans, diag) result(t)
        real(dp), intent(in) :: a(:, :)
        character, intent(in) :: uplo, trans, diag
        real(qp) :: t(size(a, 1), size(a, 1))
        integer :: i, j
        t = 0
        do j = 1, size(a, 1)
            do i = 1, size(a, 1)
                if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) t(i, j) = a(i, j)
            end do
            if (diag == 'U') t(j, j) = 1
        end do
        if (trans /= 'N') t = transpose(t)
    end function optri_q

    function zoptri_q(a, uplo, trans, diag) result(t)
        complex(dp), intent(in) :: a(:, :)
        character, intent(in) :: uplo, trans, diag
        complex(qp) :: t(size(a, 1), size(a, 1))
        integer :: i, j
        t = 0
        do j = 1, size(a, 1)
            do i = 1, size(a, 1)
                if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) t(i, j) = a(i, j)
            end do
            if (diag == 'U') t(j, j) = 1
        end do
        if (trans == 'T') t = transpose(t)
        if (trans == 'C') t = conjg(transpose(t))
    end function zoptri_q

    subroutine check_real(side, uplo, trans, diag, m, n)
        character, intent(in) :: side, uplo, trans, diag
        integer, intent(in) :: m, n
        real(dp), allocatable :: a(:, :), b0(:, :), x(:, :)
        real(qp), allocatable :: t(:, :), r(:, :), s(:, :)
        real(dp), parameter :: alpha = 0.75_dp
        integer :: na, lda, ldb
        na = merge(m, n, side == 'L')
        lda = na + 2
        ldb = m + 1
        allocate (a(lda, na), b0(ldb, n))
        call fill_tri(a(1:na, :))
        a(na + 1:, :) = 777
        call fill(b0)
        ! solve, then check op(T) X = alpha B by the backward error
        x = b0
        call fortnum_dtrsm(side, uplo, trans, diag, m, n, alpha, a, lda, x, ldb)
        t = optri_q(a(1:na, :), uplo, trans, diag)
        if (side == 'L') then
            r = matmul(t, real(x(1:m, :), qp)) - alpha*real(b0(1:m, :), qp)
            s = matmul(abs(t), abs(real(x(1:m, :), qp))) + abs(alpha*real(b0(1:m, :), qp))
        else
            r = matmul(real(x(1:m, :), qp), t) - alpha*real(b0(1:m, :), qp)
            s = matmul(abs(real(x(1:m, :), qp)), abs(t)) + abs(alpha*real(b0(1:m, :), qp))
        end if
        if (any(abs(r) > 4*(na + 2)*real(epsilon(1.0_dp), qp)*s)) then
            print '(a, 4(1x, a), 2i5)', 'dtrsm residual', side, uplo, trans, diag, m, n
            error stop 1
        end if
        if (any(x(m + 1:, :) /= b0(m + 1:, :))) error stop 'dtrsm wrote outside B'
        ! multiply and compare with the exact product
        x = b0
        call fortnum_dtrmm(side, uplo, trans, diag, m, n, alpha, a, lda, x, ldb)
        if (side == 'L') then
            r = alpha*matmul(t, real(b0(1:m, :), qp))
            s = abs(alpha)*matmul(abs(t), abs(real(b0(1:m, :), qp)))
        else
            r = alpha*matmul(real(b0(1:m, :), qp), t)
            s = abs(alpha)*matmul(abs(real(b0(1:m, :), qp)), abs(t))
        end if
        if (any(abs(real(x(1:m, :), qp) - r) > 4*(na + 2)*real(epsilon(1.0_dp), qp)*s)) then
            print '(a, 4(1x, a), 2i5)', 'dtrmm mismatch', side, uplo, trans, diag, m, n
            error stop 1
        end if
    end subroutine check_real

    subroutine check_complex(side, uplo, trans, diag, m, n)
        character, intent(in) :: side, uplo, trans, diag
        integer, intent(in) :: m, n
        complex(dp), allocatable :: a(:, :), b0(:, :), x(:, :)
        real(dp), allocatable :: re(:, :), im(:, :)
        complex(qp), allocatable :: t(:, :), r(:, :)
        real(qp), allocatable :: s(:, :)
        complex(dp), parameter :: alpha = (0.5_dp, -0.25_dp)
        integer :: na, lda, ldb
        na = merge(m, n, side == 'L')
        lda = na + 1
        ldb = m + 2
        allocate (re(na, na), im(na, na))
        call fill_tri(re)
        call fill(im)
        allocate (a(lda, na))
        a = (777.0_dp, 777.0_dp)
        a(1:na, :) = cmplx(re, 0.2_dp*im, dp)
        deallocate (re, im)
        allocate (re(ldb, n), im(ldb, n))
        call fill(re)
        call fill(im)
        b0 = cmplx(re, -im, dp)
        x = b0
        call fortnum_ztrsm(side, uplo, trans, diag, m, n, alpha, a, lda, x, ldb)
        t = zoptri_q(a(1:na, :), uplo, trans, diag)
        if (side == 'L') then
            r = matmul(t, cmplx(x(1:m, :), kind=qp)) - alpha*cmplx(b0(1:m, :), kind=qp)
            s = matmul(abs(t), abs(cmplx(x(1:m, :), kind=qp))) + abs(alpha*cmplx(b0(1:m, :), kind=qp))
        else
            r = matmul(cmplx(x(1:m, :), kind=qp), t) - alpha*cmplx(b0(1:m, :), kind=qp)
            s = matmul(abs(cmplx(x(1:m, :), kind=qp)), abs(t)) + abs(alpha*cmplx(b0(1:m, :), kind=qp))
        end if
        if (any(abs(r) > 8*(na + 2)*real(epsilon(1.0_dp), qp)*s)) then
            print '(a, 4(1x, a), 2i5)', 'ztrsm residual', side, uplo, trans, diag, m, n
            error stop 1
        end if
        x = b0
        call fortnum_ztrmm(side, uplo, trans, diag, m, n, alpha, a, lda, x, ldb)
        if (side == 'L') then
            r = alpha*matmul(t, cmplx(b0(1:m, :), kind=qp))
            s = abs(alpha)*matmul(abs(t), abs(cmplx(b0(1:m, :), kind=qp)))
        else
            r = alpha*matmul(cmplx(b0(1:m, :), kind=qp), t)
            s = abs(alpha)*matmul(abs(cmplx(b0(1:m, :), kind=qp)), abs(t))
        end if
        if (any(abs(cmplx(x(1:m, :), kind=qp) - r) > 8*(na + 2)*real(epsilon(1.0_dp), qp)*s)) then
            print '(a, 4(1x, a), 2i5)', 'ztrmm mismatch', side, uplo, trans, diag, m, n
            error stop 1
        end if
        if (any(x(m + 1:, :) /= b0(m + 1:, :))) error stop 'ztrmm wrote outside B'
    end subroutine check_complex

    subroutine check_rank_k(uplo, trans, n, k)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k
        real(dp), allocatable :: a(:, :), b(:, :), c(:, :), c0(:, :)
        real(qp), allocatable :: oa(:, :), ob(:, :), want(:, :), bound(:, :)
        real(dp), parameter :: alpha = 1.25_dp, beta = -0.5_dp
        integer :: ra, ca, i, j, pass
        ra = merge(n, k, trans == 'N')
        ca = merge(k, n, trans == 'N')
        allocate (a(max(1, ra) + 1, max(1, ca)), b(max(1, ra) + 1, max(1, ca)), c(n + 3, n))
        call fill(a)
        b = cos(a)
        call fill(c)
        c0 = c
        do pass = 1, 2
            c = c0
            if (trans == 'N') then
                oa = real(a(1:n, 1:k), qp)
                ob = real(b(1:n, 1:k), qp)
            else
                oa = transpose(real(a(1:k, 1:n), qp))
                ob = transpose(real(b(1:k, 1:n), qp))
            end if
            if (pass == 1) then
                call fortnum_dsyrk(uplo, trans, n, k, alpha, a, size(a, 1), beta, c, size(c, 1))
                want = alpha*matmul(oa, transpose(oa)) + beta*real(c0(1:n, :), qp)
                bound = abs(alpha)*matmul(abs(oa), transpose(abs(oa))) + abs(beta*real(c0(1:n, :), qp))
            else
                call fortnum_dsyr2k(uplo, trans, n, k, alpha, a, size(a, 1), b, size(b, 1), beta, c, size(c, 1))
                want = alpha*(matmul(oa, transpose(ob)) + matmul(ob, transpose(oa))) + beta*real(c0(1:n, :), qp)
                bound = abs(alpha)*(matmul(abs(oa), transpose(abs(ob))) + matmul(abs(ob), transpose(abs(oa)))) &
                    + abs(beta*real(c0(1:n, :), qp))
            end if
            do j = 1, n
                do i = 1, n
                    if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) then
                        if (abs(real(c(i, j), qp) - want(i, j)) > 4*(k + 2)*real(epsilon(1.0_dp), qp)*bound(i, j)) then
                            print '(a, i2, 2(1x, a), 2i5)', 'rank-k mismatch, pass', pass, uplo, trans, n, k
                            error stop 1
                        end if
                    else if (c(i, j) /= c0(i, j)) then
                        error stop 'rank-k update wrote the other triangle'
                    end if
                end do
            end do
            if (any(c(n + 1:, :) /= c0(n + 1:, :))) error stop 'rank-k update wrote outside C'
        end do
    end subroutine check_rank_k

    subroutine check_herk(uplo, trans, n, k)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k
        complex(dp), allocatable :: a(:, :), c(:, :), c0(:, :)
        real(dp), allocatable :: re(:, :), im(:, :)
        complex(qp), allocatable :: oa(:, :), want(:, :)
        real(qp), allocatable :: bound(:, :)
        real(dp), parameter :: alpha = 0.8_dp, beta = 1.5_dp
        integer :: ra, ca, i, j
        ra = merge(n, k, trans == 'N')
        ca = merge(k, n, trans == 'N')
        allocate (re(ra + 2, ca), im(ra + 2, ca))
        call fill(re)
        call fill(im)
        a = cmplx(re, 0.5_dp*im, dp)
        deallocate (re, im)
        allocate (re(n + 1, n), im(n + 1, n))
        call fill(re)
        call fill(im)
        c = cmplx(re, im, dp)
        c0 = c
        call fortnum_zherk(uplo, trans, n, k, alpha, a, size(a, 1), beta, c, size(c, 1))
        if (trans == 'N') then
            oa = cmplx(a(1:n, 1:k), kind=qp)
        else
            oa = conjg(transpose(cmplx(a(1:k, 1:n), kind=qp)))
        end if
        ! beta C uses only the real part of the diagonal, as in reference ZHERK
        want = alpha*matmul(oa, conjg(transpose(oa))) + beta*cmplx(c0(1:n, :), kind=qp)
        bound = abs(alpha)*matmul(abs(oa), transpose(abs(oa))) + abs(beta)*abs(cmplx(c0(1:n, :), kind=qp))
        do j = 1, n
            want(j, j) = cmplx(real(want(j, j), qp), 0, qp)
            do i = 1, n
                if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) then
                    if (abs(cmplx(c(i, j), kind=qp) - want(i, j)) > 8*(k + 2)*real(epsilon(1.0_dp), qp)*bound(i, j)) then
                        print '(2(1x, a), 2i5, 2i4)', 'zherk mismatch', uplo, n, k, i, j
                        error stop 1
                    end if
                else if (c(i, j) /= c0(i, j)) then
                    error stop 'zherk wrote the other triangle'
                end if
            end do
            if (aimag(c(j, j)) /= 0) error stop 'zherk diagonal not real'
        end do
    end subroutine check_herk

end program test_fortnum_blas3
