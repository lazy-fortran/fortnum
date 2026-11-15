program test_fortnum_lapack77
    ! The standard names of fortnum_lapack77, called through interface
    ! blocks as client codes declare them, give bitwise the results of the
    ! fortnum routines they forward to: one or two routines of every family.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm
    use fortnum_blas12, only: fortnum_dasum, fortnum_zgemv
    use fortnum_lapack_lu, only: fortnum_dgesv, fortnum_zgbsv
    use fortnum_lapack_chol, only: fortnum_dpotrf, fortnum_dpotrs
    use fortnum_lapack_qr, only: fortnum_dgels
    use fortnum_lapack_syev, only: fortnum_dsyev
    use fortnum_lapack_geev, only: fortnum_dgeev, fortnum_zgeev
    use fortnum_lapack_svd, only: fortnum_dgelss
    implicit none
    interface
        real(dp) function dasum(n, x, incx)
            import :: dp
            integer, intent(in) :: n, incx
            real(dp), intent(in) :: x(*)
        end function dasum
        subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            complex(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
            complex(dp), intent(inout) :: y(*)
        end subroutine zgemv
        subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            real(dp), intent(inout) :: c(ldc, *)
        end subroutine dgemm
        subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import :: dp
            character, intent(in) :: transa, transb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            complex(dp), intent(inout) :: c(ldc, *)
        end subroutine zgemm
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine dgesv
        subroutine zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
            complex(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
            integer, intent(out) :: ipiv(*), info
        end subroutine zgbsv
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine dgels
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine dgeev
        subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            import :: dp
            character, intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *)
            complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
            real(dp), intent(out) :: rwork(*)
            integer, intent(out) :: info
        end subroutine zgeev
        subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
            import :: dp
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: s(*), work(*)
            real(dp), intent(in) :: rcond
            integer, intent(out) :: rank, info
        end subroutine dgelss
    end interface
    integer, parameter :: n = 12
    real(dp) :: a(n, n), a1(n, n), a2(n, n), b1(n, 2), b2(n, 2), c1(n, n), c2(n, n), work(40*n)
    real(dp) :: wr1(n), wi1(n), wr2(n), wi2(n), v1(n, n), v2(n, n), u1(n, n), u2(n, n), s1(n), s2(n)
    complex(dp) :: z(n, n), z1(n, n), z2(n, n), zb1(n, 1), zb2(n, 1), zw1(n), zw2(n), zwork(4*n)
    complex(dp) :: zv1(n, n), zv2(n, n), zu1(n, n), zu2(n, n)
    complex(dp) :: ab1(7, n), ab2(7, n)
    real(dp) :: rwork(2*n)
    integer :: i, j, ipiv1(n), ipiv2(n), info1, info2, rank1, rank2

    do j = 1, n
        do i = 1, n
            a(i, j) = sin(0.37_dp*i + 1.91_dp*j + 0.11_dp*i*j)
        end do
    end do
    z = cmplx(a, cos(transpose(a)), dp)
    b1(:, 1) = a(:, 3)
    b1(:, 2) = a(5, :)

    c1 = 0.5_dp
    c2 = c1
    call dgemm('N', 'T', n, n, n, 1.5_dp, a, n, a, n, -0.5_dp, c1, n)
    call fortnum_dgemm('N', 'T', n, n, n, 1.5_dp, a, n, a, n, -0.5_dp, c2, n)
    call require(all(c1 == c2), 'dgemm')
    z1 = (0.5_dp, 1.0_dp)
    z2 = z1
    call zgemm('C', 'N', n, n, n, (1.0_dp, 2.0_dp), z, n, z, n, (0.5_dp, 0.0_dp), z1, n)
    call fortnum_zgemm('C', 'N', n, n, n, (1.0_dp, 2.0_dp), z, n, z, n, (0.5_dp, 0.0_dp), z2, n)
    call require(all(z1 == z2), 'zgemm')
    call require(dasum(n, a, 2) == fortnum_dasum(n, a, 2), 'dasum')
    zv1(:, 1) = (1.0_dp, -1.0_dp)
    zv2(:, 1) = zv1(:, 1)
    call zgemv('C', n, n, (2.0_dp, 0.5_dp), z, n, z(:, 2), 1, (0.5_dp, 0.0_dp), zv1, 1)
    call fortnum_zgemv('C', n, n, (2.0_dp, 0.5_dp), z, n, z(:, 2), 1, (0.5_dp, 0.0_dp), zv2, 1)
    call require(all(zv1(:, 1) == zv2(:, 1)), 'zgemv')

    a1 = a
    a2 = a
    v1(:, 1:2) = b1
    v2(:, 1:2) = b1
    call dgesv(n, 2, a1, n, ipiv1, v1, n, info1)
    call fortnum_dgesv(n, 2, a2, n, ipiv2, v2, n, info2)
    call require(info1 == 0 .and. info1 == info2 .and. all(a1 == a2) .and. all(ipiv1 == ipiv2) .and. &
        all(v1(:, 1:2) == v2(:, 1:2)), 'dgesv')

    ! band with kl = 2, ku = 2 stored for factorization (2 kl + ku + 1 rows)
    ab1 = 0
    do j = 1, n
        do i = max(1, j - 2), min(n, j + 2)
            ab1(5 + i - j, j) = z(i, j)
        end do
    end do
    ab2 = ab1
    zb1(:, 1) = z(:, 1)
    zb2 = zb1
    call zgbsv(n, 2, 2, 1, ab1, 7, ipiv1, zb1, n, info1)
    call fortnum_zgbsv(n, 2, 2, 1, ab2, 7, ipiv2, zb2, n, info2)
    call require(info1 == 0 .and. info1 == info2 .and. all(ab1 == ab2) .and. all(zb1 == zb2), 'zgbsv')

    a1 = matmul(transpose(a), a)
    do i = 1, n
        a1(i, i) = a1(i, i) + 1
    end do
    a2 = a1
    v1(:, 1:2) = b1
    v2(:, 1:2) = b1
    call dpotrf('U', n, a1, n, info1)
    call dpotrs('U', n, 2, a1, n, v1, n, info1)
    call fortnum_dpotrf('U', n, a2, n, info2)
    call fortnum_dpotrs('U', n, 2, a2, n, v2, n, info2)
    call require(info1 == 0 .and. info2 == 0 .and. all(a1 == a2) .and. all(v1(:, 1:2) == v2(:, 1:2)), &
        'dpotrf and dpotrs')

    a1 = a
    a2 = a
    b2 = b1
    v1(:, 1:2) = b1
    call dgels('T', n, n - 2, 2, a1, n, v1, n, work, size(work), info1)
    call fortnum_dgels('T', n, n - 2, 2, a2, n, b2, n, work, size(work), info2)
    call require(info1 == 0 .and. info2 == 0 .and. all(v1(:, 1:2) == b2), 'dgels')

    a1 = a + transpose(a)
    a2 = a1
    call dsyev('V', 'L', n, a1, n, wr1, work, size(work), info1)
    call fortnum_dsyev('V', 'L', n, a2, n, wr2, work, size(work), info2)
    call require(info1 == 0 .and. info2 == 0 .and. all(a1 == a2) .and. all(wr1 == wr2), 'dsyev')

    a1 = a
    a2 = a
    call dgeev('V', 'V', n, a1, n, wr1, wi1, u1, n, v1, n, work, size(work), info1)
    call fortnum_dgeev('V', 'V', n, a2, n, wr2, wi2, u2, n, v2, n, work, size(work), info2)
    call require(info1 == 0 .and. info2 == 0 .and. all(wr1 == wr2) .and. all(wi1 == wi2) .and. &
        all(u1 == u2) .and. all(v1 == v2), 'dgeev')

    z1 = z
    z2 = z
    call zgeev('V', 'V', n, z1, n, zw1, zu1, n, zv1, n, zwork, size(zwork), rwork, info1)
    call fortnum_zgeev('V', 'V', n, z2, n, zw2, zu2, n, zv2, n, zwork, size(zwork), rwork, info2)
    call require(info1 == 0 .and. info2 == 0 .and. all(zw1 == zw2) .and. all(zu1 == zu2) .and. &
        all(zv1 == zv2), 'zgeev')

    a1 = a
    a2 = a
    a1(:, n) = a1(:, 1)
    a2(:, n) = a2(:, 1)
    v1(:, 1:2) = b1
    v2(:, 1:2) = b1
    call dgelss(n, n, 2, a1, n, v1, n, s1, 1.0e-10_dp, rank1, work, size(work), info1)
    call fortnum_dgelss(n, n, 2, a2, n, v2, n, s2, 1.0e-10_dp, rank2, work, size(work), info2)
    call require(info1 == 0 .and. info2 == 0 .and. rank1 == n - 1 .and. rank1 == rank2 .and. &
        all(s1 == s2) .and. all(v1(:, 1:2) == v2(:, 1:2)), 'dgelss')
    print '(a)', 'standard BLAS and LAPACK names forward to fortnum'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

end program test_fortnum_lapack77
