! The BLAS and LAPACK routines of fortnum under their standard external
! names, for code that calls dgesv rather than fortnum_dgesv: link the
! library fortnum_lapack77 in place of a system BLAS and LAPACK. This is the
! only BLAS and LAPACK there is in the browser (WebAssembly). Each routine
! forwards its reference arguments unchanged.

subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    real(dp), intent(inout) :: c(ldc, *)
    call fortnum_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine dgemm

subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_zgemm
    character, intent(in) :: transa, transb
    integer, intent(in) :: m, n, k, lda, ldb, ldc
    complex(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    complex(dp), intent(inout) :: c(ldc, *)
    call fortnum_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine zgemm

subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dtrsm
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    real(dp), intent(in) :: alpha, a(lda, *)
    real(dp), intent(inout) :: b(ldb, *)
    call fortnum_dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine dtrsm

subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgesv
    integer, intent(in) :: n, nrhs, lda, ldb
    real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
end subroutine dgesv

subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgesv
    integer, intent(in) :: n, nrhs, lda, ldb
    complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
end subroutine zgesv

subroutine dgetrf(m, n, a, lda, ipiv, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgetrf
    integer, intent(in) :: m, n, lda
    real(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_dgetrf(m, n, a, lda, ipiv, info)
end subroutine dgetrf

subroutine zgetrf(m, n, a, lda, ipiv, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgetrf
    integer, intent(in) :: m, n, lda
    complex(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_zgetrf(m, n, a, lda, ipiv, info)
end subroutine zgetrf

subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgetrs
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
    real(dp), intent(in) :: a(lda, *)
    real(dp), intent(inout) :: b(ldb, *)
    integer, intent(out) :: info
    call fortnum_dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
end subroutine dgetrs

subroutine zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgetrs
    character, intent(in) :: trans
    integer, intent(in) :: n, nrhs, lda, ldb, ipiv(*)
    complex(dp), intent(in) :: a(lda, *)
    complex(dp), intent(inout) :: b(ldb, *)
    integer, intent(out) :: info
    call fortnum_zgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
end subroutine zgetrs

subroutine dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgbsv
    integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
    real(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_dgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
end subroutine dgbsv

subroutine zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgbsv
    integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
    complex(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_zgbsv(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
end subroutine zgbsv

subroutine dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgbtrf
    integer, intent(in) :: m, n, kl, ku, ldab
    real(dp), intent(inout) :: ab(ldab, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_dgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
end subroutine dgbtrf

subroutine zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgbtrf
    integer, intent(in) :: m, n, kl, ku, ldab
    complex(dp), intent(inout) :: ab(ldab, *)
    integer, intent(out) :: ipiv(*), info
    call fortnum_zgbtrf(m, n, kl, ku, ab, ldab, ipiv, info)
end subroutine zgbtrf

subroutine dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgbtrs
    character, intent(in) :: trans
    integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb, ipiv(*)
    real(dp), intent(in) :: ab(ldab, *)
    real(dp), intent(inout) :: b(ldb, *)
    integer, intent(out) :: info
    call fortnum_dgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
end subroutine dgbtrs

subroutine zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_zgbtrs
    character, intent(in) :: trans
    integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb, ipiv(*)
    complex(dp), intent(in) :: ab(ldab, *)
    complex(dp), intent(inout) :: b(ldb, *)
    integer, intent(out) :: info
    call fortnum_zgbtrs(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
end subroutine zgbtrs

subroutine dpotrf(uplo, n, a, lda, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol, only: fortnum_dpotrf
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda
    real(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: info
    call fortnum_dpotrf(uplo, n, a, lda, info)
end subroutine dpotrf

subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol, only: fortnum_dpotrs
    character, intent(in) :: uplo
    integer, intent(in) :: n, nrhs, lda, ldb
    real(dp), intent(in) :: a(lda, *)
    real(dp), intent(inout) :: b(ldb, *)
    integer, intent(out) :: info
    call fortnum_dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
end subroutine dpotrs

subroutine dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol, only: fortnum_dpstrf
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda
    real(dp), intent(inout) :: a(lda, *)
    integer, intent(out) :: piv(*), rank, info
    real(dp), intent(in) :: tol
    real(dp), intent(out) :: work(*)
    call fortnum_dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
end subroutine dpstrf

subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dgeqrf
    integer, intent(in) :: m, n, lda, lwork
    real(dp), intent(inout) :: a(lda, *)
    real(dp), intent(out) :: tau(*), work(*)
    integer, intent(out) :: info
    call fortnum_dgeqrf(m, n, a, lda, tau, work, lwork, info)
end subroutine dgeqrf

subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dormqr
    character, intent(in) :: side, trans
    integer, intent(in) :: m, n, k, lda, ldc, lwork
    real(dp), intent(in) :: a(lda, *), tau(*)
    real(dp), intent(inout) :: c(ldc, *)
    real(dp), intent(out) :: work(*)
    integer, intent(out) :: info
    call fortnum_dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
end subroutine dormqr

subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dgels
    character, intent(in) :: trans
    integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
    real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    real(dp), intent(out) :: work(*)
    integer, intent(out) :: info
    call fortnum_dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
end subroutine dgels

subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_syev, only: fortnum_dsyev
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: n, lda, lwork
    real(dp), intent(inout) :: a(lda, *)
    real(dp), intent(out) :: w(*), work(*)
    integer, intent(out) :: info
    call fortnum_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
end subroutine dsyev

subroutine dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_syev, only: fortnum_dsbevd
    character, intent(in) :: jobz, uplo
    integer, intent(in) :: n, kd, ldab, ldz, lwork, liwork
    real(dp), intent(inout) :: ab(ldab, *)
    real(dp), intent(out) :: w(*), z(ldz, *), work(*)
    integer, intent(out) :: iwork(*), info
    call fortnum_dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
end subroutine dsbevd

subroutine dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, &
        z, ldz, work, iwork, ifail, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_syev, only: fortnum_dsbevx
    character, intent(in) :: jobz, range, uplo
    integer, intent(in) :: n, kd, ldab, ldq, il, iu, ldz
    real(dp), intent(inout) :: ab(ldab, *)
    real(dp), intent(out) :: q(ldq, *), w(*), z(ldz, *), work(*)
    real(dp), intent(in) :: vl, vu, abstol
    integer, intent(out) :: m, iwork(*), ifail(*), info
    call fortnum_dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, &
        z, ldz, work, iwork, ifail, info)
end subroutine dsbevx

subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_geev, only: fortnum_dgeev
    character, intent(in) :: jobvl, jobvr
    integer, intent(in) :: n, lda, ldvl, ldvr, lwork
    real(dp), intent(inout) :: a(lda, *)
    real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
    integer, intent(out) :: info
    call fortnum_dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
end subroutine dgeev

subroutine zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_geev, only: fortnum_zgeev
    character, intent(in) :: jobvl, jobvr
    integer, intent(in) :: n, lda, ldvl, ldvr, lwork
    complex(dp), intent(inout) :: a(lda, *)
    complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
    real(dp), intent(out) :: rwork(*)
    integer, intent(out) :: info
    call fortnum_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
end subroutine zgeev

subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd, only: fortnum_dgesvd
    character, intent(in) :: jobu, jobvt
    integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
    real(dp), intent(inout) :: a(lda, *)
    real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
    integer, intent(out) :: info
    call fortnum_dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
end subroutine dgesvd

subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd, only: fortnum_dgesdd
    character, intent(in) :: jobz
    integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
    real(dp), intent(inout) :: a(lda, *)
    real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
    integer, intent(out) :: iwork(*), info
    call fortnum_dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
end subroutine dgesdd

subroutine dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd, only: fortnum_dgelss
    integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
    real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    real(dp), intent(out) :: s(*), work(*)
    real(dp), intent(in) :: rcond
    integer, intent(out) :: rank, info
    call fortnum_dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
end subroutine dgelss

subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd, only: fortnum_dgelsd
    integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
    real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    real(dp), intent(out) :: s(*), work(*)
    real(dp), intent(in) :: rcond
    integer, intent(out) :: rank, iwork(*), info
    call fortnum_dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
end subroutine dgelsd
