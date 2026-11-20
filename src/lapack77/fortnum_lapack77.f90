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

real(dp) function dasum(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dasum
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: x(*)
    dasum = fortnum_dasum(n, x, incx)
end function dasum

real(dp) function dzasum(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dzasum
    integer, intent(in) :: n, incx
    complex(dp), intent(in) :: x(*)
    dzasum = fortnum_dzasum(n, x, incx)
end function dzasum

real(dp) function ddot(n, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_ddot
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: x(*), y(*)
    ddot = fortnum_ddot(n, x, incx, y, incy)
end function ddot

real(dp) function dnrm2(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dnrm2
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: x(*)
    dnrm2 = fortnum_dnrm2(n, x, incx)
end function dnrm2

real(dp) function dznrm2(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dznrm2
    integer, intent(in) :: n, incx
    complex(dp), intent(in) :: x(*)
    dznrm2 = fortnum_dznrm2(n, x, incx)
end function dznrm2

integer function idamax(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_idamax
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: x(*)
    idamax = fortnum_idamax(n, x, incx)
end function idamax

integer function izamax(n, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_izamax
    integer, intent(in) :: n, incx
    complex(dp), intent(in) :: x(*)
    izamax = fortnum_izamax(n, x, incx)
end function izamax

subroutine daxpy(n, alpha, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_daxpy
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: alpha, x(*)
    real(dp), intent(inout) :: y(*)
    call fortnum_daxpy(n, alpha, x, incx, y, incy)
end subroutine daxpy

subroutine zaxpy(n, alpha, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zaxpy
    integer, intent(in) :: n, incx, incy
    complex(dp), intent(in) :: alpha, x(*)
    complex(dp), intent(inout) :: y(*)
    call fortnum_zaxpy(n, alpha, x, incx, y, incy)
end subroutine zaxpy

subroutine dscal(n, alpha, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dscal
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: alpha
    real(dp), intent(inout) :: x(*)
    call fortnum_dscal(n, alpha, x, incx)
end subroutine dscal

subroutine zscal(n, alpha, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zscal
    integer, intent(in) :: n, incx
    complex(dp), intent(in) :: alpha
    complex(dp), intent(inout) :: x(*)
    call fortnum_zscal(n, alpha, x, incx)
end subroutine zscal

subroutine zdscal(n, alpha, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zdscal
    integer, intent(in) :: n, incx
    real(dp), intent(in) :: alpha
    complex(dp), intent(inout) :: x(*)
    call fortnum_zdscal(n, alpha, x, incx)
end subroutine zdscal

subroutine dcopy(n, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dcopy
    integer, intent(in) :: n, incx, incy
    real(dp), intent(in) :: x(*)
    real(dp), intent(inout) :: y(*)
    call fortnum_dcopy(n, x, incx, y, incy)
end subroutine dcopy

subroutine zcopy(n, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zcopy
    integer, intent(in) :: n, incx, incy
    complex(dp), intent(in) :: x(*)
    complex(dp), intent(inout) :: y(*)
    call fortnum_zcopy(n, x, incx, y, incy)
end subroutine zcopy

subroutine dswap(n, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dswap
    integer, intent(in) :: n, incx, incy
    real(dp), intent(inout) :: x(*), y(*)
    call fortnum_dswap(n, x, incx, y, incy)
end subroutine dswap

subroutine zswap(n, x, incx, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zswap
    integer, intent(in) :: n, incx, incy
    complex(dp), intent(inout) :: x(*), y(*)
    call fortnum_zswap(n, x, incx, y, incy)
end subroutine zswap

subroutine drot(n, x, incx, y, incy, c, s)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_drot
    integer, intent(in) :: n, incx, incy
    real(dp), intent(inout) :: x(*), y(*)
    real(dp), intent(in) :: c, s
    call fortnum_drot(n, x, incx, y, incy, c, s)
end subroutine drot

subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dgemv
    character, intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    real(dp), intent(inout) :: y(*)
    call fortnum_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine dgemv

subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zgemv
    character, intent(in) :: trans
    integer, intent(in) :: m, n, lda, incx, incy
    complex(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    complex(dp), intent(inout) :: y(*)
    call fortnum_zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine zgemv

subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dsymv
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda, incx, incy
    real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    real(dp), intent(inout) :: y(*)
    call fortnum_dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine dsymv

subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zhemv
    character, intent(in) :: uplo
    integer, intent(in) :: n, lda, incx, incy
    complex(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
    complex(dp), intent(inout) :: y(*)
    call fortnum_zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
end subroutine zhemv

subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dtrsv
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: n, lda, incx
    real(dp), intent(in) :: a(lda, *)
    real(dp), intent(inout) :: x(*)
    call fortnum_dtrsv(uplo, trans, diag, n, a, lda, x, incx)
end subroutine dtrsv

subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_ztrsv
    character, intent(in) :: uplo, trans, diag
    integer, intent(in) :: n, lda, incx
    complex(dp), intent(in) :: a(lda, *)
    complex(dp), intent(inout) :: x(*)
    call fortnum_ztrsv(uplo, trans, diag, n, a, lda, x, incx)
end subroutine ztrsv

subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_dger
    integer, intent(in) :: m, n, incx, incy, lda
    real(dp), intent(in) :: alpha, x(*), y(*)
    real(dp), intent(inout) :: a(lda, *)
    call fortnum_dger(m, n, alpha, x, incx, y, incy, a, lda)
end subroutine dger

subroutine zgeru(m, n, alpha, x, incx, y, incy, a, lda)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zgeru
    integer, intent(in) :: m, n, incx, incy, lda
    complex(dp), intent(in) :: alpha, x(*), y(*)
    complex(dp), intent(inout) :: a(lda, *)
    call fortnum_zgeru(m, n, alpha, x, incx, y, incy, a, lda)
end subroutine zgeru

subroutine zgerc(m, n, alpha, x, incx, y, incy, a, lda)
    use fortnum_kinds, only: dp
    use fortnum_blas12, only: fortnum_zgerc
    integer, intent(in) :: m, n, incx, incy, lda
    complex(dp), intent(in) :: alpha, x(*), y(*)
    complex(dp), intent(inout) :: a(lda, *)
    call fortnum_zgerc(m, n, alpha, x, incx, y, incy, a, lda)
end subroutine zgerc

subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_ztrsm
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(dp), intent(in) :: alpha, a(lda, *)
    complex(dp), intent(inout) :: b(ldb, *)
    call fortnum_ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine ztrsm

subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dtrmm
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    real(dp), intent(in) :: alpha, a(lda, *)
    real(dp), intent(inout) :: b(ldb, *)
    call fortnum_dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine dtrmm

subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_ztrmm
    character, intent(in) :: side, uplo, transa, diag
    integer, intent(in) :: m, n, lda, ldb
    complex(dp), intent(in) :: alpha, a(lda, *)
    complex(dp), intent(inout) :: b(ldb, *)
    call fortnum_ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
end subroutine ztrmm

subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dsyrk
    character, intent(in) :: uplo, trans
    integer, intent(in) :: n, k, lda, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *)
    real(dp), intent(inout) :: c(ldc, *)
    call fortnum_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
end subroutine dsyrk

subroutine dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dsyr2k
    character, intent(in) :: uplo, trans
    integer, intent(in) :: n, k, lda, ldb, ldc
    real(dp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
    real(dp), intent(inout) :: c(ldc, *)
    call fortnum_dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
end subroutine dsyr2k

subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_zherk
    character, intent(in) :: uplo, trans
    integer, intent(in) :: n, k, lda, ldc
    real(dp), intent(in) :: alpha, beta
    complex(dp), intent(in) :: a(lda, *)
    complex(dp), intent(inout) :: c(ldc, *)
    call fortnum_zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
end subroutine zherk

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

subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol, only: fortnum_dposv
    character, intent(in) :: uplo
    integer, intent(in) :: n, nrhs, lda, ldb
    real(dp), intent(inout) :: a(lda, *), b(ldb, *)
    integer, intent(out) :: info
    call fortnum_dposv(uplo, n, nrhs, a, lda, b, ldb, info)
end subroutine dposv

subroutine zposv(uplo, n, nrhs, a, lda, b, ldb, info)
    use fortnum_kinds, only: dp
    use fortnum_lapack_chol, only: fortnum_zposv
    character, intent(in) :: uplo
    integer, intent(in) :: n, nrhs, lda, ldb
    complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
    integer, intent(out) :: info
    call fortnum_zposv(uplo, n, nrhs, a, lda, b, ldb, info)
end subroutine zposv

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
