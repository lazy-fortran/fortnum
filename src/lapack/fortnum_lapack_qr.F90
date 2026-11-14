module fortnum_lapack_qr
    ! Householder QR on fortnum's BLAS, with reference LAPACK arguments and
    ! info codes (prefix fortnum_): DLARFG/ZLARFG, DGEQR2/ZGEQR2, DLARFT/ZLARFT,
    ! DLARFB/ZLARFB, DGEQRF/ZGEQRF, DORMQR/ZUNMQR, DORGQR/ZUNGQR, the column
    ! pivoted DGEQP3 and the least squares driver DGELS.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm
    use fortnum_blas3, only: fortnum_dtrmm, fortnum_ztrmm, fortnum_dtrsm
    implicit none
    private

    public :: fortnum_dlarfg, fortnum_dgeqr2, fortnum_dlarft, fortnum_dlarfb, fortnum_dgeqrf
    public :: fortnum_dormqr, fortnum_dorgqr, fortnum_dgeqp3, fortnum_dgels
    public :: fortnum_zlarfg, fortnum_zgeqr2, fortnum_zlarft, fortnum_zlarfb, fortnum_zgeqrf
    public :: fortnum_zunmqr, fortnum_zungqr

    !> Panel width of the blocked QR.
    integer, parameter :: QR_NB = 32

contains

#define SCALAR real(dp)
#define ZERO 0.0_dp
#define ONE 1.0_dp
#define CONJ(x) (x)
#define IMAGPART(x) 0.0_dp
#define MAKE(re, im) (re)
#define NPARTS 1
#define TRANS_C 'T'
#define LOWER_TRANS_C 't'
#define GEMM fortnum_dgemm
#define TRMM fortnum_dtrmm
#define LARFG fortnum_dlarfg
#define VNORM dvnorm
#define GEQR2 fortnum_dgeqr2
#define APPLY_LEFT dapply_left
#define LARFT fortnum_dlarft
#define LARFB fortnum_dlarfb
#define GEQRF fortnum_dgeqrf
#define ORMQR fortnum_dormqr
#define ORGQR fortnum_dorgqr
#include "fortnum_lapack_qr.inc"
#undef SCALAR
#undef ZERO
#undef ONE
#undef CONJ
#undef IMAGPART
#undef MAKE
#undef NPARTS
#undef TRANS_C
#undef LOWER_TRANS_C
#undef GEMM
#undef TRMM
#undef LARFG
#undef VNORM
#undef GEQR2
#undef APPLY_LEFT
#undef LARFT
#undef LARFB
#undef GEQRF
#undef ORMQR
#undef ORGQR
#define SCALAR complex(dp)
#define ZERO (0.0_dp, 0.0_dp)
#define ONE (1.0_dp, 0.0_dp)
#define CONJ(x) conjg(x)
#define IMAGPART(x) aimag(x)
#define MAKE(re, im) cmplx(re, im, dp)
#define NPARTS 2
#define TRANS_C 'C'
#define LOWER_TRANS_C 'c'
#define GEMM fortnum_zgemm
#define TRMM fortnum_ztrmm
#define LARFG fortnum_zlarfg
#define VNORM zvnorm
#define GEQR2 fortnum_zgeqr2
#define APPLY_LEFT zapply_left
#define LARFT fortnum_zlarft
#define LARFB fortnum_zlarfb
#define GEQRF fortnum_zgeqrf
#define ORMQR fortnum_zunmqr
#define ORGQR fortnum_zungqr
#include "fortnum_lapack_qr.inc"

    !> QR with column pivoting (DGEQP3): A P = Q R. On entry jpvt(j) /= 0
    !> marks column j as leading; on exit jpvt(j) = k when column j of A P was
    !> column k of A. Column norms are downdated as LAPACK does and recomputed
    !> when cancellation makes the downdate unreliable.
    subroutine fortnum_dgeqp3(m, n, a, lda, jpvt, tau, work, lwork, info)
        integer, intent(in) :: m, n, lda, lwork
        real(dp), intent(inout) :: a(lda, *), work(*)
        integer, intent(inout) :: jpvt(*)
        real(dp), intent(out) :: tau(*)
        integer, intent(out) :: info
        real(dp), allocatable :: vn1(:), vn2(:), t(:)
        real(dp) :: aii, temp, temp2, tol3z
        integer :: nfxd, j, i, k, pvt, itemp
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, m)) then
            info = -4
        else if (lwork < 3*n + 1 .and. lwork /= -1) then
            info = -8
        end if
        if (info /= 0) return
        work(1) = max(1, 3*n + 1)
        if (lwork == -1) return
        k = min(m, n)
        allocate (vn1(n), vn2(n), t(m))
        ! leading columns to the front
        nfxd = 0
        do j = 1, n
            if (jpvt(j) /= 0) then
                if (j /= nfxd + 1) then
                    t = a(1:m, j)
                    a(1:m, j) = a(1:m, nfxd + 1)
                    a(1:m, nfxd + 1) = t
                    jpvt(j) = jpvt(nfxd + 1)
                    jpvt(nfxd + 1) = j
                else
                    jpvt(j) = j
                end if
                nfxd = nfxd + 1
            else
                jpvt(j) = j
            end if
        end do
        tol3z = sqrt(epsilon(1.0_dp)/2)
        ! the leading columns first, without pivoting
        do i = 1, min(nfxd, k)
            call reflect(i)
        end do
        ! then the free columns, measured below the leading block
        do j = nfxd + 1, n
            vn1(j) = norm2(a(nfxd + 1:m, j))
            vn2(j) = vn1(j)
        end do
        do i = nfxd + 1, k
            pvt = i - 1 + maxloc(vn1(i:n), 1)
            if (pvt /= i) then
                t = a(1:m, pvt)
                a(1:m, pvt) = a(1:m, i)
                a(1:m, i) = t
                itemp = jpvt(pvt)
                jpvt(pvt) = jpvt(i)
                jpvt(i) = itemp
                vn1(pvt) = vn1(i)
                vn2(pvt) = vn2(i)
            end if
            call reflect(i)
            ! downdate the remaining norms; recompute where cancellation bites
            do j = i + 1, n
                if (vn1(j) == 0) cycle
                temp = 1 - (abs(a(i, j))/vn1(j))**2
                temp = max(temp, 0.0_dp)
                temp2 = temp*(vn1(j)/vn2(j))**2
                if (temp2 <= tol3z) then
                    if (i < m) then
                        vn1(j) = norm2(a(i + 1:m, j))
                        vn2(j) = vn1(j)
                    else
                        vn1(j) = 0
                        vn2(j) = 0
                    end if
                else
                    vn1(j) = vn1(j)*sqrt(temp)
                end if
            end do
        end do
    contains
        subroutine reflect(i)
            integer, intent(in) :: i
            call fortnum_dlarfg(m - i + 1, a(i, i), a(min(i + 1, m), i), 1, tau(i))
            if (i < n) then
                aii = a(i, i)
                a(i, i) = 1
                call dapply_left(m - i + 1, n - i, a(i, i), tau(i), a(i, i + 1), lda)
                a(i, i) = aii
            end if
        end subroutine reflect
    end subroutine fortnum_dgeqp3

    !> Least squares or minimum norm solution of op(A) X = B by QR (DGELS);
    !> info = i when the triangular factor has a zero diagonal element.
    subroutine fortnum_dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
        character, intent(in) :: trans
        integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
        real(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
        integer, intent(out) :: info
        logical :: notran
        real(dp), allocatable :: tau(:), at(:, :), q(:, :)
        real(dp) :: w(1)
        integer :: mn, i, inf, rows
        info = 0
        notran = trans == 'N' .or. trans == 'n'
        mn = min(m, n)
        if (.not. (notran .or. trans == 'T' .or. trans == 't')) then
            info = -1
        else if (m < 0) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (nrhs < 0) then
            info = -4
        else if (lda < max(1, m)) then
            info = -6
        else if (ldb < max(1, m, n)) then
            info = -8
        else if (lwork < max(1, mn + max(mn, nrhs)) .and. lwork /= -1) then
            info = -10
        end if
        if (info /= 0) return
        work(1) = max(1, mn + max(mn, nrhs)*QR_NB)
        if (lwork == -1) return
        if (mn == 0 .or. nrhs == 0) then
            b(1:max(m, n), 1:nrhs) = 0
            return
        end if
        allocate (tau(mn))
        if (m >= n) then
            call fortnum_dgeqrf(m, n, a, lda, tau, w, max(1, n), inf)
            do i = 1, n
                if (a(i, i) == 0) then
                    info = i
                    return
                end if
            end do
            if (notran) then
                ! min || A x - b ||: R x = (Q^T b)(1:n)
                call fortnum_dormqr('L', 'T', m, nrhs, n, a, lda, tau, b, ldb, w, max(1, nrhs), inf)
                call fortnum_dtrsm('L', 'U', 'N', 'N', n, nrhs, 1.0_dp, a, lda, b, ldb)
            else
                ! min ||x|| with A^T x = b: R^T y = b, x = Q [y; 0]
                call fortnum_dtrsm('L', 'U', 'T', 'N', n, nrhs, 1.0_dp, a, lda, b, ldb)
                b(n + 1:m, 1:nrhs) = 0
                call fortnum_dormqr('L', 'N', m, nrhs, n, a, lda, tau, b, ldb, w, max(1, nrhs), inf)
            end if
        else
            ! m < n: factor A^T = Q R (n x m) and keep the factors in a copy
            allocate (at(n, m))
            at = transpose(a(1:m, 1:n))
            call fortnum_dgeqrf(n, m, at, n, tau, w, max(1, m), inf)
            do i = 1, m
                if (at(i, i) == 0) then
                    info = i
                    return
                end if
            end do
            rows = n
            allocate (q(rows, nrhs))
            if (notran) then
                ! min ||x|| with A x = b, A = R^T Q^T: R^T y = b, x = Q [y; 0]
                q = 0
                q(1:m, :) = b(1:m, 1:nrhs)
                call fortnum_dtrsm('L', 'U', 'T', 'N', m, nrhs, 1.0_dp, at, n, q, rows)
                call fortnum_dormqr('L', 'N', n, nrhs, m, at, n, tau, q, rows, w, max(1, nrhs), inf)
                b(1:n, 1:nrhs) = q
            else
                ! min || A^T x - b ||, A^T = Q R: R x = (Q^T b)(1:m)
                q = b(1:n, 1:nrhs)
                call fortnum_dormqr('L', 'T', n, nrhs, m, at, n, tau, q, rows, w, max(1, nrhs), inf)
                call fortnum_dtrsm('L', 'U', 'N', 'N', m, nrhs, 1.0_dp, at, n, q, rows)
                b(1:m, 1:nrhs) = q(1:m, :)
            end if
            a(1:m, 1:n) = transpose(at)
        end if
    end subroutine fortnum_dgels

end module fortnum_lapack_qr
