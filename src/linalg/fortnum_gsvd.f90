module fortnum_gsvd
    !! Generalized singular value decomposition of a real matrix pair (A, B)
    !! with A m-by-n and B p-by-n, over LAPACK DGGSVD3.
    !!
    !! On success, with r = k + l the numerical rank of [A; B]:
    !!
    !!     A = U C X^T,     B = V S X^T,     C^T C + S^T S = I_r,
    !!
    !! where U (m-by-m) and V (p-by-p) are orthogonal, X is n-by-r,
    !! C (m-by-r) and S (p-by-r) are zero except for C(i,i) = alpha(i) and
    !! S(i-k, i) = beta(i), i = 1..r, with alpha^2 + beta^2 = 1. The first k
    !! columns are the singular directions of A in the null space of B
    !! (alpha = 1, beta = 0); `gsv = alpha/beta` is the generalized singular
    !! value (huge if beta = 0, zero if alpha = 0). Values are in LAPACK order
    !! (not sorted). For B = D = diag(dscale) with D nonsingular, gsv are the
    !! ordinary singular values of A D^-1, and X^T = [0 R] Q^T is the
    !! LAPACK factor; the columns of X give the coil combination of each mode.
    !!
    !! Status contract matches fortnum_bvls (0 success): BVLS_INVALID_ARGUMENT
    !! for bad shapes or nonfinite data, BVLS_MAX_ITERATIONS when the Jacobi
    !! sweep in DGGSVD3 fails to converge. BVLS_INFEASIBLE is not used.
    use fortnum_kinds, only: dp
    use fortnum_bvls, only: BVLS_OK, BVLS_MAX_ITERATIONS, BVLS_INVALID_ARGUMENT
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: generalized_svd

    interface
        subroutine dggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, &
                beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info)
            import :: dp
            character, intent(in) :: jobu, jobv, jobq
            integer, intent(in) :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
            integer, intent(out) :: k, l, iwork(*), info
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: alpha(*), beta(*), u(ldu, *), v(ldv, *)
            real(dp), intent(out) :: q(ldq, *), work(*)
        end subroutine dggsvd3
    end interface

contains

    subroutine generalized_svd(a, b, alpha, beta, gsv, u, v, xmat, cmat, smat, &
            k, l, info)
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp), allocatable, intent(out) :: alpha(:), beta(:), gsv(:)
        real(dp), allocatable, intent(out) :: u(:, :), v(:, :), xmat(:, :)
        real(dp), allocatable, intent(out) :: cmat(:, :), smat(:, :)
        integer, intent(out) :: k, l, info
        real(dp), allocatable :: aw(:, :), bw(:, :), q(:, :), work(:), rfull(:, :)
        real(dp), allocatable :: al(:), be(:)
        integer, allocatable :: iwork(:)
        real(dp) :: wq(1)
        integer :: m, n, p, r, i

        m = size(a, 1)
        n = size(a, 2)
        p = size(b, 1)
        k = 0
        l = 0
        info = BVLS_INVALID_ARGUMENT
        if (m < 1 .or. n < 1 .or. p < 1 .or. size(b, 2) /= n) return
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return

        aw = a
        bw = b
        allocate (al(n), be(n), u(m, m), v(p, p), q(n, n), iwork(n))
        call dggsvd3('U', 'V', 'Q', m, n, p, k, l, aw, m, bw, p, al, be, u, m, v, p, &
            q, n, wq, -1, iwork, info)
        allocate (work(max(1, int(wq(1)))))
        call dggsvd3('U', 'V', 'Q', m, n, p, k, l, aw, m, bw, p, al, be, u, m, v, p, &
            q, n, work, size(work), iwork, info)
        if (info /= 0) then
            info = BVLS_MAX_ITERATIONS
            return
        end if
        info = BVLS_OK
        r = k + l

        allocate (alpha(r), beta(r), gsv(r), cmat(m, r), smat(p, r), xmat(n, r))
        alpha = al(1:r)
        beta = be(1:r)
        cmat = 0.0_dp
        smat = 0.0_dp
        do i = 1, r
            if (i <= m) cmat(i, i) = alpha(i)
            if (i > k) smat(i - k, i) = beta(i)
            if (beta(i) > 0.0_dp) then
                gsv(i) = alpha(i)/beta(i)
            else
                gsv(i) = huge(1.0_dp)
            end if
        end do

        ! X^T = [0 R] Q^T; R is stored in A, and in B when m < k + l.
        allocate (rfull(r, n))
        rfull = 0.0_dp
        if (m - r >= 0) then
            rfull(1:r, n - r + 1:n) = aw(1:r, n - r + 1:n)
        else
            rfull(1:m, n - r + 1:n) = aw(1:m, n - r + 1:n)
            rfull(m + 1:r, m + n - r + 1:n) = bw(m - k + 1:l, n + m - r + 1:n)
        end if
        xmat = transpose(matmul(rfull, transpose(q)))
    end subroutine generalized_svd

end module fortnum_gsvd
