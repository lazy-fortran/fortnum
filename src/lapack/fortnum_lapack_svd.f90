module fortnum_lapack_svd
    ! Singular value decomposition and minimum-norm least squares with
    ! reference LAPACK arguments and info codes (prefix fortnum_): DGESVD,
    ! DGESDD, DGELSS and DGELSD. The matrix is first reduced by Householder
    ! QR (of A, or of A^T when m < n) and the triangular factor is
    ! diagonalized by one-sided Jacobi rotations, which give singular values
    ! with high relative accuracy and orthogonal vectors without a
    ! bidiagonal stage. DGESDD and DGELSD return what their LAPACK
    ! namesakes return, by the same method. ZGESVD is the complex
    ! counterpart: each column pair is first given a real inner product by a
    ! unit phase, then rotated as in the real case.
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dgeqrf, fortnum_dormqr, fortnum_zgeqrf, fortnum_zunmqr
    implicit none
    private

    public :: fortnum_dgesvd, fortnum_dgesdd, fortnum_dgelss, fortnum_dgelsd, fortnum_zgesvd

    real(dp), parameter :: eps = epsilon(1.0_dp)
    !> Jacobi sweeps before giving up (convergence takes about 6 to 10).
    integer, parameter :: MAX_SWEEPS = 60

contains

    !> A = U diag(s) V^T for a real m x n matrix, s descending. jobu: 'A'
    !> all m columns of U, 'S' the first min(m, n), 'O' those written over
    !> A, 'N' none; jobvt likewise for the rows of V^T.
    subroutine fortnum_dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
        character, intent(in) :: jobu, jobvt
        integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
        integer, intent(out) :: info
        integer :: minmn, minwrk, ncu, nrv
        logical :: wntua, wntus, wntuo, wntva, wntvs, wntvo
        real(dp), allocatable :: uf(:, :), vtf(:, :)
        info = 0
        minmn = min(m, n)
        wntua = any(jobu == ['A', 'a'])
        wntus = any(jobu == ['S', 's'])
        wntuo = any(jobu == ['O', 'o'])
        wntva = any(jobvt == ['A', 'a'])
        wntvs = any(jobvt == ['S', 's'])
        wntvo = any(jobvt == ['O', 'o'])
        minwrk = max(1, 3*minmn + max(m, n), 5*minmn)
        if (.not. (wntua .or. wntus .or. wntuo .or. any(jobu == ['N', 'n']))) then
            info = -1
        else if (.not. (wntva .or. wntvs .or. wntvo .or. any(jobvt == ['N', 'n'])) .or. (wntvo .and. wntuo)) then
            info = -2
        else if (m < 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (lda < max(1, m)) then
            info = -6
        else if (ldu < 1 .or. ((wntua .or. wntus) .and. ldu < m)) then
            info = -9
        else if (ldvt < 1 .or. (wntva .and. ldvt < n) .or. (wntvs .and. ldvt < minmn)) then
            info = -11
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -13
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. minmn == 0) return
        ncu = 0
        if (wntua) ncu = m
        if (wntus .or. wntuo) ncu = minmn
        nrv = 0
        if (wntva) nrv = n
        if (wntvs .or. wntvo) nrv = minmn
        allocate (uf(m, max(1, ncu)), vtf(max(1, nrv), n))
        call svd_core(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        if (wntua .or. wntus) u(1:m, 1:ncu) = uf(:, 1:ncu)
        if (wntva .or. wntvs) vt(1:nrv, 1:n) = vtf(1:nrv, :)
        if (wntuo) a(1:m, 1:minmn) = uf(:, 1:minmn)
        if (wntvo) a(1:minmn, 1:n) = vtf(1:minmn, :)
    end subroutine fortnum_dgesvd

    !> The same decomposition with one job for both factors (LAPACK's
    !> divide-and-conquer driver): jobz 'A', 'S', 'O' (U over A when
    !> m >= n, else V^T over A) or 'N'.
    subroutine fortnum_dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
        character, intent(in) :: jobz
        integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
        integer, intent(out) :: iwork(*), info
        integer :: minmn, minwrk, ncu, nrv
        logical :: wntqa, wntqs, wntqo, wntqn
        real(dp), allocatable :: uf(:, :), vtf(:, :)
        info = 0
        minmn = min(m, n)
        wntqa = any(jobz == ['A', 'a'])
        wntqs = any(jobz == ['S', 's'])
        wntqo = any(jobz == ['O', 'o'])
        wntqn = any(jobz == ['N', 'n'])
        ! the minima LAPACK documents
        if (wntqn) then
            minwrk = 3*minmn + max(max(m, n), 7*minmn)
        else if (wntqo) then
            minwrk = 3*minmn + max(max(m, n), 5*minmn*minmn + 4*minmn)
        else if (wntqs) then
            minwrk = 4*minmn*minmn + 7*minmn
        else
            minwrk = 4*minmn*minmn + 6*minmn + max(m, n)
        end if
        minwrk = max(1, minwrk)
        if (.not. (wntqa .or. wntqs .or. wntqo .or. wntqn)) then
            info = -1
        else if (m < 0) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, m)) then
            info = -5
        else if (ldu < 1 .or. ((wntqa .or. wntqs .or. (wntqo .and. m < n)) .and. ldu < m)) then
            info = -8
        else if (ldvt < 1 .or. (wntqa .and. ldvt < n) .or. (wntqs .and. ldvt < minmn) .or. &
                (wntqo .and. m >= n .and. ldvt < n)) then
            info = -10
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -12
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. minmn == 0) return
        iwork(1) = 0
        ncu = 0
        nrv = 0
        if (wntqa) then
            ncu = m
            nrv = n
        else if (wntqs .or. wntqo) then
            ncu = minmn
            nrv = minmn
            if (wntqo .and. m >= n) nrv = n
            if (wntqo .and. m < n) ncu = m
        end if
        allocate (uf(m, max(1, ncu)), vtf(max(1, nrv), n))
        call svd_core(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        if (wntqa .or. wntqs) then
            u(1:m, 1:ncu) = uf(:, 1:ncu)
            vt(1:nrv, 1:n) = vtf(1:nrv, :)
        else if (wntqo) then
            if (m >= n) then
                a(1:m, 1:n) = uf(:, 1:n)
                vt(1:n, 1:n) = vtf(1:n, :)
            else
                u(1:m, 1:m) = uf(:, 1:m)
                a(1:m, 1:n) = vtf(1:m, :)
            end if
        end if
    end subroutine fortnum_dgesdd

    ! A (destroyed) = U diag(s) V^T: the first ncu columns of U in uf and
    ! the first nrv rows of V^T in vtf.
    subroutine svd_core(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        integer, intent(in) :: m, n, lda, ncu, nrv
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*), uf(:, :), vtf(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: at(:, :), ut(:, :), vtt(:, :)
        integer :: j
        if (m >= n) then
            call tall_svd(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        else
            ! A^T = U' S V'^T, so U = V' and V^T = U'^T
            allocate (at(n, m), ut(n, max(1, nrv)), vtt(max(1, ncu), m))
            do j = 1, m
                at(:, j) = a(j, 1:n)
            end do
            call tall_svd(n, m, at, n, s, nrv, ut, ncu, vtt, info)
            if (ncu > 0) uf(:, 1:ncu) = transpose(vtt(1:ncu, :))
            if (nrv > 0) vtf(1:nrv, :) = transpose(ut(:, 1:nrv))
        end if
    end subroutine svd_core

    ! m >= n: A = Q R, R = Ur S V^T by Jacobi, U = Q [Ur 0; 0 I]
    subroutine tall_svd(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        integer, intent(in) :: m, n, lda, ncu, nrv
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*), uf(:, :), vtf(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: tau(:), work(:), r(:, :), ur(:, :), v(:, :)
        real(dp) :: query(1)
        integer :: lw, j, ierr
        allocate (tau(n), r(n, n), ur(n, n), v(n, n))
        call fortnum_dgeqrf(m, n, a, lda, tau, query, -1, ierr)
        lw = max(int(query(1)), n, max(ncu, 1)*64)
        allocate (work(lw))
        call fortnum_dgeqrf(m, n, a, lda, tau, work, lw, ierr)
        r = 0
        do j = 1, n
            r(1:j, j) = a(1:j, j)
        end do
        call jacobi_svd(n, r, s, ur, v, info)
        if (nrv > 0) vtf(1:nrv, :) = transpose(v(:, 1:nrv))
        if (ncu > 0) then
            uf = 0
            uf(1:n, 1:min(ncu, n)) = ur(:, 1:min(ncu, n))
            do j = n + 1, ncu
                uf(j, j) = 1
            end do
            call fortnum_dormqr('L', 'N', m, ncu, n, a, lda, tau, uf, m, work, lw, ierr)
        end if
    end subroutine tall_svd

    ! One-sided Jacobi on the square r: r = ur diag(s) v^T with s
    ! descending, ur and v orthogonal; info > 0 if the sweeps do not
    ! converge.
    subroutine jacobi_svd(n, r, s, ur, v, info)
        integer, intent(in) :: n
        real(dp), intent(inout) :: r(n, n)
        real(dp), intent(out) :: s(*), ur(n, n), v(n, n)
        integer, intent(out) :: info
        real(dp) :: scl, alpha, beta, gamma, zeta, t, c, sn, tol, small
        real(dp) :: wp(n), sl(n)
        integer :: p, q, sweep, j, order(n)
        logical :: rotated
        info = 0
        v = 0
        do j = 1, n
            v(j, j) = 1
        end do
        ! a power of 2 near the largest entry: exact, and keeps the squared
        ! column norms in range
        scl = maxval(abs(r))
        if (scl > 0) then
            scl = 2.0_dp**exponent(scl)
            r = r/scl
        else
            scl = 1
        end if
        tol = eps*sqrt(real(n, dp))
        do sweep = 1, MAX_SWEEPS
            rotated = .false.
            ! a column below rounding of the largest is numerically zero: the
            ! cosine of its pairs is rounding noise and never falls below tol
            small = 0
            do j = 1, n
                small = max(small, dot_product(r(:, j), r(:, j)))
            end do
            small = (n*eps)**2*small
            do p = 1, n - 1
                do q = p + 1, n
                    alpha = dot_product(r(:, p), r(:, p))
                    beta = dot_product(r(:, q), r(:, q))
                    if (alpha <= small .or. beta <= small) cycle
                    gamma = dot_product(r(:, p), r(:, q))
                    if (abs(gamma) <= tol*sqrt(alpha)*sqrt(beta)) cycle
                    rotated = .true.
                    zeta = (beta - alpha)/(2*gamma)
                    t = sign(1.0_dp, zeta)/(abs(zeta) + sqrt(1 + zeta*zeta))
                    c = 1/sqrt(1 + t*t)
                    sn = c*t
                    wp = r(:, p)
                    r(:, p) = c*wp - sn*r(:, q)
                    r(:, q) = sn*wp + c*r(:, q)
                    wp = v(:, p)
                    v(:, p) = c*wp - sn*v(:, q)
                    v(:, q) = sn*wp + c*v(:, q)
                end do
            end do
            if (.not. rotated) exit
        end do
        if (rotated) info = n
        do j = 1, n
            sl(j) = norm2(r(:, j))
            ! never rotated against the others, so not orthogonal to them
            if (sl(j)**2 <= small) sl(j) = 0
        end do
        ! descending order
        do j = 1, n
            order(j) = j
        end do
        do j = 2, n
            p = order(j)
            q = j - 1
            do while (q >= 1)
                if (sl(order(q)) >= sl(p)) exit
                order(q + 1) = order(q)
                q = q - 1
            end do
            order(q + 1) = p
        end do
        sl = sl(order)
        s(1:n) = sl
        r = r(:, order)
        v = v(:, order)
        do j = 1, n
            if (s(j) > 0) then
                ur(:, j) = r(:, j)/s(j)
            else
                ur(:, j) = 0
            end if
        end do
        call complete_basis(n, s, ur)
        s(1:n) = s(1:n)*scl
    end subroutine jacobi_svd

    ! columns of ur for zero singular values: unit vectors orthogonalized
    ! (twice) against the columns already there
    subroutine complete_basis(n, s, ur)
        integer, intent(in) :: n
        real(dp), intent(in) :: s(*)
        real(dp), intent(inout) :: ur(n, n)
        real(dp) :: x(n), w(n), nx
        logical :: basis(n)
        integer :: j, k, i, pass
        ! w(k): squared norm of e_k after projecting out the basis so far. It
        ! sums to the number of missing columns, so its largest entry is at
        ! least 1/n: the unit vector with the largest w is always usable (a
        ! fixed threshold on e_1, e_2, ... in turn can run out of candidates
        ! when the rank is high)
        basis = s(1:n) > 0
        w = 1
        do i = 1, n
            if (basis(i)) w = w - ur(:, i)**2
        end do
        do j = 1, n
            if (basis(j)) cycle
            k = maxloc(w, 1)
            x = 0
            x(k) = 1
            do pass = 1, 2
                do i = 1, n
                    if (basis(i)) x = x - dot_product(ur(:, i), x)*ur(:, i)
                end do
            end do
            nx = norm2(x)
            ur(:, j) = x/nx
            basis(j) = .true.
            w = w - ur(:, j)**2
        end do
    end subroutine complete_basis

    !> A = U diag(s) V^H for a complex m x n matrix, s descending, with the
    !> jobs of fortnum_dgesvd. rwork (5 min(m, n)) is not used.
    subroutine fortnum_zgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
        character, intent(in) :: jobu, jobvt
        integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
        complex(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*), rwork(*)
        complex(dp), intent(out) :: u(ldu, *), vt(ldvt, *), work(*)
        integer, intent(out) :: info
        integer :: minmn, minwrk, ncu, nrv
        logical :: wntua, wntus, wntuo, wntva, wntvs, wntvo
        complex(dp), allocatable :: uf(:, :), vtf(:, :)
        info = 0
        minmn = min(m, n)
        wntua = any(jobu == ['A', 'a'])
        wntus = any(jobu == ['S', 's'])
        wntuo = any(jobu == ['O', 'o'])
        wntva = any(jobvt == ['A', 'a'])
        wntvs = any(jobvt == ['S', 's'])
        wntvo = any(jobvt == ['O', 'o'])
        minwrk = max(1, 2*minmn + max(m, n))
        if (.not. (wntua .or. wntus .or. wntuo .or. any(jobu == ['N', 'n']))) then
            info = -1
        else if (.not. (wntva .or. wntvs .or. wntvo .or. any(jobvt == ['N', 'n'])) .or. (wntvo .and. wntuo)) then
            info = -2
        else if (m < 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (lda < max(1, m)) then
            info = -6
        else if (ldu < 1 .or. ((wntua .or. wntus) .and. ldu < m)) then
            info = -9
        else if (ldvt < 1 .or. (wntva .and. ldvt < n) .or. (wntvs .and. ldvt < minmn)) then
            info = -11
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -13
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. minmn == 0) return
        ncu = 0
        if (wntua) ncu = m
        if (wntus .or. wntuo) ncu = minmn
        nrv = 0
        if (wntva) nrv = n
        if (wntvs .or. wntvo) nrv = minmn
        allocate (uf(m, max(1, ncu)), vtf(max(1, nrv), n))
        call zsvd_core(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        if (wntua .or. wntus) u(1:m, 1:ncu) = uf(:, 1:ncu)
        if (wntva .or. wntvs) vt(1:nrv, 1:n) = vtf(1:nrv, :)
        if (wntuo) a(1:m, 1:minmn) = uf(:, 1:minmn)
        if (wntvo) a(1:minmn, 1:n) = vtf(1:minmn, :)
    end subroutine fortnum_zgesvd

    ! A (destroyed) = U diag(s) V^H: the first ncu columns of U in uf and
    ! the first nrv rows of V^H in vtf.
    subroutine zsvd_core(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        integer, intent(in) :: m, n, lda, ncu, nrv
        complex(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*)
        complex(dp), intent(out) :: uf(:, :), vtf(:, :)
        integer, intent(out) :: info
        complex(dp), allocatable :: ah(:, :), ut(:, :), vtt(:, :)
        integer :: j
        if (m >= n) then
            call ztall_svd(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        else
            ! A^H = U' S V'^H, so U = V' and V^H = U'^H
            allocate (ah(n, m), ut(n, max(1, nrv)), vtt(max(1, ncu), m))
            do j = 1, m
                ah(:, j) = conjg(a(j, 1:n))
            end do
            call ztall_svd(n, m, ah, n, s, nrv, ut, ncu, vtt, info)
            if (ncu > 0) uf(:, 1:ncu) = conjg(transpose(vtt(1:ncu, :)))
            if (nrv > 0) vtf(1:nrv, :) = conjg(transpose(ut(:, 1:nrv)))
        end if
    end subroutine zsvd_core

    ! m >= n: A = Q R, R = Ur S V^H by Jacobi, U = Q [Ur 0; 0 I]
    subroutine ztall_svd(m, n, a, lda, s, ncu, uf, nrv, vtf, info)
        integer, intent(in) :: m, n, lda, ncu, nrv
        complex(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: s(*)
        complex(dp), intent(out) :: uf(:, :), vtf(:, :)
        integer, intent(out) :: info
        complex(dp), allocatable :: tau(:), work(:), r(:, :), ur(:, :), v(:, :)
        complex(dp) :: query(1)
        integer :: lw, j, ierr
        allocate (tau(n), r(n, n), ur(n, n), v(n, n))
        call fortnum_zgeqrf(m, n, a, lda, tau, query, -1, ierr)
        lw = max(int(real(query(1))), n, max(ncu, 1)*64)
        allocate (work(lw))
        call fortnum_zgeqrf(m, n, a, lda, tau, work, lw, ierr)
        r = 0
        do j = 1, n
            r(1:j, j) = a(1:j, j)
        end do
        call zjacobi_svd(n, r, s, ur, v, info)
        if (nrv > 0) vtf(1:nrv, :) = conjg(transpose(v(:, 1:nrv)))
        if (ncu > 0) then
            uf = 0
            uf(1:n, 1:min(ncu, n)) = ur(:, 1:min(ncu, n))
            do j = n + 1, ncu
                uf(j, j) = 1
            end do
            call fortnum_zunmqr('L', 'N', m, ncu, n, a, lda, tau, uf, m, work, lw, ierr)
        end if
    end subroutine ztall_svd

    ! One-sided Jacobi on the square complex r: r = ur diag(s) v^H.
    ! Multiplying column q of r and v by the unit phase conj(gamma)/|gamma|
    ! makes the pair's inner product real; the real rotation follows.
    subroutine zjacobi_svd(n, r, s, ur, v, info)
        integer, intent(in) :: n
        complex(dp), intent(inout) :: r(n, n)
        real(dp), intent(out) :: s(*)
        complex(dp), intent(out) :: ur(n, n), v(n, n)
        integer, intent(out) :: info
        real(dp) :: scl, alpha, beta, g, zeta, t, c, sn, tol
        real(dp) :: sl(n)
        complex(dp) :: gamma, ph, wp(n)
        integer :: p, q, sweep, j, order(n)
        logical :: rotated
        info = 0
        v = 0
        do j = 1, n
            v(j, j) = 1
        end do
        scl = maxval(abs(r))
        if (scl > 0) then
            scl = 2.0_dp**exponent(scl)
            r = r/scl
        else
            scl = 1
        end if
        tol = eps*sqrt(real(n, dp))
        do sweep = 1, MAX_SWEEPS
            rotated = .false.
            do p = 1, n - 1
                do q = p + 1, n
                    alpha = real(dot_product(r(:, p), r(:, p)))
                    beta = real(dot_product(r(:, q), r(:, q)))
                    if (alpha == 0 .or. beta == 0) cycle
                    gamma = dot_product(r(:, p), r(:, q))
                    g = abs(gamma)
                    if (g <= tol*sqrt(alpha)*sqrt(beta)) cycle
                    rotated = .true.
                    ph = conjg(gamma)/g
                    r(:, q) = ph*r(:, q)
                    v(:, q) = ph*v(:, q)
                    zeta = (beta - alpha)/(2*g)
                    t = sign(1.0_dp, zeta)/(abs(zeta) + sqrt(1 + zeta*zeta))
                    c = 1/sqrt(1 + t*t)
                    sn = c*t
                    wp = r(:, p)
                    r(:, p) = c*wp - sn*r(:, q)
                    r(:, q) = sn*wp + c*r(:, q)
                    wp = v(:, p)
                    v(:, p) = c*wp - sn*v(:, q)
                    v(:, q) = sn*wp + c*v(:, q)
                end do
            end do
            if (.not. rotated) exit
        end do
        if (rotated) info = n
        do j = 1, n
            sl(j) = sqrt(real(dot_product(r(:, j), r(:, j))))
        end do
        do j = 1, n
            order(j) = j
        end do
        do j = 2, n
            p = order(j)
            q = j - 1
            do while (q >= 1)
                if (sl(order(q)) >= sl(p)) exit
                order(q + 1) = order(q)
                q = q - 1
            end do
            order(q + 1) = p
        end do
        sl = sl(order)
        s(1:n) = sl
        r = r(:, order)
        v = v(:, order)
        do j = 1, n
            if (s(j) > 0) then
                ur(:, j) = r(:, j)/s(j)
            else
                ur(:, j) = 0
            end if
        end do
        call zcomplete_basis(n, s, ur)
        s(1:n) = s(1:n)*scl
    end subroutine zjacobi_svd

    subroutine zcomplete_basis(n, s, ur)
        integer, intent(in) :: n
        real(dp), intent(in) :: s(*)
        complex(dp), intent(inout) :: ur(n, n)
        complex(dp) :: x(n)
        real(dp) :: w(n), nx
        logical :: basis(n)
        integer :: j, k, i, pass
        ! w(k): squared norm of e_k after projecting out the basis so far. It
        ! sums to the number of missing columns, so its largest entry is at
        ! least 1/n: the unit vector with the largest w is always usable (a
        ! fixed threshold on e_1, e_2, ... in turn can run out of candidates
        ! when the rank is high)
        basis = s(1:n) > 0
        w = 1
        do i = 1, n
            if (basis(i)) w = w - abs(ur(:, i))**2
        end do
        do j = 1, n
            if (basis(j)) cycle
            k = maxloc(w, 1)
            x = 0
            x(k) = 1
            do pass = 1, 2
                do i = 1, n
                    if (basis(i)) x = x - dot_product(ur(:, i), x)*ur(:, i)
                end do
            end do
            nx = sqrt(real(dot_product(x, x)))
            ur(:, j) = x/nx
            basis(j) = .true.
            w = w - abs(ur(:, j))**2
        end do
    end subroutine zcomplete_basis

    !> Minimum-norm least squares X = A^+ B through the SVD: singular values
    !> below rcond times the largest count as zero (rcond < 0: machine
    !> precision); rank is how many remain. B (ldb >= max(m, n)) holds the
    !> solution in rows 1:n; for m > n and full rank its rows n+1:m hold the
    !> components of the residual. A is overwritten by the first min(m, n)
    !> rows of V^T.
    subroutine fortnum_dgelss(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, info)
        integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
        real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        real(dp), intent(out) :: s(*), work(*)
        real(dp), intent(in) :: rcond
        integer, intent(out) :: rank, info
        integer :: minwrk
        minwrk = max(1, 3*min(m, n) + max(2*min(m, n), max(m, n), nrhs))
        call lss_check(m, n, nrhs, lda, ldb, lwork, minwrk, info)
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1) return
        call lss_core(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, .true., info)
    end subroutine fortnum_dgelss

    !> The same solution as DGELSS (LAPACK's divide-and-conquer driver);
    !> A is destroyed.
    subroutine fortnum_dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, iwork, info)
        integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
        real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        real(dp), intent(out) :: s(*), work(*)
        real(dp), intent(in) :: rcond
        integer, intent(out) :: rank, iwork(*), info
        integer, parameter :: smlsiz = 25
        integer :: minwrk, minmn, nlvl
        ! the minimum LAPACK documents for its divide-and-conquer tree
        minmn = max(1, min(m, n))
        nlvl = max(int(log(real(minmn, dp)/(smlsiz + 1))/log(2.0_dp)) + 1, 0)
        minwrk = max(1, 12*minmn + 2*minmn*smlsiz + 8*minmn*nlvl + minmn*nrhs + (smlsiz + 1)**2)
        call lss_check(m, n, nrhs, lda, ldb, lwork, minwrk, info)
        if (info /= 0) return
        work(1) = minwrk
        iwork(1) = max(1, 3*minmn*nlvl + 11*minmn)
        if (lwork == -1) return
        call lss_core(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, .false., info)
    end subroutine fortnum_dgelsd

    subroutine lss_check(m, n, nrhs, lda, ldb, lwork, minwrk, info)
        integer, intent(in) :: m, n, nrhs, lda, ldb, lwork, minwrk
        integer, intent(out) :: info
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (nrhs < 0) then
            info = -3
        else if (lda < max(1, m)) then
            info = -5
        else if (ldb < max(1, m, n)) then
            info = -7
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -12
        end if
    end subroutine lss_check

    subroutine lss_core(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, keep_vt, info)
        integer, intent(in) :: m, n, nrhs, lda, ldb
        real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        real(dp), intent(out) :: s(*)
        real(dp), intent(in) :: rcond
        integer, intent(out) :: rank, info
        logical, intent(in) :: keep_vt
        real(dp), allocatable :: tau(:), work(:), r(:, :), ur(:, :), v(:, :), c(:, :), at(:, :)
        real(dp) :: thresh, query(1)
        integer :: minmn, lw, j, ierr
        info = 0
        rank = 0
        minmn = min(m, n)
        if (minmn == 0) then
            if (nrhs > 0) b(1:max(m, n), 1:nrhs) = 0
            return
        end if
        allocate (tau(minmn), r(minmn, minmn), ur(minmn, minmn), v(minmn, minmn))
        if (m >= n) then
            ! A = Q R, c = Q^T b, R = Ur S V^T: x = V S^+ Ur^T c(1:n)
            call fortnum_dgeqrf(m, n, a, lda, tau, query, -1, ierr)
            lw = max(int(query(1)), max(1, nrhs)*64, n)
            allocate (work(lw))
            call fortnum_dgeqrf(m, n, a, lda, tau, work, lw, ierr)
            if (nrhs > 0) call fortnum_dormqr('L', 'T', m, nrhs, n, a, lda, tau, b, ldb, work, lw, ierr)
            r = 0
            do j = 1, n
                r(1:j, j) = a(1:j, j)
            end do
            call jacobi_svd(n, r, s, ur, v, info)
            if (info /= 0) return
            thresh = merge(rcond, eps, rcond >= 0)*s(1)
            rank = count(s(1:n) > thresh .and. s(1:n) > 0)
            if (nrhs > 0) then
                allocate (c(n, nrhs))
                c = 0
                c(1:rank, :) = matmul(transpose(ur(:, 1:rank)), b(1:n, 1:nrhs))
                do j = 1, rank
                    c(j, :) = c(j, :)/s(j)
                end do
                b(1:n, 1:nrhs) = matmul(v(:, 1:rank), c(1:rank, :))
            end if
            if (keep_vt) a(1:n, 1:n) = transpose(v)
        else
            ! A^T = Q R, so A = R^T Q^T; R^T = Ur S W^T gives A = Ur S (Q W)^T
            ! and x = Q [W S^+ Ur^T b; 0]
            allocate (at(n, m))
            do j = 1, m
                at(:, j) = a(j, 1:n)
            end do
            call fortnum_dgeqrf(n, m, at, n, tau, query, -1, ierr)
            lw = max(int(query(1)), max(1, nrhs, m)*64)
            allocate (work(lw))
            call fortnum_dgeqrf(n, m, at, n, tau, work, lw, ierr)
            ! r = R^T, lower triangular
            r = 0
            do j = 1, m
                r(j, 1:j) = at(1:j, j)
            end do
            call jacobi_svd(m, r, s, ur, v, info)
            if (info /= 0) return
            thresh = merge(rcond, eps, rcond >= 0)*s(1)
            rank = count(s(1:m) > thresh .and. s(1:m) > 0)
            if (nrhs > 0) then
                allocate (c(n, nrhs))
                c = 0
                c(1:rank, :) = matmul(transpose(ur(:, 1:rank)), b(1:m, 1:nrhs))
                do j = 1, rank
                    c(j, :) = c(j, :)/s(j)
                end do
                c(1:m, :) = matmul(v(:, 1:rank), c(1:rank, :))
                call fortnum_dormqr('L', 'N', n, nrhs, m, at, n, tau, c, n, work, lw, ierr)
                b(1:n, 1:nrhs) = c
            end if
            if (keep_vt) then
                ! rows of (Q W)^T
                if (allocated(c)) deallocate (c)
                allocate (c(n, m))
                c = 0
                c(1:m, :) = v
                call fortnum_dormqr('L', 'N', n, m, m, at, n, tau, c, n, work, lw, ierr)
                a(1:m, 1:n) = transpose(c)
            end if
        end if
    end subroutine lss_core

end module fortnum_lapack_svd
