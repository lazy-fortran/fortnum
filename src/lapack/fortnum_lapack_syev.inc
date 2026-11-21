! Symmetric / Hermitian eigenproblems, real or complex through the macros of
! fortnum_lapack_syev.F90: tridiagonal reduction, Q, implicit QL, drivers.

    !> Reduction to real symmetric tridiagonal form Q^H A Q = T (xSYTD2 /
    !> xHETD2): d the diagonal, e the off-diagonal, reflectors in A and tau.
    subroutine SYTRD(uplo, n, a, lda, d, e, tau, work, lwork, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, lwork
        SCALAR, intent(inout) :: a(lda, *), work(*)
        real(dp), intent(out) :: d(*), e(*)
        SCALAR, intent(out) :: tau(*)
        integer, intent(out) :: info
        logical :: upper
        integer :: i, j
        SCALAR :: taui, alpha, x(max(1, n))
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        else if (lwork < 1 .and. lwork /= -1) then
            info = -9
        end if
        if (info /= 0) return
        work(1) = max(1, 32*n)
        if (lwork == -1 .or. n == 0) return
        if (upper) then
            a(n, n) = real(a(n, n), dp)
            do i = n - 1, 1, -1
                alpha = a(i, i + 1)
                call LARFG(i, alpha, a(1, i + 1), 1, taui)
                e(i) = real(alpha, dp)
                if (taui /= ZERO) then
                    a(i, i + 1) = ONE
                    ! x = taui A(1:i, 1:i) v, then x += alpha v, A -= v x^H + x v^H
                    x(1:i) = taui*hemv_upper(i, a(1:i, i + 1))
                    alpha = -0.5_dp*taui*sum(CONJ(x(1:i))*a(1:i, i + 1))
                    x(1:i) = x(1:i) + alpha*a(1:i, i + 1)
                    do j = 1, i
                        a(1:j, j) = a(1:j, j) - a(1:j, i + 1)*CONJ(x(j)) - x(1:j)*CONJ(a(j, i + 1))
                        a(j, j) = real(a(j, j), dp)
                    end do
                else
                    a(i, i) = real(a(i, i), dp)
                end if
                a(i, i + 1) = e(i)
                d(i + 1) = real(a(i + 1, i + 1), dp)
                tau(i) = taui
            end do
            d(1) = real(a(1, 1), dp)
        else
            a(1, 1) = real(a(1, 1), dp)
            do i = 1, n - 1
                alpha = a(i + 1, i)
                call LARFG(n - i, alpha, a(min(i + 2, n), i), 1, taui)
                e(i) = real(alpha, dp)
                if (taui /= ZERO) then
                    a(i + 1, i) = ONE
                    x(i + 1:n) = taui*hemv_lower(i + 1, n, a(i + 1:n, i))
                    alpha = -0.5_dp*taui*sum(CONJ(x(i + 1:n))*a(i + 1:n, i))
                    x(i + 1:n) = x(i + 1:n) + alpha*a(i + 1:n, i)
                    do j = i + 1, n
                        a(j:n, j) = a(j:n, j) - a(j:n, i)*CONJ(x(j)) - x(j:n)*CONJ(a(j, i))
                        a(j, j) = real(a(j, j), dp)
                    end do
                else
                    a(i + 1, i + 1) = real(a(i + 1, i + 1), dp)
                end if
                a(i + 1, i) = e(i)
                d(i) = real(a(i, i), dp)
                tau(i) = taui
            end do
            d(n) = real(a(n, n), dp)
        end if
    contains
        ! A(1:m, 1:m) v from the upper triangle
        function hemv_upper(m, v) result(y)
            integer, intent(in) :: m
            SCALAR, intent(in) :: v(:)
            SCALAR :: y(m)
            integer :: r, c
            y = ZERO
            do c = 1, m
                do r = 1, c - 1
                    y(r) = y(r) + a(r, c)*v(c)
                    y(c) = y(c) + CONJ(a(r, c))*v(r)
                end do
                y(c) = y(c) + real(a(c, c), dp)*v(c)
            end do
        end function hemv_upper
        ! A(lo:n, lo:n) v from the lower triangle
        function hemv_lower(lo, nn, v) result(y)
            integer, intent(in) :: lo, nn
            SCALAR, intent(in) :: v(:)
            SCALAR :: y(nn - lo + 1)
            integer :: r, c
            y = ZERO
            do c = lo, nn
                y(c - lo + 1) = y(c - lo + 1) + real(a(c, c), dp)*v(c - lo + 1)
                do r = c + 1, nn
                    y(r - lo + 1) = y(r - lo + 1) + a(r, c)*v(c - lo + 1)
                    y(c - lo + 1) = y(c - lo + 1) + CONJ(a(r, c))*v(r - lo + 1)
                end do
            end do
        end function hemv_lower
    end subroutine SYTRD

    !> Q of SYTRD as an explicit n x n matrix, in place (xORGTR / xUNGTR).
    subroutine ORGTR(uplo, n, a, lda, tau, work, lwork, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, lwork
        SCALAR, intent(inout) :: a(lda, *), work(*)
        SCALAR, intent(in) :: tau(*)
        integer, intent(out) :: info
        logical :: upper
        SCALAR, allocatable :: q(:, :), v(:, :), w(:)
        integer :: i, j, c, iinfo
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        else if (lwork < max(1, n - 1) .and. lwork /= -1) then
            info = -7
        end if
        if (info /= 0) return
        work(1) = max(1, 32*(n - 1))
        if (lwork == -1 .or. n == 0) return
        ! As reference LAPACK dorgtr: the reflectors of the lower case are
        ! those of a QR factorization of order n - 1 (unit in row i + 1), so
        ! Q(2:n, 2:n) comes from the blocked ORGQR. The upper case
        ! (Q = H(n-1) ... H(1), unit in row i, entries above) is the same
        ! with rows and columns reversed: J Q J = H'(1) ... H'(n-1) with
        ! H'(j) = J H(n-j) J.
        allocate (q(n, n))
        q = ZERO
        q(1, 1) = ONE
        if (n > 1) then
            allocate (v(n - 1, n - 1), w(n - 1))
            v = ZERO
            do j = 1, n - 1
                if (upper) then
                    i = n - j
                    ! row r of the reversed order n - 1 problem is row n - r of a
                    do c = 1, i - 1
                        v(n - c, j) = a(c, i + 1)
                    end do
                else
                    i = j
                    v(j + 1:n - 1, j) = a(j + 2:n, j)
                end if
                w(j) = tau(i)
            end do
            call ORGQR(n - 1, n - 1, n - 1, v, n - 1, w, work, -1, iinfo)
            call ORGQR(n - 1, n - 1, n - 1, v, n - 1, w, work, lwork, iinfo)
            q(2:n, 2:n) = v
        end if
        if (upper) q = q(n:1:-1, n:1:-1)
        a(1:n, 1:n) = q
    end subroutine ORGTR

    !> Eigenvalues (ascending) and, when wantz, eigenvectors of the symmetric
    !> tridiagonal matrix (d, e) by implicit QL or QR with Wilkinson shifts,
    !> as LAPACK's DSTEQR: z holds Q on entry and Q times the eigenvectors of
    !> T on exit. A subdiagonal entry is neglected when e(m)^2 <= eps^2
    !> |d(m) d(m+1)|, and each unreduced block is chased from its larger end,
    !> so graded matrices keep their small eigenvalues to high relative
    !> accuracy. info > 0 counts the subdiagonal entries left after 30 n
    !> iterations.
    subroutine STEQR(wantz, n, d, e, z, ldz, info)
        logical, intent(in) :: wantz
        integer, intent(in) :: n, ldz
        real(dp), intent(inout) :: d(*), e(*)
        SCALAR, intent(inout) :: z(ldz, *)
        integer, intent(out) :: info
        real(dp) :: eps2, safmin, ssfmin, ssfmax, anorm, scl, tst, p, g, r, s, c, f, b, rt1, rt2
        real(dp) :: cs(max(1, n)), sn(max(1, n))
        integer :: l1, l, m, lsv, lend, lendsv, i, jtot, nmaxit, ii, k
        SCALAR :: zt(max(1, n))
        info = 0
        if (n <= 1) return
        eps2 = epsilon(1.0_dp)**2
        safmin = tiny(1.0_dp)
        ! each block is scaled into [ssfmin, ssfmax] so that the squares in the
        ! deflation test neither underflow nor overflow
        ssfmin = sqrt(safmin)/eps2
        ssfmax = sqrt(1/safmin)/3
        nmaxit = 30*n
        jtot = 0
        l1 = 1
        blocks: do while (l1 <= n)
            if (l1 > 1) e(l1 - 1) = 0
            m = n
            do i = l1, n - 1
                tst = abs(e(i))
                if (tst == 0) then
                    m = i
                    exit
                end if
                if (tst <= sqrt(abs(d(i)))*sqrt(abs(d(i + 1)))*epsilon(1.0_dp)) then
                    e(i) = 0
                    m = i
                    exit
                end if
            end do
            l = l1
            lsv = l
            lend = m
            lendsv = lend
            l1 = m + 1
            if (lend == l) cycle blocks
            anorm = max(maxval(abs(d(l:lend))), maxval(abs(e(l:lend - 1))))
            if (anorm == 0) cycle blocks
            scl = 1
            if (anorm > ssfmax) scl = ssfmax/anorm
            if (anorm < ssfmin) scl = ssfmin/anorm
            if (scl /= 1) then
                d(l:lend) = scl*d(l:lend)
                e(l:lend - 1) = scl*e(l:lend - 1)
            end if
            ! chase from the end with the larger diagonal entry
            if (abs(d(lend)) < abs(d(l))) then
                lend = lsv
                l = lendsv
            end if
            if (lend > l) then
                ql: do
                    m = lend
                    do i = l, lend - 1
                        if (abs(e(i))**2 <= (eps2*abs(d(i)))*abs(d(i + 1)) + safmin) then
                            m = i
                            exit
                        end if
                    end do
                    if (m < lend) e(m) = 0
                    p = d(l)
                    if (m == l) then
                        d(l) = p
                        l = l + 1
                        if (l <= lend) cycle ql
                        exit ql
                    end if
                    if (m == l + 1) then
                        call sym2x2(d(l), e(l), d(l + 1), rt1, rt2, c, s)
                        if (wantz) call rotate_back([c], [s], z(1:n, l:l + 1))
                        d(l) = rt1
                        d(l + 1) = rt2
                        e(l) = 0
                        l = l + 2
                        if (l <= lend) cycle ql
                        exit ql
                    end if
                    if (jtot == nmaxit) exit ql
                    jtot = jtot + 1
                    g = (d(l + 1) - p)/(2*e(l))
                    r = hypot(g, 1.0_dp)
                    g = d(m) - p + (e(l)/(g + sign(r, g)))
                    s = 1
                    c = 1
                    p = 0
                    do i = m - 1, l, -1
                        f = s*e(i)
                        b = c*e(i)
                        call givens(g, f, c, s, r)
                        if (i /= m - 1) e(i + 1) = r
                        g = d(i + 1) - p
                        r = (d(i) - g)*s + 2*c*b
                        p = s*r
                        d(i + 1) = g + p
                        g = c*r - b
                        cs(i) = c
                        sn(i) = -s
                    end do
                    if (wantz) call rotate_back(cs(l:m - 1), sn(l:m - 1), z(1:n, l:m))
                    d(l) = d(l) - p
                    e(l) = g
                end do ql
            else
                qr: do
                    m = lend
                    do i = l, lend + 1, -1
                        if (abs(e(i - 1))**2 <= (eps2*abs(d(i)))*abs(d(i - 1)) + safmin) then
                            m = i
                            exit
                        end if
                    end do
                    if (m > lend) e(m - 1) = 0
                    p = d(l)
                    if (m == l) then
                        d(l) = p
                        l = l - 1
                        if (l >= lend) cycle qr
                        exit qr
                    end if
                    if (m == l - 1) then
                        call sym2x2(d(l - 1), e(l - 1), d(l), rt1, rt2, c, s)
                        if (wantz) call rotate_forward([c], [s], z(1:n, l - 1:l))
                        d(l - 1) = rt1
                        d(l) = rt2
                        e(l - 1) = 0
                        l = l - 2
                        if (l >= lend) cycle qr
                        exit qr
                    end if
                    if (jtot == nmaxit) exit qr
                    jtot = jtot + 1
                    g = (d(l - 1) - p)/(2*e(l - 1))
                    r = hypot(g, 1.0_dp)
                    g = d(m) - p + (e(l - 1)/(g + sign(r, g)))
                    s = 1
                    c = 1
                    p = 0
                    do i = m, l - 1
                        f = s*e(i)
                        b = c*e(i)
                        call givens(g, f, c, s, r)
                        if (i /= m) e(i - 1) = r
                        g = d(i) - p
                        r = (d(i + 1) - g)*s + 2*c*b
                        p = s*r
                        d(i) = g + p
                        g = c*r - b
                        cs(i) = c
                        sn(i) = s
                    end do
                    if (wantz) call rotate_forward(cs(m:l - 1), sn(m:l - 1), z(1:n, m:l))
                    d(l) = d(l) - p
                    e(l - 1) = g
                end do qr
            end if
            if (scl /= 1) then
                d(lsv:lendsv) = d(lsv:lendsv)/scl
                e(lsv:lendsv - 1) = e(lsv:lendsv - 1)/scl
            end if
            if (jtot >= nmaxit) then
                info = count(e(1:n - 1) /= 0)
                return
            end if
        end do blocks
        ! ascending order, vectors with them
        do ii = 2, n
            i = ii - 1
            k = i
            p = d(i)
            do l = ii, n
                if (d(l) < p) then
                    k = l
                    p = d(l)
                end if
            end do
            if (k /= i) then
                d(k) = d(i)
                d(i) = p
                if (wantz) then
                    zt(1:n) = z(1:n, i)
                    z(1:n, i) = z(1:n, k)
                    z(1:n, k) = zt(1:n)
                end if
            end if
        end do
        e(1:n - 1) = 0
    contains
        ! a(:, j:j+1) times the plane rotation (c(j), s(j)), j descending
        ! (LAPACK DLASR 'R', 'V', 'B')
        subroutine rotate_back(c, s, a)
            real(dp), intent(in) :: c(:), s(:)
            SCALAR, intent(inout) :: a(:, :)
            SCALAR :: t(size(a, 1))
            integer :: j
            do j = size(c), 1, -1
                t = a(:, j + 1)
                a(:, j + 1) = c(j)*t - s(j)*a(:, j)
                a(:, j) = s(j)*t + c(j)*a(:, j)
            end do
        end subroutine rotate_back

        ! the same with j ascending (DLASR 'R', 'V', 'F')
        subroutine rotate_forward(c, s, a)
            real(dp), intent(in) :: c(:), s(:)
            SCALAR, intent(inout) :: a(:, :)
            SCALAR :: t(size(a, 1))
            integer :: j
            do j = 1, size(c)
                t = a(:, j + 1)
                a(:, j + 1) = c(j)*t - s(j)*a(:, j)
                a(:, j) = s(j)*t + c(j)*a(:, j)
            end do
        end subroutine rotate_forward
    end subroutine STEQR

    !> Eigenvalues w (ascending) and, jobz V, orthonormal eigenvectors in A
    !> of a symmetric (Hermitian) matrix (xSYEV / xHEEV core), with LAPACK's
    !> scaling of matrices near the overflow or underflow threshold.
    subroutine SYEV_CORE(jobz, uplo, n, a, lda, w, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        real(dp), intent(out) :: w(*)
        integer, intent(out) :: info
        logical :: wantz, upper
        real(dp) :: anrm, rmin, rmax, sigma, safmin
        real(dp), allocatable :: e(:)
        SCALAR, allocatable :: tau(:), work(:)
        SCALAR :: q(1, 1)
        integer :: j, iinfo
        wantz = jobz == 'V' .or. jobz == 'v'
        upper = uplo == 'U' .or. uplo == 'u'
        info = 0
        if (n == 0) return
        if (n == 1) then
            w(1) = real(a(1, 1), dp)
            if (wantz) a(1, 1) = ONE
            return
        end if
        safmin = tiny(1.0_dp)
        rmin = sqrt(safmin/epsilon(1.0_dp))
        rmax = sqrt(1/(safmin/epsilon(1.0_dp)))
        anrm = 0
        do j = 1, n
            if (upper) then
                anrm = max(anrm, maxval(abs(a(1:j, j))))
            else
                anrm = max(anrm, maxval(abs(a(j:n, j))))
            end if
        end do
        sigma = 1
        if (anrm > 0 .and. anrm < rmin) then
            sigma = rmin/anrm
        else if (anrm > rmax) then
            sigma = rmax/anrm
        end if
        if (sigma /= 1) then
            do j = 1, n
                if (upper) then
                    a(1:j, j) = sigma*a(1:j, j)
                else
                    a(j:n, j) = sigma*a(j:n, j)
                end if
            end do
        end if
        allocate (e(n), tau(n), work(1))
        call SYTRD_DRIVER(uplo, n, a, lda, w, e, tau, work, 1, iinfo)
        if (wantz) then
            call TRIDIAG_VECTORS(uplo, n, a, lda, tau, w, e, info)
        else
            call STEQR(.false., n, w, e, q, 1, info)
        end if
        if (sigma /= 1) w(1:n) = w(1:n)/sigma
    end subroutine SYEV_CORE

    !> A := inv(L) A inv(L)^H or inv(U)^H A inv(U) (itype 1), or L^H A L or
    !> U A U^H (itype 2, 3), B = L L^H or U^H U from POTRF (xSYGST / xHEGST).
    subroutine SYGST(itype, uplo, n, a, lda, b, ldb, info)
        integer, intent(in) :: itype, n, lda, ldb
        character, intent(in) :: uplo
        SCALAR, intent(inout) :: a(lda, *)
        SCALAR, intent(in) :: b(ldb, *)
        integer, intent(out) :: info
        logical :: upper
        SCALAR, allocatable :: full(:, :)
        integer :: i, j
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (itype < 1 .or. itype > 3) then
            info = -1
        else if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldb < max(1, n)) then
            info = -7
        end if
        if (info /= 0 .or. n == 0) return
        allocate (full(n, n))
        do j = 1, n
            do i = 1, n
                if ((i <= j) .eqv. upper) then
                    full(i, j) = a(i, j)
                else
                    full(i, j) = CONJ(a(j, i))
                end if
            end do
            full(j, j) = real(a(j, j), dp)
        end do
        if (itype == 1) then
            if (upper) then
                call TRSM('L', 'U', TRANS_C, 'N', n, n, ONE, b, ldb, full, n)
                call TRSM('R', 'U', 'N', 'N', n, n, ONE, b, ldb, full, n)
            else
                call TRSM('L', 'L', 'N', 'N', n, n, ONE, b, ldb, full, n)
                call TRSM('R', 'L', TRANS_C, 'N', n, n, ONE, b, ldb, full, n)
            end if
        else
            if (upper) then
                call TRMM('L', 'U', 'N', 'N', n, n, ONE, b, ldb, full, n)
                call TRMM('R', 'U', TRANS_C, 'N', n, n, ONE, b, ldb, full, n)
            else
                call TRMM('L', 'L', TRANS_C, 'N', n, n, ONE, b, ldb, full, n)
                call TRMM('R', 'L', 'N', 'N', n, n, ONE, b, ldb, full, n)
            end if
        end if
        do j = 1, n
            if (upper) then
                a(1:j, j) = full(1:j, j)
            else
                a(j:n, j) = full(j:n, j)
            end if
            a(j, j) = real(full(j, j), dp)
        end do
    end subroutine SYGST

    !> Generalized problem A x = lambda B x (itype 1), A B x = lambda x (2) or
    !> B A x = lambda x (3), B positive definite (xSYGV / xHEGV core); info =
    !> n + i when B is not positive definite at order i.
    subroutine SYGV_CORE(itype, jobz, uplo, n, a, lda, b, ldb, w, info)
        integer, intent(in) :: itype, n, lda, ldb
        character, intent(in) :: jobz, uplo
        SCALAR, intent(inout) :: a(lda, *), b(ldb, *)
        real(dp), intent(out) :: w(*)
        integer, intent(out) :: info
        logical :: upper, wantz
        character :: tr
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        wantz = jobz == 'V' .or. jobz == 'v'
        call POTRF(uplo, n, b, ldb, info)
        if (info /= 0) then
            info = n + info
            return
        end if
        call SYGST(itype, uplo, n, a, lda, b, ldb, info)
        call SYEV_CORE(jobz, uplo, n, a, lda, w, info)
        if (.not. wantz .or. info /= 0) return
        if (itype == 1 .or. itype == 2) then
            ! x = inv(L)^H y or inv(U) y
            tr = merge('N', TRANS_C, upper)
            call TRSM('L', uplo, tr, 'N', n, n, ONE, b, ldb, a, lda)
        else
            ! x = L y or U^H y
            tr = merge(TRANS_C, 'N', upper)
            call TRMM('L', uplo, tr, 'N', n, n, ONE, b, ldb, a, lda)
        end if
    end subroutine SYGV_CORE
