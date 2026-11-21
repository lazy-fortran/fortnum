! LU factorization with partial pivoting, dense and banded, real or complex
! through the macros of fortnum_lapack_lu.F90. Arguments and info codes are
! those of reference LAPACK; an illegal argument returns info = -position.

    !> A = P L U with partial pivoting through the recursive GETRF2, so that
    !> also the panel work is TRSM and GEMM.
    subroutine GETRF(m, n, a, lda, ipiv, info)
        integer, intent(in) :: m, n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, m)) then
            info = -4
        end if
        if (info /= 0 .or. m == 0 .or. n == 0) return
        call GETRF2(m, n, a, lda, ipiv, info)
    end subroutine GETRF

    !> Recursive LU with partial pivoting (xGETRF2): the left half of the
    !> columns recursively, then LASWP, TRSM and one GEMM update the right
    !> half, which is factored recursively in turn; panels of at most
    !> GETRF_LEAF columns are factored unblocked.
    recursive subroutine GETRF2(m, n, a, lda, ipiv, info)
        integer, intent(in) :: m, n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        integer :: n1, n2, iinfo, mn, i
        info = 0
        if (m == 0 .or. n == 0) return
        mn = min(m, n)
        if (m == 1) then
            ipiv(1) = 1
            if (a(1, 1) == ZERO) info = 1
            return
        end if
        if (n <= GETRF_LEAF) then
            call GETF2(m, n, a, lda, ipiv, info)
            return
        end if
        n1 = mn/2
        n2 = n - n1
        call GETRF2(m, n1, a, lda, ipiv, iinfo)
        if (info == 0 .and. iinfo > 0) info = iinfo
        call LASWP(n2, a(1, n1 + 1), lda, 1, n1, ipiv, 1)
        call TRSM('L', 'L', 'N', 'U', n1, n2, ONE, a, lda, a(1, n1 + 1), lda)
        call GEMM('N', 'N', m - n1, n2, n1, -ONE, a(n1 + 1, 1), lda, a(1, n1 + 1), lda, ONE, &
            a(n1 + 1, n1 + 1), lda)
        call GETRF2(m - n1, n2, a(n1 + 1, n1 + 1), lda, ipiv(n1 + 1), iinfo)
        if (info == 0 .and. iinfo > 0) info = iinfo + n1
        do i = n1 + 1, mn
            ipiv(i) = ipiv(i) + n1
        end do
        call LASWP(n1, a, lda, n1 + 1, mn, ipiv, 1)
    end subroutine GETRF2

    !> Unblocked LU of an m x n panel (LAPACK xGETF2): pivot by the largest
    !> ABS1, scale by the reciprocal unless the pivot is tiny.
    subroutine GETF2(m, n, a, lda, ipiv, info)
        integer, intent(in) :: m, n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        integer :: j, p, i
        SCALAR :: t
        real(dp) :: big
        info = 0
        do j = 1, min(m, n)
            p = j
            big = ABS1(a(j, j))
            do i = j + 1, m
                if (ABS1(a(i, j)) > big) then
                    big = ABS1(a(i, j))
                    p = i
                end if
            end do
            ipiv(j) = p
            if (a(p, j) /= ZERO) then
                if (p /= j) then
                    do i = 1, n
                        t = a(j, i)
                        a(j, i) = a(p, i)
                        a(p, i) = t
                    end do
                end if
                if (j < m) then
                    if (abs(a(j, j)) >= tiny(1.0_dp)) then
                        a(j + 1:m, j) = (ONE/a(j, j))*a(j + 1:m, j)
                    else
                        a(j + 1:m, j) = a(j + 1:m, j)/a(j, j)
                    end if
                end if
            else if (info == 0) then
                info = j
            end if
            ! rank-one update of the rest of the panel
            do i = j + 1, n
                if (a(j, i) /= ZERO) a(j + 1:m, i) = a(j + 1:m, i) - a(j, i)*a(j + 1:m, j)
            end do
        end do
    end subroutine GETF2

    !> Row interchanges ipiv(k1..k2) applied to the n columns of a (xLASWP);
    !> incx < 0 applies them in reverse order.
    subroutine LASWP(n, a, lda, k1, k2, ipiv, incx)
        integer, intent(in) :: n, lda, k1, k2, incx
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(in) :: ipiv(*)
        integer :: k, p, j, j0, j1, first, last, step, ip, ip0, nth
        SCALAR :: t
        if (incx == 0 .or. n <= 0) return
        ! as reference LAPACK: positions k1.. of ipiv, walked backwards for incx < 0
        if (incx > 0) then
            first = k1
            last = k2
            step = 1
            ip0 = k1
        else
            first = k2
            last = k1
            step = -1
            ip0 = k1 + (k1 - k2)*incx
        end if
        ! all interchanges on LASWP_NB columns at a time (one interchange
        ! over all columns touches a page per column); blocks are
        ! independent, so they go to threads
        nth = thread_count(real(n, dp)*abs(k2 - k1 + 1))
        nth = min(nth, (n + LASWP_NB - 1)/LASWP_NB)
        !$omp parallel do num_threads(nth) if(nth > 1) schedule(static) private(j1, ip, k, p, j, t)
        do j0 = 1, n, LASWP_NB
            j1 = min(n, j0 + LASWP_NB - 1)
            ip = ip0
            do k = first, last, step
                p = ipiv(ip)
                if (p /= k) then
                    do j = j0, j1
                        t = a(k, j)
                        a(k, j) = a(p, j)
                        a(p, j) = t
                    end do
                end if
                ip = ip + incx
            end do
        end do
        !$omp end parallel do
    end subroutine LASWP

    !> Solve op(A) X = B with the factors of GETRF; trans N, T or C.
    subroutine GETRS(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
        character, intent(in) :: trans
        integer, intent(in) :: n, nrhs, lda, ldb
        SCALAR, intent(in) :: a(lda, *)
        integer, intent(in) :: ipiv(*)
        SCALAR, intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        info = 0
        if (index('NnTtCc', trans) == 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (nrhs < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldb < max(1, n)) then
            info = -8
        end if
        if (info /= 0 .or. n == 0 .or. nrhs == 0) return
        if (trans == 'N' .or. trans == 'n') then
            call LASWP(nrhs, b, ldb, 1, n, ipiv, 1)
            call TRSM('L', 'L', 'N', 'U', n, nrhs, ONE, a, lda, b, ldb)
            call TRSM('L', 'U', 'N', 'N', n, nrhs, ONE, a, lda, b, ldb)
        else
            call TRSM('L', 'U', trans, 'N', n, nrhs, ONE, a, lda, b, ldb)
            call TRSM('L', 'L', trans, 'U', n, nrhs, ONE, a, lda, b, ldb)
            call LASWP(nrhs, b, ldb, 1, n, ipiv, -1)
        end if
    end subroutine GETRS

    !> Solve A X = B by LU with partial pivoting (xGESV).
    subroutine GESV(n, nrhs, a, lda, ipiv, b, ldb, info)
        integer, intent(in) :: n, nrhs, lda, ldb
        SCALAR, intent(inout) :: a(lda, *), b(ldb, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        info = 0
        if (n < 0) then
            info = -1
        else if (nrhs < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        else if (ldb < max(1, n)) then
            info = -7
        end if
        if (info /= 0) return
        call GETRF(n, n, a, lda, ipiv, info)
        if (info == 0) call GETRS('N', n, nrhs, a, lda, ipiv, b, ldb, info)
    end subroutine GESV

    !> Inverse from the factors of GETRF (xGETRI); lwork = -1 returns the
    !> optimal workspace in work(1).
    subroutine GETRI(n, a, lda, ipiv, work, lwork, info)
        integer, intent(in) :: n, lda, lwork
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(in) :: ipiv(*)
        SCALAR, intent(inout) :: work(*)
        integer, intent(out) :: info
        SCALAR, allocatable :: x(:, :)
        integer :: j
        info = 0
        if (n < 0) then
            info = -1
        else if (lda < max(1, n)) then
            info = -3
        else if (lwork < max(1, n) .and. lwork /= -1) then
            info = -6
        end if
        if (info /= 0) return
        work(1) = max(1, n)
        if (lwork == -1 .or. n == 0) return
        do j = 1, n
            if (a(j, j) == ZERO) then
                info = j
                return
            end if
        end do
        allocate (x(n, n))
        x = ZERO
        do j = 1, n
            x(j, j) = ONE
        end do
        call GETRS('N', n, n, a, lda, ipiv, x, n, info)
        do j = 1, n
            a(1:n, j) = x(:, j)
        end do
    end subroutine GETRI

    !> Band LU with partial pivoting (xGBTRF/xGBTF2): A(i, j) is
    !> ab(kl + ku + 1 + i - j, j); rows 1..kl of ab take the fill-in.
    subroutine GBTRF(m, n, kl, ku, ab, ldab, ipiv, info)
        integer, intent(in) :: m, n, kl, ku, ldab
        SCALAR, intent(inout) :: ab(ldab, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        integer :: kv, j, jp, ju, km, i, c
        real(dp) :: big
        SCALAR :: t, r
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (kl < 0) then
            info = -3
        else if (ku < 0) then
            info = -4
        else if (ldab < 2*kl + ku + 1) then
            info = -6
        end if
        if (info /= 0 .or. m == 0 .or. n == 0) return
        kv = ku + kl
        ! the fill-in rows of columns ku+2 .. kv start as zero
        do j = ku + 2, min(kv, n)
            ab(kv - j + 2:kl, j) = ZERO
        end do
        ju = 1
        do j = 1, min(m, n)
            if (j + kv <= n) ab(1:kl, j + kv) = ZERO
            km = min(kl, m - j)
            jp = 1
            big = ABS1(ab(kv + 1, j))
            do i = 2, km + 1
                if (ABS1(ab(kv + i, j)) > big) then
                    big = ABS1(ab(kv + i, j))
                    jp = i
                end if
            end do
            ipiv(j) = jp + j - 1
            if (ab(kv + jp, j) /= ZERO) then
                ju = max(ju, min(j + ku + jp - 1, n))
                ! swap rows j and j+jp-1 in columns j..ju (stride ldab - 1)
                if (jp /= 1) then
                    do c = 0, ju - j
                        t = ab(kv + jp - c, j + c)
                        ab(kv + jp - c, j + c) = ab(kv + 1 - c, j + c)
                        ab(kv + 1 - c, j + c) = t
                    end do
                end if
                if (km > 0) then
                    r = ONE/ab(kv + 1, j)
                    ab(kv + 2:kv + km + 1, j) = r*ab(kv + 2:kv + km + 1, j)
                    ! rank-one update of the km x (ju - j) block to the right
                    do c = 1, ju - j
                        t = ab(kv + 1 - c, j + c)
                        if (t /= ZERO) ab(kv + 2 - c:kv + km + 1 - c, j + c) = &
                            ab(kv + 2 - c:kv + km + 1 - c, j + c) - t*ab(kv + 2:kv + km + 1, j)
                    end do
                end if
            else if (info == 0) then
                info = j
            end if
        end do
    end subroutine GBTRF

    !> Solve op(A) X = B with the band factors of GBTRF; trans N, T or C.
    subroutine GBTRS(trans, n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
        character, intent(in) :: trans
        integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
        SCALAR, intent(in) :: ab(ldab, *)
        integer, intent(in) :: ipiv(*)
        SCALAR, intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        logical :: notrans, conj
        integer :: kv, j, lm, l, r, i, kd
        SCALAR :: t, u
        info = 0
        notrans = trans == 'N' .or. trans == 'n'
        conj = trans == 'C' .or. trans == 'c'
        if (.not. (notrans .or. conj .or. trans == 'T' .or. trans == 't')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (kl < 0) then
            info = -3
        else if (ku < 0) then
            info = -4
        else if (nrhs < 0) then
            info = -5
        else if (ldab < 2*kl + ku + 1) then
            info = -7
        else if (ldb < max(1, n)) then
            info = -10
        end if
        if (info /= 0 .or. n == 0 .or. nrhs == 0) return
        kv = ku + kl
        kd = kv + 1
        if (notrans) then
            ! L: interchanges and unit lower band, column by column
            if (kl > 0) then
                do j = 1, n - 1
                    lm = min(kl, n - j)
                    l = ipiv(j)
                    if (l /= j) then
                        do r = 1, nrhs
                            t = b(l, r)
                            b(l, r) = b(j, r)
                            b(j, r) = t
                        end do
                    end if
                    do r = 1, nrhs
                        if (b(j, r) /= ZERO) b(j + 1:j + lm, r) = b(j + 1:j + lm, r) - b(j, r)*ab(kd + 1:kd + lm, j)
                    end do
                end do
            end if
            ! U: upper band with kl + ku superdiagonals
            do r = 1, nrhs
                do j = n, 1, -1
                    if (b(j, r) /= ZERO) then
                        b(j, r) = b(j, r)/ab(kd, j)
                        t = b(j, r)
                        do i = max(1, j - kv), j - 1
                            b(i, r) = b(i, r) - t*ab(kd + i - j, j)
                        end do
                    end if
                end do
            end do
        else
            ! U^T (or U^H) forward
            do r = 1, nrhs
                do j = 1, n
                    t = b(j, r)
                    do i = max(1, j - kv), j - 1
                        u = ab(kd + i - j, j)
                        if (conj) u = CONJ(u)
                        t = t - u*b(i, r)
                    end do
                    u = ab(kd, j)
                    if (conj) u = CONJ(u)
                    b(j, r) = t/u
                end do
            end do
            ! L^T (or L^H) backward with the interchanges in reverse
            if (kl > 0) then
                do j = n - 1, 1, -1
                    lm = min(kl, n - j)
                    do r = 1, nrhs
                        t = ZERO
                        do i = 1, lm
                            u = ab(kd + i, j)
                            if (conj) u = CONJ(u)
                            t = t + u*b(j + i, r)
                        end do
                        b(j, r) = b(j, r) - t
                    end do
                    l = ipiv(j)
                    if (l /= j) then
                        do r = 1, nrhs
                            t = b(l, r)
                            b(l, r) = b(j, r)
                            b(j, r) = t
                        end do
                    end if
                end do
            end if
        end if
    end subroutine GBTRS

    !> Solve A X = B for a band matrix (xGBSV).
    subroutine GBSV(n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
        integer, intent(in) :: n, kl, ku, nrhs, ldab, ldb
        SCALAR, intent(inout) :: ab(ldab, *), b(ldb, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        info = 0
        if (n < 0) then
            info = -1
        else if (kl < 0) then
            info = -2
        else if (ku < 0) then
            info = -3
        else if (nrhs < 0) then
            info = -4
        else if (ldab < 2*kl + ku + 1) then
            info = -6
        else if (ldb < max(1, n)) then
            info = -9
        end if
        if (info /= 0) return
        call GBTRF(n, n, kl, ku, ab, ldab, ipiv, info)
        if (info == 0) call GBTRS('N', n, kl, ku, nrhs, ab, ldab, ipiv, b, ldb, info)
    end subroutine GBSV
