! Householder QR, real or complex through the macros of fortnum_lapack_qr.F90.
! Reflectors H = I - tau v v^H with v(1) = 1 stored below the diagonal, as
! in reference LAPACK; Q = H(1) H(2) ... H(k).

    !> Elementary reflector (xLARFG): H^H [alpha; x] = [beta; 0] with beta
    !> real; on exit alpha = beta, x = v(2:n), tau. tau = 0 when nothing to do.
    subroutine LARFG(n, alpha, x, incx, tau)
        integer, intent(in) :: n, incx
        SCALAR, intent(inout) :: alpha, x(*)
        SCALAR, intent(out) :: tau
        real(dp) :: xnorm, beta, alphr, safmin, rsafmn
        integer :: i, knt
        SCALAR :: scal
        tau = ZERO
        if (n <= 0) return
        xnorm = 0
        if (n > 1) xnorm = VNORM(n - 1, x, incx)
        alphr = real(alpha, dp)
        if (xnorm == 0 .and. IMAGPART(alpha) == 0) return
        beta = -sign(hypot3(alpha, xnorm), alphr)
        safmin = tiny(1.0_dp)/epsilon(1.0_dp)
        rsafmn = 1/safmin
        knt = 0
        ! rescale while beta is tiny, as LAPACK does
        do while (abs(beta) < safmin .and. knt < 20)
            knt = knt + 1
            do i = 1, n - 1
                x(1 + (i - 1)*incx) = rsafmn*x(1 + (i - 1)*incx)
            end do
            beta = beta*rsafmn
            alpha = alpha*rsafmn
            alphr = real(alpha, dp)
            xnorm = VNORM(n - 1, x, incx)
            beta = -sign(hypot3(alpha, xnorm), alphr)
        end do
        tau = MAKE((beta - alphr)/beta, -IMAGPART(alpha)/beta)
        scal = ONE/(alpha - beta)
        do i = 1, n - 1
            x(1 + (i - 1)*incx) = scal*x(1 + (i - 1)*incx)
        end do
        do i = 1, knt
            beta = beta*safmin
        end do
        alpha = beta
    contains
        real(dp) function hypot3(z, r)
            SCALAR, intent(in) :: z
            real(dp), intent(in) :: r
            hypot3 = hypot(hypot(real(z, dp), IMAGPART(z)), r)
        end function hypot3
    end subroutine LARFG

    ! 2-norm of a vector with stride, without overflow
    real(dp) function VNORM(n, x, incx)
        integer, intent(in) :: n, incx
        SCALAR, intent(in) :: x(*)
        real(dp) :: scale, ssq, v
        integer :: i, part
        scale = 0
        ssq = 1
        do i = 1, n
            do part = 1, NPARTS
                if (part == 1) then
                    v = abs(real(x(1 + (i - 1)*incx), dp))
                else
                    v = abs(IMAGPART(x(1 + (i - 1)*incx)))
                end if
                if (v /= 0) then
                    if (scale < v) then
                        ssq = 1 + ssq*(scale/v)**2
                        scale = v
                    else
                        ssq = ssq + (v/scale)**2
                    end if
                end if
            end do
        end do
        VNORM = scale*sqrt(ssq)
    end function VNORM

    !> Unblocked QR (xGEQR2).
    subroutine GEQR2(m, n, a, lda, tau, info)
        integer, intent(in) :: m, n, lda
        SCALAR, intent(inout) :: a(lda, *)
        SCALAR, intent(out) :: tau(*)
        integer, intent(out) :: info
        integer :: i, k
        SCALAR :: aii
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, m)) then
            info = -4
        end if
        if (info /= 0) return
        k = min(m, n)
        do i = 1, k
            call LARFG(m - i + 1, a(i, i), a(min(i + 1, m), i), 1, tau(i))
            if (i < n) then
                aii = a(i, i)
                a(i, i) = ONE
                call APPLY_LEFT(m - i + 1, n - i, a(i, i), CONJ(tau(i)), a(i, i + 1), lda)
                a(i, i) = aii
            end if
        end do
    end subroutine GEQR2

    ! C := (I - tau v v^H) C for an m x n block C
    subroutine APPLY_LEFT(m, n, v, tau, c, ldc)
        integer, intent(in) :: m, n, ldc
        SCALAR, intent(in) :: v(*), tau
        SCALAR, intent(inout) :: c(ldc, *)
        SCALAR :: w
        integer :: j
        if (tau == ZERO) return
        do j = 1, n
            w = sum(CONJ(v(1:m))*c(1:m, j))
            c(1:m, j) = c(1:m, j) - tau*w*v(1:m)
        end do
    end subroutine APPLY_LEFT

    !> Triangular factor T of the block reflector H = I - V T V^H of k
    !> forward column reflectors (xLARFT, direct F, storev C).
    subroutine LARFT(n, k, v, ldv, tau, t, ldt)
        integer, intent(in) :: n, k, ldv, ldt
        SCALAR, intent(in) :: v(ldv, *), tau(*)
        SCALAR, intent(out) :: t(ldt, *)
        integer :: i, j
        SCALAR :: vij
        do i = 1, k
            t(1:k, i) = ZERO
            if (tau(i) == ZERO) cycle
            ! t(1:i-1, i) = -tau(i) V(i:n, 1:i-1)^H V(i:n, i), V unit lower
            do j = 1, i - 1
                vij = CONJ(v(i, j))
                if (i < n) vij = vij + sum(CONJ(v(i + 1:n, j))*v(i + 1:n, i))
                t(j, i) = -tau(i)*vij
            end do
            ! t(1:i-1, i) = T(1:i-1, 1:i-1) t(1:i-1, i)
            do j = 1, i - 1
                t(j, i) = sum(t(j, j:i - 1)*t(j:i - 1, i))
            end do
            t(i, i) = tau(i)
        end do
    end subroutine LARFT

    !> C := H C, H^H C, C H or C H^H for the block reflector of LARFT
    !> (xLARFB, forward, columnwise): GEMM and TRMM on the unit lower V.
    subroutine LARFB(side, trans, m, n, k, v, ldv, t, ldt, c, ldc)
        character, intent(in) :: side, trans
        integer, intent(in) :: m, n, k, ldv, ldt, ldc
        SCALAR, intent(in) :: v(ldv, *), t(ldt, *)
        SCALAR, intent(inout) :: c(ldc, *)
        SCALAR, allocatable :: w(:, :), vf(:, :)
        character :: opt
        integer :: j, nv
        if (m <= 0 .or. n <= 0 .or. k <= 0) return
        ! H = I - V T V^H, H^H = I - V T^H V^H
        opt = merge('N', TRANS_C, trans == 'N' .or. trans == 'n')
        nv = merge(m, n, side == 'L' .or. side == 'l')
        ! V as a dense nv x k matrix with its unit diagonal and zeros above
        allocate (vf(nv, k))
        do j = 1, k
            vf(1:j - 1, j) = ZERO
            vf(j, j) = ONE
            vf(j + 1:nv, j) = v(j + 1:nv, j)
        end do
        if (side == 'L' .or. side == 'l') then
            ! W = V^H C, W := op(T) W, C := C - V W
            allocate (w(k, n))
            call GEMM(TRANS_C, 'N', k, n, m, ONE, vf, nv, c, ldc, ZERO, w, k)
            call TRMM('L', 'U', opt, 'N', k, n, ONE, t, ldt, w, k)
            call GEMM('N', 'N', m, n, k, -ONE, vf, nv, w, k, ONE, c, ldc)
        else
            ! W = C V, W := W op(T), C := C - W V^H
            allocate (w(m, k))
            call GEMM('N', 'N', m, k, n, ONE, c, ldc, vf, nv, ZERO, w, m)
            call TRMM('R', 'U', opt, 'N', m, k, ONE, t, ldt, w, m)
            call GEMM('N', TRANS_C, m, n, k, -ONE, w, m, vf, nv, ONE, c, ldc)
        end if
    end subroutine LARFB

    !> Blocked QR (xGEQRF): panels of QR_NB columns by GEQR2, trailing matrix
    !> by the block reflector. lwork = -1 returns the workspace query.
    subroutine GEQRF(m, n, a, lda, tau, work, lwork, info)
        integer, intent(in) :: m, n, lda, lwork
        SCALAR, intent(inout) :: a(lda, *)
        SCALAR, intent(out) :: tau(*)
        SCALAR, intent(inout) :: work(*)
        integer, intent(out) :: info
        SCALAR :: t(QR_NB, QR_NB)
        integer :: i, ib, k, iinfo
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, m)) then
            info = -4
        else if (lwork < max(1, n) .and. lwork /= -1) then
            info = -7
        end if
        if (info /= 0) return
        work(1) = max(1, n*QR_NB)
        if (lwork == -1) return
        k = min(m, n)
        do i = 1, k, QR_NB
            ib = min(QR_NB, k - i + 1)
            call GEQR2(m - i + 1, ib, a(i, i), lda, tau(i), iinfo)
            if (i + ib <= n) then
                call LARFT(m - i + 1, ib, a(i, i), lda, tau(i), t, QR_NB)
                call LARFB('L', TRANS_C, m - i + 1, n - i - ib + 1, ib, a(i, i), lda, t, QR_NB, &
                    a(i, i + ib), lda)
            end if
        end do
    end subroutine GEQRF

    !> C := op(Q) C or C op(Q), Q from GEQRF (DORMQR / ZUNMQR).
    subroutine ORMQR(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
        character, intent(in) :: side, trans
        integer, intent(in) :: m, n, k, lda, ldc, lwork
        SCALAR, intent(in) :: a(lda, *), tau(*)
        SCALAR, intent(inout) :: c(ldc, *), work(*)
        integer, intent(out) :: info
        logical :: left, notran, forward
        integer :: nq, i, ib, first, last, step
        SCALAR :: t(QR_NB, QR_NB)
        info = 0
        left = side == 'L' .or. side == 'l'
        notran = trans == 'N' .or. trans == 'n'
        nq = merge(m, n, left)
        if (.not. (left .or. side == 'R' .or. side == 'r')) then
            info = -1
        else if (.not. (notran .or. trans == TRANS_C .or. trans == LOWER_TRANS_C)) then
            info = -2
        else if (m < 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (k < 0 .or. k > nq) then
            info = -5
        else if (lda < max(1, nq)) then
            info = -7
        else if (ldc < max(1, m)) then
            info = -10
        else if (lwork < max(1, merge(n, m, left)) .and. lwork /= -1) then
            info = -12
        end if
        if (info /= 0) return
        work(1) = max(1, merge(n, m, left)*QR_NB)
        if (lwork == -1 .or. m == 0 .or. n == 0 .or. k == 0) return
        ! Q = H1 ... Hk: Q C applies Hk first; Q^H C applies H1 first
        forward = left .neqv. notran
        if (forward) then
            first = 1
            last = k
            step = QR_NB
        else
            first = ((k - 1)/QR_NB)*QR_NB + 1
            last = 1
            step = -QR_NB
        end if
        do i = first, last, step
            ib = min(QR_NB, k - i + 1)
            call LARFT(nq - i + 1, ib, a(i, i), lda, tau(i), t, QR_NB)
            if (left) then
                call LARFB('L', trans, m - i + 1, n, ib, a(i, i), lda, t, QR_NB, c(i, 1), ldc)
            else
                call LARFB('R', trans, m, n - i + 1, ib, a(i, i), lda, t, QR_NB, c(1, i), ldc)
            end if
        end do
    end subroutine ORMQR

    !> The first n columns of Q = H1 ... Hk from GEQRF, in place
    !> (DORGQR / ZUNGQR).
    subroutine ORGQR(m, n, k, a, lda, tau, work, lwork, info)
        integer, intent(in) :: m, n, k, lda, lwork
        SCALAR, intent(inout) :: a(lda, *), work(*)
        SCALAR, intent(in) :: tau(*)
        integer, intent(out) :: info
        SCALAR, allocatable :: q(:, :), v(:, :)
        SCALAR :: t(QR_NB, QR_NB)
        integer :: i, ib, j
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0 .or. n > m) then
            info = -2
        else if (k < 0 .or. k > n) then
            info = -3
        else if (lda < max(1, m)) then
            info = -5
        else if (lwork < max(1, n) .and. lwork /= -1) then
            info = -8
        end if
        if (info /= 0) return
        work(1) = max(1, n*QR_NB)
        if (lwork == -1 .or. n == 0) return
        allocate (q(m, n), v(m, max(k, 1)))
        v(:, 1:k) = a(1:m, 1:k)
        q = ZERO
        do j = 1, n
            q(j, j) = ONE
        end do
        ! Q = H1 ... Hk applied to the identity from the last block on, as
        ! reference LAPACK dorgqr: columns before block i are still unit
        ! vectors with zeros in rows i..m, so block i acts on q(i:m, i:n)
        ! only (4/3 n^3 instead of 4 n^3 for a square Q)
        do i = ((k - 1)/QR_NB)*QR_NB + 1, 1, -QR_NB
            ib = min(QR_NB, k - i + 1)
            call LARFT(m - i + 1, ib, v(i, i), m, tau(i), t, QR_NB)
            call LARFB('L', 'N', m - i + 1, n - i + 1, ib, v(i, i), m, t, QR_NB, q(i, i), m)
        end do
        a(1:m, 1:n) = q
    end subroutine ORGQR
