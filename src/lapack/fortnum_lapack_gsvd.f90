module fortnum_lapack_gsvd
    ! Generalized singular value decomposition of a real matrix pair with
    ! reference LAPACK arguments and info codes (prefix fortnum_): the driver
    ! DGGSVD3, U^T A Q = D1 (0 R), V^T B Q = D2 (0 R). The pair is first
    ! reduced to upper triangular form by QR with column pivoting and RQ
    ! factorizations (LAPACK's DGGSVP3), then the 2 x 2 Jacobi-Kogbetliantz
    ! iteration of Paige and of Bai and Demmel (LAPACK's DTGSJA) makes the
    ! rows of the triangular blocks parallel. Also the building blocks DLAGS2
    ! (2 x 2 orthogonal triplet that triangularizes a 2 x 2 pair) and DLASV2
    ! (SVD of a 2 x 2 upper triangular matrix).
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dlarfg, fortnum_dgeqr2, fortnum_dgeqp3, &
        fortnum_dormqr, fortnum_dorgqr
    implicit none
    private

    public :: fortnum_dggsvd3, fortnum_dlags2, fortnum_dlasv2

    real(dp), parameter :: ulp = epsilon(1.0_dp)
    real(dp), parameter :: safmin = tiny(1.0_dp)

contains

    ! [x; y] := [c s; -s c] [x; y]
    pure subroutine rot(x, y, c, s)
        real(dp), intent(inout) :: x(:), y(:)
        real(dp), intent(in) :: c, s
        real(dp) :: t(size(x))
        t = c*x + s*y
        y = c*y - s*x
        x = t
    end subroutine rot

    ! [c s; -s c] [f; g] = [r; 0], c >= 0 (LAPACK 3.10's DLARTG)
    pure subroutine lartg(f, g, c, s, r)
        real(dp), intent(in) :: f, g
        real(dp), intent(out) :: c, s, r
        real(dp) :: d
        if (g == 0) then
            c = 1
            s = 0
            r = f
        else if (f == 0) then
            c = 0
            s = sign(1.0_dp, g)
            r = abs(g)
        else
            d = hypot(f, g)
            c = abs(f)/d
            s = sign(1.0_dp, f)*g/d
            r = sign(d, f)
        end if
    end subroutine lartg

    !> SVD of the 2 x 2 upper triangular [f g; 0 h]:
    !> [csl snl; -snl csl] [f g; 0 h] [csr -snr; snr csr] = diag(ssmax, ssmin)
    !> with abs(ssmax) >= abs(ssmin).
    pure subroutine fortnum_dlasv2(f, g, h, ssmin, ssmax, snr, csr, snl, csl)
        real(dp), intent(in) :: f, g, h
        real(dp), intent(out) :: ssmin, ssmax, snr, csr, snl, csl
        real(dp) :: ft, fa, ht, ha, gt, ga, temp, d, l, m, t, mm, tt, s, r, a
        real(dp) :: clt, crt, slt, srt, tsign
        integer :: pmax
        logical :: swap, gasmal
        ft = f
        fa = abs(ft)
        ht = h
        ha = abs(h)
        ! pmax: the largest element is f (1), g (2) or h (3)
        pmax = 1
        swap = ha > fa
        if (swap) then
            pmax = 3
            temp = ft
            ft = ht
            ht = temp
            temp = fa
            fa = ha
            ha = temp
        end if
        gt = g
        ga = abs(gt)
        if (ga == 0) then
            ssmin = ha
            ssmax = fa
            clt = 1
            crt = 1
            slt = 0
            srt = 0
        else
            gasmal = .true.
            if (ga > fa) then
                pmax = 2
                if (fa/ga < ulp) then
                    ! g dominates
                    gasmal = .false.
                    ssmax = ga
                    if (ha > 1) then
                        ssmin = fa/(ga/ha)
                    else
                        ssmin = (fa/ga)*ha
                    end if
                    clt = 1
                    slt = ht/gt
                    srt = 1
                    crt = ft/gt
                end if
            end if
            if (gasmal) then
                d = fa - ha
                if (d == fa) then
                    l = 1
                else
                    l = d/fa
                end if
                m = gt/ft
                t = 2 - l
                mm = m*m
                tt = t*t
                s = sqrt(tt + mm)
                if (l == 0) then
                    r = abs(m)
                else
                    r = sqrt(l*l + mm)
                end if
                a = (s + r)/2
                ssmin = ha/a
                ssmax = fa*a
                if (mm == 0) then
                    if (l == 0) then
                        t = sign(2.0_dp, ft)*sign(1.0_dp, gt)
                    else
                        t = gt/sign(d, ft) + m/t
                    end if
                else
                    t = (m/(s + t) + m/(r + l))*(1 + a)
                end if
                l = sqrt(t*t + 4)
                crt = 2/l
                srt = t/l
                clt = (crt + srt*m)/a
                slt = (ht/ft)*srt/a
            end if
        end if
        if (swap) then
            csl = srt
            snl = crt
            csr = slt
            snr = clt
        else
            csl = clt
            snl = slt
            csr = crt
            snr = srt
        end if
        if (pmax == 1) tsign = sign(1.0_dp, csr)*sign(1.0_dp, csl)*sign(1.0_dp, f)
        if (pmax == 2) tsign = sign(1.0_dp, snr)*sign(1.0_dp, csl)*sign(1.0_dp, g)
        if (pmax == 3) tsign = sign(1.0_dp, snr)*sign(1.0_dp, snl)*sign(1.0_dp, h)
        ssmax = sign(ssmax, tsign)
        ssmin = sign(ssmin, tsign*sign(1.0_dp, f)*sign(1.0_dp, h))
    end subroutine fortnum_dlasv2

    !> Orthogonal U = [csu snu; -snu csu], V and Q such that for upper
    !> triangular A = [a1 a2; 0 a3], B = [b1 b2; 0 b3] (upper) the (1,2)
    !> entries of U^T A Q and V^T B Q are zero, and for lower triangular
    !> A = [a1 0; a2 a3], B = [b1 0; b2 b3] (.not. upper) the (2,1) entries.
    !> The two remaining rows are parallel to those of the other matrix.
    pure subroutine fortnum_dlags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, &
            csq, snq)
        logical, intent(in) :: upper
        real(dp), intent(in) :: a1, a2, a3, b1, b2, b3
        real(dp), intent(out) :: csu, snu, csv, snv, csq, snq
        real(dp) :: a, b, c, d, s1, s2, snr, csr, snl, csl, r
        real(dp) :: ua11r, ua12, vb11r, vb12, aua12, avb12, ua21, ua22, vb21, vb22
        real(dp) :: aua22, avb22, ua22r, vb22r, aua21, avb21, ua11, vb11, aua11, avb11
        if (upper) then
            ! C = A adj(B) = [a b; 0 d] and its SVD
            a = a1*b3
            d = a3*b1
            b = a2*b1 - a1*b2
            call fortnum_dlasv2(a, b, d, s1, s2, snr, csr, snl, csl)
            if (abs(csl) >= abs(snl) .or. abs(csr) >= abs(snr)) then
                ! the (1,1) and (1,2) entries of U^T A and V^T B
                ua11r = csl*a1
                ua12 = csl*a2 + snl*a3
                vb11r = csr*b1
                vb12 = csr*b2 + snr*b3
                aua12 = abs(csl)*abs(a2) + abs(snl)*abs(a3)
                avb12 = abs(csr)*abs(b2) + abs(snr)*abs(b3)
                if (abs(ua11r) + abs(ua12) /= 0) then
                    if (aua12/(abs(ua11r) + abs(ua12)) <= avb12/(abs(vb11r) + abs(vb12))) then
                        call lartg(-ua11r, ua12, csq, snq, r)
                    else
                        call lartg(-vb11r, vb12, csq, snq, r)
                    end if
                else
                    call lartg(-vb11r, vb12, csq, snq, r)
                end if
                csu = csl
                snu = -snl
                csv = csr
                snv = -snr
            else
                ! the (2,1) and (2,2) entries, zeroed and swapped
                ua21 = -snl*a1
                ua22 = -snl*a2 + csl*a3
                vb21 = -snr*b1
                vb22 = -snr*b2 + csr*b3
                aua22 = abs(snl)*abs(a2) + abs(csl)*abs(a3)
                avb22 = abs(snr)*abs(b2) + abs(csr)*abs(b3)
                if (abs(ua21) + abs(ua22) /= 0) then
                    if (aua22/(abs(ua21) + abs(ua22)) <= avb22/(abs(vb21) + abs(vb22))) then
                        call lartg(-ua21, ua22, csq, snq, r)
                    else
                        call lartg(-vb21, vb22, csq, snq, r)
                    end if
                else
                    call lartg(-vb21, vb22, csq, snq, r)
                end if
                csu = snl
                snu = csl
                csv = snr
                snv = csr
            end if
        else
            ! C = A adj(B) = [a 0; c d] and its SVD
            a = a1*b3
            d = a3*b1
            c = a2*b3 - a3*b2
            call fortnum_dlasv2(a, c, d, s1, s2, snr, csr, snl, csl)
            if (abs(csr) >= abs(snr) .or. abs(csl) >= abs(snl)) then
                ! the (2,1) and (2,2) entries of U^T A and V^T B
                ua21 = -snr*a1 + csr*a2
                ua22r = csr*a3
                vb21 = -snl*b1 + csl*b2
                vb22r = csl*b3
                aua21 = abs(snr)*abs(a1) + abs(csr)*abs(a2)
                avb21 = abs(snl)*abs(b1) + abs(csl)*abs(b2)
                if (abs(ua21) + abs(ua22r) /= 0) then
                    if (aua21/(abs(ua21) + abs(ua22r)) <= avb21/(abs(vb21) + abs(vb22r))) then
                        call lartg(ua22r, ua21, csq, snq, r)
                    else
                        call lartg(vb22r, vb21, csq, snq, r)
                    end if
                else
                    call lartg(vb22r, vb21, csq, snq, r)
                end if
                csu = csr
                snu = -snr
                csv = csl
                snv = -snl
            else
                ! the (1,1) and (1,2) entries, zeroed and swapped
                ua11 = csr*a1 + snr*a2
                ua12 = snr*a3
                vb11 = csl*b1 + snl*b2
                vb12 = snl*b3
                aua11 = abs(csr)*abs(a1) + abs(snr)*abs(a2)
                avb11 = abs(csl)*abs(b1) + abs(snl)*abs(b2)
                if (abs(ua11) + abs(ua12) /= 0) then
                    if (aua11/(abs(ua11) + abs(ua12)) <= avb11/(abs(vb11) + abs(vb12))) then
                        call lartg(ua12, ua11, csq, snq, r)
                    else
                        call lartg(vb12, vb11, csq, snq, r)
                    end if
                else
                    call lartg(vb12, vb11, csq, snq, r)
                end if
                csu = snr
                snu = csr
                csv = snl
                snv = csl
            end if
        end if
    end subroutine fortnum_dlags2

    ! smallest singular value of the n x 2 matrix [x y]
    real(dp) function lapll(x, y) result(ssmin)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: xv(size(x)), yv(size(y)), tau, a11, a12, a22, smax, sr, cr, sl, cl
        integer :: n
        n = size(x)
        ssmin = 0
        if (n <= 1) return
        xv = x
        yv = y
        call fortnum_dlarfg(n, xv(1), xv(2:), 1, tau)
        a11 = xv(1)
        xv(1) = 1
        yv = yv - tau*dot_product(xv, yv)*xv
        call fortnum_dlarfg(n - 1, yv(2), yv(3:), 1, tau)
        a12 = yv(1)
        a22 = yv(2)
        call fortnum_dlasv2(a11, a12, a22, ssmin, smax, sr, cr, sl, cl)
        ssmin = abs(ssmin)
    end function lapll

    ! A := A P, column j of the result is column perm(j) of A
    subroutine permute_columns(m, n, a, lda, perm)
        integer, intent(in) :: m, n, lda, perm(*)
        real(dp), intent(inout) :: a(lda, *)
        real(dp), allocatable :: t(:, :)
        integer :: j
        allocate (t(m, n))
        do j = 1, n
            t(:, j) = a(1:m, perm(j))
        end do
        a(1:m, 1:n) = t
    end subroutine permute_columns

    ! RQ factorization A = R Q of the m x n matrix A, m <= n (LAPACK's
    ! DGERQ2): Q = H(1) ... H(m), H(i) = I - tau(i) v v^T with v(n-m+i) = 1
    ! and v(1:n-m+i-1) stored in row i of a.
    subroutine gerq2(m, n, a, lda, tau)
        integer, intent(in) :: m, n, lda
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: tau(*)
        real(dp), allocatable :: v(:), w(:)
        real(dp) :: alpha
        integer :: i, nv
        allocate (v(n), w(m))
        do i = m, 1, -1
            nv = n - m + i
            alpha = a(i, nv)
            v(1:nv - 1) = a(i, 1:nv - 1)
            call fortnum_dlarfg(nv, alpha, v, 1, tau(i))
            a(i, 1:nv - 1) = v(1:nv - 1)
            a(i, nv) = alpha
            ! rows 1..i-1 times H(i)
            if (i > 1 .and. tau(i) /= 0) then
                v(nv) = 1
                w(1:i - 1) = matmul(a(1:i - 1, 1:nv), v(1:nv))
                a(1:i - 1, 1:nv) = a(1:i - 1, 1:nv) - tau(i)*spread(w(1:i - 1), 2, nv) &
                    *spread(v(1:nv), 1, i - 1)
            end if
        end do
    end subroutine gerq2

    ! C := C Q^T for the m x n matrix C and Q = H(1) ... H(k) of gerq2
    ! applied to vectors of length n (LAPACK's DORMR2, side 'R', trans 'T')
    subroutine ormr2_rt(m, n, k, a, lda, tau, c, ldc)
        integer, intent(in) :: m, n, k, lda, ldc
        real(dp), intent(in) :: a(lda, *), tau(*)
        real(dp), intent(inout) :: c(ldc, *)
        real(dp), allocatable :: v(:), w(:)
        integer :: i, nv
        allocate (v(n), w(m))
        do i = k, 1, -1
            if (tau(i) == 0) cycle
            nv = n - k + i
            v(1:nv - 1) = a(i, 1:nv - 1)
            v(nv) = 1
            w = matmul(c(1:m, 1:nv), v(1:nv))
            c(1:m, 1:nv) = c(1:m, 1:nv) - tau(i)*spread(w, 2, nv)*spread(v(1:nv), 1, m)
        end do
    end subroutine ormr2_rt

    !> Generalized singular value decomposition of the m x n matrix A and the
    !> p x n matrix B: U^T A Q = D1 (0 R), V^T B Q = D2 (0 R) with U, V, Q
    !> orthogonal (jobu 'U', jobv 'V', jobq 'Q' compute them), k + l the
    !> effective rank of [A; B], R (k+l) x (k+l) upper triangular, returned
    !> in a (and in b when m < k + l) as LAPACK's DGGSVD3 does, and
    !> alpha(i)^2 + beta(i)^2 = 1. iwork records the sort of
    !> alpha(k+1:min(m,k+l)) into decreasing order. info = 1: the Jacobi
    !> iteration did not converge.
    subroutine fortnum_dggsvd3(jobu, jobv, jobq, m, n, p, k, l, a, lda, b, ldb, alpha, &
            beta, u, ldu, v, ldv, q, ldq, work, lwork, iwork, info)
        character, intent(in) :: jobu, jobv, jobq
        integer, intent(in) :: m, n, p, lda, ldb, ldu, ldv, ldq, lwork
        integer, intent(out) :: k, l, iwork(*), info
        real(dp), intent(inout) :: a(lda, *), b(ldb, *)
        real(dp), intent(out) :: alpha(*), beta(*), u(ldu, *), v(ldv, *), q(ldq, *), work(*)
        logical :: wantu, wantv, wantq
        real(dp) :: anorm, bnorm, tola, tolb, smax, temp
        real(dp), allocatable :: sorted(:)
        integer :: i, j, ibnd, isub, lwmin
        info = 0
        wantu = any(jobu == ['U', 'u'])
        wantv = any(jobv == ['V', 'v'])
        wantq = any(jobq == ['Q', 'q'])
        lwmin = max(1, 4*n + m + p, 3*n + 1 + n)
        if (.not. wantu .and. .not. any(jobu == ['N', 'n'])) then
            info = -1
        else if (.not. wantv .and. .not. any(jobv == ['N', 'n'])) then
            info = -2
        else if (.not. wantq .and. .not. any(jobq == ['N', 'n'])) then
            info = -3
        else if (m < 0) then
            info = -4
        else if (n < 0) then
            info = -5
        else if (p < 0) then
            info = -6
        else if (lda < max(1, m)) then
            info = -10
        else if (ldb < max(1, p)) then
            info = -12
        else if (ldu < 1 .or. (wantu .and. ldu < m)) then
            info = -16
        else if (ldv < 1 .or. (wantv .and. ldv < p)) then
            info = -18
        else if (ldq < 1 .or. (wantq .and. ldq < n)) then
            info = -20
        else if (lwork < 1 .and. lwork /= -1) then
            info = -24
        end if
        if (info /= 0) return
        work(1) = lwmin
        if (lwork == -1) return
        k = 0
        l = 0
        ! 1-norms decide the numerical ranks
        anorm = 0
        bnorm = 0
        do j = 1, n
            if (m > 0) anorm = max(anorm, sum(abs(a(1:m, j))))
            if (p > 0) bnorm = max(bnorm, sum(abs(b(1:p, j))))
        end do
        tola = max(m, n)*max(anorm, safmin)*ulp
        tolb = max(p, n)*max(bnorm, safmin)*ulp
        call ggsvp3(wantu, wantv, wantq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, u, ldu, &
            v, ldv, q, ldq)
        call tgsja(wantu, wantv, wantq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, alpha, &
            beta, u, ldu, v, ldv, q, ldq, info)
        ! sort alpha(k+1:k+ibnd) decreasingly in a copy, recording the swaps
        allocate (sorted(max(1, n)))
        sorted(1:n) = alpha(1:n)
        ibnd = min(l, m - k)
        do i = 1, ibnd
            isub = i
            smax = sorted(k + i)
            do j = i + 1, ibnd
                temp = sorted(k + j)
                if (temp > smax) then
                    isub = j
                    smax = temp
                end if
            end do
            if (isub /= i) then
                sorted(k + isub) = sorted(k + i)
                sorted(k + i) = smax
                iwork(k + i) = k + isub
            else
                iwork(k + i) = k + i
            end if
        end do
    end subroutine fortnum_dggsvd3

    ! LAPACK's DGGSVP3: orthogonal U, V, Q with
    !   U^T A Q = [0 A12 A13; 0 0 A23; 0 0 0] (rows k, l, m-k-l; columns
    !   n-k-l, k, l), V^T B Q = [0 0 B13; 0 0 0] (rows l, p-l),
    ! A12 and B13 upper triangular and nonsingular, A23 upper trapezoidal.
    subroutine ggsvp3(wantu, wantv, wantq, m, p, n, a, lda, b, ldb, tola, tolb, k, l, &
            u, ldu, v, ldv, q, ldq)
        logical, intent(in) :: wantu, wantv, wantq
        integer, intent(in) :: m, p, n, lda, ldb, ldu, ldv, ldq
        real(dp), intent(inout) :: a(lda, *), b(ldb, *), u(ldu, *), v(ldv, *), q(ldq, *)
        real(dp), intent(in) :: tola, tolb
        integer, intent(out) :: k, l
        integer, allocatable :: piv(:)
        real(dp), allocatable :: tau(:), work(:)
        integer :: i, j, info
        allocate (piv(max(1, n)), tau(max(1, n, m, p)), work(max(1, 64*(n + m + p))))
        ! B P = V [S11 S12; 0 0] by QR with column pivoting
        piv = 0
        call fortnum_dgeqp3(p, n, b, ldb, piv, tau, work, size(work), info)
        call permute_columns(m, n, a, lda, piv)
        l = 0
        do i = 1, min(p, n)
            if (abs(b(i, i)) > tolb) l = l + 1
        end do
        if (wantv) then
            v(1:p, 1:p) = 0
            do j = 1, min(p, n)
                v(j + 1:p, j) = b(j + 1:p, j)
            end do
            call fortnum_dorgqr(p, p, min(p, n), v, ldv, tau, work, size(work), info)
        end if
        do j = 1, l - 1
            b(j + 1:l, j) = 0
        end do
        if (p > l) b(l + 1:p, 1:n) = 0
        if (wantq) then
            q(1:n, 1:n) = 0
            do j = 1, n
                q(j, j) = 1
            end do
            call permute_columns(n, n, q, ldq, piv)
        end if
        if (p >= l .and. n /= l) then
            ! [S11 S12] = [0 S12'] Z by RQ, A := A Z^T, Q := Q Z^T
            call gerq2(l, n, b, ldb, tau)
            call ormr2_rt(m, n, l, b, ldb, tau, a, lda)
            if (wantq) call ormr2_rt(n, n, l, b, ldb, tau, q, ldq)
            b(1:l, 1:n - l) = 0
            do j = n - l + 1, n
                b(j - n + l + 1:l, j) = 0
            end do
        end if
        ! A11 = U [0 T12; 0 0] P1^T: QR with column pivoting of A(:, 1:n-l)
        piv = 0
        call fortnum_dgeqp3(m, n - l, a, lda, piv, tau, work, size(work), info)
        k = 0
        do i = 1, min(m, n - l)
            if (abs(a(i, i)) > tola) k = k + 1
        end do
        if (l > 0) call fortnum_dormqr('L', 'T', m, l, min(m, n - l), a, lda, tau, &
            a(1, n - l + 1), lda, work, size(work), info)
        if (wantu) then
            u(1:m, 1:m) = 0
            do j = 1, min(m, n - l)
                u(j + 1:m, j) = a(j + 1:m, j)
            end do
            call fortnum_dorgqr(m, m, min(m, n - l), u, ldu, tau, work, size(work), info)
        end if
        if (wantq) call permute_columns(n, n - l, q, ldq, piv)
        do j = 1, k - 1
            a(j + 1:k, j) = 0
        end do
        if (m > k) a(k + 1:m, 1:n - l) = 0
        if (n - l > k) then
            ! [T11 T12] = [0 T12'] Z1 by RQ, Q(:, 1:n-l) := Q(:, 1:n-l) Z1^T
            call gerq2(k, n - l, a, lda, tau)
            if (wantq) call ormr2_rt(n, n - l, k, a, lda, tau, q, ldq)
            a(1:k, 1:n - l - k) = 0
            do j = n - l - k + 1, n - l
                a(j - n + l + k + 1:k, j) = 0
            end do
        end if
        if (m > k .and. l > 0) then
            ! QR of A(k+1:m, n-l+1:n), U(:, k+1:m) := U(:, k+1:m) U1
            call fortnum_dgeqr2(m - k, l, a(k + 1, n - l + 1), lda, tau, info)
            if (wantu) call fortnum_dormqr('R', 'N', m, m - k, min(m - k, l), &
                a(k + 1, n - l + 1), lda, tau, u(1, k + 1), ldu, work, size(work), info)
            do j = n - l + 1, n
                a(j - n + k + l + 1:m, j) = 0
            end do
        end if
    end subroutine ggsvp3

    ! LAPACK's DTGSJA: the GSVD of the upper triangular A13 = A(k+1:k+l or m,
    ! n-l+1:n) and B13 = B(1:l, n-l+1:n) by 2 x 2 rotations until their rows
    ! are parallel; U, V and Q accumulate the rotations.
    subroutine tgsja(wantu, wantv, wantq, m, p, n, k, l, a, lda, b, ldb, tola, tolb, &
            alpha, beta, u, ldu, v, ldv, q, ldq, info)
        logical, intent(in) :: wantu, wantv, wantq
        integer, intent(in) :: m, p, n, k, l, lda, ldb, ldu, ldv, ldq
        real(dp), intent(inout) :: a(lda, *), b(ldb, *), u(ldu, *), v(ldv, *), q(ldq, *)
        real(dp), intent(in) :: tola, tolb
        real(dp), intent(out) :: alpha(*), beta(*)
        integer, intent(out) :: info
        integer, parameter :: maxit = 40
        logical :: upper, converged
        real(dp) :: a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, csq, snq, error, ssmin, gamma
        real(dp) :: rwk
        integer :: kcycle, i, j, nl
        info = 0
        upper = .false.
        converged = .false.
        nl = n - l
        do kcycle = 1, maxit
            upper = .not. upper
            do i = 1, l - 1
                do j = i + 1, l
                    a1 = 0
                    a2 = 0
                    a3 = 0
                    if (k + i <= m) a1 = a(k + i, nl + i)
                    if (k + j <= m) a3 = a(k + j, nl + j)
                    b1 = b(i, nl + i)
                    b3 = b(j, nl + j)
                    if (upper) then
                        if (k + i <= m) a2 = a(k + i, nl + j)
                        b2 = b(i, nl + j)
                    else
                        if (k + j <= m) a2 = a(k + j, nl + i)
                        b2 = b(j, nl + i)
                    end if
                    call fortnum_dlags2(upper, a1, a2, a3, b1, b2, b3, csu, snu, csv, snv, &
                        csq, snq)
                    ! rows k+j, k+i of A, rows j, i of B, columns nl+j, nl+i
                    if (k + j <= m) call rot(a(k + j, nl + 1:n), a(k + i, nl + 1:n), csu, snu)
                    call rot(b(j, nl + 1:n), b(i, nl + 1:n), csv, snv)
                    call rot(a(1:min(k + l, m), nl + j), a(1:min(k + l, m), nl + i), csq, snq)
                    call rot(b(1:l, nl + j), b(1:l, nl + i), csq, snq)
                    if (upper) then
                        if (k + i <= m) a(k + i, nl + j) = 0
                        b(i, nl + j) = 0
                    else
                        if (k + j <= m) a(k + j, nl + i) = 0
                        b(j, nl + i) = 0
                    end if
                    if (wantu .and. k + j <= m) call rot(u(1:m, k + j), u(1:m, k + i), csu, snu)
                    if (wantv) call rot(v(1:p, j), v(1:p, i), csv, snv)
                    if (wantq) call rot(q(1:n, nl + j), q(1:n, nl + i), csq, snq)
                end do
            end do
            if (.not. upper) then
                ! converged when the rows i of A13 and B13 are parallel
                error = 0
                do i = 1, min(l, m - k)
                    ssmin = lapll(a(k + i, nl + i:n), b(i, nl + i:n))
                    error = max(error, ssmin)
                end do
                if (abs(error) <= min(tola, tolb)) then
                    converged = .true.
                    exit
                end if
            end if
        end do
        if (.not. converged) then
            info = 1
            return
        end if
        alpha(1:k) = 1
        beta(1:k) = 0
        do i = 1, min(l, m - k)
            a1 = a(k + i, nl + i)
            b1 = b(i, nl + i)
            gamma = huge(1.0_dp)
            if (a1 /= 0) gamma = b1/a1
            if (abs(gamma) < huge(1.0_dp)) then
                if (gamma < 0) then
                    b(i, nl + i:n) = -b(i, nl + i:n)
                    if (wantv) v(1:p, i) = -v(1:p, i)
                end if
                call lartg(abs(gamma), 1.0_dp, beta(k + i), alpha(k + i), rwk)
                if (alpha(k + i) >= beta(k + i)) then
                    a(k + i, nl + i:n) = a(k + i, nl + i:n)/alpha(k + i)
                else
                    b(i, nl + i:n) = b(i, nl + i:n)/beta(k + i)
                    a(k + i, nl + i:n) = b(i, nl + i:n)
                end if
            else
                alpha(k + i) = 0
                beta(k + i) = 1
                a(k + i, nl + i:n) = b(i, nl + i:n)
            end if
        end do
        do i = m + 1, k + l
            alpha(i) = 0
            beta(i) = 1
        end do
        do i = k + l + 1, n
            alpha(i) = 0
            beta(i) = 0
        end do
    end subroutine tgsja

end module fortnum_lapack_gsvd
