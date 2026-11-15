module fortnum_lapack_geev
    ! Nonsymmetric eigenproblems with reference LAPACK arguments and info
    ! codes (prefix fortnum_): balancing DGEBAL/ZGEBAL and its inverse
    ! DGEBAK/ZGEBAK, Hessenberg reduction DGEHRD/ZGEHRD with DORGHR/ZUNGHR,
    ! the Hessenberg QR algorithm DHSEQR/ZHSEQR (LAPACK's xLAHQR: Francis
    ! double shifts for real matrices with standardized 2 x 2 blocks, single
    ! Wilkinson shifts for complex ones, Ahues-Kressner deflation), and the
    ! drivers DGEEV/ZGEEV. Eigenvectors come from back substitution on the
    ! Schur form in complex arithmetic and are normalized as LAPACK does.
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_dlarfg, fortnum_zlarfg
    implicit none
    private

    public :: fortnum_dgebal, fortnum_dgebak, fortnum_dgehrd, fortnum_dorghr, fortnum_dhseqr, fortnum_dgeev
    public :: fortnum_zgebal, fortnum_zgebak, fortnum_zgehrd, fortnum_zunghr, fortnum_zhseqr, fortnum_zgeev

    real(dp), parameter :: ulp = epsilon(1.0_dp)
    real(dp), parameter :: safmin = tiny(1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp), cone = (1.0_dp, 0.0_dp)

contains

    ! ------------------------------------------------------------ balancing

    !> Balance A (job 'N', 'P', 'S' or 'B'): permute to isolate eigenvalues
    !> in A(1:ilo-1, :) and A(ihi+1:n, :), then scale rows and columns of
    !> A(ilo:ihi, ilo:ihi) by powers of 2 so that their norms are close.
    subroutine fortnum_dgebal(job, n, a, lda, ilo, ihi, scale, info)
        character, intent(in) :: job
        integer, intent(in) :: n, lda
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: ilo, ihi, info
        real(dp), intent(out) :: scale(*)
        complex(dp), allocatable :: z(:, :)
        integer :: i
        info = 0
        if (.not. any(job == ['N', 'n', 'P', 'p', 'S', 's', 'B', 'b'])) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        end if
        if (info /= 0) return
        ! the complex balancing on the real matrix: same permutations, same
        ! powers of 2, since only magnitudes enter
        allocate (z(max(1, n), max(1, n)))
        do i = 1, n
            z(1:n, i) = a(1:n, i)
        end do
        call zgebal_core(job, n, z, max(1, n), ilo, ihi, scale, info)
        do i = 1, n
            a(1:n, i) = real(z(1:n, i), dp)
        end do
    end subroutine fortnum_dgebal

    subroutine fortnum_zgebal(job, n, a, lda, ilo, ihi, scale, info)
        character, intent(in) :: job
        integer, intent(in) :: n, lda
        complex(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: ilo, ihi, info
        real(dp), intent(out) :: scale(*)
        info = 0
        if (.not. any(job == ['N', 'n', 'P', 'p', 'S', 's', 'B', 'b'])) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        end if
        if (info /= 0) return
        call zgebal_core(job, n, a, lda, ilo, ihi, scale, info)
    end subroutine fortnum_zgebal

    subroutine zgebal_core(job, n, a, lda, ilo, ihi, scale, info)
        character, intent(in) :: job
        integer, intent(in) :: n, lda
        complex(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: ilo, ihi, info
        real(dp), intent(out) :: scale(*)
        real(dp), parameter :: sclfac = 2, factor = 0.95_dp
        real(dp) :: c, r, ca, ra, g, f, s, sfmin1, sfmax1, sfmin2, sfmax2
        logical :: noconv, canswap, permute, scaling
        integer :: i, j, k, l, ica, ira
        info = 0
        ilo = 1
        ihi = n
        if (n == 0) return
        permute = any(job == ['P', 'p', 'B', 'b'])
        scaling = any(job == ['S', 's', 'B', 'b'])
        if (.not. permute .and. .not. scaling) then
            scale(1:n) = 1
            return
        end if
        k = 1
        l = n
        if (permute) then
            ! rows isolating an eigenvalue go to the bottom
            noconv = .true.
            do while (noconv)
                noconv = .false.
                do i = l, 1, -1
                    canswap = .true.
                    do j = 1, l
                        if (i /= j .and. a(i, j) /= czero) then
                            canswap = .false.
                            exit
                        end if
                    end do
                    if (canswap) then
                        scale(l) = i
                        if (i /= l) then
                            call zswap_cols(l, a, lda, i, l)
                            call zswap_rows(n - k + 1, a(1, k), lda, i, l)
                        end if
                        noconv = .true.
                        if (l == 1) then
                            ilo = 1
                            ihi = 1
                            return
                        end if
                        l = l - 1
                    end if
                end do
            end do
            ! columns isolating an eigenvalue go to the left
            noconv = .true.
            do while (noconv)
                noconv = .false.
                do j = k, l
                    canswap = .true.
                    do i = k, l
                        if (i /= j .and. a(i, j) /= czero) then
                            canswap = .false.
                            exit
                        end if
                    end do
                    if (canswap) then
                        scale(k) = j
                        if (j /= k) then
                            call zswap_cols(l, a, lda, j, k)
                            call zswap_rows(n - k + 1, a(1, k), lda, j, k)
                        end if
                        noconv = .true.
                        k = k + 1
                    end if
                end do
            end do
        end if
        scale(k:l) = 1
        ilo = k
        ihi = l
        if (.not. scaling) return
        sfmin1 = safmin/ulp
        sfmax1 = 1/sfmin1
        sfmin2 = sfmin1*sclfac
        sfmax2 = 1/sfmin2
        noconv = .true.
        do while (noconv)
            noconv = .false.
            do i = k, l
                c = norm2_c(a(k:l, i))
                r = norm2_c(a(i, k:l))
                ica = maxloc(abs(a(1:l, i)), 1)
                ca = abs(a(ica, i))
                ira = maxloc(abs(a(i, k:n)), 1)
                ra = abs(a(i, ira + k - 1))
                if (c == 0 .or. r == 0) cycle
                if (c + ca + r + ra /= c + ca + r + ra) then
                    info = -3
                    return
                end if
                g = r/sclfac
                f = 1
                s = c + r
                do while (c < g .and. max(f, c, ca) < sfmax2 .and. min(r, g, ra) > sfmin2)
                    f = f*sclfac
                    c = c*sclfac
                    ca = ca*sclfac
                    r = r/sclfac
                    g = g/sclfac
                    ra = ra/sclfac
                end do
                g = c/sclfac
                do while (g >= r .and. max(r, ra) < sfmax2 .and. min(f, c, g, ca) > sfmin2)
                    f = f/sclfac
                    c = c/sclfac
                    g = g/sclfac
                    ca = ca/sclfac
                    r = r*sclfac
                    ra = ra*sclfac
                end do
                if (c + r >= factor*s) cycle
                if (f < 1 .and. scale(i) < 1) then
                    if (f*scale(i) <= sfmin1) cycle
                end if
                if (f > 1 .and. scale(i) > 1) then
                    if (scale(i) >= sfmax1/f) cycle
                end if
                scale(i) = scale(i)*f
                noconv = .true.
                a(i, k:n) = a(i, k:n)/f
                a(1:l, i) = a(1:l, i)*f
            end do
        end do
    end subroutine zgebal_core

    real(dp) function norm2_c(x)
        complex(dp), intent(in) :: x(:)
        real(dp) :: scl, ssq, v
        integer :: i, part
        scl = 0
        ssq = 1
        do i = 1, size(x)
            do part = 1, 2
                if (part == 1) then
                    v = abs(real(x(i), dp))
                else
                    v = abs(aimag(x(i)))
                end if
                if (v /= 0) then
                    if (scl < v) then
                        ssq = 1 + ssq*(scl/v)**2
                        scl = v
                    else
                        ssq = ssq + (v/scl)**2
                    end if
                end if
            end do
        end do
        norm2_c = scl*sqrt(ssq)
    end function norm2_c

    subroutine zswap_cols(m, a, lda, i, j)
        integer, intent(in) :: m, lda, i, j
        complex(dp), intent(inout) :: a(lda, *)
        complex(dp) :: t(m)
        t = a(1:m, i)
        a(1:m, i) = a(1:m, j)
        a(1:m, j) = t
    end subroutine zswap_cols

    ! rows i and j of the m columns starting at a
    subroutine zswap_rows(m, a, lda, i, j)
        integer, intent(in) :: m, lda, i, j
        complex(dp), intent(inout) :: a(lda, *)
        complex(dp) :: t(m)
        t = a(i, 1:m)
        a(i, 1:m) = a(j, 1:m)
        a(j, 1:m) = t
    end subroutine zswap_rows

    !> Back-transform the m eigenvectors in V (side 'R' right, 'L' left) of
    !> a matrix balanced by xGEBAL with job.
    subroutine fortnum_dgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
        character, intent(in) :: job, side
        integer, intent(in) :: n, ilo, ihi, m, ldv
        real(dp), intent(in) :: scale(*)
        real(dp), intent(inout) :: v(ldv, *)
        integer, intent(out) :: info
        integer :: i, ii, k
        call gebak_check(job, side, n, ilo, ihi, m, ldv, info)
        if (info /= 0 .or. n == 0 .or. m == 0 .or. any(job == ['N', 'n'])) return
        if (ilo /= ihi .and. any(job == ['S', 's', 'B', 'b'])) then
            do i = ilo, ihi
                if (any(side == ['R', 'r'])) then
                    v(i, 1:m) = v(i, 1:m)*scale(i)
                else
                    v(i, 1:m) = v(i, 1:m)/scale(i)
                end if
            end do
        end if
        if (any(job == ['P', 'p', 'B', 'b'])) then
            do ii = 1, n
                i = ii
                if (i >= ilo .and. i <= ihi) cycle
                if (i < ilo) i = ilo - ii
                k = int(scale(i))
                if (k == i) cycle
                call dswap_rows(m, v, ldv, i, k)
            end do
        end if
    end subroutine fortnum_dgebak

    subroutine fortnum_zgebak(job, side, n, ilo, ihi, scale, m, v, ldv, info)
        character, intent(in) :: job, side
        integer, intent(in) :: n, ilo, ihi, m, ldv
        real(dp), intent(in) :: scale(*)
        complex(dp), intent(inout) :: v(ldv, *)
        integer, intent(out) :: info
        integer :: i, ii, k
        call gebak_check(job, side, n, ilo, ihi, m, ldv, info)
        if (info /= 0 .or. n == 0 .or. m == 0 .or. any(job == ['N', 'n'])) return
        if (ilo /= ihi .and. any(job == ['S', 's', 'B', 'b'])) then
            do i = ilo, ihi
                if (any(side == ['R', 'r'])) then
                    v(i, 1:m) = v(i, 1:m)*scale(i)
                else
                    v(i, 1:m) = v(i, 1:m)/scale(i)
                end if
            end do
        end if
        if (any(job == ['P', 'p', 'B', 'b'])) then
            do ii = 1, n
                i = ii
                if (i >= ilo .and. i <= ihi) cycle
                if (i < ilo) i = ilo - ii
                k = int(scale(i))
                if (k == i) cycle
                call zswap_rows(m, v, ldv, i, k)
            end do
        end if
    end subroutine fortnum_zgebak

    subroutine gebak_check(job, side, n, ilo, ihi, m, ldv, info)
        character, intent(in) :: job, side
        integer, intent(in) :: n, ilo, ihi, m, ldv
        integer, intent(out) :: info
        info = 0
        if (.not. any(job == ['N', 'n', 'P', 'p', 'S', 's', 'B', 'b'])) then
            info = -1
        else if (.not. any(side == ['R', 'r', 'L', 'l'])) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (ilo < 1 .or. ilo > max(1, n)) then
            info = -4
        else if (ihi < min(ilo, n) .or. ihi > n) then
            info = -5
        else if (m < 0) then
            info = -7
        else if (ldv < max(1, n)) then
            info = -9
        end if
    end subroutine gebak_check

    subroutine dswap_rows(m, a, lda, i, j)
        integer, intent(in) :: m, lda, i, j
        real(dp), intent(inout) :: a(lda, *)
        real(dp) :: t(m)
        t = a(i, 1:m)
        a(i, 1:m) = a(j, 1:m)
        a(j, 1:m) = t
    end subroutine dswap_rows

    ! ------------------------------------------------- Hessenberg reduction

    !> Reduce A(ilo:ihi, ilo:ihi) to upper Hessenberg form H = Q^T A Q,
    !> Q = H(ilo) ... H(ihi-1) stored below the subdiagonal with tau.
    subroutine fortnum_dgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: tau(*), work(*)
        integer, intent(out) :: info
        real(dp) :: alpha
        integer :: i
        call gehrd_check(n, ilo, ihi, lda, lwork, info)
        if (info /= 0) return
        work(1) = max(1, n)
        if (lwork == -1) return
        tau(1:ilo - 1) = 0
        tau(max(1, ihi):n - 1) = 0
        do i = ilo, ihi - 1
            alpha = a(i + 1, i)
            call fortnum_dlarfg(ihi - i, alpha, a(min(i + 2, n):, i), 1, tau(i))
            a(i + 1, i) = 1
            call dapply_right(ihi, ihi - i, a(i + 1:ihi, i), tau(i), a(1:ihi, i + 1:ihi))
            call dapply_left(ihi - i, n - i, a(i + 1:ihi, i), tau(i), a(i + 1:ihi, i + 1:n))
            a(i + 1, i) = alpha
        end do
    end subroutine fortnum_dgehrd

    subroutine fortnum_zgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        complex(dp), intent(inout) :: a(lda, *)
        complex(dp), intent(out) :: tau(*), work(*)
        integer, intent(out) :: info
        complex(dp) :: alpha
        integer :: i
        call gehrd_check(n, ilo, ihi, lda, lwork, info)
        if (info /= 0) return
        work(1) = max(1, n)
        if (lwork == -1) return
        tau(1:ilo - 1) = czero
        tau(max(1, ihi):n - 1) = czero
        do i = ilo, ihi - 1
            alpha = a(i + 1, i)
            call fortnum_zlarfg(ihi - i, alpha, a(min(i + 2, n):, i), 1, tau(i))
            a(i + 1, i) = cone
            call zapply_right(ihi, ihi - i, a(i + 1:ihi, i), tau(i), a(1:ihi, i + 1:ihi))
            call zapply_left(ihi - i, n - i, a(i + 1:ihi, i), conjg(tau(i)), a(i + 1:ihi, i + 1:n))
            a(i + 1, i) = alpha
        end do
    end subroutine fortnum_zgehrd

    subroutine gehrd_check(n, ilo, ihi, lda, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        integer, intent(out) :: info
        info = 0
        if (n < 0) then
            info = -1
        else if (ilo < 1 .or. ilo > max(1, n)) then
            info = -2
        else if (ihi < min(ilo, n) .or. ihi > n) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (lwork < max(1, n) .and. lwork /= -1) then
            info = -8
        end if
    end subroutine gehrd_check

    ! c := (I - tau v v^T) c
    subroutine dapply_left(m, n, v, tau, c)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: v(:), tau
        real(dp), intent(inout) :: c(:, :)
        real(dp) :: w
        integer :: j
        if (tau == 0 .or. m == 0) return
        do j = 1, n
            w = tau*dot_product(v(1:m), c(1:m, j))
            c(1:m, j) = c(1:m, j) - w*v(1:m)
        end do
    end subroutine dapply_left

    ! c := c (I - tau v v^T)
    subroutine dapply_right(m, n, v, tau, c)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: v(:), tau
        real(dp), intent(inout) :: c(:, :)
        real(dp) :: w(m)
        integer :: j
        if (tau == 0 .or. n == 0) return
        w = 0
        do j = 1, n
            w = w + v(j)*c(1:m, j)
        end do
        w = tau*w
        do j = 1, n
            c(1:m, j) = c(1:m, j) - w*v(j)
        end do
    end subroutine dapply_right

    ! c := (I - tau v v^H) c
    subroutine zapply_left(m, n, v, tau, c)
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: v(:), tau
        complex(dp), intent(inout) :: c(:, :)
        complex(dp) :: w
        integer :: j
        if (tau == czero .or. m == 0) return
        do j = 1, n
            w = tau*dot_product(v(1:m), c(1:m, j))
            c(1:m, j) = c(1:m, j) - w*v(1:m)
        end do
    end subroutine zapply_left

    ! c := c (I - tau v v^H)
    subroutine zapply_right(m, n, v, tau, c)
        integer, intent(in) :: m, n
        complex(dp), intent(in) :: v(:), tau
        complex(dp), intent(inout) :: c(:, :)
        complex(dp) :: w(m)
        integer :: j
        if (tau == czero .or. n == 0) return
        w = czero
        do j = 1, n
            w = w + v(j)*c(1:m, j)
        end do
        w = tau*w
        do j = 1, n
            c(1:m, j) = c(1:m, j) - w*conjg(v(j))
        end do
    end subroutine zapply_right

    !> Q of xGEHRD explicitly, overwriting A (n x n).
    subroutine fortnum_dorghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(in) :: tau(*)
        real(dp), intent(out) :: work(*)
        integer, intent(out) :: info
        real(dp), allocatable :: v(:, :)
        integer :: i, j
        call orghr_check(n, ilo, ihi, lda, lwork, info)
        if (info /= 0) return
        work(1) = max(1, ihi - ilo)
        if (lwork == -1 .or. n == 0) return
        allocate (v(n, n))
        v = 0
        do i = ilo, ihi - 1
            v(i + 1, i) = 1
            v(i + 2:ihi, i) = a(i + 2:ihi, i)
        end do
        do j = 1, n
            a(1:n, j) = 0
            a(j, j) = 1
        end do
        ! backward accumulation: H(i) touches rows and columns i+1:ihi
        do i = ihi - 1, ilo, -1
            call dapply_left(ihi - i, ihi - i, v(i + 1:ihi, i), tau(i), a(i + 1:ihi, i + 1:ihi))
        end do
    end subroutine fortnum_dorghr

    subroutine fortnum_zunghr(n, ilo, ihi, a, lda, tau, work, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        complex(dp), intent(inout) :: a(lda, *)
        complex(dp), intent(in) :: tau(*)
        complex(dp), intent(out) :: work(*)
        integer, intent(out) :: info
        complex(dp), allocatable :: v(:, :)
        integer :: i, j
        call orghr_check(n, ilo, ihi, lda, lwork, info)
        if (info /= 0) return
        work(1) = max(1, ihi - ilo)
        if (lwork == -1 .or. n == 0) return
        allocate (v(n, n))
        v = czero
        do i = ilo, ihi - 1
            v(i + 1, i) = cone
            v(i + 2:ihi, i) = a(i + 2:ihi, i)
        end do
        do j = 1, n
            a(1:n, j) = czero
            a(j, j) = cone
        end do
        do i = ihi - 1, ilo, -1
            call zapply_left(ihi - i, ihi - i, v(i + 1:ihi, i), tau(i), a(i + 1:ihi, i + 1:ihi))
        end do
    end subroutine fortnum_zunghr

    subroutine orghr_check(n, ilo, ihi, lda, lwork, info)
        integer, intent(in) :: n, ilo, ihi, lda, lwork
        integer, intent(out) :: info
        info = 0
        if (n < 0) then
            info = -1
        else if (ilo < 1 .or. ilo > max(1, n)) then
            info = -2
        else if (ihi < min(ilo, n) .or. ihi > n) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (lwork < max(1, ihi - ilo) .and. lwork /= -1) then
            info = -8
        end if
    end subroutine orghr_check

    ! ----------------------------------------------- Hessenberg QR, real

    !> Eigenvalues (job 'E') or Schur form (job 'S') of the Hessenberg
    !> matrix H, with Schur vectors not computed (compz 'N'), of H ('I') or
    !> multiplied into Z ('V').
    subroutine fortnum_dhseqr(job, compz, n, ilo, ihi, h, ldh, wr, wi, z, ldz, work, lwork, info)
        character, intent(in) :: job, compz
        integer, intent(in) :: n, ilo, ihi, ldh, ldz, lwork
        real(dp), intent(inout) :: h(ldh, *), z(ldz, *)
        real(dp), intent(out) :: wr(*), wi(*), work(*)
        integer, intent(out) :: info
        logical :: wantt, wantz
        integer :: i, j
        call hseqr_check(job, compz, n, ilo, ihi, ldh, ldz, lwork, info)
        if (info /= 0) return
        work(1) = max(1, n)
        if (lwork == -1 .or. n == 0) return
        wantt = any(job == ['S', 's'])
        wantz = .not. any(compz == ['N', 'n'])
        if (any(compz == ['I', 'i'])) then
            do j = 1, n
                z(1:n, j) = 0
                z(j, j) = 1
            end do
        end if
        ! eigenvalues isolated by balancing
        do i = 1, ilo - 1
            wr(i) = h(i, i)
            wi(i) = 0
        end do
        do i = ihi + 1, n
            wr(i) = h(i, i)
            wi(i) = 0
        end do
        call dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, 1, n, z, ldz, info)
        if (wantt .and. info == 0) then
            do j = 1, n - 2
                h(j + 2:n, j) = 0
            end do
        end if
    end subroutine fortnum_dhseqr

    subroutine hseqr_check(job, compz, n, ilo, ihi, ldh, ldz, lwork, info)
        character, intent(in) :: job, compz
        integer, intent(in) :: n, ilo, ihi, ldh, ldz, lwork
        integer, intent(out) :: info
        info = 0
        if (.not. any(job == ['E', 'e', 'S', 's'])) then
            info = -1
        else if (.not. any(compz == ['N', 'n', 'I', 'i', 'V', 'v'])) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (ilo < 1 .or. ilo > max(1, n)) then
            info = -4
        else if (ihi < min(ilo, n) .or. ihi > n) then
            info = -5
        else if (ldh < max(1, n)) then
            info = -7
        else if (ldz < 1 .or. (.not. any(compz == ['N', 'n']) .and. ldz < max(1, n))) then
            info = -11
        else if (lwork < max(1, n) .and. lwork /= -1) then
            info = -13
        end if
    end subroutine hseqr_check

    ! LAPACK's DLAHQR: double-shift QR on H(ilo:ihi, ilo:ihi)
    subroutine dlahqr(wantt, wantz, n, ilo, ihi, h, ldh, wr, wi, iloz, ihiz, z, ldz, info)
        logical, intent(in) :: wantt, wantz
        integer, intent(in) :: n, ilo, ihi, ldh, iloz, ihiz, ldz
        real(dp), intent(inout) :: h(ldh, *), z(ldz, *)
        real(dp), intent(inout) :: wr(*), wi(*)
        integer, intent(out) :: info
        real(dp), parameter :: dat1 = 0.75_dp, dat2 = -0.4375_dp
        integer, parameter :: kexsh = 10
        real(dp) :: smlnum, tst, ab, ba, aa, bb, s, h11, h12, h21, h22, tr, det, rtdisc
        real(dp) :: rt1r, rt1i, rt2r, rt2i, h21s, h00, h01, t1, t2, t3, v2, v3, sm, cs, sn
        real(dp) :: v(3)
        integer :: i, j, k, l, m, its, itmax, kdefl, nr, i1, i2, nh, nz
        logical :: found
        info = 0
        if (n == 0) return
        if (ilo == ihi) then
            wr(ilo) = h(ilo, ilo)
            wi(ilo) = 0
            return
        end if
        do j = ilo, ihi - 3
            h(j + 2, j) = 0
            h(j + 3, j) = 0
        end do
        if (ilo <= ihi - 2) h(ihi, ihi - 2) = 0
        nh = ihi - ilo + 1
        nz = ihiz - iloz + 1
        smlnum = safmin*(real(nh, dp)/ulp)
        i1 = 1
        i2 = n
        itmax = 30*max(10, nh)
        kdefl = 0
        i = ihi
        do while (i >= ilo)
            l = ilo
            found = .false.
            do its = 0, itmax
                ! a single small subdiagonal element
                do k = i, l + 1, -1
                    if (abs(h(k, k - 1)) <= smlnum) exit
                    tst = abs(h(k - 1, k - 1)) + abs(h(k, k))
                    if (tst == 0) then
                        if (k - 2 >= ilo) tst = tst + abs(h(k - 1, k - 2))
                        if (k + 1 <= ihi) tst = tst + abs(h(k + 1, k))
                    end if
                    if (abs(h(k, k - 1)) <= ulp*tst) then
                        ab = max(abs(h(k, k - 1)), abs(h(k - 1, k)))
                        ba = min(abs(h(k, k - 1)), abs(h(k - 1, k)))
                        aa = max(abs(h(k, k)), abs(h(k - 1, k - 1) - h(k, k)))
                        bb = min(abs(h(k, k)), abs(h(k - 1, k - 1) - h(k, k)))
                        s = aa + ab
                        if (ba*(ab/s) <= max(smlnum, ulp*(bb*(aa/s)))) exit
                    end if
                end do
                l = max(k, l)
                if (l > ilo) h(l, l - 1) = 0
                if (l >= i - 1) then
                    found = .true.
                    exit
                end if
                kdefl = kdefl + 1
                if (.not. wantt) then
                    i1 = l
                    i2 = i
                end if
                if (mod(kdefl, 2*kexsh) == 0) then
                    s = abs(h(i, i - 1)) + abs(h(i - 1, i - 2))
                    h11 = dat1*s + h(i, i)
                    h12 = dat2*s
                    h21 = s
                    h22 = h11
                else if (mod(kdefl, kexsh) == 0) then
                    s = abs(h(l + 1, l)) + abs(h(l + 2, l + 1))
                    h11 = dat1*s + h(l, l)
                    h12 = dat2*s
                    h21 = s
                    h22 = h11
                else
                    h11 = h(i - 1, i - 1)
                    h21 = h(i, i - 1)
                    h12 = h(i - 1, i)
                    h22 = h(i, i)
                end if
                s = abs(h11) + abs(h12) + abs(h21) + abs(h22)
                if (s == 0) then
                    rt1r = 0
                    rt1i = 0
                    rt2r = 0
                    rt2i = 0
                else
                    h11 = h11/s
                    h21 = h21/s
                    h12 = h12/s
                    h22 = h22/s
                    tr = (h11 + h22)/2
                    det = (h11 - tr)*(h22 - tr) - h12*h21
                    rtdisc = sqrt(abs(det))
                    if (det >= 0) then
                        rt1r = tr*s
                        rt2r = rt1r
                        rt1i = rtdisc*s
                        rt2i = -rt1i
                    else
                        ! real shifts: use the one closer to h22 twice
                        rt1r = tr + rtdisc
                        rt2r = tr - rtdisc
                        if (abs(rt1r - h22) <= abs(rt2r - h22)) then
                            rt1r = rt1r*s
                            rt2r = rt1r
                        else
                            rt2r = rt2r*s
                            rt1r = rt2r
                        end if
                        rt1i = 0
                        rt2i = 0
                    end if
                end if
                ! two consecutive small subdiagonal elements
                do m = i - 2, l, -1
                    h21s = h(m + 1, m)
                    s = abs(h(m, m) - rt2r) + abs(rt2i) + abs(h21s)
                    h21s = h(m + 1, m)/s
                    v(1) = h21s*h(m, m + 1) + (h(m, m) - rt1r)*((h(m, m) - rt2r)/s) - rt1i*(rt2i/s)
                    v(2) = h21s*(h(m, m) + h(m + 1, m + 1) - rt1r - rt2r)
                    v(3) = h21s*h(m + 2, m + 1)
                    s = abs(v(1)) + abs(v(2)) + abs(v(3))
                    v = v/s
                    if (m == l) exit
                    h00 = abs(h(m, m - 1))*(abs(v(2)) + abs(v(3)))
                    h01 = ulp*abs(v(1))*(abs(h(m - 1, m - 1)) + abs(h(m, m)) + abs(h(m + 1, m + 1)))
                    if (h00 <= h01) exit
                end do
                ! double-shift QR sweep
                do k = m, i - 1
                    nr = min(3, i - k + 1)
                    if (k > m) v(1:nr) = h(k:k + nr - 1, k - 1)
                    call fortnum_dlarfg(nr, v(1), v(2:), 1, t1)
                    if (k > m) then
                        h(k, k - 1) = v(1)
                        h(k + 1, k - 1) = 0
                        if (k < i - 1) h(k + 2, k - 1) = 0
                    else if (m > l) then
                        ! not -h(k, k-1): v(2) and v(3) may have underflowed
                        h(k, k - 1) = h(k, k - 1)*(1 - t1)
                    end if
                    v2 = v(2)
                    t2 = t1*v2
                    if (nr == 3) then
                        v3 = v(3)
                        t3 = t1*v3
                        do j = k, i2
                            sm = h(k, j) + v2*h(k + 1, j) + v3*h(k + 2, j)
                            h(k, j) = h(k, j) - sm*t1
                            h(k + 1, j) = h(k + 1, j) - sm*t2
                            h(k + 2, j) = h(k + 2, j) - sm*t3
                        end do
                        do j = i1, min(k + 3, i)
                            sm = h(j, k) + v2*h(j, k + 1) + v3*h(j, k + 2)
                            h(j, k) = h(j, k) - sm*t1
                            h(j, k + 1) = h(j, k + 1) - sm*t2
                            h(j, k + 2) = h(j, k + 2) - sm*t3
                        end do
                        if (wantz) then
                            do j = iloz, ihiz
                                sm = z(j, k) + v2*z(j, k + 1) + v3*z(j, k + 2)
                                z(j, k) = z(j, k) - sm*t1
                                z(j, k + 1) = z(j, k + 1) - sm*t2
                                z(j, k + 2) = z(j, k + 2) - sm*t3
                            end do
                        end if
                    else if (nr == 2) then
                        do j = k, i2
                            sm = h(k, j) + v2*h(k + 1, j)
                            h(k, j) = h(k, j) - sm*t1
                            h(k + 1, j) = h(k + 1, j) - sm*t2
                        end do
                        do j = i1, i
                            sm = h(j, k) + v2*h(j, k + 1)
                            h(j, k) = h(j, k) - sm*t1
                            h(j, k + 1) = h(j, k + 1) - sm*t2
                        end do
                        if (wantz) then
                            do j = iloz, ihiz
                                sm = z(j, k) + v2*z(j, k + 1)
                                z(j, k) = z(j, k) - sm*t1
                                z(j, k + 1) = z(j, k + 1) - sm*t2
                            end do
                        end if
                    end if
                end do
            end do
            if (.not. found) then
                info = i
                return
            end if
            if (l == i) then
                wr(i) = h(i, i)
                wi(i) = 0
            else if (l == i - 1) then
                ! a converged 2 x 2 block, in standard form
                call dlanv2(h(i - 1, i - 1), h(i - 1, i), h(i, i - 1), h(i, i), &
                    wr(i - 1), wi(i - 1), wr(i), wi(i), cs, sn)
                if (wantt) then
                    if (i2 > i) call drot_rows(i2 - i, h, ldh, i - 1, i, i + 1, cs, sn)
                    call drot_cols(i - i1 - 1, h, ldh, i1, i - 1, i, cs, sn)
                end if
                if (wantz) call drot_cols(nz, z, ldz, iloz, i - 1, i, cs, sn)
            end if
            kdefl = 0
            i = l - 1
        end do
    end subroutine dlahqr

    ! rows p and q of h over columns c0 .. c0+m-1: (x, y) := (c x + s y, c y - s x)
    subroutine drot_rows(m, h, ldh, p, q, c0, c, s)
        integer, intent(in) :: m, ldh, p, q, c0
        real(dp), intent(inout) :: h(ldh, *)
        real(dp), intent(in) :: c, s
        real(dp) :: x, y
        integer :: j
        do j = c0, c0 + m - 1
            x = h(p, j)
            y = h(q, j)
            h(p, j) = c*x + s*y
            h(q, j) = c*y - s*x
        end do
    end subroutine drot_rows

    ! columns p and q of h over rows r0 .. r0+m-1
    subroutine drot_cols(m, h, ldh, r0, p, q, c, s)
        integer, intent(in) :: m, ldh, r0, p, q
        real(dp), intent(inout) :: h(ldh, *)
        real(dp), intent(in) :: c, s
        real(dp) :: x, y
        integer :: j
        do j = r0, r0 + m - 1
            x = h(j, p)
            y = h(j, q)
            h(j, p) = c*x + s*y
            h(j, q) = c*y - s*x
        end do
    end subroutine drot_cols

    ! LAPACK's DLANV2: Schur factorization of a real 2 x 2 block in standard
    ! form (a = d and b c < 0 for complex eigenvalues, else c = 0)
    subroutine dlanv2(a, b, c, d, rt1r, rt1i, rt2r, rt2i, cs, sn)
        real(dp), intent(inout) :: a, b, c, d
        real(dp), intent(out) :: rt1r, rt1i, rt2r, rt2i, cs, sn
        real(dp), parameter :: multpl = 4
        real(dp) :: eps, safmn2, safmx2, temp, p, bcmax, bcmis, scl, z, tau, sigma
        real(dp) :: aa, bb, cc, dd, sab, sac, cs1, sn1
        integer :: count
        eps = ulp
        safmn2 = 2.0_dp**int(log(safmin/eps)/log(2.0_dp)/2)
        safmx2 = 1/safmn2
        if (c == 0) then
            cs = 1
            sn = 0
        else if (b == 0) then
            cs = 0
            sn = 1
            temp = d
            d = a
            a = temp
            b = -c
            c = 0
        else if (a - d == 0 .and. sign(1.0_dp, b) /= sign(1.0_dp, c)) then
            cs = 1
            sn = 0
        else
            temp = a - d
            p = temp/2
            bcmax = max(abs(b), abs(c))
            bcmis = min(abs(b), abs(c))*sign(1.0_dp, b)*sign(1.0_dp, c)
            scl = max(abs(p), bcmax)
            z = (p/scl)*p + (bcmax/scl)*bcmis
            if (z >= multpl*eps) then
                ! real eigenvalues
                z = p + sign(sqrt(scl)*sqrt(z), p)
                a = d + z
                d = d - (bcmax/z)*bcmis
                tau = hypot(c, z)
                cs = z/tau
                sn = c/tau
                b = b - c
                c = 0
            else
                ! complex or nearly equal real eigenvalues: equal diagonal
                count = 0
                sigma = b + c
                do
                    count = count + 1
                    scl = max(abs(temp), abs(sigma))
                    if (scl >= safmx2) then
                        sigma = sigma*safmn2
                        temp = temp*safmn2
                        if (count <= 20) cycle
                    end if
                    if (scl <= safmn2) then
                        sigma = sigma*safmx2
                        temp = temp*safmx2
                        if (count <= 20) cycle
                    end if
                    exit
                end do
                p = temp/2
                tau = hypot(sigma, temp)
                cs = sqrt((1 + abs(sigma)/tau)/2)
                sn = -(p/(tau*cs))*sign(1.0_dp, sigma)
                aa = a*cs + b*sn
                bb = -a*sn + b*cs
                cc = c*cs + d*sn
                dd = -c*sn + d*cs
                a = aa*cs + cc*sn
                b = bb*cs + dd*sn
                c = -aa*sn + cc*cs
                d = -bb*sn + dd*cs
                temp = (a + d)/2
                a = temp
                d = temp
                if (c /= 0) then
                    if (b /= 0) then
                        if (sign(1.0_dp, b) == sign(1.0_dp, c)) then
                            ! real eigenvalues: upper triangular
                            sab = sqrt(abs(b))
                            sac = sqrt(abs(c))
                            p = sign(sab*sac, c)
                            tau = 1/sqrt(abs(b + c))
                            a = temp + p
                            d = temp - p
                            b = b - c
                            c = 0
                            cs1 = sab*tau
                            sn1 = sac*tau
                            temp = cs*cs1 - sn*sn1
                            sn = cs*sn1 + sn*cs1
                            cs = temp
                        end if
                    else
                        b = -c
                        c = 0
                        temp = cs
                        cs = -sn
                        sn = temp
                    end if
                end if
            end if
        end if
        rt1r = a
        rt2r = d
        if (c == 0) then
            rt1i = 0
            rt2i = 0
        else
            rt1i = sqrt(abs(b))*sqrt(abs(c))
            rt2i = -rt1i
        end if
    end subroutine dlanv2

    ! -------------------------------------------- Hessenberg QR, complex

    subroutine fortnum_zhseqr(job, compz, n, ilo, ihi, h, ldh, w, z, ldz, work, lwork, info)
        character, intent(in) :: job, compz
        integer, intent(in) :: n, ilo, ihi, ldh, ldz, lwork
        complex(dp), intent(inout) :: h(ldh, *), z(ldz, *)
        complex(dp), intent(out) :: w(*), work(*)
        integer, intent(out) :: info
        logical :: wantt, wantz
        integer :: i, j
        call hseqr_check(job, compz, n, ilo, ihi, ldh, ldz, lwork, info)
        if (info == -11) info = -10
        if (info == -13) info = -12
        if (info /= 0) return
        work(1) = max(1, n)
        if (lwork == -1 .or. n == 0) return
        wantt = any(job == ['S', 's'])
        wantz = .not. any(compz == ['N', 'n'])
        if (any(compz == ['I', 'i'])) then
            do j = 1, n
                z(1:n, j) = czero
                z(j, j) = cone
            end do
        end if
        do i = 1, ilo - 1
            w(i) = h(i, i)
        end do
        do i = ihi + 1, n
            w(i) = h(i, i)
        end do
        call zlahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, 1, n, z, ldz, info)
        if (wantt .and. info == 0) then
            do j = 1, n - 2
                h(j + 2:n, j) = czero
            end do
        end if
    end subroutine fortnum_zhseqr

    real(dp) elemental function cabs1(x)
        complex(dp), intent(in) :: x
        cabs1 = abs(real(x, dp)) + abs(aimag(x))
    end function cabs1

    ! LAPACK's ZLAHQR: single-shift QR on H(ilo:ihi, ilo:ihi), real subdiagonal
    subroutine zlahqr(wantt, wantz, n, ilo, ihi, h, ldh, w, iloz, ihiz, z, ldz, info)
        logical, intent(in) :: wantt, wantz
        integer, intent(in) :: n, ilo, ihi, ldh, iloz, ihiz, ldz
        complex(dp), intent(inout) :: h(ldh, *), z(ldz, *), w(*)
        integer, intent(out) :: info
        real(dp), parameter :: dat1 = 0.75_dp
        integer, parameter :: kexsh = 10
        complex(dp) :: sc, t, u, x, y, h11, h22, h11s, t1, v2, sm, temp
        complex(dp) :: v(2)
        real(dp) :: smlnum, tst, ab, ba, aa, bb, s, sx, h21, h10, rtemp
        integer :: i, j, k, l, m, its, itmax, kdefl, i1, i2, jlo, jhi, nh, nz
        logical :: found
        info = 0
        if (n == 0) return
        if (ilo == ihi) then
            w(ilo) = h(ilo, ilo)
            return
        end if
        do j = ilo, ihi - 3
            h(j + 2, j) = czero
            h(j + 3, j) = czero
        end do
        if (ilo <= ihi - 2) h(ihi, ihi - 2) = czero
        if (wantt) then
            jlo = 1
            jhi = n
        else
            jlo = ilo
            jhi = ihi
        end if
        ! make the subdiagonal real
        do i = ilo + 1, ihi
            if (aimag(h(i, i - 1)) /= 0) then
                sc = h(i, i - 1)/cabs1(h(i, i - 1))
                sc = conjg(sc)/abs(sc)
                h(i, i - 1) = abs(h(i, i - 1))
                h(i, i:jhi) = h(i, i:jhi)*sc
                h(jlo:min(jhi, i + 1), i) = h(jlo:min(jhi, i + 1), i)*conjg(sc)
                if (wantz) z(iloz:ihiz, i) = z(iloz:ihiz, i)*conjg(sc)
            end if
        end do
        nh = ihi - ilo + 1
        nz = ihiz - iloz + 1
        smlnum = safmin*(real(nh, dp)/ulp)
        i1 = 1
        i2 = n
        itmax = 30*max(10, nh)
        kdefl = 0
        i = ihi
        do while (i >= ilo)
            l = ilo
            found = .false.
            do its = 0, itmax
                do k = i, l + 1, -1
                    if (cabs1(h(k, k - 1)) <= smlnum) exit
                    tst = cabs1(h(k - 1, k - 1)) + cabs1(h(k, k))
                    if (tst == 0) then
                        if (k - 2 >= ilo) tst = tst + abs(real(h(k - 1, k - 2), dp))
                        if (k + 1 <= ihi) tst = tst + abs(real(h(k + 1, k), dp))
                    end if
                    if (abs(real(h(k, k - 1), dp)) <= ulp*tst) then
                        ab = max(cabs1(h(k, k - 1)), cabs1(h(k - 1, k)))
                        ba = min(cabs1(h(k, k - 1)), cabs1(h(k - 1, k)))
                        aa = max(cabs1(h(k, k)), cabs1(h(k - 1, k - 1) - h(k, k)))
                        bb = min(cabs1(h(k, k)), cabs1(h(k - 1, k - 1) - h(k, k)))
                        s = aa + ab
                        if (ba*(ab/s) <= max(smlnum, ulp*(bb*(aa/s)))) exit
                    end if
                end do
                l = max(k, l)
                if (l > ilo) h(l, l - 1) = czero
                if (l >= i) then
                    found = .true.
                    exit
                end if
                kdefl = kdefl + 1
                if (.not. wantt) then
                    i1 = l
                    i2 = i
                end if
                if (mod(kdefl, 2*kexsh) == 0) then
                    s = dat1*abs(real(h(i, i - 1), dp))
                    t = s + h(i, i)
                else if (mod(kdefl, kexsh) == 0) then
                    s = dat1*abs(real(h(l + 1, l), dp))
                    t = s + h(l, l)
                else
                    ! Wilkinson shift
                    t = h(i, i)
                    u = sqrt(h(i - 1, i))*sqrt(h(i, i - 1))
                    s = cabs1(u)
                    if (s /= 0) then
                        x = (h(i - 1, i - 1) - t)/2
                        sx = cabs1(x)
                        s = max(s, cabs1(x))
                        y = s*sqrt((x/s)**2 + (u/s)**2)
                        if (sx > 0) then
                            if (real(x/sx, dp)*real(y, dp) + aimag(x/sx)*aimag(y) < 0) y = -y
                        end if
                        t = t - u*(u/(x + y))
                    end if
                end if
                ! two consecutive small subdiagonal elements
                m = l
                do j = i - 1, l + 1, -1
                    h11 = h(j, j)
                    h22 = h(j + 1, j + 1)
                    h11s = h11 - t
                    h21 = real(h(j + 1, j), dp)
                    s = cabs1(h11s) + abs(h21)
                    h11s = h11s/s
                    h21 = h21/s
                    h10 = real(h(j, j - 1), dp)
                    if (abs(h10)*abs(h21) <= ulp*(cabs1(h11s)*(cabs1(h11) + cabs1(h22)))) then
                        m = j
                        exit
                    end if
                end do
                h11 = h(m, m)
                h11s = h11 - t
                h21 = real(h(m + 1, m), dp)
                s = cabs1(h11s) + abs(h21)
                v(1) = h11s/s
                v(2) = h21/s
                ! single-shift QR sweep
                do k = m, i - 1
                    if (k > m) v = h(k:k + 1, k - 1)
                    call fortnum_zlarfg(2, v(1), v(2:), 1, t1)
                    if (k > m) then
                        h(k, k - 1) = v(1)
                        h(k + 1, k - 1) = czero
                    end if
                    v2 = v(2)
                    do j = k, i2
                        sm = conjg(t1)*(h(k, j) + conjg(v2)*h(k + 1, j))
                        h(k, j) = h(k, j) - sm
                        h(k + 1, j) = h(k + 1, j) - sm*v2
                    end do
                    do j = i1, min(k + 2, i)
                        sm = t1*(h(j, k) + v2*h(j, k + 1))
                        h(j, k) = h(j, k) - sm
                        h(j, k + 1) = h(j, k + 1) - sm*conjg(v2)
                    end do
                    if (wantz) then
                        do j = iloz, ihiz
                            sm = t1*(z(j, k) + v2*z(j, k + 1))
                            z(j, k) = z(j, k) - sm
                            z(j, k + 1) = z(j, k + 1) - sm*conjg(v2)
                        end do
                    end if
                    if (k == m .and. m > l) then
                        ! keep h(m, m-1) real after a sweep started below l
                        temp = cone - t1
                        temp = temp/abs(temp)
                        h(m + 1, m) = h(m + 1, m)*conjg(temp)
                        if (m + 2 <= i) h(m + 2, m + 1) = h(m + 2, m + 1)*temp
                        do j = m, i
                            if (j /= m + 1) then
                                if (i2 > j) h(j, j + 1:i2) = h(j, j + 1:i2)*temp
                                h(i1:j - 1, j) = h(i1:j - 1, j)*conjg(temp)
                                if (wantz) z(iloz:ihiz, j) = z(iloz:ihiz, j)*conjg(temp)
                            end if
                        end do
                    end if
                end do
                ! keep h(i, i-1) real
                temp = h(i, i - 1)
                if (aimag(temp) /= 0) then
                    rtemp = abs(temp)
                    h(i, i - 1) = rtemp
                    temp = temp/rtemp
                    if (i2 > i) h(i, i + 1:i2) = h(i, i + 1:i2)*conjg(temp)
                    h(i1:i - 1, i) = h(i1:i - 1, i)*temp
                    if (wantz) z(iloz:ihiz, i) = z(iloz:ihiz, i)*temp
                end if
            end do
            if (.not. found) then
                info = i
                return
            end if
            w(i) = h(i, i)
            kdefl = 0
            i = l - 1
        end do
    end subroutine zlahqr

    ! ---------------------------------------------------------- drivers

    !> Eigenvalues wr + i wi and, with jobvl/jobvr 'V', left and right
    !> eigenvectors of a real general matrix, as LAPACK's DGEEV: a complex
    !> pair has wi > 0 first and its vectors in two columns (real and
    !> imaginary part); each vector has Euclidean norm 1 and its largest
    !> component real.
    subroutine fortnum_dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: n, lda, ldvl, ldvr, lwork
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
        integer, intent(out) :: info
        logical :: wantvl, wantvr, scalea
        real(dp), allocatable :: tau(:), q(:, :), scl(:), wk(:)
        complex(dp), allocatable :: x(:)
        real(dp) :: anrm, smlnum, bignum, cscale, nrm
        integer :: ilo, ihi, minwrk, j, k, ierr
        info = 0
        wantvl = any(jobvl == ['V', 'v'])
        wantvr = any(jobvr == ['V', 'v'])
        minwrk = max(1, merge(4*n, 3*n, wantvl .or. wantvr))
        if (.not. wantvl .and. .not. any(jobvl == ['N', 'n'])) then
            info = -1
        else if (.not. wantvr .and. .not. any(jobvr == ['N', 'n'])) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
            info = -9
        else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
            info = -11
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -13
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. n == 0) return
        ! scale A into [smlnum, bignum] as LAPACK does
        smlnum = sqrt(safmin)/ulp
        bignum = 1/smlnum
        anrm = 0
        do j = 1, n
            anrm = max(anrm, maxval(abs(a(1:n, j))))
        end do
        if (anrm /= anrm) then
            info = -5
            return
        end if
        scalea = .false.
        if (anrm > 0 .and. anrm < smlnum) then
            scalea = .true.
            cscale = smlnum
        else if (anrm > bignum) then
            scalea = .true.
            cscale = bignum
        end if
        if (scalea) then
            do j = 1, n
                a(1:n, j) = a(1:n, j)*(cscale/anrm)
            end do
        end if
        allocate (tau(max(1, n - 1)), scl(n), wk(max(1, n)))
        call fortnum_dgebal('B', n, a, lda, ilo, ihi, scl, ierr)
        call fortnum_dgehrd(n, ilo, ihi, a, lda, tau, wk, size(wk), ierr)
        if (wantvl .or. wantvr) then
            allocate (q(n, n))
            do j = 1, n
                q(1:n, j) = a(1:n, j)
            end do
            call fortnum_dorghr(n, ilo, ihi, q, n, tau, wk, size(wk), ierr)
            call fortnum_dhseqr('S', 'V', n, ilo, ihi, a, lda, wr, wi, q, n, wk, size(wk), info)
        else
            allocate (q(1, 1))
            call fortnum_dhseqr('E', 'N', n, ilo, ihi, a, lda, wr, wi, q, 1, wk, size(wk), info)
        end if
        if (info == 0 .and. (wantvl .or. wantvr)) then
            allocate (x(n))
            j = 1
            do while (j <= n)
                k = merge(2, 1, wi(j) /= 0)
                if (wantvr) then
                    call dtrevc_right(n, a, lda, j, k, cmplx(wr(j), wi(j), dp), x)
                    call store_real_vector(n, q, x, vr(:, j:j + k - 1), k)
                end if
                if (wantvl) then
                    call dtrevc_left(n, a, lda, j, k, cmplx(wr(j), -wi(j), dp), x)
                    call store_real_vector(n, q, x, vl(:, j:j + k - 1), k)
                end if
                j = j + k
            end do
            if (wantvr) call fortnum_dgebak('B', 'R', n, ilo, ihi, scl, n, vr, ldvr, ierr)
            if (wantvl) call fortnum_dgebak('B', 'L', n, ilo, ihi, scl, n, vl, ldvl, ierr)
            j = 1
            do while (j <= n)
                k = merge(2, 1, wi(j) /= 0)
                if (wantvr) call normalize_real(n, vr(:, j:j + k - 1), k)
                if (wantvl) call normalize_real(n, vl(:, j:j + k - 1), k)
                j = j + k
            end do
        end if
        if (scalea) then
            nrm = anrm/cscale
            wr(info + 1:n) = wr(info + 1:n)*nrm
            wi(info + 1:n) = wi(info + 1:n)*nrm
            if (info > 0) then
                wr(1:ilo - 1) = wr(1:ilo - 1)*nrm
                wi(1:ilo - 1) = wi(1:ilo - 1)*nrm
            end if
        end if
    contains
        subroutine store_real_vector(n, q, x, v, k)
            integer, intent(in) :: n, k
            real(dp), intent(in) :: q(:, :)
            complex(dp), intent(in) :: x(:)
            real(dp), intent(out) :: v(:, :)
            v(1:n, 1) = matmul(q, real(x, dp))
            if (k == 2) v(1:n, 2) = matmul(q, aimag(x))
        end subroutine store_real_vector
    end subroutine fortnum_dgeev

    ! Euclidean norm 1 and, for a pair, the largest component real
    subroutine normalize_real(n, v, k)
        integer, intent(in) :: n, k
        real(dp), intent(inout) :: v(:, :)
        real(dp) :: s, c, sn, r, x, y
        integer :: i, imax
        if (k == 1) then
            s = norm2(v(1:n, 1))
            if (s > 0) v(1:n, 1) = v(1:n, 1)/s
            return
        end if
        s = hypot(norm2(v(1:n, 1)), norm2(v(1:n, 2)))
        if (s > 0) v(1:n, 1:2) = v(1:n, 1:2)/s
        imax = maxloc(v(1:n, 1)**2 + v(1:n, 2)**2, 1)
        r = hypot(v(imax, 1), v(imax, 2))
        if (r == 0) return
        c = v(imax, 1)/r
        sn = v(imax, 2)/r
        do i = 1, n
            x = v(i, 1)
            y = v(i, 2)
            v(i, 1) = c*x + sn*y
            v(i, 2) = c*y - sn*x
        end do
        v(imax, 2) = 0
    end subroutine normalize_real

    ! diagonal block starting at j, of size k (2 for a complex pair), in T
    pure integer function block_size(n, t, ldt, j) result(k)
        integer, intent(in) :: n, ldt, j
        real(dp), intent(in) :: t(ldt, *)
        k = 1
        if (j < n) then
            if (t(j + 1, j) /= 0) k = 2
        end if
    end function block_size

    ! Right eigenvector of the quasi-triangular T for eigenvalue lam of the
    ! diagonal block (j, k): back substitution (T - lam I) x = 0 in complex
    ! arithmetic, small pivots perturbed and growth rescaled as xTREVC does.
    subroutine dtrevc_right(n, t, ldt, j, k, lam, x)
        integer, intent(in) :: n, ldt, j, k
        real(dp), intent(in) :: t(ldt, *)
        complex(dp), intent(in) :: lam
        complex(dp), intent(out) :: x(:)
        real(dp) :: smin, bignum
        complex(dp) :: rhs(2), m(2, 2)
        integer :: r, b, c
        bignum = 1/(safmin*n/ulp)
        smin = max(ulp*abs(lam), safmin*n/ulp)
        x = czero
        if (k == 1) then
            x(j) = cone
        else
            ! null vector of the standardized block [a b; c a]
            if (abs(t(j, j + 1)) >= abs(t(j + 1, j))) then
                x(j) = t(j, j + 1)
                x(j + 1) = lam - t(j, j)
            else
                x(j) = lam - t(j + 1, j + 1)
                x(j + 1) = t(j + 1, j)
            end if
            x(j:j + 1) = x(j:j + 1)/maxval(abs(x(j:j + 1)))
        end if
        r = j - 1
        do while (r >= 1)
            b = 1
            if (r > 1) then
                if (t(r, r - 1) /= 0) b = 2
            end if
            c = r - b + 1
            rhs(1:b) = czero
            rhs(1:b) = -matmul(t(c:r, r + 1:j + k - 1), x(r + 1:j + k - 1))
            m(1:b, 1:b) = t(c:r, c:r)
            m(1, 1) = m(1, 1) - lam
            if (b == 2) m(2, 2) = m(2, 2) - lam
            call solve_small(b, m, rhs, smin)
            x(c:r) = rhs(1:b)
            call rescale(x, bignum)
            r = c - 1
        end do
    end subroutine dtrevc_right

    ! Left eigenvector: (T^T - mu I) w = 0 by forward substitution
    subroutine dtrevc_left(n, t, ldt, j, k, mu, x)
        integer, intent(in) :: n, ldt, j, k
        real(dp), intent(in) :: t(ldt, *)
        complex(dp), intent(in) :: mu
        complex(dp), intent(out) :: x(:)
        real(dp) :: smin, bignum
        complex(dp) :: rhs(2), m(2, 2)
        integer :: r, b, e
        bignum = 1/(safmin*n/ulp)
        smin = max(ulp*abs(mu), safmin*n/ulp)
        x = czero
        if (k == 1) then
            x(j) = cone
        else
            ! null vector of the transposed block [a c; b a]
            if (abs(t(j + 1, j)) >= abs(t(j, j + 1))) then
                x(j) = t(j + 1, j)
                x(j + 1) = mu - t(j, j)
            else
                x(j) = mu - t(j + 1, j + 1)
                x(j + 1) = t(j, j + 1)
            end if
            x(j:j + 1) = x(j:j + 1)/maxval(abs(x(j:j + 1)))
        end if
        r = j + k
        do while (r <= n)
            b = block_size(n, t, ldt, r)
            e = r + b - 1
            rhs(1:b) = -matmul(transpose(t(j:r - 1, r:e)), x(j:r - 1))
            m(1:b, 1:b) = transpose(t(r:e, r:e))
            m(1, 1) = m(1, 1) - mu
            if (b == 2) m(2, 2) = m(2, 2) - mu
            call solve_small(b, m, rhs, smin)
            x(r:e) = rhs(1:b)
            call rescale(x, bignum)
            r = e + 1
        end do
    end subroutine dtrevc_left

    ! m x = rhs for b = 1 or 2 by Gaussian elimination with pivoting;
    ! pivots below smin are raised to smin
    subroutine solve_small(b, m, rhs, smin)
        integer, intent(in) :: b
        complex(dp), intent(inout) :: m(2, 2), rhs(2)
        real(dp), intent(in) :: smin
        complex(dp) :: f, tmp(2)
        if (b == 1) then
            if (abs(m(1, 1)) < smin) m(1, 1) = smin
            rhs(1) = rhs(1)/m(1, 1)
            return
        end if
        if (abs(m(2, 1)) > abs(m(1, 1))) then
            tmp = m(1, :)
            m(1, :) = m(2, :)
            m(2, :) = tmp
            f = rhs(1)
            rhs(1) = rhs(2)
            rhs(2) = f
        end if
        if (abs(m(1, 1)) < smin) m(1, 1) = smin
        f = m(2, 1)/m(1, 1)
        m(2, 2) = m(2, 2) - f*m(1, 2)
        rhs(2) = rhs(2) - f*rhs(1)
        if (abs(m(2, 2)) < smin) m(2, 2) = smin
        rhs(2) = rhs(2)/m(2, 2)
        rhs(1) = (rhs(1) - m(1, 2)*rhs(2))/m(1, 1)
    end subroutine solve_small

    ! keep the vector far from overflow; only its direction matters
    subroutine rescale(x, bignum)
        complex(dp), intent(inout) :: x(:)
        real(dp), intent(in) :: bignum
        real(dp) :: s
        s = maxval(cabs1(x))
        if (s > bignum) x = x/s
    end subroutine rescale

    !> Eigenvalues w and, with jobvl/jobvr 'V', left and right eigenvectors
    !> of a complex general matrix, as LAPACK's ZGEEV: each vector has
    !> Euclidean norm 1 and its largest component real.
    subroutine fortnum_zgeev(jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: n, lda, ldvl, ldvr, lwork
        complex(dp), intent(inout) :: a(lda, *)
        complex(dp), intent(out) :: w(*), vl(ldvl, *), vr(ldvr, *), work(*)
        real(dp), intent(out) :: rwork(*)
        integer, intent(out) :: info
        logical :: wantvl, wantvr, scalea
        complex(dp), allocatable :: tau(:), q(:, :), wk(:), x(:)
        real(dp), allocatable :: scl(:)
        real(dp) :: anrm, smlnum, bignum, cscale
        integer :: ilo, ihi, minwrk, j, ierr
        info = 0
        wantvl = any(jobvl == ['V', 'v'])
        wantvr = any(jobvr == ['V', 'v'])
        minwrk = max(1, 2*n)
        if (.not. wantvl .and. .not. any(jobvl == ['N', 'n'])) then
            info = -1
        else if (.not. wantvr .and. .not. any(jobvr == ['N', 'n'])) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldvl < 1 .or. (wantvl .and. ldvl < n)) then
            info = -8
        else if (ldvr < 1 .or. (wantvr .and. ldvr < n)) then
            info = -10
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -12
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. n == 0) return
        rwork(1) = 0
        smlnum = sqrt(safmin)/ulp
        bignum = 1/smlnum
        anrm = 0
        do j = 1, n
            anrm = max(anrm, maxval(abs(a(1:n, j))))
        end do
        if (anrm /= anrm) then
            info = -5
            return
        end if
        scalea = .false.
        if (anrm > 0 .and. anrm < smlnum) then
            scalea = .true.
            cscale = smlnum
        else if (anrm > bignum) then
            scalea = .true.
            cscale = bignum
        end if
        if (scalea) then
            do j = 1, n
                a(1:n, j) = a(1:n, j)*(cscale/anrm)
            end do
        end if
        allocate (tau(max(1, n - 1)), scl(n), wk(max(1, n)))
        call fortnum_zgebal('B', n, a, lda, ilo, ihi, scl, ierr)
        call fortnum_zgehrd(n, ilo, ihi, a, lda, tau, wk, size(wk), ierr)
        if (wantvl .or. wantvr) then
            allocate (q(n, n))
            do j = 1, n
                q(1:n, j) = a(1:n, j)
            end do
            call fortnum_zunghr(n, ilo, ihi, q, n, tau, wk, size(wk), ierr)
            call fortnum_zhseqr('S', 'V', n, ilo, ihi, a, lda, w, q, n, wk, size(wk), info)
        else
            allocate (q(1, 1))
            call fortnum_zhseqr('E', 'N', n, ilo, ihi, a, lda, w, q, 1, wk, size(wk), info)
        end if
        if (info == 0 .and. (wantvl .or. wantvr)) then
            allocate (x(n))
            do j = 1, n
                if (wantvr) then
                    call ztrevc_right(n, a, lda, j, x)
                    vr(1:n, j) = matmul(q, x)
                end if
                if (wantvl) then
                    call ztrevc_left(n, a, lda, j, x)
                    vl(1:n, j) = matmul(q, x)
                end if
            end do
            if (wantvr) call fortnum_zgebak('B', 'R', n, ilo, ihi, scl, n, vr, ldvr, ierr)
            if (wantvl) call fortnum_zgebak('B', 'L', n, ilo, ihi, scl, n, vl, ldvl, ierr)
            do j = 1, n
                if (wantvr) call normalize_complex(vr(1:n, j))
                if (wantvl) call normalize_complex(vl(1:n, j))
            end do
        end if
        if (scalea) then
            w(info + 1:n) = w(info + 1:n)*(anrm/cscale)
            if (info > 0) w(1:ilo - 1) = w(1:ilo - 1)*(anrm/cscale)
        end if
    end subroutine fortnum_zgeev

    subroutine normalize_complex(v)
        complex(dp), intent(inout) :: v(:)
        real(dp) :: s
        integer :: imax
        s = norm2_c(v)
        if (s == 0) return
        v = v/s
        imax = maxloc(real(v, dp)**2 + aimag(v)**2, 1)
        v = v*(conjg(v(imax))/abs(v(imax)))
        v(imax) = cmplx(real(v(imax), dp), 0.0_dp, dp)
    end subroutine normalize_complex

    subroutine ztrevc_right(n, t, ldt, j, x)
        integer, intent(in) :: n, ldt, j
        complex(dp), intent(in) :: t(ldt, *)
        complex(dp), intent(out) :: x(:)
        real(dp) :: smin, bignum
        complex(dp) :: d
        integer :: r
        bignum = 1/(safmin*n/ulp)
        smin = max(ulp*cabs1(t(j, j)), safmin*n/ulp)
        x = czero
        x(j) = cone
        do r = j - 1, 1, -1
            d = t(r, r) - t(j, j)
            if (cabs1(d) < smin) d = smin
            x(r) = -sum(t(r, r + 1:j)*x(r + 1:j))/d
            call rescale(x, bignum)
        end do
    end subroutine ztrevc_right

    ! (T^H - conj(t(j,j)) I) w = 0
    subroutine ztrevc_left(n, t, ldt, j, x)
        integer, intent(in) :: n, ldt, j
        complex(dp), intent(in) :: t(ldt, *)
        complex(dp), intent(out) :: x(:)
        real(dp) :: smin, bignum
        complex(dp) :: d
        integer :: r
        bignum = 1/(safmin*n/ulp)
        smin = max(ulp*cabs1(t(j, j)), safmin*n/ulp)
        x = czero
        x(j) = cone
        do r = j + 1, n
            d = conjg(t(r, r) - t(j, j))
            if (cabs1(d) < smin) d = smin
            x(r) = -sum(conjg(t(j:r - 1, r))*x(j:r - 1))/d
            call rescale(x, bignum)
        end do
    end subroutine ztrevc_left

end module fortnum_lapack_geev
