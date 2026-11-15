module fortnum_lapack_gees
    ! Complex Schur factorization with reference LAPACK arguments and info
    ! codes (prefix fortnum_): the driver ZGEES (Schur form A = Z T Z^H,
    ! optionally with selected eigenvalues moved to the top left of T),
    ! ZTREXC (moves one diagonal entry of a Schur form by unitary swaps of
    ! neighbours), ZTREVC (eigenvectors of an upper triangular matrix by back
    ! substitution, optionally multiplied by the Schur vectors) and the
    ! complex plane rotation ZLARTG.
    use fortnum_kinds, only: dp
    use fortnum_lapack_geev, only: fortnum_zgebal, fortnum_zgebak, fortnum_zgehrd, &
        fortnum_zunghr, fortnum_zhseqr
    implicit none
    private

    public :: fortnum_zgees, fortnum_ztrexc, fortnum_ztrevc, fortnum_zlartg
    public :: fortnum_zgees_select

    abstract interface
        !> SELECT of ZGEES: true for an eigenvalue to be moved to the top left
        logical function fortnum_zgees_select(w)
            import :: dp
            complex(dp), intent(in) :: w
        end function fortnum_zgees_select
    end interface

    real(dp), parameter :: ulp = epsilon(1.0_dp)
    real(dp), parameter :: safmin = tiny(1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp), cone = (1.0_dp, 0.0_dp)

contains

    !> The rotation [c s; -conj(s) c] [f; g] = [r; 0] with real c >= 0.
    pure subroutine fortnum_zlartg(f, g, c, s, r)
        complex(dp), intent(in) :: f, g
        real(dp), intent(out) :: c
        complex(dp), intent(out) :: s, r
        real(dp) :: af, d
        if (g == czero) then
            c = 1
            s = czero
            r = f
        else if (f == czero) then
            c = 0
            d = abs(g)
            s = conjg(g)/d
            r = d
        else
            af = abs(f)
            d = hypot(af, abs(g))
            c = af/d
            s = (f/af)*(conjg(g)/d)
            r = (f/af)*d
        end if
    end subroutine fortnum_zlartg

    ! [x; y] := [c s; -conj(s) c] [x; y], LAPACK's ZROT
    pure subroutine rot(x, y, c, s)
        complex(dp), intent(inout) :: x(:), y(:)
        real(dp), intent(in) :: c
        complex(dp), intent(in) :: s
        complex(dp) :: t(size(x))
        t = c*x + s*y
        y = c*y - conjg(s)*x
        x = t
    end subroutine rot

    real(dp) elemental function cabs1(x)
        complex(dp), intent(in) :: x
        cabs1 = abs(real(x, dp)) + abs(aimag(x))
    end function cabs1

    !> Moves the diagonal entry ifst of the upper triangular T to position
    !> ilst by unitary similarity, T := Q_s^H T Q_s; with compq 'V' the Schur
    !> vectors are updated, Q := Q Q_s.
    subroutine fortnum_ztrexc(compq, n, t, ldt, q, ldq, ifst, ilst, info)
        character, intent(in) :: compq
        integer, intent(in) :: n, ldt, ldq, ifst, ilst
        complex(dp), intent(inout) :: t(ldt, *), q(ldq, *)
        integer, intent(out) :: info
        logical :: wantq
        integer :: k, k1, k2, step
        complex(dp) :: t11, t22, s, r
        real(dp) :: c
        info = 0
        wantq = any(compq == ['V', 'v'])
        if (.not. wantq .and. .not. any(compq == ['N', 'n'])) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (ldt < max(1, n)) then
            info = -4
        else if (ldq < 1 .or. (wantq .and. ldq < max(1, n))) then
            info = -6
        else if ((ifst < 1 .or. ifst > n) .and. n > 0) then
            info = -7
        else if ((ilst < 1 .or. ilst > n) .and. n > 0) then
            info = -8
        end if
        if (info /= 0) return
        if (n <= 1 .or. ifst == ilst) return
        if (ifst < ilst) then
            k1 = ifst
            k2 = ilst - 1
            step = 1
        else
            k1 = ifst - 1
            k2 = ilst
            step = -1
        end if
        do k = k1, k2, step
            ! swap the neighbours t(k,k) and t(k+1,k+1)
            t11 = t(k, k)
            t22 = t(k + 1, k + 1)
            call fortnum_zlartg(t(k, k + 1), t22 - t11, c, s, r)
            if (k + 2 <= n) call rot(t(k, k + 2:n), t(k + 1, k + 2:n), c, s)
            call rot(t(1:k - 1, k), t(1:k - 1, k + 1), c, conjg(s))
            t(k, k) = t22
            t(k + 1, k + 1) = t11
            if (wantq) call rot(q(1:n, k), q(1:n, k + 1), c, conjg(s))
        end do
    end subroutine fortnum_ztrexc

    !> Right (side 'R') and/or left ('L', 'B' for both) eigenvectors of the
    !> upper triangular T: T x = lambda x, y^H T = lambda y^H. howmny 'A'
    !> returns those of T, 'S' those selected, in order, and 'B' multiplies
    !> the matrices given in vr and vl (the Schur vectors Q of A = Q T Q^H)
    !> by them, giving eigenvectors of A. Each vector is scaled so that its
    !> largest component has |re| + |im| = 1. m is the number of columns used.
    subroutine fortnum_ztrevc(side, howmny, select, n, t, ldt, vl, ldvl, vr, ldvr, &
            mm, m, work, rwork, info)
        character, intent(in) :: side, howmny
        logical, intent(in) :: select(*)
        integer, intent(in) :: n, ldt, ldvl, ldvr, mm
        complex(dp), intent(in) :: t(ldt, *)
        complex(dp), intent(inout) :: vl(ldvl, *), vr(ldvr, *)
        integer, intent(out) :: m, info
        complex(dp), intent(out) :: work(*)
        real(dp), intent(out) :: rwork(*)
        logical :: rightv, leftv, allv, over, somev
        complex(dp), allocatable :: x(:)
        real(dp) :: smlnum, bignum, smin
        complex(dp) :: d
        integer :: ki, j, is
        info = 0
        rightv = any(side == ['R', 'r', 'B', 'b'])
        leftv = any(side == ['L', 'l', 'B', 'b'])
        allv = any(howmny == ['A', 'a'])
        over = any(howmny == ['B', 'b'])
        somev = any(howmny == ['S', 's'])
        m = n
        if (somev .and. n > 0) m = count(select(1:n))
        if (.not. rightv .and. .not. leftv) then
            info = -1
        else if (.not. allv .and. .not. over .and. .not. somev) then
            info = -2
        else if (n < 0) then
            info = -4
        else if (ldt < max(1, n)) then
            info = -6
        else if (ldvl < 1 .or. (leftv .and. ldvl < n)) then
            info = -8
        else if (ldvr < 1 .or. (rightv .and. ldvr < n)) then
            info = -10
        else if (mm < m) then
            info = -11
        end if
        if (info /= 0) return
        if (n == 0) return
        work(1) = czero
        rwork(1) = 0
        smlnum = safmin*(n/ulp)
        bignum = (1 - ulp)/smlnum
        allocate (x(n))
        if (rightv) then
            is = m
            do ki = n, 1, -1
                if (somev) then
                    if (.not. select(ki)) cycle
                end if
                smin = max(ulp*cabs1(t(ki, ki)), smlnum)
                ! (T(1:ki-1,1:ki-1) - t(ki,ki)) x = -T(1:ki-1,ki), x(ki) = 1
                x = czero
                x(ki) = cone
                x(1:ki - 1) = -t(1:ki - 1, ki)
                do j = ki - 1, 1, -1
                    d = t(j, j) - t(ki, ki)
                    if (cabs1(d) < smin) d = smin
                    x(j) = (x(j) - sum(t(j, j + 1:ki - 1)*x(j + 1:ki - 1)))/d
                    if (cabs1(x(j)) > bignum) then
                        x(1:ki) = x(1:ki)/cabs1(x(j))
                    end if
                end do
                if (over) then
                    vr(1:n, ki) = x(ki)*vr(1:n, ki)
                    if (ki > 1) vr(1:n, ki) = vr(1:n, ki) + matmul(vr(1:n, 1:ki - 1), x(1:ki - 1))
                    vr(1:n, ki) = vr(1:n, ki)/maxval(cabs1(vr(1:n, ki)))
                else
                    vr(1:ki, is) = x(1:ki)/maxval(cabs1(x(1:ki)))
                    vr(ki + 1:n, is) = czero
                    is = is - 1
                end if
            end do
        end if
        if (leftv) then
            is = 1
            do ki = 1, n
                if (somev) then
                    if (.not. select(ki)) cycle
                end if
                smin = max(ulp*cabs1(t(ki, ki)), smlnum)
                ! (T(ki+1:n,ki+1:n) - t(ki,ki))^H x = -T(ki,ki+1:n)^H
                x = czero
                x(ki) = cone
                x(ki + 1:n) = -conjg(t(ki, ki + 1:n))
                do j = ki + 1, n
                    d = conjg(t(j, j) - t(ki, ki))
                    if (cabs1(d) < smin) d = smin
                    x(j) = (x(j) - sum(conjg(t(ki + 1:j - 1, j))*x(ki + 1:j - 1)))/d
                    if (cabs1(x(j)) > bignum) then
                        x(ki:n) = x(ki:n)/cabs1(x(j))
                    end if
                end do
                if (over) then
                    vl(1:n, ki) = x(ki)*vl(1:n, ki)
                    if (ki < n) vl(1:n, ki) = vl(1:n, ki) + matmul(vl(1:n, ki + 1:n), x(ki + 1:n))
                    vl(1:n, ki) = vl(1:n, ki)/maxval(cabs1(vl(1:n, ki)))
                else
                    vl(ki:n, is) = x(ki:n)/maxval(cabs1(x(ki:n)))
                    vl(1:ki - 1, is) = czero
                    is = is + 1
                end if
            end do
        end if
    end subroutine fortnum_ztrevc

    !> Schur factorization A = Z T Z^H of a complex general matrix: T (in a)
    !> upper triangular with the eigenvalues w on its diagonal, Z (in vs,
    !> jobvs 'V') unitary. With sort 'S' the sdim eigenvalues for which
    !> select is true lead the diagonal of T.
    subroutine fortnum_zgees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, &
            work, lwork, rwork, bwork, info)
        character, intent(in) :: jobvs, sort
        procedure(fortnum_zgees_select) :: select
        integer, intent(in) :: n, lda, ldvs, lwork
        complex(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: sdim, info
        complex(dp), intent(out) :: w(*), vs(ldvs, *), work(*)
        real(dp), intent(out) :: rwork(*)
        logical, intent(out) :: bwork(*)
        logical :: wantvs, wantst, scalea
        complex(dp), allocatable :: tau(:), wk(:)
        real(dp), allocatable :: scl(:)
        real(dp) :: anrm, smlnum, bignum, cscale
        integer :: ilo, ihi, minwrk, i, j, k, ierr
        character :: compz
        info = 0
        wantvs = any(jobvs == ['V', 'v'])
        wantst = any(sort == ['S', 's'])
        minwrk = max(1, 2*n)
        if (.not. wantvs .and. .not. any(jobvs == ['N', 'n'])) then
            info = -1
        else if (.not. wantst .and. .not. any(sort == ['N', 'n'])) then
            info = -2
        else if (n < 0) then
            info = -4
        else if (lda < max(1, n)) then
            info = -6
        else if (ldvs < 1 .or. (wantvs .and. ldvs < n)) then
            info = -10
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -12
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1) return
        sdim = 0
        if (n == 0) return
        rwork(1) = 0
        smlnum = sqrt(safmin)/ulp
        bignum = 1/smlnum
        anrm = 0
        do j = 1, n
            anrm = max(anrm, maxval(abs(a(1:n, j))))
        end do
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
        ! permute only: scaling would make the Schur vectors non-unitary
        allocate (tau(max(1, n - 1)), scl(n), wk(max(1, n)))
        call fortnum_zgebal('P', n, a, lda, ilo, ihi, scl, ierr)
        call fortnum_zgehrd(n, ilo, ihi, a, lda, tau, wk, size(wk), ierr)
        compz = 'N'
        if (wantvs) then
            compz = 'V'
            do j = 1, n
                vs(1:n, j) = a(1:n, j)
            end do
            call fortnum_zunghr(n, ilo, ihi, vs, ldvs, tau, wk, size(wk), ierr)
        end if
        call fortnum_zhseqr('S', compz, n, ilo, ihi, a, lda, w, vs, ldvs, wk, size(wk), info)
        if (wantst .and. info == 0) then
            if (scalea) w(1:n) = w(1:n)*(anrm/cscale)
            do i = 1, n
                bwork(i) = select(w(i))
            end do
            k = 0
            do i = 1, n
                if (.not. bwork(i)) cycle
                k = k + 1
                if (i /= k) call fortnum_ztrexc(compz, n, a, lda, vs, ldvs, i, k, ierr)
            end do
            sdim = k
            do i = 1, n
                w(i) = a(i, i)
            end do
        end if
        if (wantvs) call fortnum_zgebak('P', 'R', n, ilo, ihi, scl, n, vs, ldvs, ierr)
        if (scalea) then
            do j = 1, n
                a(1:j, j) = a(1:j, j)*(anrm/cscale)
                w(j) = a(j, j)
            end do
        end if
    end subroutine fortnum_zgees

end module fortnum_lapack_gees
