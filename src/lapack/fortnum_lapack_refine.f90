module fortnum_lapack_refine
    ! Matrix norms, condition estimates, iterative refinement with error
    ! bounds and the expert LU driver, with reference LAPACK arguments and
    ! info codes (prefix fortnum_): DLANGE/ZLANGE, DGECON/ZGECON, DGERFS,
    ! DSYRFS, DGEEQU, DLAQGE, DGESVX.
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgetrf, fortnum_dgetrs, fortnum_zgetrs
    use fortnum_lapack_ldlt, only: fortnum_dsytrs
    use fortnum_lapack_chol, only: norm1_inverse_estimate, znorm1_inverse_estimate
    implicit none
    private

    public :: fortnum_dlange, fortnum_zlange, fortnum_dgecon, fortnum_zgecon
    public :: fortnum_dgerfs, fortnum_dsyrfs, fortnum_dgeequ, fortnum_dlaqge, fortnum_dgesvx

    !> LAPACK's machine constants: eps = epsilon/2 (rounding unit), sfmin.
    real(dp), parameter :: UNIT_ROUNDOFF = epsilon(1.0_dp)/2
    real(dp), parameter :: SFMIN = tiny(1.0_dp)
    !> Refinement steps at most, as LAPACK.
    integer, parameter :: ITMAX = 5
    !> Kinds of op(A) in the refinement: general, transposed, symmetric.
    integer, parameter :: GE_N = 1, GE_T = 2, SY = 3

contains

    !> max |a|, 1-norm ('1', 'O'), infinity norm ('I') or Frobenius norm
    !> ('F', 'E') of an m x n matrix (DLANGE).
    real(dp) function fortnum_dlange(norm, m, n, a, lda, work) result(v)
        character, intent(in) :: norm
        integer, intent(in) :: m, n, lda
        real(dp), intent(in) :: a(lda, *)
        real(dp), intent(inout) :: work(*)
        integer :: j
        v = 0
        if (min(m, n) == 0) return
        select case (norm)
        case ('M', 'm')
            do j = 1, n
                v = max(v, maxval(abs(a(1:m, j))))
            end do
        case ('1', 'O', 'o')
            do j = 1, n
                v = max(v, sum(abs(a(1:m, j))))
            end do
        case ('I', 'i')
            work(1:m) = 0
            do j = 1, n
                work(1:m) = work(1:m) + abs(a(1:m, j))
            end do
            v = maxval(work(1:m))
        case ('F', 'f', 'E', 'e')
            v = norm2(a(1:m, 1:n))
        end select
    end function fortnum_dlange

    real(dp) function fortnum_zlange(norm, m, n, a, lda, work) result(v)
        character, intent(in) :: norm
        integer, intent(in) :: m, n, lda
        complex(dp), intent(in) :: a(lda, *)
        real(dp), intent(inout) :: work(*)
        integer :: j
        v = 0
        if (min(m, n) == 0) return
        select case (norm)
        case ('M', 'm')
            do j = 1, n
                v = max(v, maxval(abs(a(1:m, j))))
            end do
        case ('1', 'O', 'o')
            do j = 1, n
                v = max(v, sum(abs(a(1:m, j))))
            end do
        case ('I', 'i')
            work(1:m) = 0
            do j = 1, n
                work(1:m) = work(1:m) + abs(a(1:m, j))
            end do
            v = maxval(work(1:m))
        case ('F', 'f', 'E', 'e')
            v = sqrt(sum(real(a(1:m, 1:n), dp)**2 + aimag(a(1:m, 1:n))**2))
        end select
    end function fortnum_zlange

    !> Reciprocal condition number in the 1-norm or infinity norm from the LU
    !> factors of GETRF (DGECON); anorm is the norm of A.
    subroutine fortnum_dgecon(norm, n, a, lda, anorm, rcond, work, iwork, info)
        character, intent(in) :: norm
        integer, intent(in) :: n, lda
        real(dp), intent(in) :: a(lda, *), anorm
        real(dp), intent(out) :: rcond
        real(dp), intent(inout) :: work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: onenrm
        real(dp) :: ainvnm, v(max(1, n)), x(max(1, n))
        integer :: kase, isave(3), inf, ipiv(max(1, n)), i
        info = 0
        onenrm = norm == '1' .or. norm == 'O' .or. norm == 'o'
        if (.not. (onenrm .or. norm == 'I' .or. norm == 'i')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        else if (anorm < 0) then
            info = -5
        end if
        rcond = 0
        if (info /= 0) return
        if (n == 0) then
            rcond = 1
            return
        end if
        if (anorm == 0) return
        ! the factors without interchanges: the norm of inv(A) is that of inv(L U)
        ipiv(1:n) = [(i, i = 1, n)]
        kase = 0
        ainvnm = 0
        do
            call norm1_inverse_estimate(n, v, x, ainvnm, kase, isave)
            if (kase == 0) exit
            if ((kase == 1) .eqv. onenrm) then
                call fortnum_dgetrs('N', n, 1, a, lda, ipiv, x, n, inf)
            else
                call fortnum_dgetrs('T', n, 1, a, lda, ipiv, x, n, inf)
            end if
        end do
        if (ainvnm > 0) rcond = (1/ainvnm)/anorm
    end subroutine fortnum_dgecon

    subroutine fortnum_zgecon(norm, n, a, lda, anorm, rcond, work, rwork, info)
        character, intent(in) :: norm
        integer, intent(in) :: n, lda
        complex(dp), intent(in) :: a(lda, *)
        real(dp), intent(in) :: anorm
        real(dp), intent(out) :: rcond
        complex(dp), intent(inout) :: work(*)
        real(dp), intent(inout) :: rwork(*)
        integer, intent(out) :: info
        logical :: onenrm
        real(dp) :: ainvnm
        complex(dp) :: v(max(1, n)), x(max(1, n))
        integer :: kase, isave(3), inf, ipiv(max(1, n)), i
        info = 0
        onenrm = norm == '1' .or. norm == 'O' .or. norm == 'o'
        if (.not. (onenrm .or. norm == 'I' .or. norm == 'i')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        else if (anorm < 0) then
            info = -5
        end if
        rcond = 0
        if (info /= 0) return
        if (n == 0) then
            rcond = 1
            return
        end if
        if (anorm == 0) return
        ipiv(1:n) = [(i, i = 1, n)]
        kase = 0
        ainvnm = 0
        do
            call znorm1_inverse_estimate(n, v, x, ainvnm, kase, isave)
            if (kase == 0) exit
            if ((kase == 1) .eqv. onenrm) then
                call fortnum_zgetrs('N', n, 1, a, lda, ipiv, x, n, inf)
            else
                call fortnum_zgetrs('C', n, 1, a, lda, ipiv, x, n, inf)
            end if
        end do
        if (ainvnm > 0) rcond = (1/ainvnm)/anorm
    end subroutine fortnum_zgecon

    !> Iterative refinement and error bounds for op(A) X = B with the LU
    !> factors af, ipiv of GETRF (DGERFS).
    subroutine fortnum_dgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, &
            work, iwork, info)
        character, intent(in) :: trans
        integer, intent(in) :: n, nrhs, lda, ldaf, ldb, ldx
        real(dp), intent(in) :: a(lda, *), af(ldaf, *), b(ldb, *)
        integer, intent(in) :: ipiv(*)
        real(dp), intent(inout) :: x(ldx, *)
        real(dp), intent(out) :: ferr(*), berr(*)
        real(dp), intent(inout) :: work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: notran
        info = 0
        notran = trans == 'N' .or. trans == 'n'
        if (.not. (notran .or. index('TtCc', trans) > 0)) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (nrhs < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldaf < max(1, n)) then
            info = -7
        else if (ldb < max(1, n)) then
            info = -10
        else if (ldx < max(1, n)) then
            info = -12
        end if
        if (info /= 0) return
        call refine(merge(GE_N, GE_T, notran), .false., n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, &
            ferr, berr)
    end subroutine fortnum_dgerfs

    !> Iterative refinement and error bounds for A X = B, A symmetric, with
    !> the factors of DSYTRF (DSYRFS).
    subroutine fortnum_dsyrfs(uplo, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, &
            work, iwork, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, nrhs, lda, ldaf, ldb, ldx
        real(dp), intent(in) :: a(lda, *), af(ldaf, *), b(ldb, *)
        integer, intent(in) :: ipiv(*)
        real(dp), intent(inout) :: x(ldx, *)
        real(dp), intent(out) :: ferr(*), berr(*)
        real(dp), intent(inout) :: work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: upper
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (nrhs < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldaf < max(1, n)) then
            info = -7
        else if (ldb < max(1, n)) then
            info = -10
        else if (ldx < max(1, n)) then
            info = -12
        end if
        if (info /= 0) return
        call refine(SY, upper, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr)
    end subroutine fortnum_dsyrfs

    ! LAPACK's refinement loop and forward error bound (xGERFS/xSYRFS) for
    ! op(A) X = B: kind GE_N (A), GE_T (A^T) with LU factors, or SY (A
    ! symmetric, uplo upper) with Bunch-Kaufman factors.
    subroutine refine(kind, upper, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr)
        integer, intent(in) :: kind, n, nrhs, lda, ldaf, ldb, ldx
        logical, intent(in) :: upper
        real(dp), intent(in) :: a(lda, *), af(ldaf, *), b(ldb, *)
        integer, intent(in) :: ipiv(*)
        real(dp), intent(inout) :: x(ldx, *)
        real(dp), intent(out) :: ferr(*), berr(*)
        real(dp) :: r(n), s(n), w(n), v(n), y(n), lstres, safe1, safe2, est, xnorm
        integer :: j, i, count, kase, isave(3), nz
        if (n == 0 .or. nrhs == 0) then
            ferr(1:max(0, nrhs)) = 0
            berr(1:max(0, nrhs)) = 0
            return
        end if
        nz = n + 1
        safe1 = nz*SFMIN
        safe2 = safe1/UNIT_ROUNDOFF
        do j = 1, nrhs
            count = 1
            lstres = 3
            do
                call residual(kind, upper, n, a, lda, b(1:n, j), x(1:n, j), r, s)
                berr(j) = 0
                do i = 1, n
                    if (s(i) > safe2) then
                        berr(j) = max(berr(j), abs(r(i))/s(i))
                    else
                        berr(j) = max(berr(j), (abs(r(i)) + safe1)/(s(i) + safe1))
                    end if
                end do
                if (berr(j) > UNIT_ROUNDOFF .and. 2*berr(j) <= lstres .and. count <= ITMAX) then
                    y = r
                    call solve(kind, upper, n, af, ldaf, ipiv, y, .false.)
                    x(1:n, j) = x(1:n, j) + y
                    lstres = berr(j)
                    count = count + 1
                else
                    exit
                end if
            end do
            ! bound || inv(op(A)) diag(w) ||_inf with w = |r| + nz eps (|A||x| + |b|)
            do i = 1, n
                if (s(i) > safe2) then
                    w(i) = abs(r(i)) + nz*UNIT_ROUNDOFF*s(i)
                else
                    w(i) = abs(r(i)) + nz*UNIT_ROUNDOFF*s(i) + safe1
                end if
            end do
            kase = 0
            est = 0
            do
                call norm1_inverse_estimate(n, v, y, est, kase, isave)
                if (kase == 0) exit
                if (kase == 1) then
                    call solve(kind, upper, n, af, ldaf, ipiv, y, .true.)
                    y = w*y
                else
                    y = w*y
                    call solve(kind, upper, n, af, ldaf, ipiv, y, .false.)
                end if
            end do
            xnorm = maxval(abs(x(1:n, j)))
            ferr(j) = est
            if (xnorm /= 0) ferr(j) = ferr(j)/xnorm
        end do
    end subroutine refine

    ! r := b - op(A) x and s := |b| + |op(A)| |x|
    subroutine residual(kind, upper, n, a, lda, bj, xj, r, s)
        integer, intent(in) :: kind, n, lda
        logical, intent(in) :: upper
        real(dp), intent(in) :: a(lda, *), bj(:), xj(:)
        real(dp), intent(out) :: r(:), s(:)
        integer :: i, k
        real(dp) :: aik
        select case (kind)
        case (GE_N)
            r = bj - matmul(a(1:n, 1:n), xj)
            s = abs(bj) + matmul(abs(a(1:n, 1:n)), abs(xj))
        case (GE_T)
            do i = 1, n
                r(i) = bj(i) - dot_product(a(1:n, i), xj)
                s(i) = abs(bj(i)) + dot_product(abs(a(1:n, i)), abs(xj))
            end do
        case default
            do i = 1, n
                r(i) = bj(i)
                s(i) = abs(bj(i))
                do k = 1, n
                    if ((i <= k) .eqv. upper) then
                        aik = a(i, k)
                    else
                        aik = a(k, i)
                    end if
                    r(i) = r(i) - aik*xj(k)
                    s(i) = s(i) + abs(aik*xj(k))
                end do
            end do
        end select
    end subroutine residual

    ! y := inv(op(A)) y, or inv(op(A))^T y when transposed
    subroutine solve(kind, upper, n, af, ldaf, ipiv, y, transposed)
        integer, intent(in) :: kind, n, ldaf
        logical, intent(in) :: upper, transposed
        real(dp), intent(in) :: af(ldaf, *)
        integer, intent(in) :: ipiv(*)
        real(dp), intent(inout) :: y(:)
        integer :: inf
        select case (kind)
        case (GE_N)
            call fortnum_dgetrs(merge('T', 'N', transposed), n, 1, af, ldaf, ipiv, y, n, inf)
        case (GE_T)
            call fortnum_dgetrs(merge('N', 'T', transposed), n, 1, af, ldaf, ipiv, y, n, inf)
        case default
            call fortnum_dsytrs(merge('U', 'L', upper), n, 1, af, ldaf, ipiv, y, n, inf)
        end select
    end subroutine solve

    !> Row and column scalings that equilibrate A (DGEEQU): r(i) = 1 / max_j
    !> |a(i, j)|, c(j) = 1 / max_i |a(i, j)| r(i), clamped to the safe range;
    !> info = i (or m + j) for an exactly zero row (column).
    subroutine fortnum_dgeequ(m, n, a, lda, r, c, rowcnd, colcnd, amax, info)
        integer, intent(in) :: m, n, lda
        real(dp), intent(in) :: a(lda, *)
        real(dp), intent(out) :: r(*), c(*), rowcnd, colcnd, amax
        integer, intent(out) :: info
        real(dp) :: smlnum, bignum, rcmin, rcmax
        integer :: i, j
        info = 0
        if (m < 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, m)) then
            info = -4
        end if
        rowcnd = 1
        colcnd = 1
        amax = 0
        if (info /= 0 .or. m == 0 .or. n == 0) return
        smlnum = SFMIN
        bignum = 1/smlnum
        r(1:m) = 0
        do j = 1, n
            r(1:m) = max(r(1:m), abs(a(1:m, j)))
        end do
        rcmin = minval(r(1:m))
        rcmax = maxval(r(1:m))
        amax = rcmax
        if (rcmin == 0) then
            do i = 1, m
                if (r(i) == 0) then
                    info = i
                    return
                end if
            end do
        end if
        r(1:m) = 1/min(max(r(1:m), smlnum), bignum)
        rowcnd = max(rcmin, smlnum)/min(rcmax, bignum)
        do j = 1, n
            c(j) = maxval(abs(a(1:m, j))*r(1:m))
        end do
        rcmin = minval(c(1:n))
        rcmax = maxval(c(1:n))
        if (rcmin == 0) then
            do j = 1, n
                if (c(j) == 0) then
                    info = m + j
                    return
                end if
            end do
        end if
        c(1:n) = 1/min(max(c(1:n), smlnum), bignum)
        colcnd = max(rcmin, smlnum)/min(rcmax, bignum)
    end subroutine fortnum_dgeequ

    !> Apply the scalings of DGEEQU when they are needed (DLAQGE): rows when
    !> rowcnd < 0.1 or amax is near the limits of the range, columns when
    !> colcnd < 0.1; equed reports N, R, C or B.
    subroutine fortnum_dlaqge(m, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
        integer, intent(in) :: m, n, lda
        real(dp), intent(inout) :: a(lda, *)
        real(dp), intent(in) :: r(*), c(*), rowcnd, colcnd, amax
        character, intent(out) :: equed
        real(dp), parameter :: thresh = 0.1_dp
        real(dp) :: small, large
        integer :: j
        equed = 'N'
        if (m <= 0 .or. n <= 0) return
        small = SFMIN/UNIT_ROUNDOFF
        large = 1/small
        if (rowcnd >= thresh .and. amax >= small .and. amax <= large) then
            if (colcnd < thresh) then
                do j = 1, n
                    a(1:m, j) = c(j)*a(1:m, j)
                end do
                equed = 'C'
            end if
        else if (colcnd >= thresh) then
            do j = 1, n
                a(1:m, j) = r(1:m)*a(1:m, j)
            end do
            equed = 'R'
        else
            do j = 1, n
                a(1:m, j) = c(j)*r(1:m)*a(1:m, j)
            end do
            equed = 'B'
        end if
    end subroutine fortnum_dlaqge

    !> Expert driver for A X = B or A^T X = B by LU (DGESVX): optional
    !> equilibration (fact E) or given factors (fact F), condition estimate,
    !> refinement with error bounds; work(1) returns the reciprocal pivot
    !> growth, info = n + 1 when rcond < eps.
    subroutine fortnum_dgesvx(fact, trans, n, nrhs, a, lda, af, ldaf, ipiv, equed, r, c, b, ldb, x, ldx, &
            rcond, ferr, berr, work, iwork, info)
        character, intent(in) :: fact, trans
        character, intent(inout) :: equed
        integer, intent(in) :: n, nrhs, lda, ldaf, ldb, ldx
        real(dp), intent(inout) :: a(lda, *), af(ldaf, *), r(*), c(*), b(ldb, *)
        integer, intent(inout) :: ipiv(*)
        real(dp), intent(out) :: x(ldx, *), rcond, ferr(*), berr(*)
        real(dp), intent(inout) :: work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: nofact, equil, notran, rowequ, colequ
        real(dp) :: rowcnd, colcnd, amax, anorm, rpvgrw, umax
        integer :: infequ, j
        character :: norm
        info = 0
        nofact = fact == 'N' .or. fact == 'n'
        equil = fact == 'E' .or. fact == 'e'
        notran = trans == 'N' .or. trans == 'n'
        rowequ = .false.
        colequ = .false.
        if (nofact .or. equil) then
            equed = 'N'
        else
            rowequ = index('RrBb', equed) > 0
            colequ = index('CcBb', equed) > 0
        end if
        if (.not. (nofact .or. equil .or. fact == 'F' .or. fact == 'f')) then
            info = -1
        else if (.not. (notran .or. index('TtCc', trans) > 0)) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (nrhs < 0) then
            info = -4
        else if (lda < max(1, n)) then
            info = -6
        else if (ldaf < max(1, n)) then
            info = -8
        else if (.not. (nofact .or. equil) .and. index('NnRrCcBb', equed) == 0) then
            info = -10
        else if (ldb < max(1, n)) then
            info = -14
        else if (ldx < max(1, n)) then
            info = -16
        end if
        if (info /= 0) return
        if (equil) then
            call fortnum_dgeequ(n, n, a, lda, r, c, rowcnd, colcnd, amax, infequ)
            if (infequ == 0) then
                call fortnum_dlaqge(n, n, a, lda, r, c, rowcnd, colcnd, amax, equed)
                rowequ = index('RrBb', equed) > 0
                colequ = index('CcBb', equed) > 0
            end if
        end if
        if (notran) then
            if (rowequ) then
                do j = 1, nrhs
                    b(1:n, j) = r(1:n)*b(1:n, j)
                end do
            end if
        else if (colequ) then
            do j = 1, nrhs
                b(1:n, j) = c(1:n)*b(1:n, j)
            end do
        end if
        if (nofact .or. equil) then
            af(1:n, 1:n) = a(1:n, 1:n)
            call fortnum_dgetrf(n, n, af, ldaf, ipiv, info)
            if (info > 0) then
                work(1) = pivot_growth(info)
                rcond = 0
                return
            end if
        end if
        norm = merge('1', 'I', notran)
        anorm = fortnum_dlange(norm, n, n, a, lda, work)
        rpvgrw = pivot_growth(n)
        call fortnum_dgecon(norm, n, af, ldaf, anorm, rcond, work, iwork, infequ)
        x(1:n, 1:nrhs) = b(1:n, 1:nrhs)
        call fortnum_dgetrs(trans, n, nrhs, af, ldaf, ipiv, x, ldx, infequ)
        call fortnum_dgerfs(trans, n, nrhs, a, lda, af, ldaf, ipiv, b, ldb, x, ldx, ferr, berr, work, iwork, infequ)
        if (notran) then
            if (colequ) then
                do j = 1, nrhs
                    x(1:n, j) = c(1:n)*x(1:n, j)
                end do
                ferr(1:nrhs) = ferr(1:nrhs)/colcnd_of()
            end if
        else if (rowequ) then
            do j = 1, nrhs
                x(1:n, j) = r(1:n)*x(1:n, j)
            end do
            ferr(1:nrhs) = ferr(1:nrhs)/rowcnd_of()
        end if
        work(1) = rpvgrw
        if (rcond < UNIT_ROUNDOFF) info = n + 1
    contains
        ! max |A(:, j)| / max |U(:, j)| over the first ncols columns, minimised
        real(dp) function pivot_growth(ncols)
            integer, intent(in) :: ncols
            integer :: jj
            real(dp) :: amaxj
            pivot_growth = 1
            do jj = 1, ncols
                amaxj = maxval(abs(a(1:n, jj)))
                umax = maxval(abs(af(1:min(jj, n), jj)))
                if (umax /= 0) pivot_growth = min(amaxj/umax, pivot_growth)
            end do
        end function pivot_growth
        real(dp) function colcnd_of()
            colcnd_of = minval(c(1:n))/maxval(c(1:n))
        end function colcnd_of
        real(dp) function rowcnd_of()
            rowcnd_of = minval(r(1:n))/maxval(r(1:n))
        end function rowcnd_of
    end subroutine fortnum_dgesvx

end module fortnum_lapack_refine
