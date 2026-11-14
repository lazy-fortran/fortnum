module fortnum_lapack_ldlt
    ! Symmetric indefinite factorization A = L D L^T or U D U^T with
    ! Bunch-Kaufman pivoting (1 x 1 and 2 x 2 pivots): DSYTF2, the partial
    ! panel DLASYF, DSYTRF and DSYTRS, with the arguments, info codes and
    ! pivot encoding of reference LAPACK (prefix fortnum_): ipiv(k) > 0 is a
    ! 1 x 1 pivot after interchanging k and ipiv(k); a 2 x 2 pivot has
    ! ipiv(k) = ipiv(k+1) = -p (lower; ipiv(k-1) = ipiv(k) = -p upper) after
    ! interchanging k+1 (k-1) and p. Interchanges touch only the part not yet
    ! factored, so the stored L is that of LAPACK's xSYTF2.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm
    use fortnum_blas12, only: fortnum_dgemv
    implicit none
    private

    public :: fortnum_dsytf2, fortnum_dlasyf, fortnum_dsytrf, fortnum_dsytrs

    !> Bunch-Kaufman threshold (1 + sqrt(17))/8.
    real(dp), parameter :: BK_ALPHA = (1.0_dp + sqrt(17.0_dp))/8.0_dp
    !> Panel width reported by the workspace query of DSYTRF.
    integer, parameter :: LDLT_NB = 64

contains

    !> Unblocked Bunch-Kaufman factorization (DSYTF2).
    subroutine fortnum_dsytf2(uplo, n, a, lda, ipiv, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        logical :: upper
        integer :: k, kstep
        call check(uplo, n, lda, upper, info)
        if (info /= 0) return
        if (upper) then
            k = n
            do while (k >= 1)
                call step_upper(n, a, lda, k, ipiv, info, kstep)
                k = k - kstep
            end do
        else
            k = 1
            do while (k <= n)
                call step_lower(n, a, lda, k, ipiv, info, kstep)
                k = k + kstep
            end do
        end if
    end subroutine fortnum_dsytf2

    !> Partial factorization of a panel (DLASYF): factors kb columns, nb - 1
    !> or nb (all n when n <= nb), the first ones for uplo L and the last ones
    !> for uplo U, by LAPACK's blocked algorithm: W = L D is built column by
    !> column with GEMV, the not yet factored block is then updated with GEMM,
    !> and the interchanges on earlier columns are partially undone so the
    !> stored factor has the standard form of DSYTF2. ldw >= n.
    subroutine fortnum_dlasyf(uplo, n, nb, kb, a, lda, ipiv, w, ldw, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, nb, lda, ldw
        integer, intent(out) :: kb
        real(dp), intent(inout) :: a(lda, *), w(ldw, *)
        integer, intent(out) :: ipiv(*)
        integer, intent(out) :: info
        info = 0
        kb = 0
        if (n == 0) return
        if (uplo == 'U' .or. uplo == 'u') then
            call lasyf_upper(n, nb, kb, a, lda, ipiv, w, ldw, info)
        else
            call lasyf_lower(n, nb, kb, a, lda, ipiv, w, ldw, info)
        end if
    end subroutine fortnum_dlasyf

    subroutine lasyf_lower(n, nb, kb, a, lda, ipiv, w, ldw, info)
        integer, intent(in) :: n, nb, lda, ldw
        integer, intent(out) :: kb
        real(dp), intent(inout) :: a(lda, *), w(ldw, *)
        integer, intent(inout) :: ipiv(*), info
        real(dp) :: absakk, colmax, rowmax, t, d11, d22, d21, r1
        integer :: k, kstep, imax, kp, kk, j, jj, jb, jp
        k = 1
        do while (.not. ((k >= nb .and. nb < n) .or. k > n))
            kstep = 1
            ! column k of the trailing matrix, updated: W(k:n, k)
            w(k:n, k) = a(k:n, k)
            if (k > 1) call fortnum_dgemv('N', n - k + 1, k - 1, -1.0_dp, a(k, 1), lda, w(k, 1), ldw, &
                1.0_dp, w(k, k), 1)
            absakk = abs(w(k, k))
            imax = k
            colmax = 0
            if (k < n) then
                imax = k + maxloc(abs(w(k + 1:n, k)), 1)
                colmax = abs(w(imax, k))
            end if
            if (max(absakk, colmax) == 0 .or. absakk /= absakk) then
                if (info == 0) info = k
                kp = k
                a(k:n, k) = w(k:n, k)
            else
                if (absakk >= BK_ALPHA*colmax) then
                    kp = k
                else
                    ! column imax, updated, into W(:, k+1)
                    w(k:imax - 1, k + 1) = a(imax, k:imax - 1)
                    w(imax:n, k + 1) = a(imax:n, imax)
                    if (k > 1) call fortnum_dgemv('N', n - k + 1, k - 1, -1.0_dp, a(k, 1), lda, w(imax, 1), ldw, &
                        1.0_dp, w(k, k + 1), 1)
                    rowmax = 0
                    if (imax > k) rowmax = maxval(abs(w(k:imax - 1, k + 1)))
                    if (imax < n) rowmax = max(rowmax, maxval(abs(w(imax + 1:n, k + 1))))
                    if (absakk >= BK_ALPHA*colmax*(colmax/rowmax)) then
                        kp = k
                    else if (abs(w(imax, k + 1)) >= BK_ALPHA*rowmax) then
                        kp = imax
                        w(k:n, k) = w(k:n, k + 1)
                    else
                        kp = imax
                        kstep = 2
                    end if
                end if
                kk = k + kstep - 1
                if (kp /= kk) then
                    ! the not yet updated column kk into column kp
                    a(kp, kp) = a(kk, kk)
                    do j = kk + 1, kp - 1
                        a(kp, j) = a(j, kk)
                    end do
                    if (kp < n) a(kp + 1:n, kp) = a(kp + 1:n, kk)
                    ! rows kk and kp of the factored columns of A and of W
                    if (k > 1) call swap_rows(a, lda, kk, kp, 1, k - 1)
                    call swap_rows(w, ldw, kk, kp, 1, kk)
                end if
                if (kstep == 1) then
                    a(k:n, k) = w(k:n, k)
                    if (k < n) then
                        r1 = 1/a(k, k)
                        a(k + 1:n, k) = r1*a(k + 1:n, k)
                    end if
                else
                    if (k < n - 1) then
                        d21 = w(k + 1, k)
                        d11 = w(k + 1, k + 1)/d21
                        d22 = w(k, k)/d21
                        t = 1/(d11*d22 - 1)
                        d21 = t/d21
                        do j = k + 2, n
                            a(j, k) = d21*(d11*w(j, k) - w(j, k + 1))
                            a(j, k + 1) = d21*(d22*w(j, k + 1) - w(j, k))
                        end do
                    end if
                    a(k, k) = w(k, k)
                    a(k + 1, k) = w(k + 1, k)
                    a(k + 1, k + 1) = w(k + 1, k + 1)
                end if
            end if
            if (kstep == 1) then
                ipiv(k) = kp
            else
                ipiv(k) = -kp
                ipiv(k + 1) = -kp
            end if
            k = k + kstep
        end do
        ! A22 := A22 - L21 W21^T, the diagonal blocks column by column
        if (k > 1) then
            do j = k, n, nb
                jb = min(nb, n - j + 1)
                do jj = j, j + jb - 1
                    call fortnum_dgemv('N', j + jb - jj, k - 1, -1.0_dp, a(jj, 1), lda, w(jj, 1), ldw, &
                        1.0_dp, a(jj, jj), 1)
                end do
                if (j + jb <= n) call fortnum_dgemm('N', 'T', n - j - jb + 1, jb, k - 1, -1.0_dp, a(j + jb, 1), lda, &
                    w(j, 1), ldw, 1.0_dp, a(j + jb, j), lda)
            end do
        end if
        ! standard form: undo the interchanges in the columns before each pivot
        j = k - 1
        do while (j >= 1)
            jj = j
            jp = ipiv(j)
            if (jp < 0) then
                jp = -jp
                j = j - 1
            end if
            j = j - 1
            if (jp /= jj .and. j >= 1) call swap_rows(a, lda, jp, jj, 1, j)
        end do
        kb = k - 1
    end subroutine lasyf_lower

    subroutine lasyf_upper(n, nb, kb, a, lda, ipiv, w, ldw, info)
        integer, intent(in) :: n, nb, lda, ldw
        integer, intent(out) :: kb
        real(dp), intent(inout) :: a(lda, *), w(ldw, *)
        integer, intent(inout) :: ipiv(*), info
        real(dp) :: absakk, colmax, rowmax, t, d11, d22, d21, r1
        integer :: k, kw, kstep, imax, kp, kk, kkw, j, jj, jb, jp
        k = n
        do while (.not. ((k <= n - nb + 1 .and. nb < n) .or. k < 1))
            kw = nb + k - n
            kstep = 1
            w(1:k, kw) = a(1:k, k)
            if (k < n) call fortnum_dgemv('N', k, n - k, -1.0_dp, a(1, k + 1), lda, w(k, kw + 1), ldw, &
                1.0_dp, w(1, kw), 1)
            absakk = abs(w(k, kw))
            imax = k
            colmax = 0
            if (k > 1) then
                imax = maxloc(abs(w(1:k - 1, kw)), 1)
                colmax = abs(w(imax, kw))
            end if
            if (max(absakk, colmax) == 0 .or. absakk /= absakk) then
                if (info == 0) info = k
                kp = k
                a(1:k, k) = w(1:k, kw)
            else
                if (absakk >= BK_ALPHA*colmax) then
                    kp = k
                else
                    w(1:imax, kw - 1) = a(1:imax, imax)
                    w(imax + 1:k, kw - 1) = a(imax, imax + 1:k)
                    if (k < n) call fortnum_dgemv('N', k, n - k, -1.0_dp, a(1, k + 1), lda, w(imax, kw + 1), ldw, &
                        1.0_dp, w(1, kw - 1), 1)
                    rowmax = 0
                    if (imax < k) rowmax = maxval(abs(w(imax + 1:k, kw - 1)))
                    if (imax > 1) rowmax = max(rowmax, maxval(abs(w(1:imax - 1, kw - 1))))
                    if (absakk >= BK_ALPHA*colmax*(colmax/rowmax)) then
                        kp = k
                    else if (abs(w(imax, kw - 1)) >= BK_ALPHA*rowmax) then
                        kp = imax
                        w(1:k, kw) = w(1:k, kw - 1)
                    else
                        kp = imax
                        kstep = 2
                    end if
                end if
                kk = k - kstep + 1
                kkw = nb + kk - n
                if (kp /= kk) then
                    a(kp, kp) = a(kk, kk)
                    do j = kp + 1, kk - 1
                        a(kp, j) = a(j, kk)
                    end do
                    if (kp > 1) a(1:kp - 1, kp) = a(1:kp - 1, kk)
                    if (k < n) call swap_rows(a, lda, kk, kp, k + 1, n)
                    call swap_rows(w, ldw, kk, kp, kkw, nb)
                end if
                if (kstep == 1) then
                    a(1:k, k) = w(1:k, kw)
                    r1 = 1/a(k, k)
                    a(1:k - 1, k) = r1*a(1:k - 1, k)
                else
                    if (k > 2) then
                        d21 = w(k - 1, kw)
                        d11 = w(k, kw)/d21
                        d22 = w(k - 1, kw - 1)/d21
                        t = 1/(d11*d22 - 1)
                        d21 = t/d21
                        do j = 1, k - 2
                            a(j, k - 1) = d21*(d11*w(j, kw - 1) - w(j, kw))
                            a(j, k) = d21*(d22*w(j, kw) - w(j, kw - 1))
                        end do
                    end if
                    a(k - 1, k - 1) = w(k - 1, kw - 1)
                    a(k - 1, k) = w(k - 1, kw)
                    a(k, k) = w(k, kw)
                end if
            end if
            if (kstep == 1) then
                ipiv(k) = kp
            else
                ipiv(k) = -kp
                ipiv(k - 1) = -kp
            end if
            k = k - kstep
        end do
        ! A11 := A11 - U12 W12^T, the diagonal blocks column by column
        if (k < n) then
            kw = nb + k - n
            do j = ((k - 1)/nb)*nb + 1, 1, -nb
                jb = min(nb, k - j + 1)
                do jj = j, j + jb - 1
                    call fortnum_dgemv('N', jj - j + 1, n - k, -1.0_dp, a(j, k + 1), lda, w(jj, kw + 1), ldw, &
                        1.0_dp, a(j, jj), 1)
                end do
                if (j > 1) call fortnum_dgemm('N', 'T', j - 1, jb, n - k, -1.0_dp, a(1, k + 1), lda, &
                    w(j, kw + 1), ldw, 1.0_dp, a(1, j), lda)
            end do
        end if
        ! standard form: undo the interchanges in the columns after each pivot
        j = k + 1
        do while (j <= n)
            jj = j
            jp = ipiv(j)
            if (jp < 0) then
                jp = -jp
                j = j + 1
            end if
            j = j + 1
            if (jp /= jj .and. j <= n) call swap_rows(a, lda, jp, jj, j, n)
        end do
        kb = n - k
    end subroutine lasyf_upper

    ! swap rows r1 and r2 of x in columns c1..c2
    subroutine swap_rows(x, ldx, r1, r2, c1, c2)
        integer, intent(in) :: ldx, r1, r2, c1, c2
        real(dp), intent(inout) :: x(ldx, *)
        real(dp) :: t
        integer :: c
        do c = c1, c2
            t = x(r1, c)
            x(r1, c) = x(r2, c)
            x(r2, c) = t
        end do
    end subroutine swap_rows

    !> Bunch-Kaufman factorization (DSYTRF): panels of LDLT_NB columns by
    !> DLASYF, the rest by DSYTF2, as LAPACK. lwork = -1 returns the
    !> workspace query (n LDLT_NB); W is allocated internally.
    subroutine fortnum_dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, lwork
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: ipiv(*)
        real(dp), intent(inout) :: work(*)
        integer, intent(out) :: info
        logical :: upper
        real(dp), allocatable :: w(:, :)
        integer :: k, kb, iinfo, j
        call check(uplo, n, lda, upper, info)
        if (info == 0 .and. lwork < 1 .and. lwork /= -1) info = -7
        if (info /= 0) return
        work(1) = max(1, n*LDLT_NB)
        if (lwork == -1 .or. n == 0) return
        allocate (w(n, LDLT_NB))
        if (upper) then
            k = n
            do while (k >= 1)
                if (k > LDLT_NB) then
                    call fortnum_dlasyf('U', k, LDLT_NB, kb, a, lda, ipiv, w, n, iinfo)
                else
                    call fortnum_dsytf2('U', k, a, lda, ipiv, iinfo)
                    kb = k
                end if
                if (info == 0 .and. iinfo > 0) info = iinfo
                k = k - kb
            end do
        else
            k = 1
            do while (k <= n)
                if (k <= n - LDLT_NB) then
                    call fortnum_dlasyf('L', n - k + 1, LDLT_NB, kb, a(k, k), lda, ipiv(k), w, n, iinfo)
                else
                    call fortnum_dsytf2('L', n - k + 1, a(k, k), lda, ipiv(k), iinfo)
                    kb = n - k + 1
                end if
                if (info == 0 .and. iinfo > 0) info = iinfo + k - 1
                do j = k, k + kb - 1
                    if (ipiv(j) > 0) then
                        ipiv(j) = ipiv(j) + k - 1
                    else
                        ipiv(j) = ipiv(j) - k + 1
                    end if
                end do
                k = k + kb
            end do
        end if
    end subroutine fortnum_dsytrf

    !> Solve A X = B with the factors of DSYTRF (DSYTRS).
    subroutine fortnum_dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, nrhs, lda, ldb
        real(dp), intent(in) :: a(lda, *)
        integer, intent(in) :: ipiv(*)
        real(dp), intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        logical :: upper
        integer :: k, kp, j
        call check(uplo, n, lda, upper, info)
        if (info == 0 .and. nrhs < 0) info = -3
        if (info == 0 .and. ldb < max(1, n)) info = -8
        if (info /= 0 .or. n == 0 .or. nrhs == 0) return
        if (upper) then
            ! U D Y = B, from the last column back
            k = n
            do while (k >= 1)
                if (ipiv(k) > 0) then
                    call swap_rows(k, ipiv(k))
                    do j = 1, nrhs
                        b(1:k - 1, j) = b(1:k - 1, j) - b(k, j)*a(1:k - 1, k)
                    end do
                    b(k, 1:nrhs) = b(k, 1:nrhs)/a(k, k)
                    k = k - 1
                else
                    kp = -ipiv(k)
                    call swap_rows(k - 1, kp)
                    do j = 1, nrhs
                        b(1:k - 2, j) = b(1:k - 2, j) - b(k, j)*a(1:k - 2, k) - b(k - 1, j)*a(1:k - 2, k - 1)
                    end do
                    call solve2(k - 1, k, a(k - 1, k), a(k - 1, k - 1), a(k, k))
                    k = k - 2
                end if
            end do
            ! U^T X = Y, from the first column on
            k = 1
            do while (k <= n)
                if (ipiv(k) > 0) then
                    do j = 1, nrhs
                        b(k, j) = b(k, j) - dot_product(a(1:k - 1, k), b(1:k - 1, j))
                    end do
                    call swap_rows(k, ipiv(k))
                    k = k + 1
                else
                    do j = 1, nrhs
                        b(k, j) = b(k, j) - dot_product(a(1:k - 1, k), b(1:k - 1, j))
                        b(k + 1, j) = b(k + 1, j) - dot_product(a(1:k - 1, k + 1), b(1:k - 1, j))
                    end do
                    call swap_rows(k, -ipiv(k))
                    k = k + 2
                end if
            end do
        else
            ! L D Y = B, from the first column on
            k = 1
            do while (k <= n)
                if (ipiv(k) > 0) then
                    call swap_rows(k, ipiv(k))
                    do j = 1, nrhs
                        b(k + 1:n, j) = b(k + 1:n, j) - b(k, j)*a(k + 1:n, k)
                    end do
                    b(k, 1:nrhs) = b(k, 1:nrhs)/a(k, k)
                    k = k + 1
                else
                    kp = -ipiv(k)
                    call swap_rows(k + 1, kp)
                    do j = 1, nrhs
                        b(k + 2:n, j) = b(k + 2:n, j) - b(k, j)*a(k + 2:n, k) - b(k + 1, j)*a(k + 2:n, k + 1)
                    end do
                    call solve2(k, k + 1, a(k + 1, k), a(k, k), a(k + 1, k + 1))
                    k = k + 2
                end if
            end do
            ! L^T X = Y, from the last column back
            k = n
            do while (k >= 1)
                if (ipiv(k) > 0) then
                    do j = 1, nrhs
                        b(k, j) = b(k, j) - dot_product(a(k + 1:n, k), b(k + 1:n, j))
                    end do
                    call swap_rows(k, ipiv(k))
                    k = k - 1
                else
                    do j = 1, nrhs
                        b(k, j) = b(k, j) - dot_product(a(k + 1:n, k), b(k + 1:n, j))
                        b(k - 1, j) = b(k - 1, j) - dot_product(a(k + 1:n, k - 1), b(k + 1:n, j))
                    end do
                    call swap_rows(k, -ipiv(k))
                    k = k - 2
                end if
            end do
        end if

    contains

        subroutine swap_rows(i, p)
            integer, intent(in) :: i, p
            real(dp) :: t
            integer :: jj
            if (i == p) return
            do jj = 1, nrhs
                t = b(i, jj)
                b(i, jj) = b(p, jj)
                b(p, jj) = t
            end do
        end subroutine swap_rows

        ! the 2 x 2 diagonal block [[d1, e], [e, d2]] applied inversely to rows i1, i2
        subroutine solve2(i1, i2, e, d1, d2)
            integer, intent(in) :: i1, i2
            real(dp), intent(in) :: e, d1, d2
            real(dp) :: akm1, ak, denom, bkm1, bk
            integer :: jj
            akm1 = d1/e
            ak = d2/e
            denom = akm1*ak - 1
            do jj = 1, nrhs
                bkm1 = b(i1, jj)/e
                bk = b(i2, jj)/e
                b(i1, jj) = (ak*bkm1 - bk)/denom
                b(i2, jj) = (akm1*bk - bkm1)/denom
            end do
        end subroutine solve2

    end subroutine fortnum_dsytrs

    ! One Bunch-Kaufman step at column k of the lower triangle: choose the
    ! pivot, interchange within A(k:n, k:n), and update the trailing block.
    subroutine step_lower(n, a, lda, k, ipiv, info, kstep)
        integer, intent(in) :: n, lda, k
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(inout) :: ipiv(*), info
        integer, intent(out) :: kstep
        real(dp) :: absakk, colmax, rowmax, t, d11, d22, d21, wk, wkp1, r1
        integer :: imax, kp, kk, i, j
        kstep = 1
        absakk = abs(a(k, k))
        imax = k
        colmax = 0
        if (k < n) then
            imax = k + maxloc(abs(a(k + 1:n, k)), 1)
            colmax = abs(a(imax, k))
        end if
        if (max(absakk, colmax) == 0 .or. absakk /= absakk) then
            if (info == 0) info = k
            kp = k
        else
            if (absakk >= BK_ALPHA*colmax) then
                kp = k
            else
                rowmax = maxval(abs(a(imax, k:imax - 1)))
                if (imax < n) rowmax = max(rowmax, maxval(abs(a(imax + 1:n, imax))))
                if (absakk >= BK_ALPHA*colmax*(colmax/rowmax)) then
                    kp = k
                else if (abs(a(imax, imax)) >= BK_ALPHA*rowmax) then
                    kp = imax
                else
                    kp = imax
                    kstep = 2
                end if
            end if
            kk = k + kstep - 1
            if (kp /= kk) then
                do i = kp + 1, n
                    t = a(i, kk)
                    a(i, kk) = a(i, kp)
                    a(i, kp) = t
                end do
                do j = kk + 1, kp - 1
                    t = a(j, kk)
                    a(j, kk) = a(kp, j)
                    a(kp, j) = t
                end do
                t = a(kk, kk)
                a(kk, kk) = a(kp, kp)
                a(kp, kp) = t
                if (kstep == 2) then
                    t = a(k + 1, k)
                    a(k + 1, k) = a(kp, k)
                    a(kp, k) = t
                end if
            end if
            if (kstep == 1) then
                if (k < n) then
                    r1 = 1/a(k, k)
                    do j = k + 1, n
                        a(j:n, j) = a(j:n, j) - r1*a(j, k)*a(j:n, k)
                    end do
                    a(k + 1:n, k) = r1*a(k + 1:n, k)
                end if
            else if (k < n - 1) then
                d21 = a(k + 1, k)
                d11 = a(k + 1, k + 1)/d21
                d22 = a(k, k)/d21
                t = 1/(d11*d22 - 1)
                d21 = t/d21
                do j = k + 2, n
                    wk = d21*(d11*a(j, k) - a(j, k + 1))
                    wkp1 = d21*(d22*a(j, k + 1) - a(j, k))
                    a(j:n, j) = a(j:n, j) - a(j:n, k)*wk - a(j:n, k + 1)*wkp1
                    a(j, k) = wk
                    a(j, k + 1) = wkp1
                end do
            end if
        end if
        if (kstep == 1) then
            ipiv(k) = kp
        else
            ipiv(k) = -kp
            ipiv(k + 1) = -kp
        end if
    end subroutine step_lower

    ! One Bunch-Kaufman step at column k of the upper triangle.
    subroutine step_upper(n, a, lda, k, ipiv, info, kstep)
        integer, intent(in) :: n, lda, k
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(inout) :: ipiv(*), info
        integer, intent(out) :: kstep
        real(dp) :: absakk, colmax, rowmax, t, d11, d22, d12, wk, wkm1, r1
        integer :: imax, kp, kk, i, j
        kstep = 1
        absakk = abs(a(k, k))
        imax = k
        colmax = 0
        if (k > 1) then
            imax = maxloc(abs(a(1:k - 1, k)), 1)
            colmax = abs(a(imax, k))
        end if
        if (max(absakk, colmax) == 0 .or. absakk /= absakk) then
            if (info == 0) info = k
            kp = k
        else
            if (absakk >= BK_ALPHA*colmax) then
                kp = k
            else
                rowmax = maxval(abs(a(imax, imax + 1:k)))
                if (imax > 1) rowmax = max(rowmax, maxval(abs(a(1:imax - 1, imax))))
                if (absakk >= BK_ALPHA*colmax*(colmax/rowmax)) then
                    kp = k
                else if (abs(a(imax, imax)) >= BK_ALPHA*rowmax) then
                    kp = imax
                else
                    kp = imax
                    kstep = 2
                end if
            end if
            kk = k - kstep + 1
            if (kp /= kk) then
                do i = 1, kp - 1
                    t = a(i, kk)
                    a(i, kk) = a(i, kp)
                    a(i, kp) = t
                end do
                do j = kp + 1, kk - 1
                    t = a(j, kk)
                    a(j, kk) = a(kp, j)
                    a(kp, j) = t
                end do
                t = a(kk, kk)
                a(kk, kk) = a(kp, kp)
                a(kp, kp) = t
                if (kstep == 2) then
                    t = a(k - 1, k)
                    a(k - 1, k) = a(kp, k)
                    a(kp, k) = t
                end if
            end if
            if (kstep == 1) then
                if (k > 1) then
                    r1 = 1/a(k, k)
                    do j = 1, k - 1
                        a(1:j, j) = a(1:j, j) - r1*a(j, k)*a(1:j, k)
                    end do
                    a(1:k - 1, k) = r1*a(1:k - 1, k)
                end if
            else if (k > 2) then
                d12 = a(k - 1, k)
                d22 = a(k - 1, k - 1)/d12
                d11 = a(k, k)/d12
                t = 1/(d11*d22 - 1)
                d12 = t/d12
                do j = k - 2, 1, -1
                    wkm1 = d12*(d11*a(j, k - 1) - a(j, k))
                    wk = d12*(d22*a(j, k) - a(j, k - 1))
                    a(1:j, j) = a(1:j, j) - a(1:j, k)*wk - a(1:j, k - 1)*wkm1
                    a(j, k) = wk
                    a(j, k - 1) = wkm1
                end do
            end if
        end if
        if (kstep == 1) then
            ipiv(k) = kp
        else
            ipiv(k) = -kp
            ipiv(k - 1) = -kp
        end if
    end subroutine step_upper

    subroutine check(uplo, n, lda, upper, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        logical, intent(out) :: upper
        integer, intent(out) :: info
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -4
        end if
    end subroutine check

end module fortnum_lapack_ldlt
