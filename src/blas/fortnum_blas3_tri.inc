! Triangular solve and multiply, real or complex through the macros SCALAR,
! GEMM, TRSM, TRMM, CONJ, ONE and ZERO (fortnum_blas3.F90). Blocked: the
! diagonal blocks are solved or multiplied directly, everything else is GEMM
! on the remaining blocks, so the cost is that of GEMM for large problems.

    !> B := alpha inv(op(A)) B (side L) or alpha B inv(op(A)) (side R), A
    !> triangular, op(A) = A, A^T or A^H; reference BLAS xTRSM.
    subroutine TRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        character, intent(in) :: side, uplo, transa, diag
        integer, intent(in) :: m, n, lda, ldb
        SCALAR, intent(in) :: alpha
        SCALAR, intent(in) :: a(lda, *)
        SCALAR, intent(inout) :: b(ldb, *)
        real(dp), allocatable :: work(:)
        integer :: lwork, nth, nfree, chunk, c0, cn, j
        logical :: left, upper, notrans, conj, unit
        character :: tr
        left = side == 'L' .or. side == 'l'
        ! small problems: substitution in place, without packing, workspace
        ! or threads, whose set-up would cost more than the solve
        if (real(merge(m, n, left), dp)**2*merge(n, m, left) <= TRSM_SMALL_WORK) then
            call tri_flags('TRSM', side, uplo, transa, diag, m, n, lda, ldb, left, upper, &
                notrans, conj, unit, tr)
            if (m == 0 .or. n == 0) return
            if (alpha == ZERO) then
                do j = 1, n
                    b(1:m, j) = ZERO
                end do
                return
            end if
            if (alpha /= ONE) then
                do j = 1, n
                    b(1:m, j) = alpha*b(1:m, j)
                end do
            end if
            call small_solve()
            return
        end if
        ! the columns (left) or rows (right) of B are independent: on threads,
        ! each solves its own chunk with serial GEMMs, without the barriers of
        ! threaded GEMMs on thin blocks
        nfree = merge(n, m, left)
        nth = thread_count(real(SCALAR_WORDS, dp)*nfree*real(merge(m, n, left), dp)**2)
        nth = min(nth, nfree/TRSM_CHUNK)
        if (nth > 1) then
            chunk = (nfree + nth - 1)/nth
            !$omp parallel do num_threads(nth) schedule(static) private(c0, cn, lwork, work)
            do c0 = 1, nfree, chunk
                cn = min(chunk, nfree - c0 + 1)
                if (left) then
                    lwork = TRSM_LWORK(side, m, cn, 1)
                    allocate (work(lwork))
                    call TRSM_WS(side, uplo, transa, diag, m, cn, alpha, a, lda, b(1, c0), ldb, work, lwork)
                else
                    lwork = TRSM_LWORK(side, cn, n, 1)
                    allocate (work(lwork))
                    call TRSM_WS(side, uplo, transa, diag, cn, n, alpha, a, lda, b(c0, 1), ldb, work, lwork)
                end if
                deallocate (work)
            end do
            !$omp end parallel do
            return
        end if
        lwork = 0
        if (m > 0 .and. n > 0) lwork = TRSM_LWORK(side, m, n)
        allocate (work(lwork))
        call TRSM_WS(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, work, lwork)

    contains

        ! The column-oriented substitutions of reference BLAS xTRSM, alpha
        ! already applied.
        subroutine small_solve()
            integer :: i, k
            SCALAR :: t, akk
            if (left .and. notrans) then
                do j = 1, n
                    if (upper) then
                        do k = m, 1, -1
                            if (b(k, j) == ZERO) cycle
                            if (.not. unit) b(k, j) = b(k, j)/a(k, k)
                            b(1:k - 1, j) = b(1:k - 1, j) - b(k, j)*a(1:k - 1, k)
                        end do
                    else
                        do k = 1, m
                            if (b(k, j) == ZERO) cycle
                            if (.not. unit) b(k, j) = b(k, j)/a(k, k)
                            b(k + 1:m, j) = b(k + 1:m, j) - b(k, j)*a(k + 1:m, k)
                        end do
                    end if
                end do
            else if (left) then
                do j = 1, n
                    if (upper) then
                        do i = 1, m
                            t = b(i, j)
                            do k = 1, i - 1
                                t = t - op(a(k, i))*b(k, j)
                            end do
                            if (.not. unit) t = t/op(a(i, i))
                            b(i, j) = t
                        end do
                    else
                        do i = m, 1, -1
                            t = b(i, j)
                            do k = i + 1, m
                                t = t - op(a(k, i))*b(k, j)
                            end do
                            if (.not. unit) t = t/op(a(i, i))
                            b(i, j) = t
                        end do
                    end if
                end do
            else if (notrans) then
                if (upper) then
                    do j = 1, n
                        do k = 1, j - 1
                            if (a(k, j) /= ZERO) b(1:m, j) = b(1:m, j) - a(k, j)*b(1:m, k)
                        end do
                        if (.not. unit) b(1:m, j) = b(1:m, j)/a(j, j)
                    end do
                else
                    do j = n, 1, -1
                        do k = j + 1, n
                            if (a(k, j) /= ZERO) b(1:m, j) = b(1:m, j) - a(k, j)*b(1:m, k)
                        end do
                        if (.not. unit) b(1:m, j) = b(1:m, j)/a(j, j)
                    end do
                end if
            else
                if (upper) then
                    do k = n, 1, -1
                        if (.not. unit) then
                            akk = op(a(k, k))
                            b(1:m, k) = b(1:m, k)/akk
                        end if
                        do j = 1, k - 1
                            if (a(j, k) /= ZERO) b(1:m, j) = b(1:m, j) - op(a(j, k))*b(1:m, k)
                        end do
                    end do
                else
                    do k = 1, n
                        if (.not. unit) then
                            akk = op(a(k, k))
                            b(1:m, k) = b(1:m, k)/akk
                        end if
                        do j = k + 1, n
                            if (a(j, k) /= ZERO) b(1:m, j) = b(1:m, j) - op(a(j, k))*b(1:m, k)
                        end do
                    end do
                end if
            end if
        end subroutine small_solve

        elemental SCALAR function op(x)
            SCALAR, intent(in) :: x
            if (conj) then
                op = CONJ(x)
            else
                op = x
            end if
        end function op
    end subroutine TRSM

    !> Workspace in real(dp) elements that the workspace variant of TRSM
    !> needs for B of m x n on up to `threads` OpenMP threads (default as
    !> GEMM's). Nondecreasing in m, n and threads.
    integer function TRSM_LWORK(side, m, n, threads) result(lwork)
        character, intent(in) :: side
        integer, intent(in) :: m, n
        integer, intent(in), optional :: threads
        if (side == 'L' .or. side == 'l') then
            lwork = max(GEMM_LWORK(m, n, m, threads), SCALAR_WORDS*max(0, n)*TRSM_LEAF)
        else
            lwork = GEMM_LWORK(m, n, n, threads)
        end if
#ifdef TRSM_PACKED
        if (side == 'L' .or. side == 'l') then
            lwork = max(lwork, TRSM_PACKED_LWORK(m, n))
        else
            lwork = max(lwork, TRSM_PACKED_LWORK(n, m))
        end if
#endif
    end function TRSM_LWORK

    !> TRSM with caller-provided workspace work(lwork), lwork at least
    !> TRSM_LWORK(side, m, n); no memory is allocated.
    subroutine TRSM_WS(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb, work, lwork)
        character, intent(in) :: side, uplo, transa, diag
        integer, intent(in) :: m, n, lda, ldb, lwork
        SCALAR, intent(in) :: alpha
        SCALAR, intent(in) :: a(lda, *)
        SCALAR, intent(inout) :: b(ldb, *)
        real(dp), intent(inout), target :: work(lwork)
        logical :: left, upper, notrans, conj, unit
        character :: tr
        integer :: j

        call tri_flags('TRSM', side, uplo, transa, diag, m, n, lda, ldb, left, upper, notrans, &
            conj, unit, tr)
        if (m == 0 .or. n == 0) return
        if (alpha == ZERO) then
            do j = 1, n
                b(1:m, j) = ZERO
            end do
            return
        end if
        if (alpha /= ONE) then
            do j = 1, n
                b(1:m, j) = alpha*b(1:m, j)
            end do
        end if
#ifdef TRSM_PACKED
        ! real data: packed panels on the GEMM register tile, no recursion
        call TRSM_PACKED(left, upper, notrans, unit, m, n, a, lda, b, ldb, work, lwork)
        return
#endif
        if (left) then
            call solve(1, m)
        else
            call solve(1, n)
        end if

    contains

        ! The triangle op(A)(i0:i0+ni-1, same) against rows (left) or columns
        ! (right) i0..i0+ni-1 of B, by halves: one half solved, the other
        ! updated by GEMM, then solved. No inverse is formed and no flop is
        ! spent on the zero triangle; inverting 64 x 64 diagonal blocks and
        ! multiplying by them in full cost DPOTRF a third of its time at
        ! order 200.
        recursive subroutine solve(i0, ni)
            integer, intent(in) :: i0, ni
            integer :: n1, n2, i1
            if (ni <= TRSM_LEAF) then
                call substitute(i0, ni)
                return
            end if
            n1 = ni/2
            n2 = ni - n1
            i1 = i0 + n1
            if (left) then
                if (upper .neqv. notrans) then
                    ! op(A) lower: forward
                    call solve(i0, n1)
                    if (notrans) then
                        call GEMM_WS('N', 'N', n2, n, n1, -ONE, a(i1, i0), lda, b(i0, 1), ldb, &
                            ONE, b(i1, 1), ldb, work, lwork)
                    else
                        call GEMM_WS(tr, 'N', n2, n, n1, -ONE, a(i0, i1), lda, b(i0, 1), ldb, &
                            ONE, b(i1, 1), ldb, work, lwork)
                    end if
                    call solve(i1, n2)
                else
                    ! op(A) upper: backward
                    call solve(i1, n2)
                    if (notrans) then
                        call GEMM_WS('N', 'N', n1, n, n2, -ONE, a(i0, i1), lda, b(i1, 1), ldb, &
                            ONE, b(i0, 1), ldb, work, lwork)
                    else
                        call GEMM_WS(tr, 'N', n1, n, n2, -ONE, a(i1, i0), lda, b(i1, 1), ldb, &
                            ONE, b(i0, 1), ldb, work, lwork)
                    end if
                    call solve(i0, n1)
                end if
            else
                if (upper .eqv. notrans) then
                    ! op(A) upper: forward over columns
                    call solve(i0, n1)
                    if (notrans) then
                        call GEMM_WS('N', 'N', m, n2, n1, -ONE, b(1, i0), ldb, a(i0, i1), lda, &
                            ONE, b(1, i1), ldb, work, lwork)
                    else
                        call GEMM_WS('N', tr, m, n2, n1, -ONE, b(1, i0), ldb, a(i1, i0), lda, &
                            ONE, b(1, i1), ldb, work, lwork)
                    end if
                    call solve(i1, n2)
                else
                    ! op(A) lower: backward over columns
                    call solve(i1, n2)
                    if (notrans) then
                        call GEMM_WS('N', 'N', m, n1, n2, -ONE, b(1, i1), ldb, a(i1, i0), lda, &
                            ONE, b(1, i0), ldb, work, lwork)
                    else
                        call GEMM_WS('N', tr, m, n1, n2, -ONE, b(1, i1), ldb, a(i0, i1), lda, &
                            ONE, b(1, i0), ldb, work, lwork)
                    end if
                    call solve(i0, n1)
                end if
            end if
        end subroutine solve

        ! substitution with at most TRSM_LEAF rows or columns: op(A) copied
        ! once, each column of B (left) or the columns of B in order (right)
        subroutine substitute(i0, ni)
            integer, intent(in) :: i0, ni
            SCALAR :: t(TRSM_LEAF, TRSM_LEAF)
            SCALAR, pointer, contiguous :: w(:, :)
            integer :: r, c, r0, r1, nth
            do c = 1, ni
                do r = 1, ni
                    if (notrans) then
                        t(r, c) = a(i0 + r - 1, i0 + c - 1)
                    else if (conj) then
                        t(r, c) = CONJ(a(i0 + c - 1, i0 + r - 1))
                    else
                        t(r, c) = a(i0 + c - 1, i0 + r - 1)
                    end if
                end do
            end do
            ! the rows (right) or columns (left) of B are independent: chunks of
            ! them go to threads. Vector updates run at about a tenth of the
            ! GEMM rate that thread_count assumes, hence the weight.
            nth = thread_count(10*real(SCALAR_WORDS, dp)*merge(n, m, left)*real(ni, dp)*ni)
            if (left) then
                ! the rows of B transposed into w, so that each step is a
                ! vector operation along the columns of B
                if (lwork < SCALAR_WORDS*n*ni) call argument_error('TRSM', 13)
                call c_f_pointer(c_loc(work), w, [n, ni])
                w = transpose(b(i0:i0 + ni - 1, 1:n))
                !$omp parallel do num_threads(nth) if(nth > 1) schedule(static) private(r0, r1, c, r)
                do r0 = 1, n, SUBST_CHUNK
                    r1 = min(n, r0 + SUBST_CHUNK - 1)
                    if (upper .neqv. notrans) then
                        do c = 1, ni
                            if (.not. unit) w(r0:r1, c) = w(r0:r1, c)/t(c, c)
                            do r = c + 1, ni
                                w(r0:r1, r) = w(r0:r1, r) - w(r0:r1, c)*t(r, c)
                            end do
                        end do
                    else
                        do c = ni, 1, -1
                            if (.not. unit) w(r0:r1, c) = w(r0:r1, c)/t(c, c)
                            do r = 1, c - 1
                                w(r0:r1, r) = w(r0:r1, r) - w(r0:r1, c)*t(r, c)
                            end do
                        end do
                    end if
                end do
                !$omp end parallel do
                b(i0:i0 + ni - 1, 1:n) = transpose(w)
            else
                !$omp parallel do num_threads(nth) if(nth > 1) schedule(static) private(r0, r1, c, r)
                do r0 = 1, m, SUBST_CHUNK
                    r1 = min(m, r0 + SUBST_CHUNK - 1)
                    if (upper .eqv. notrans) then
                        do c = 1, ni
                            if (.not. unit) b(r0:r1, i0 + c - 1) = b(r0:r1, i0 + c - 1)/t(c, c)
                            do r = c + 1, ni
                                b(r0:r1, i0 + r - 1) = b(r0:r1, i0 + r - 1) - b(r0:r1, i0 + c - 1)*t(c, r)
                            end do
                        end do
                    else
                        do c = ni, 1, -1
                            if (.not. unit) b(r0:r1, i0 + c - 1) = b(r0:r1, i0 + c - 1)/t(c, c)
                            do r = 1, c - 1
                                b(r0:r1, i0 + r - 1) = b(r0:r1, i0 + r - 1) - b(r0:r1, i0 + c - 1)*t(c, r)
                            end do
                        end do
                    end if
                end do
                !$omp end parallel do
            end if
        end subroutine substitute

    end subroutine TRSM_WS

    !> B := alpha op(A) B (side L) or alpha B op(A) (side R), A triangular;
    !> reference BLAS xTRMM.
    subroutine TRMM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
        character, intent(in) :: side, uplo, transa, diag
        integer, intent(in) :: m, n, lda, ldb
        SCALAR, intent(in) :: alpha
        SCALAR, intent(in) :: a(lda, *)
        SCALAR, intent(inout) :: b(ldb, *)
        logical :: left, upper, notrans, conj, unit
        character :: tr
        integer :: i0, ib, j0, jb, rest, j
        SCALAR, allocatable :: t(:, :), w(:, :)

        call tri_flags('TRMM', side, uplo, transa, diag, m, n, lda, ldb, left, upper, notrans, &
            conj, unit, tr)
        if (m == 0 .or. n == 0) return
        if (alpha == ZERO) then
            do j = 1, n
                b(1:m, j) = ZERO
            end do
            return
        end if
        allocate (t(TRI_NB, TRI_NB))
        if (left) then
            allocate (w(TRI_NB, n))
            if (upper .eqv. notrans) then
                ! op(A) upper: top-down, rows below are still the old B
                do i0 = 1, m, TRI_NB
                    ib = min(TRI_NB, m - i0 + 1)
                    call diag_left(i0, ib, .true.)
                    rest = m - i0 - ib + 1
                    if (rest > 0) then
                        if (notrans) then
                            call GEMM('N', 'N', ib, n, rest, ONE, a(i0, i0 + ib), lda, b(i0 + ib, 1), ldb, &
                                ONE, b(i0, 1), ldb)
                        else
                            call GEMM(tr, 'N', ib, n, rest, ONE, a(i0 + ib, i0), lda, b(i0 + ib, 1), ldb, &
                                ONE, b(i0, 1), ldb)
                        end if
                    end if
                end do
            else
                ! op(A) lower: bottom-up, rows above are still the old B
                do i0 = ((m - 1)/TRI_NB)*TRI_NB + 1, 1, -TRI_NB
                    ib = min(TRI_NB, m - i0 + 1)
                    call diag_left(i0, ib, .false.)
                    if (i0 > 1) then
                        if (notrans) then
                            call GEMM('N', 'N', ib, n, i0 - 1, ONE, a(i0, 1), lda, b, ldb, &
                                ONE, b(i0, 1), ldb)
                        else
                            call GEMM(tr, 'N', ib, n, i0 - 1, ONE, a(1, i0), lda, b, ldb, &
                                ONE, b(i0, 1), ldb)
                        end if
                    end if
                end do
            end if
        else
            allocate (w(m, TRI_NB))
            if (upper .eqv. notrans) then
                ! op(A) upper: right to left, columns to the left are still the old B
                do j0 = ((n - 1)/TRI_NB)*TRI_NB + 1, 1, -TRI_NB
                    jb = min(TRI_NB, n - j0 + 1)
                    call diag_right(j0, jb, .true.)
                    if (j0 > 1) then
                        if (notrans) then
                            call GEMM('N', 'N', m, jb, j0 - 1, ONE, b, ldb, a(1, j0), lda, &
                                ONE, b(1, j0), ldb)
                        else
                            call GEMM('N', tr, m, jb, j0 - 1, ONE, b, ldb, a(j0, 1), lda, &
                                ONE, b(1, j0), ldb)
                        end if
                    end if
                end do
            else
                ! op(A) lower: left to right, columns to the right are still the old B
                do j0 = 1, n, TRI_NB
                    jb = min(TRI_NB, n - j0 + 1)
                    call diag_right(j0, jb, .false.)
                    rest = n - j0 - jb + 1
                    if (rest > 0) then
                        if (notrans) then
                            call GEMM('N', 'N', m, jb, rest, ONE, b(1, j0 + jb), ldb, a(j0 + jb, j0), lda, &
                                ONE, b(1, j0), ldb)
                        else
                            call GEMM('N', tr, m, jb, rest, ONE, b(1, j0 + jb), ldb, a(j0, j0 + jb), lda, &
                                ONE, b(1, j0), ldb)
                        end if
                    end if
                end do
            end if
        end if
        if (alpha /= ONE) then
            do j = 1, n
                b(1:m, j) = alpha*b(1:m, j)
            end do
        end if

    contains

        ! t(1:nb, 1:nb) := the triangle op(A)(k0:k0+nb-1, same), zeros elsewhere.
        subroutine diag_block(k0, nb, upper_op)
            integer, intent(in) :: k0, nb
            logical, intent(in) :: upper_op
            integer :: r, c
            do c = 1, nb
                do r = 1, nb
                    if ((upper_op .and. r > c) .or. (.not. upper_op .and. r < c)) then
                        t(r, c) = ZERO
                    else if (r == c .and. unit) then
                        t(r, c) = ONE
                    else if (notrans) then
                        t(r, c) = a(k0 + r - 1, k0 + c - 1)
                    else if (conj) then
                        t(r, c) = CONJ(a(k0 + c - 1, k0 + r - 1))
                    else
                        t(r, c) = a(k0 + c - 1, k0 + r - 1)
                    end if
                end do
            end do
        end subroutine diag_block

        subroutine diag_left(i0, ib, upper_op)
            integer, intent(in) :: i0, ib
            logical, intent(in) :: upper_op
            call diag_block(i0, ib, upper_op)
            w(1:ib, 1:n) = b(i0:i0 + ib - 1, 1:n)
            call GEMM('N', 'N', ib, n, ib, ONE, t, TRI_NB, w, TRI_NB, ZERO, b(i0, 1), ldb)
        end subroutine diag_left

        subroutine diag_right(j0, jb, upper_op)
            integer, intent(in) :: j0, jb
            logical, intent(in) :: upper_op
            call diag_block(j0, jb, upper_op)
            w(1:m, 1:jb) = b(1:m, j0:j0 + jb - 1)
            call GEMM('N', 'N', m, jb, jb, ONE, w, m, t, TRI_NB, ZERO, b(1, j0), ldb)
        end subroutine diag_right

    end subroutine TRMM
