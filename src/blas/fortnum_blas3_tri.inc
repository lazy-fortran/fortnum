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
        logical :: left, upper, notrans, conj, unit
        character :: tr
        integer :: i0, ib, j0, jb, rest, j

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
        if (left) then
            if (upper .neqv. notrans) then
                ! op(A) lower: forward over row blocks
                do i0 = 1, m, TRI_NB
                    ib = min(TRI_NB, m - i0 + 1)
                    call solve_left(i0, ib, .true.)
                    rest = m - i0 - ib + 1
                    if (rest > 0) then
                        if (notrans) then
                            call GEMM('N', 'N', rest, n, ib, -ONE, a(i0 + ib, i0), lda, b(i0, 1), ldb, &
                                ONE, b(i0 + ib, 1), ldb)
                        else
                            call GEMM(tr, 'N', rest, n, ib, -ONE, a(i0, i0 + ib), lda, b(i0, 1), ldb, &
                                ONE, b(i0 + ib, 1), ldb)
                        end if
                    end if
                end do
            else
                ! op(A) upper: backward over row blocks
                do i0 = ((m - 1)/TRI_NB)*TRI_NB + 1, 1, -TRI_NB
                    ib = min(TRI_NB, m - i0 + 1)
                    call solve_left(i0, ib, .false.)
                    if (i0 > 1) then
                        if (notrans) then
                            call GEMM('N', 'N', i0 - 1, n, ib, -ONE, a(1, i0), lda, b(i0, 1), ldb, &
                                ONE, b, ldb)
                        else
                            call GEMM(tr, 'N', i0 - 1, n, ib, -ONE, a(i0, 1), lda, b(i0, 1), ldb, &
                                ONE, b, ldb)
                        end if
                    end if
                end do
            end if
        else
            if (upper .eqv. notrans) then
                ! op(A) upper: forward over column blocks
                do j0 = 1, n, TRI_NB
                    jb = min(TRI_NB, n - j0 + 1)
                    call solve_right(j0, jb, .true.)
                    rest = n - j0 - jb + 1
                    if (rest > 0) then
                        if (notrans) then
                            call GEMM('N', 'N', m, rest, jb, -ONE, b(1, j0), ldb, a(j0, j0 + jb), lda, &
                                ONE, b(1, j0 + jb), ldb)
                        else
                            call GEMM('N', tr, m, rest, jb, -ONE, b(1, j0), ldb, a(j0 + jb, j0), lda, &
                                ONE, b(1, j0 + jb), ldb)
                        end if
                    end if
                end do
            else
                ! op(A) lower: backward over column blocks
                do j0 = ((n - 1)/TRI_NB)*TRI_NB + 1, 1, -TRI_NB
                    jb = min(TRI_NB, n - j0 + 1)
                    call solve_right(j0, jb, .false.)
                    if (j0 > 1) then
                        if (notrans) then
                            call GEMM('N', 'N', m, j0 - 1, jb, -ONE, b(1, j0), ldb, a(j0, 1), lda, &
                                ONE, b, ldb)
                        else
                            call GEMM('N', tr, m, j0 - 1, jb, -ONE, b(1, j0), ldb, a(1, j0), lda, &
                                ONE, b, ldb)
                        end if
                    end if
                end do
            end if
        end if

    contains

        ! op(A)(i0:i0+ib-1, same) X = B(i0:i0+ib-1, :): the inverse of the
        ! diagonal block, formed once, applied by GEMM (as optimized BLAS do)
        subroutine solve_left(i0, ib, lower)
            integer, intent(in) :: i0, ib
            logical, intent(in) :: lower
            SCALAR :: tinv(TRI_NB, TRI_NB), w(TRI_NB, n)
            call inverse_block(i0, ib, lower, tinv)
            w(1:ib, 1:n) = b(i0:i0 + ib - 1, 1:n)
            call GEMM('N', 'N', ib, n, ib, ONE, tinv, TRI_NB, w, TRI_NB, ZERO, b(i0, 1), ldb)
        end subroutine solve_left

        ! X op(A)(j0:j0+jb-1, same) = B(:, j0:j0+jb-1)
        subroutine solve_right(j0, jb, upper_op)
            integer, intent(in) :: j0, jb
            logical, intent(in) :: upper_op
            SCALAR :: tinv(TRI_NB, TRI_NB)
            SCALAR, allocatable :: w(:, :)
            call inverse_block(j0, jb, .not. upper_op, tinv)
            allocate (w(m, jb))
            w = b(1:m, j0:j0 + jb - 1)
            call GEMM('N', 'N', m, jb, jb, ONE, w, m, tinv, TRI_NB, ZERO, b(1, j0), ldb)
        end subroutine solve_right

        ! tinv := inverse of op(A)(k0:k0+nb-1, same), a lower (or upper)
        ! triangle, by column-oriented substitution on the identity
        subroutine inverse_block(k0, nb, lower, tinv)
            integer, intent(in) :: k0, nb
            logical, intent(in) :: lower
            SCALAR, intent(out) :: tinv(TRI_NB, TRI_NB)
            SCALAR :: t(TRI_NB, TRI_NB)
            integer :: r, c, j
            do c = 1, nb
                do r = 1, nb
                    if (notrans) then
                        t(r, c) = a(k0 + r - 1, k0 + c - 1)
                    else if (conj) then
                        t(r, c) = CONJ(a(k0 + c - 1, k0 + r - 1))
                    else
                        t(r, c) = a(k0 + c - 1, k0 + r - 1)
                    end if
                end do
                if (unit) t(c, c) = ONE
            end do
            tinv(1:nb, 1:nb) = ZERO
            do j = 1, nb
                tinv(j, j) = ONE
                if (lower) then
                    ! forward: column j of inv(T) has zeros above j
                    do c = j, nb
                        tinv(c, j) = tinv(c, j)/t(c, c)
                        tinv(c + 1:nb, j) = tinv(c + 1:nb, j) - tinv(c, j)*t(c + 1:nb, c)
                    end do
                else
                    do c = j, 1, -1
                        tinv(c, j) = tinv(c, j)/t(c, c)
                        tinv(1:c - 1, j) = tinv(1:c - 1, j) - tinv(c, j)*t(1:c - 1, c)
                    end do
                end if
            end do
        end subroutine inverse_block

    end subroutine TRSM

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
