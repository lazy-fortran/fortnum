! Cholesky factorization and triangular inverses, real or complex through
! the macros of fortnum_lapack_chol.F90. Arguments and info codes are those
! of reference LAPACK; an illegal argument returns info = -position.

    !> Unblocked Cholesky (xPOTF2): A = U^H U (uplo U) or L L^H (uplo L).
    subroutine POTF2(uplo, n, a, lda, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        logical :: upper
        integer :: j, i
        real(dp) :: ajj
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        if (info /= 0) return
        do j = 1, n
            if (upper) then
                ajj = real(a(j, j), dp) - real(sum(CONJ(a(1:j - 1, j))*a(1:j - 1, j)), dp)
            else
                ajj = real(a(j, j), dp) - real(sum(a(j, 1:j - 1)*CONJ(a(j, 1:j - 1))), dp)
            end if
            if (ajj <= 0 .or. ajj /= ajj) then
                a(j, j) = ajj
                info = j
                return
            end if
            ajj = sqrt(ajj)
            a(j, j) = ajj
            if (upper) then
                do i = j + 1, n
                    a(j, i) = (a(j, i) - sum(CONJ(a(1:j - 1, j))*a(1:j - 1, i)))/ajj
                end do
            else
                do i = j + 1, n
                    a(i, j) = (a(i, j) - sum(a(i, 1:j - 1)*CONJ(a(j, 1:j - 1))))/ajj
                end do
            end if
        end do
    end subroutine POTF2

    !> Blocked Cholesky (xPOTRF): per block of columns, the diagonal block is
    !> updated by a rank-k update and factored directly, the block beside it
    !> by GEMM and a triangular solve.
    subroutine POTRF(uplo, n, a, lda, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        logical :: upper
        integer :: j, jb, rest
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        if (info /= 0 .or. n == 0) return
        do j = 1, n, CHOL_NB
            jb = min(CHOL_NB, n - j + 1)
            rest = n - j - jb + 1
            if (upper) then
                call RANKK('U', TRANS_C, jb, j - 1, -1.0_dp, a(1, j), lda, 1.0_dp, a(j, j), lda)
                call POTF2('U', jb, a(j, j), lda, info)
                if (info /= 0) then
                    info = info + j - 1
                    return
                end if
                if (rest > 0) then
                    if (j > 1) call GEMM(TRANS_C, 'N', jb, rest, j - 1, -ONE, a(1, j), lda, a(1, j + jb), lda, &
                        ONE, a(j, j + jb), lda)
                    call TRSM('L', 'U', TRANS_C, 'N', jb, rest, ONE, a(j, j), lda, a(j, j + jb), lda)
                end if
            else
                call RANKK('L', 'N', jb, j - 1, -1.0_dp, a(j, 1), lda, 1.0_dp, a(j, j), lda)
                call POTF2('L', jb, a(j, j), lda, info)
                if (info /= 0) then
                    info = info + j - 1
                    return
                end if
                if (rest > 0) then
                    if (j > 1) call GEMM('N', TRANS_C, rest, jb, j - 1, -ONE, a(j + jb, 1), lda, a(j, 1), lda, &
                        ONE, a(j + jb, j), lda)
                    call TRSM('R', 'L', TRANS_C, 'N', rest, jb, ONE, a(j, j), lda, a(j + jb, j), lda)
                end if
            end if
        end do
    end subroutine POTRF

    !> Solve A X = B with the Cholesky factor of POTRF (xPOTRS).
    subroutine POTRS(uplo, n, nrhs, a, lda, b, ldb, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, nrhs, lda, ldb
        SCALAR, intent(in) :: a(lda, *)
        SCALAR, intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        logical :: upper
        call check_uplo_n(uplo, n, lda, 5, upper, info)
        if (info == 0 .and. nrhs < 0) info = -3
        if (info == 0 .and. ldb < max(1, n)) info = -7
        if (info /= 0 .or. n == 0 .or. nrhs == 0) return
        if (upper) then
            call TRSM('L', 'U', TRANS_C, 'N', n, nrhs, ONE, a, lda, b, ldb)
            call TRSM('L', 'U', 'N', 'N', n, nrhs, ONE, a, lda, b, ldb)
        else
            call TRSM('L', 'L', 'N', 'N', n, nrhs, ONE, a, lda, b, ldb)
            call TRSM('L', 'L', TRANS_C, 'N', n, nrhs, ONE, a, lda, b, ldb)
        end if
    end subroutine POTRS

    !> Solve A X = B for A Hermitian positive definite (xPOSV).
    subroutine POSV(uplo, n, nrhs, a, lda, b, ldb, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, nrhs, lda, ldb
        SCALAR, intent(inout) :: a(lda, *), b(ldb, *)
        integer, intent(out) :: info
        logical :: upper
        call check_uplo_n(uplo, n, lda, 5, upper, info)
        if (info == 0 .and. nrhs < 0) info = -3
        if (info == 0 .and. ldb < max(1, n)) info = -7
        if (info /= 0) return
        call POTRF(uplo, n, a, lda, info)
        if (info == 0) call POTRS(uplo, n, nrhs, a, lda, b, ldb, info)
    end subroutine POSV

    !> Inverse of a triangular matrix in place (xTRTRI); info = j for a
    !> zero diagonal element of a non-unit matrix.
    subroutine TRTRI(uplo, diag, n, a, lda, info)
        character, intent(in) :: uplo, diag
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        logical :: upper, unit
        integer :: j
        SCALAR, allocatable :: x(:, :)
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        unit = diag == 'U' .or. diag == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (.not. (unit .or. diag == 'N' .or. diag == 'n')) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        end if
        if (info /= 0 .or. n == 0) return
        if (.not. unit) then
            do j = 1, n
                if (a(j, j) == ZERO) then
                    info = j
                    return
                end if
            end do
        end if
        allocate (x(n, n))
        x = ZERO
        do j = 1, n
            x(j, j) = ONE
        end do
        call TRSM('L', uplo, 'N', diag, n, n, ONE, a, lda, x, n)
        ! a unit triangle keeps its stored diagonal, as in LAPACK
        do j = 1, n
            if (upper) then
                a(1:j - 1, j) = x(1:j - 1, j)
            else
                a(j + 1:n, j) = x(j + 1:n, j)
            end if
            if (.not. unit) a(j, j) = x(j, j)
        end do
    end subroutine TRTRI

    !> Solve op(A) X = B for triangular A after checking for singularity
    !> (xTRTRS).
    subroutine TRTRS(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
        character, intent(in) :: uplo, trans, diag
        integer, intent(in) :: n, nrhs, lda, ldb
        SCALAR, intent(in) :: a(lda, *)
        SCALAR, intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        logical :: unit
        integer :: j
        info = 0
        unit = diag == 'U' .or. diag == 'u'
        if (index('UuLl', uplo) == 0) then
            info = -1
        else if (index('NnTtCc', trans) == 0) then
            info = -2
        else if (.not. (unit .or. diag == 'N' .or. diag == 'n')) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (nrhs < 0) then
            info = -5
        else if (lda < max(1, n)) then
            info = -7
        else if (ldb < max(1, n)) then
            info = -9
        end if
        if (info /= 0 .or. n == 0) return
        if (.not. unit) then
            do j = 1, n
                if (a(j, j) == ZERO) then
                    info = j
                    return
                end if
            end do
        end if
        call TRSM('L', uplo, trans, diag, n, nrhs, ONE, a, lda, b, ldb)
    end subroutine TRTRS

    !> Inverse of a Hermitian positive definite matrix from its Cholesky
    !> factor (xPOTRI): inv(A) = inv(U) inv(U)^H or inv(L)^H inv(L).
    subroutine POTRI(uplo, n, a, lda, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        logical :: upper
        integer :: j
        SCALAR, allocatable :: w(:, :)
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        if (info /= 0 .or. n == 0) return
        allocate (w(n, n))
        w = ZERO
        do j = 1, n
            if (upper) then
                w(1:j, j) = a(1:j, j)
            else
                w(j:n, j) = a(j:n, j)
            end if
        end do
        call TRTRI(uplo, 'N', n, w, n, info)
        if (info /= 0) return
        if (upper) then
            call RANKK('U', 'N', n, n, 1.0_dp, w, n, 0.0_dp, a, lda)
        else
            call RANKK('L', TRANS_C, n, n, 1.0_dp, w, n, 0.0_dp, a, lda)
        end if
    end subroutine POTRI
