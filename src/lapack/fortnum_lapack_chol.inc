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
                ! below the diagonal by columns 1..j-1 (contiguous), not by
                ! dot products along the rows of a
                do i = 1, j - 1
                    a(j + 1:n, j) = a(j + 1:n, j) - a(j + 1:n, i)*CONJ(a(j, i))
                end do
                a(j + 1:n, j) = a(j + 1:n, j)/ajj
            end if
        end do
    end subroutine POTF2

    !> Cholesky (xPOTRF) through the recursive POTRF2.
    subroutine POTRF(uplo, n, a, lda, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        logical :: upper
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        if (info /= 0 .or. n == 0) return
        call POTRF2(uplo, n, a, lda, info)
    end subroutine POTRF

    !> Recursive Cholesky (xPOTRF2): the leading half, a triangular solve for
    !> the off-diagonal block, a rank-k update of the trailing half, and the
    !> trailing half recursively.
    recursive subroutine POTRF2(uplo, n, a, lda, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        SCALAR, intent(inout) :: a(lda, *)
        integer, intent(out) :: info
        integer :: n1, n2, iinfo
        info = 0
        if (n == 0) return
        if (n <= 32) then
            ! small blocks directly: the recursion would be all call overhead
            call POTF2(uplo, n, a, lda, info)
            return
        end if
        n1 = n/2
        n2 = n - n1
        call POTRF2(uplo, n1, a, lda, iinfo)
        if (iinfo /= 0) then
            info = iinfo
            return
        end if
        if (uplo == 'U' .or. uplo == 'u') then
            call TRSM('L', 'U', TRANS_C, 'N', n1, n2, ONE, a, lda, a(1, n1 + 1), lda)
            call RANKK('U', TRANS_C, n2, n1, -1.0_dp, a(1, n1 + 1), lda, 1.0_dp, a(n1 + 1, n1 + 1), lda)
        else
            call TRSM('R', 'L', TRANS_C, 'N', n2, n1, ONE, a, lda, a(n1 + 1, 1), lda)
            call RANKK('L', 'N', n2, n1, -1.0_dp, a(n1 + 1, 1), lda, 1.0_dp, a(n1 + 1, n1 + 1), lda)
        end if
        call POTRF2(uplo, n2, a(n1 + 1, n1 + 1), lda, iinfo)
        if (iinfo /= 0) info = iinfo + n1
    end subroutine POTRF2

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
