module fortnum_blas3
    ! Level-3 BLAS on fortnum_dgemm/zgemm: triangular solve and multiply
    ! (fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm) and the
    ! symmetric and Hermitian rank-k updates (fortnum_dsyrk, fortnum_dsyr2k,
    ! fortnum_zherk). Arguments and semantics are those of reference BLAS.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm
    implicit none
    private

    public :: fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm
    public :: fortnum_dsyrk, fortnum_dsyr2k, fortnum_zherk

    !> Order of the diagonal blocks handled directly; the rest is GEMM.
    integer, parameter :: TRI_NB = 64
    !> Column blocks of the rank-k updates.
    integer, parameter :: SYRK_NB = 128

contains

#define SCALAR real(dp)
#define GEMM fortnum_dgemm
#define TRSM fortnum_dtrsm
#define TRMM fortnum_dtrmm
#define CONJ(x) (x)
#define ONE 1.0_dp
#define ZERO 0.0_dp
#include "fortnum_blas3_tri.inc"
#undef SCALAR
#undef GEMM
#undef TRSM
#undef TRMM
#undef CONJ
#undef ONE
#undef ZERO
#define SCALAR complex(dp)
#define GEMM fortnum_zgemm
#define TRSM fortnum_ztrsm
#define TRMM fortnum_ztrmm
#define CONJ(x) conjg(x)
#define ONE (1.0_dp, 0.0_dp)
#define ZERO (0.0_dp, 0.0_dp)
#include "fortnum_blas3_tri.inc"

    subroutine tri_flags(name, side, uplo, transa, diag, m, n, lda, ldb, left, upper, notrans, &
            conj, unit, tr)
        character(*), intent(in) :: name
        character, intent(in) :: side, uplo, transa, diag
        integer, intent(in) :: m, n, lda, ldb
        logical, intent(out) :: left, upper, notrans, conj, unit
        character, intent(out) :: tr
        select case (side)
        case ('L', 'l')
            left = .true.
        case ('R', 'r')
            left = .false.
        case default
            call argument_error(name, 1)
        end select
        upper = flag(uplo, 'U', 'L', name, 2)
        select case (transa)
        case ('N', 'n')
            notrans = .true.
            conj = .false.
            tr = 'N'
        case ('T', 't')
            notrans = .false.
            conj = .false.
            tr = 'T'
        case ('C', 'c')
            notrans = .false.
            conj = .true.
            tr = 'C'
        case default
            call argument_error(name, 3)
        end select
        unit = flag(diag, 'U', 'N', name, 4)
        if (m < 0) call argument_error(name, 5)
        if (n < 0) call argument_error(name, 6)
        if (lda < max(1, merge(m, n, left))) call argument_error(name, 9)
        if (ldb < max(1, m)) call argument_error(name, 11)
    end subroutine tri_flags

    !> C := alpha op(A) op(A)^T + beta C on the uplo triangle, op(A) = A (n x k)
    !> or A^T; reference BLAS DSYRK.
    subroutine fortnum_dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k, lda, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *)
        real(dp), intent(inout) :: c(ldc, *)
        call fortnum_dsyr2k_impl('DSYRK', uplo, trans, n, k, alpha, a, lda, a, lda, beta, c, ldc, .false.)
    end subroutine fortnum_dsyrk

    !> C := alpha (op(A) op(B)^T + op(B) op(A)^T) + beta C on the uplo
    !> triangle; reference BLAS DSYR2K.
    subroutine fortnum_dsyr2k(uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        call fortnum_dsyr2k_impl('DSYR2K', uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, .true.)
    end subroutine fortnum_dsyr2k

    subroutine fortnum_dsyr2k_impl(name, uplo, trans, n, k, alpha, a, lda, b, ldb, beta, c, ldc, two)
        character(*), intent(in) :: name
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        logical, intent(in) :: two
        logical :: upper, notrans
        integer :: j0, jb, i, j
        real(dp), allocatable :: t(:, :)

        upper = flag(uplo, 'U', 'L', name, 1)
        select case (trans)
        case ('N', 'n')
            notrans = .true.
        case ('T', 't', 'C', 'c')
            notrans = .false.
        case default
            call argument_error(name, 2)
        end select
        if (n < 0) call argument_error(name, 3)
        if (k < 0) call argument_error(name, 4)
        if (lda < max(1, merge(n, k, notrans))) call argument_error(name, 7)
        if (two .and. ldb < max(1, merge(n, k, notrans))) call argument_error(name, 9)
        if (ldc < max(1, n)) call argument_error(name, merge(12, 10, two))
        if (n == 0) return
        if ((alpha == 0 .or. k == 0) .and. beta == 1) return
        allocate (t(SYRK_NB, SYRK_NB))
        do j0 = 1, n, SYRK_NB
            jb = min(SYRK_NB, n - j0 + 1)
            ! the diagonal block through a full product, then its triangle
            t(1:jb, 1:jb) = 0
            if (alpha /= 0 .and. k > 0) then
                call dprod(notrans, k, a, lda, j0, b, ldb, j0, jb, jb, alpha, 0.0_dp, t, SYRK_NB)
                if (two) call dprod(notrans, k, b, ldb, j0, a, lda, j0, jb, jb, alpha, 1.0_dp, t, SYRK_NB)
            end if
            do j = 1, jb
                do i = merge(1, j, upper), merge(j, jb, upper)
                    if (beta == 0) then
                        c(j0 + i - 1, j0 + j - 1) = t(i, j)
                    else
                        c(j0 + i - 1, j0 + j - 1) = beta*c(j0 + i - 1, j0 + j - 1) + t(i, j)
                    end if
                end do
            end do
            ! the block of the triangle beside it
            if (upper .and. j0 > 1) then
                call rect(1, j0 - 1)
            else if (.not. upper .and. j0 + jb <= n) then
                call rect(j0 + jb, n - j0 - jb + 1)
            end if
        end do

    contains

        subroutine rect(r0, nr)
            integer, intent(in) :: r0, nr
            integer :: jj
            if (alpha == 0 .or. k == 0) then
                do jj = j0, j0 + jb - 1
                    if (beta == 0) then
                        c(r0:r0 + nr - 1, jj) = 0
                    else
                        c(r0:r0 + nr - 1, jj) = beta*c(r0:r0 + nr - 1, jj)
                    end if
                end do
                return
            end if
            call dprod(notrans, k, a, lda, r0, b, ldb, j0, nr, jb, alpha, beta, c(r0, j0), ldc)
            if (two) call dprod(notrans, k, b, ldb, r0, a, lda, j0, nr, jb, alpha, 1.0_dp, c(r0, j0), ldc)
        end subroutine rect

    end subroutine fortnum_dsyr2k_impl

    !> C := alpha op(A) op(A)^H + beta C on the uplo triangle, alpha and beta
    !> real, op(A) = A (n x k) or A^H; the diagonal stays real. Reference
    !> BLAS ZHERK.
    subroutine fortnum_zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k, lda, ldc
        real(dp), intent(in) :: alpha, beta
        complex(dp), intent(in) :: a(lda, *)
        complex(dp), intent(inout) :: c(ldc, *)
        logical :: upper, notrans
        integer :: j0, jb, i, j, r0, nr
        complex(dp), allocatable :: t(:, :)

        upper = flag(uplo, 'U', 'L', 'ZHERK', 1)
        select case (trans)
        case ('N', 'n')
            notrans = .true.
        case ('C', 'c')
            notrans = .false.
        case default
            call argument_error('ZHERK', 2)
        end select
        if (n < 0) call argument_error('ZHERK', 3)
        if (k < 0) call argument_error('ZHERK', 4)
        if (lda < max(1, merge(n, k, notrans))) call argument_error('ZHERK', 7)
        if (ldc < max(1, n)) call argument_error('ZHERK', 10)
        if (n == 0) return
        if ((alpha == 0 .or. k == 0) .and. beta == 1) then
            do j = 1, n
                c(j, j) = cmplx(real(c(j, j), dp), 0.0_dp, dp)
            end do
            return
        end if
        allocate (t(SYRK_NB, SYRK_NB))
        do j0 = 1, n, SYRK_NB
            jb = min(SYRK_NB, n - j0 + 1)
            t(1:jb, 1:jb) = 0
            if (alpha /= 0 .and. k > 0) &
                call zprod(notrans, k, a, lda, j0, j0, jb, jb, cmplx(alpha, 0.0_dp, dp), &
                    (0.0_dp, 0.0_dp), t, SYRK_NB)
            do j = 1, jb
                do i = merge(1, j, upper), merge(j, jb, upper)
                    if (beta == 0) then
                        c(j0 + i - 1, j0 + j - 1) = t(i, j)
                    else
                        c(j0 + i - 1, j0 + j - 1) = beta*c(j0 + i - 1, j0 + j - 1) + t(i, j)
                    end if
                end do
                c(j0 + j - 1, j0 + j - 1) = cmplx(real(c(j0 + j - 1, j0 + j - 1), dp), 0.0_dp, dp)
            end do
            if (upper) then
                r0 = 1
                nr = j0 - 1
            else
                r0 = j0 + jb
                nr = n - j0 - jb + 1
            end if
            if (nr <= 0) cycle
            if (alpha == 0 .or. k == 0) then
                do j = j0, j0 + jb - 1
                    if (beta == 0) then
                        c(r0:r0 + nr - 1, j) = 0
                    else
                        c(r0:r0 + nr - 1, j) = beta*c(r0:r0 + nr - 1, j)
                    end if
                end do
            else
                call zprod(notrans, k, a, lda, r0, j0, nr, jb, cmplx(alpha, 0.0_dp, dp), &
                    cmplx(beta, 0.0_dp, dp), c(r0, j0), ldc)
            end if
        end do

    end subroutine fortnum_zherk

    ! c := alpha op(x)(rx:rx+rows-1, :) op(y)(ry:ry+cols-1, :)^T + beta c with
    ! op(x) = x (rows of x) or x^T (columns of x), depth k.
    subroutine dprod(notrans, k, x, ldx, rx, y, ldy, ry, rows, cols, alpha, beta, c, ldc)
        logical, intent(in) :: notrans
        integer, intent(in) :: k, ldx, rx, ldy, ry, rows, cols, ldc
        real(dp), intent(in) :: x(ldx, *), y(ldy, *), alpha, beta
        real(dp), intent(inout) :: c(ldc, *)
        if (notrans) then
            call fortnum_dgemm('N', 'T', rows, cols, k, alpha, x(rx, 1), ldx, y(ry, 1), ldy, beta, c, ldc)
        else
            call fortnum_dgemm('T', 'N', rows, cols, k, alpha, x(1, rx), ldx, y(1, ry), ldy, beta, c, ldc)
        end if
    end subroutine dprod

    ! c := alpha op(x)(rx.., :) op(x)(ry.., :)^H + beta c, op(x) = x or x^H.
    subroutine zprod(notrans, k, x, ldx, rx, ry, rows, cols, alpha, beta, c, ldc)
        logical, intent(in) :: notrans
        integer, intent(in) :: k, ldx, rx, ry, rows, cols, ldc
        complex(dp), intent(in) :: x(ldx, *), alpha, beta
        complex(dp), intent(inout) :: c(ldc, *)
        if (notrans) then
            call fortnum_zgemm('N', 'C', rows, cols, k, alpha, x(rx, 1), ldx, x(ry, 1), ldx, beta, c, ldc)
        else
            call fortnum_zgemm('C', 'N', rows, cols, k, alpha, x(1, rx), ldx, x(1, ry), ldx, beta, c, ldc)
        end if
    end subroutine zprod

    logical function flag(value, yes, no, name, position)
        character, intent(in) :: value, yes, no
        character(*), intent(in) :: name
        integer, intent(in) :: position
        if (value == yes .or. value == achar(iachar(yes) + 32)) then
            flag = .true.
        else if (value == no .or. value == achar(iachar(no) + 32)) then
            flag = .false.
        else
            flag = .false.
            call argument_error(name, position)
        end if
    end function flag

    subroutine argument_error(routine, position)
        character(*), intent(in) :: routine
        integer, intent(in) :: position
        character(len=16) :: number
        write (number, '(i0)') position
        error stop routine//': argument '//trim(number)//' has an illegal value'
    end subroutine argument_error

end module fortnum_blas3
