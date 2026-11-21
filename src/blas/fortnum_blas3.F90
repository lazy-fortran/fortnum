module fortnum_blas3
    ! Level-3 BLAS on fortnum_dgemm/zgemm: triangular solve and multiply
    ! (fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm) and the
    ! symmetric and Hermitian rank-k updates (fortnum_dsyrk, fortnum_dsyr2k,
    ! fortnum_zherk). Arguments and semantics are those of reference BLAS.
    use, intrinsic :: iso_c_binding, only: c_loc, c_f_pointer
    use fortnum_kinds, only: dp
    use fortnum_blas, only: thread_count, fortnum_dgemm, fortnum_dgemm_triangle, fortnum_zgemm, fortnum_dgemm_ws, fortnum_zgemm_ws, &
        fortnum_dgemm_lwork, fortnum_zgemm_lwork, fortnum_dtrsm_packed, fortnum_dtrsm_packed_lwork
    implicit none
    private

    public :: fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm
    public :: fortnum_dtrsm_ws, fortnum_ztrsm_ws, fortnum_dtrsm_lwork, fortnum_ztrsm_lwork
    public :: fortnum_dsyrk, fortnum_dsyr2k, fortnum_zherk

    !> Order of the diagonal blocks TRMM handles directly; the rest is GEMM.
    integer, parameter :: TRI_NB = 64
    !> Rows or columns TRSM solves by substitution; it halves larger triangles.
    integer, parameter :: TRSM_LEAF = 16
    !> Rows or columns of B per thread chunk in the TRSM substitution.
    integer, parameter :: SUBST_CHUNK = 256
    !> Fewest rows or columns of B per thread when TRSM splits B over threads.
    integer, parameter :: TRSM_CHUNK = 32
    !> Order below which the rank-k updates compute diagonal blocks in full.
    integer, parameter :: SYRK_LEAF = 32

contains

#define SCALAR real(dp)
#define SCALAR_WORDS 1
#define GEMM fortnum_dgemm
#define GEMM_WS fortnum_dgemm_ws
#define GEMM_LWORK fortnum_dgemm_lwork
#define TRSM fortnum_dtrsm
#define TRSM_WS fortnum_dtrsm_ws
#define TRSM_LWORK fortnum_dtrsm_lwork
#define TRMM fortnum_dtrmm
#define CONJ(x) (x)
#define ONE 1.0_dp
#define ZERO 0.0_dp
#define TRSM_PACKED fortnum_dtrsm_packed
#define TRSM_PACKED_LWORK fortnum_dtrsm_packed_lwork
#include "fortnum_blas3_tri.inc"
#undef TRSM_PACKED
#undef TRSM_PACKED_LWORK
#undef SCALAR
#undef SCALAR_WORDS
#undef GEMM
#undef GEMM_WS
#undef GEMM_LWORK
#undef TRSM
#undef TRSM_WS
#undef TRSM_LWORK
#undef TRMM
#undef CONJ
#undef ONE
#undef ZERO
#define SCALAR complex(dp)
#define SCALAR_WORDS 2
#define GEMM fortnum_zgemm
#define GEMM_WS fortnum_zgemm_ws
#define GEMM_LWORK fortnum_zgemm_lwork
#define TRSM fortnum_ztrsm
#define TRSM_WS fortnum_ztrsm_ws
#define TRSM_LWORK fortnum_ztrsm_lwork
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
        if (notrans) then
            call fortnum_dgemm_triangle(uplo, 'N', 'T', n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            if (two) call fortnum_dgemm_triangle(uplo, 'N', 'T', n, k, alpha, b, ldb, a, lda, 1.0_dp, c, ldc)
        else
            call fortnum_dgemm_triangle(uplo, 'T', 'N', n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            if (two) call fortnum_dgemm_triangle(uplo, 'T', 'N', n, k, alpha, b, ldb, a, lda, 1.0_dp, c, ldc)
        end if
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
        integer :: j

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
        call tri(1, n)

    contains

        recursive subroutine tri(r0, nb)
            integer, intent(in) :: r0, nb
            integer :: n1, i, j
            complex(dp) :: t(SYRK_LEAF, SYRK_LEAF)
            if (nb <= SYRK_LEAF) then
                t(1:nb, 1:nb) = 0
                if (alpha /= 0 .and. k > 0) &
                    call zprod(notrans, k, a, lda, r0, r0, nb, nb, cmplx(alpha, 0.0_dp, dp), &
                        (0.0_dp, 0.0_dp), t, SYRK_LEAF)
                do j = 1, nb
                    do i = merge(1, j, upper), merge(j, nb, upper)
                        if (beta == 0) then
                            c(r0 + i - 1, r0 + j - 1) = t(i, j)
                        else
                            c(r0 + i - 1, r0 + j - 1) = beta*c(r0 + i - 1, r0 + j - 1) + t(i, j)
                        end if
                    end do
                    c(r0 + j - 1, r0 + j - 1) = cmplx(real(c(r0 + j - 1, r0 + j - 1), dp), 0.0_dp, dp)
                end do
                return
            end if
            n1 = nb/2
            call tri(r0, n1)
            if (upper) then
                call rect(r0, n1, r0 + n1, nb - n1)
            else
                call rect(r0 + n1, nb - n1, r0, n1)
            end if
            call tri(r0 + n1, nb - n1)
        end subroutine tri

        subroutine rect(r0, nr, c0, nc)
            integer, intent(in) :: r0, nr, c0, nc
            integer :: jj
            if (alpha == 0 .or. k == 0) then
                do jj = c0, c0 + nc - 1
                    if (beta == 0) then
                        c(r0:r0 + nr - 1, jj) = 0
                    else
                        c(r0:r0 + nr - 1, jj) = beta*c(r0:r0 + nr - 1, jj)
                    end if
                end do
                return
            end if
            call zprod(notrans, k, a, lda, r0, c0, nr, nc, cmplx(alpha, 0.0_dp, dp), &
                cmplx(beta, 0.0_dp, dp), c(r0, c0), ldc)
        end subroutine rect

    end subroutine fortnum_zherk

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
