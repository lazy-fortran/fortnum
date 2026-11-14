module fortnum_lapack_chol
    ! Cholesky factorization and triangular inverses on fortnum's BLAS:
    ! DPOTRF/ZPOTRF, DPOTF2/ZPOTF2, DPOTRS/ZPOTRS, DPOSV/ZPOSV, DPOTRI/ZPOTRI,
    ! DTRTRI/ZTRTRI, DTRTRS/ZTRTRS, the condition estimate DPOCON and the
    ! pivoted Cholesky DPSTRF, with reference LAPACK arguments and info codes
    ! (prefix fortnum_). Also the 1-norm estimator the condition numbers use.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm
    use fortnum_blas3, only: fortnum_dtrsm, fortnum_ztrsm, fortnum_dsyrk, fortnum_zherk
    implicit none
    private

    public :: fortnum_dpotrf, fortnum_dpotf2, fortnum_dpotrs, fortnum_dposv, fortnum_dpotri
    public :: fortnum_dtrtri, fortnum_dtrtrs, fortnum_dpocon, fortnum_dpstrf
    public :: fortnum_zpotrf, fortnum_zpotf2, fortnum_zpotrs, fortnum_zposv, fortnum_zpotri
    public :: fortnum_ztrtri, fortnum_ztrtrs
    public :: norm1_inverse_estimate, znorm1_inverse_estimate

    !> Column block of the blocked Cholesky.
    integer, parameter :: CHOL_NB = 64


contains

#define SCALAR real(dp)
#define CONJ(x) (x)
#define ONE 1.0_dp
#define ZERO 0.0_dp
#define TRANS_C 'T'
#define GEMM fortnum_dgemm
#define TRSM fortnum_dtrsm
#define RANKK fortnum_dsyrk
#define POTF2 fortnum_dpotf2
#define POTRF fortnum_dpotrf
#define POTRS fortnum_dpotrs
#define POSV fortnum_dposv
#define POTRI fortnum_dpotri
#define TRTRI fortnum_dtrtri
#define TRTRS fortnum_dtrtrs
#include "fortnum_lapack_chol.inc"
#undef SCALAR
#undef CONJ
#undef ONE
#undef ZERO
#undef TRANS_C
#undef GEMM
#undef TRSM
#undef RANKK
#undef POTF2
#undef POTRF
#undef POTRS
#undef POSV
#undef POTRI
#undef TRTRI
#undef TRTRS
#define SCALAR complex(dp)
#define CONJ(x) conjg(x)
#define ONE (1.0_dp, 0.0_dp)
#define ZERO (0.0_dp, 0.0_dp)
#define TRANS_C 'C'
#define GEMM fortnum_zgemm
#define TRSM fortnum_ztrsm
#define RANKK fortnum_zherk
#define POTF2 fortnum_zpotf2
#define POTRF fortnum_zpotrf
#define POTRS fortnum_zpotrs
#define POSV fortnum_zposv
#define POTRI fortnum_zpotri
#define TRTRI fortnum_ztrtri
#define TRTRS fortnum_ztrtrs
#include "fortnum_lapack_chol.inc"

    subroutine check_uplo_n(uplo, n, lda, lda_position, upper, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, lda_position
        logical, intent(out) :: upper
        integer, intent(out) :: info
        info = 0
        upper = uplo == 'U' .or. uplo == 'u'
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (lda < max(1, n)) then
            info = -lda_position
        end if
    end subroutine check_uplo_n

    !> 1-norm of inv(A) of order n estimated by reverse communication, as
    !> LAPACK xLACN2 (Hager's method with Higham's refinements). Start with
    !> kase = 0; on return kase = 1 asks for x := inv(A) x, kase = 2 for
    !> x := inv(A)^T x (inv(A)^H x for complex), then call again; kase = 0 on
    !> return means est is final. isave keeps the state between calls.
    subroutine norm1_inverse_estimate(n, v, x, est, kase, isave)
        integer, intent(in) :: n
        real(dp), intent(inout) :: v(n), x(n), est
        integer, intent(inout) :: kase, isave(3)
        integer :: i, j, jlast
        real(dp) :: estold, temp
        if (kase == 0) then
            x = 1.0_dp/n
            kase = 1
            isave(1) = 1
            return
        end if
        select case (isave(1))
        case (1)
            if (n == 1) then
                v(1) = x(1)
                est = abs(v(1))
                kase = 0
                return
            end if
            est = sum(abs(x))
            x = sign(1.0_dp, x)
            kase = 2
            isave(1) = 2
        case (2)
            isave(2) = maxloc(abs(x), 1)
            isave(3) = 2
            call unit_vector(isave(2))
        case (3)
            v = x
            estold = est
            est = sum(abs(v))
            if (est <= estold) then
                call alternating()
                return
            end if
            x = sign(1.0_dp, x)
            kase = 2
            isave(1) = 4
        case (4)
            jlast = isave(2)
            j = maxloc(abs(x), 1)
            if (abs(x(jlast)) /= abs(x(j)) .and. isave(3) < 5) then
                isave(2) = j
                isave(3) = isave(3) + 1
                call unit_vector(j)
            else
                call alternating()
            end if
        case (5)
            temp = 2*sum(abs(x))/(3*n)
            if (temp > est) then
                v = x
                est = temp
            end if
            kase = 0
        end select
    contains
        subroutine unit_vector(k)
            integer, intent(in) :: k
            x = 0
            x(k) = 1
            kase = 1
            isave(1) = 3
        end subroutine unit_vector
        subroutine alternating()
            do i = 1, n
                x(i) = merge(1.0_dp, -1.0_dp, mod(i, 2) == 1)*(1 + real(i - 1, dp)/max(1, n - 1))
            end do
            kase = 1
            isave(1) = 5
        end subroutine alternating
    end subroutine norm1_inverse_estimate

    !> Complex version of norm1_inverse_estimate (xLACN2 for complex).
    subroutine znorm1_inverse_estimate(n, v, x, est, kase, isave)
        integer, intent(in) :: n
        complex(dp), intent(inout) :: v(n), x(n)
        real(dp), intent(inout) :: est
        integer, intent(inout) :: kase, isave(3)
        integer :: i, j, jlast
        real(dp) :: estold, temp
        if (kase == 0) then
            x = cmplx(1.0_dp/n, 0.0_dp, dp)
            kase = 1
            isave(1) = 1
            return
        end if
        select case (isave(1))
        case (1)
            if (n == 1) then
                v(1) = x(1)
                est = abs(v(1))
                kase = 0
                return
            end if
            est = sum(abs(x))
            call unit_phase()
            kase = 2
            isave(1) = 2
        case (2)
            isave(2) = maxloc(real(x, dp), 1)
            isave(3) = 2
            call unit_vector(isave(2))
        case (3)
            v = x
            estold = est
            est = sum(abs(v))
            if (est <= estold) then
                call alternating()
                return
            end if
            call unit_phase()
            kase = 2
            isave(1) = 4
        case (4)
            jlast = isave(2)
            j = maxloc(abs(x), 1)
            if (abs(x(jlast)) /= abs(x(j)) .and. isave(3) < 5) then
                isave(2) = j
                isave(3) = isave(3) + 1
                call unit_vector(j)
            else
                call alternating()
            end if
        case (5)
            temp = 2*sum(abs(x))/(3*n)
            if (temp > est) then
                v = x
                est = temp
            end if
            kase = 0
        end select
    contains
        subroutine unit_phase()
            do i = 1, n
                if (abs(x(i)) > 0) then
                    x(i) = x(i)/abs(x(i))
                else
                    x(i) = 1
                end if
            end do
        end subroutine unit_phase
        subroutine unit_vector(k)
            integer, intent(in) :: k
            x = 0
            x(k) = 1
            kase = 1
            isave(1) = 3
        end subroutine unit_vector
        subroutine alternating()
            do i = 1, n
                x(i) = merge(1.0_dp, -1.0_dp, mod(i, 2) == 1)*(1 + real(i - 1, dp)/max(1, n - 1))
            end do
            kase = 1
            isave(1) = 5
        end subroutine alternating
    end subroutine znorm1_inverse_estimate

    !> Reciprocal 1-norm condition number of a symmetric positive definite
    !> matrix from its Cholesky factor (DPOCON); anorm is the 1-norm of A.
    subroutine fortnum_dpocon(uplo, n, a, lda, anorm, rcond, work, iwork, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        real(dp), intent(in) :: a(lda, *), anorm
        real(dp), intent(out) :: rcond
        real(dp), intent(inout) :: work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: upper
        real(dp) :: ainvnm, v(max(1, n)), x(max(1, n))
        integer :: kase, isave(3), inf
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        if (info == 0 .and. anorm < 0) info = -5
        rcond = 0
        if (info /= 0) return
        if (n == 0) then
            rcond = 1
            return
        end if
        if (anorm == 0) return
        kase = 0
        ainvnm = 0
        do
            call norm1_inverse_estimate(n, v, x, ainvnm, kase, isave)
            if (kase == 0) exit
            ! A is symmetric: inv(A)^T = inv(A)
            call fortnum_dpotrs(uplo, n, 1, a, lda, x, n, inf)
        end do
        if (ainvnm > 0) rcond = (1/ainvnm)/anorm
    end subroutine fortnum_dpocon

    !> Cholesky with complete (diagonal) pivoting of a positive semidefinite
    !> matrix (DPSTRF): P^T A P = L L^T or U^T U up to the computed rank, which
    !> stops at the first pivot <= tol (tol < 0: n eps max diag(A)); info = 1
    !> when the rank is below n.
    subroutine fortnum_dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda
        real(dp), intent(inout) :: a(lda, *)
        integer, intent(out) :: piv(*), rank, info
        real(dp), intent(in) :: tol
        real(dp), intent(inout) :: work(*)
        logical :: upper
        real(dp), allocatable :: s(:, :), d(:), t(:)
        real(dp) :: dstop, tt
        integer :: i, j, p, ti
        call check_uplo_n(uplo, n, lda, 4, upper, info)
        rank = 0
        if (info /= 0 .or. n == 0) return
        ! a full symmetric copy makes the symmetric interchanges simple
        allocate (s(n, n), d(n), t(n))
        do j = 1, n
            do i = 1, n
                if ((i <= j) .eqv. upper .or. i == j) then
                    s(i, j) = a(i, j)
                else
                    s(i, j) = a(j, i)
                end if
            end do
        end do
        piv(1:n) = [(i, i = 1, n)]
        d = [(s(i, i), i = 1, n)]
        p = maxloc(d, 1)
        if (d(p) <= 0 .or. d(p) /= d(p)) then
            info = 1
            return
        end if
        if (tol < 0) then
            dstop = n*epsilon(1.0_dp)*d(p)
        else
            dstop = tol
        end if
        rank = n
        do j = 1, n
            p = j - 1 + maxloc(d(j:n), 1)
            if (d(p) <= dstop .or. d(p) /= d(p)) then
                rank = j - 1
                info = 1
                exit
            end if
            if (p /= j) then
                t = s(j, :)
                s(j, :) = s(p, :)
                s(p, :) = t
                t = s(:, j)
                s(:, j) = s(:, p)
                s(:, p) = t
                ti = piv(j)
                piv(j) = piv(p)
                piv(p) = ti
                tt = d(j)
                d(j) = d(p)
                d(p) = tt
            end if
            s(j, j) = sqrt(d(j))
            do i = j + 1, n
                s(i, j) = (s(i, j) - dot_product(s(i, 1:j - 1), s(j, 1:j - 1)))/s(j, j)
                d(i) = d(i) - s(i, j)**2
            end do
        end do
        do j = 1, n
            do i = 1, n
                if (upper .and. i <= j) then
                    a(i, j) = merge(s(j, i), s(i, j), i <= rank)
                else if (.not. upper .and. i >= j) then
                    a(i, j) = s(i, j)
                end if
            end do
        end do
    end subroutine fortnum_dpstrf

end module fortnum_lapack_chol
