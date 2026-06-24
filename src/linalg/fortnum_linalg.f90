module fortnum_linalg
    ! Small dense linear algebra for orbit hot paths and device kernels:
    ! closed-form 2x2/3x3 determinant and inverse with a near-singular guard,
    ! plus a fixed-MAX_N dense LU solve with partial pivoting.  Every public
    ! routine is pure and carries !$acc routine seq, so one source serves the
    ! CPU integrators and the OpenACC device kernels with no branching.  No
    ! automatic or allocatable arrays appear anywhere, so the routines are
    ! valid leaf routines inside an OpenACC compute region.
    !
    ! Singularity is a return status, never an error stop or print: device code
    ! cannot abort and host hot loops must not do I/O.  The status is a bare
    ! integer (not fortnum_status_t): a derived type with a character(120)
    ! component is not portable to assign inside !$acc routine seq.  By design
    ! LINALG_OK == FORTNUM_OK == 0, so host callers map the flag trivially.
    !
    ! DERIVATIVE POLICY (ad.md S1): det2/det3/inv2/inv3/jacobian_ok3 are
    !   primal_only as standalone algebra (their derivatives, when needed, are
    !   produced analytically by the consumer, not by AD through the cofactor
    !   arithmetic).  lu_solve realises the linear map x = A^{-1} b; the
    !   implicit/linear-solve rule (d x = A^{-1}(d b - dA x)) belongs to the
    !   consumer that owns A and b, not to the in-place factorisation here.
    !   Do not differentiate through the elimination.
    !
    ! The pivot/singularity logic mirrors solve_linear in fortnum_multiroot
    ! (sing_tol = eps*maxval(|A|)*n); the closed-form guard mirrors the FO Boris
    ! jacobian_ok (|det| > 1e-8 * scale**rank, scale = sqrt(sum(A**2))).
    use fortnum_kinds, only: dp
    implicit none
    private

    ! Largest system lu_solve accepts; documented upper bound for callers and
    ! the test driver.  lu_solve works in place, so it needs no MAX_N scratch.
    integer, parameter, public :: LINALG_MAX_N = 16

    ! Return status.  Zero is success (== FORTNUM_OK); positive is failure.
    ! For inv2/inv3 the failure value is LINALG_SINGULAR; for lu_solve it is the
    ! failing pivot column k (so the value localises the failure, like LAPACK
    ! info).
    integer, parameter, public :: LINALG_OK = 0
    integer, parameter, public :: LINALG_SINGULAR = 1

    ! Relative near-singular threshold for the closed-form inverses, matching
    ! the FO Boris guard so the chartmap inversion behaviour is preserved.
    real(dp), parameter :: SING_TOL_REL = 1.0e-8_dp

    public :: det2, det3, inv2, inv3, jacobian_ok3, lu_solve

contains

    ! Determinant of a 2x2 matrix.
    pure function det2(a) result(d)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2)
        real(dp) :: d
        d = a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1)
    end function det2

    ! Determinant of a 3x3 matrix by cofactor expansion along the first row.
    pure function det3(a) result(d)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3)
        real(dp) :: d
        d = a(1, 1)*(a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2)) &
            - a(1, 2)*(a(2, 1)*a(3, 3) - a(2, 3)*a(3, 1)) &
            + a(1, 3)*(a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))
    end function det3

    ! Near-singular predicate for a 3x3 Jacobian, byte-identical to the FO
    ! Boris guard: reject NaN, then require |det| above a scale-cubed floor.
    pure function jacobian_ok3(a) result(ok)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3)
        logical :: ok
        real(dp) :: d, scale
        d = det3(a)
        scale = sqrt(sum(a**2))
        ok = (d == d) .and. (abs(d) > SING_TOL_REL*max(scale, 1.0e-30_dp)**3)
    end function jacobian_ok3

    ! Closed-form 2x2 inverse with a near-singular guard.  On LINALG_SINGULAR
    ! the output is zeroed and no division is performed.
    pure subroutine inv2(a, ainv, info)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2)
        real(dp), intent(out) :: ainv(2, 2)
        integer, intent(out) :: info
        real(dp) :: d, scale
        d = det2(a)
        scale = sqrt(sum(a**2))
        if ((d /= d) .or. (abs(d) <= SING_TOL_REL*max(scale, 1.0e-30_dp)**2)) then
            ainv = 0.0_dp
            info = LINALG_SINGULAR
            return
        end if
        info = LINALG_OK
        ainv(1, 1) = a(2, 2)/d
        ainv(1, 2) = -a(1, 2)/d
        ainv(2, 1) = -a(2, 1)/d
        ainv(2, 2) = a(1, 1)/d
    end subroutine inv2

    ! Closed-form 3x3 inverse (cofactor / adjugate) with the jacobian_ok3 guard
    ! folded in.  On LINALG_SINGULAR the output is zeroed.
    pure subroutine inv3(a, ainv, info)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3)
        real(dp), intent(out) :: ainv(3, 3)
        integer, intent(out) :: info
        real(dp) :: d
        if (.not. jacobian_ok3(a)) then
            ainv = 0.0_dp
            info = LINALG_SINGULAR
            return
        end if
        info = LINALG_OK
        d = det3(a)
        ainv(1, 1) = (a(2, 2)*a(3, 3) - a(2, 3)*a(3, 2))/d
        ainv(1, 2) = (a(1, 3)*a(3, 2) - a(1, 2)*a(3, 3))/d
        ainv(1, 3) = (a(1, 2)*a(2, 3) - a(1, 3)*a(2, 2))/d
        ainv(2, 1) = (a(2, 3)*a(3, 1) - a(2, 1)*a(3, 3))/d
        ainv(2, 2) = (a(1, 1)*a(3, 3) - a(1, 3)*a(3, 1))/d
        ainv(2, 3) = (a(1, 3)*a(2, 1) - a(1, 1)*a(2, 3))/d
        ainv(3, 1) = (a(2, 1)*a(3, 2) - a(2, 2)*a(3, 1))/d
        ainv(3, 2) = (a(1, 2)*a(3, 1) - a(1, 1)*a(3, 2))/d
        ainv(3, 3) = (a(1, 1)*a(2, 2) - a(1, 2)*a(2, 1))/d
    end subroutine inv3

    ! Solve A x = b (n x n, 1 <= n <= LINALG_MAX_N) by Gaussian elimination with
    ! row partial pivoting, in place: A is overwritten with the elimination
    ! factors and b carries the solution on return.  info = 0 on success; info =
    ! k > 0 when column k has a pivot at or below the scaled singularity
    ! threshold (eps*maxval(|A|)*n), the same relative test as the multiroot
    ! solver.  No scratch arrays, so the routine is valid under routine seq.
    pure subroutine lu_solve(n, a, b, info)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n, n)
        real(dp), intent(inout) :: b(n)
        integer, intent(out) :: info
        real(dp) :: factor, pivmax, s, sing_tol, amax, tmp
        integer :: i, j, k, p

        info = LINALG_OK
        amax = 0.0_dp
        do j = 1, n
            do i = 1, n
                if (abs(a(i, j)) > amax) amax = abs(a(i, j))
            end do
        end do
        sing_tol = epsilon(1.0_dp)*max(amax, tiny(1.0_dp))*real(n, dp)

        do k = 1, n - 1
            ! Partial pivot: largest |a(i,k)| for i >= k.
            p = k
            pivmax = abs(a(k, k))
            do i = k + 1, n
                if (abs(a(i, k)) > pivmax) then
                    pivmax = abs(a(i, k))
                    p = i
                end if
            end do
            if (pivmax <= sing_tol) then
                info = k
                return
            end if
            if (p /= k) then
                do j = 1, n
                    tmp = a(k, j)
                    a(k, j) = a(p, j)
                    a(p, j) = tmp
                end do
                s = b(k); b(k) = b(p); b(p) = s
            end if
            do i = k + 1, n
                factor = a(i, k)/a(k, k)
                do j = k, n
                    a(i, j) = a(i, j) - factor*a(k, j)
                end do
                b(i) = b(i) - factor*b(k)
            end do
        end do

        if (abs(a(n, n)) <= sing_tol) then
            info = n
            return
        end if

        ! Back substitution; b is overwritten with the solution.
        do i = n, 1, -1
            s = b(i)
            do j = i + 1, n
                s = s - a(i, j)*b(j)
            end do
            b(i) = s/a(i, i)
        end do
    end subroutine lu_solve

end module fortnum_linalg
