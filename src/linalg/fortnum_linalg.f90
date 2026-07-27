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
    ! DERIVATIVE POLICY (ad.md S1): det2/det3 have fortsym-generated analytical
    !   JVP and VJP candidates. inv2/inv3 have fused analytical JVP and VJP
    !   candidates; jacobian_ok3 remains primal-only. lu_solve realises
    !   A^{-1} b; the
    !   implicit/linear-solve rule (d x = A^{-1}(d b - dA x)) belongs to the
    !   consumer that owns A and b, not to the in-place factorisation here.
    !   Do not differentiate through the elimination.
    !
    ! The pivot/singularity logic mirrors solve_linear in fortnum_multiroot
    ! (sing_tol = eps*maxval(|A|)*n); the closed-form guard mirrors the FO Boris
    ! jacobian_ok (|det| > 1e-8 * scale**rank, scale = sqrt(sum(A**2))).
    use fortnum_kinds, only: dp
    use fortnum_generated_det2_jvp, only: fortnum_det2_jvp_kernel
    use fortnum_generated_det2_vjp, only: fortnum_det2_vjp_kernel
    use fortnum_generated_det3_jvp, only: fortnum_det3_jvp_kernel
    use fortnum_generated_det3_vjp, only: fortnum_det3_vjp_kernel
    use fortnum_generated_inv2_jvp, only: fortnum_inv2_jvp_kernel
    use fortnum_generated_inv2_vjp, only: fortnum_inv2_vjp_kernel
    use fortnum_generated_inv3_jvp, only: fortnum_inv3_jvp_kernel
    use fortnum_generated_inv3_vjp, only: fortnum_inv3_vjp_kernel
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

    public :: det2, det3, det2_jvp, det3_jvp, det2_vjp, det3_vjp
    public :: inv2, inv3, inv2_jvp, inv3_jvp, inv2_vjp, inv3_vjp, jacobian_ok3
    public :: lu_factor, lu_solve_factored, lu_solve
    public :: linear_solve_jvp, linear_solve_jvp_factored
    public :: linear_solve_vjp, linear_solve_vjp_factored

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

    ! Fortsym-generated analytical directional derivative of det2.
    pure subroutine det2_jvp(a, va, jv)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2), va(2, 2)
        real(dp), intent(out) :: jv

        call fortnum_det2_jvp_kernel(a(1, 1), a(2, 1), a(1, 2), a(2, 2), &
            va(1, 1), va(2, 1), va(1, 2), va(2, 2), jv)
    end subroutine det2_jvp

    ! Fortsym-generated analytical directional derivative of det3.
    pure subroutine det3_jvp(a, va, jv)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3), va(3, 3)
        real(dp), intent(out) :: jv

        call fortnum_det3_jvp_kernel(a(1, 1), a(2, 1), a(3, 1), &
            a(1, 2), a(2, 2), a(3, 2), a(1, 3), a(2, 3), a(3, 3), &
            va(1, 1), va(2, 1), va(3, 1), va(1, 2), va(2, 2), va(3, 2), &
            va(1, 3), va(2, 3), va(3, 3), jv)
    end subroutine det3_jvp

    ! Fortsym-generated analytical transpose product for det2.
    pure subroutine det2_vjp(a, u, abar)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2), u
        real(dp), intent(out) :: abar(2, 2)

        call fortnum_det2_vjp_kernel(a(1, 1), a(2, 1), a(1, 2), a(2, 2), u, &
            abar(1, 1), abar(2, 1), abar(1, 2), abar(2, 2))
    end subroutine det2_vjp

    ! Fortsym-generated analytical transpose product for det3.
    pure subroutine det3_vjp(a, u, abar)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3), u
        real(dp), intent(out) :: abar(3, 3)

        call fortnum_det3_vjp_kernel(a(1, 1), a(2, 1), a(3, 1), &
            a(1, 2), a(2, 2), a(3, 2), a(1, 3), a(2, 3), a(3, 3), u, &
            abar(1, 1), abar(2, 1), abar(3, 1), abar(1, 2), abar(2, 2), &
            abar(3, 2), abar(1, 3), abar(2, 3), abar(3, 3))
    end subroutine det3_vjp

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

    ! Guarded inverse value and analytical JVP, sharing the primal inverse.
    pure subroutine inv2_jvp(a, va, ainv, vainv, info)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2), va(2, 2)
        real(dp), intent(out) :: ainv(2, 2), vainv(2, 2)
        integer, intent(out) :: info

        call inv2(a, ainv, info)
        if (info /= LINALG_OK) then
            vainv = 0.0_dp
            return
        end if
        call fortnum_inv2_jvp_kernel( &
            ainv(1, 1), ainv(2, 1), ainv(1, 2), ainv(2, 2), &
            va(1, 1), va(2, 1), va(1, 2), va(2, 2), &
            vainv(1, 1), vainv(2, 1), vainv(1, 2), vainv(2, 2))
    end subroutine inv2_jvp

    ! Guarded inverse value and analytical JVP, sharing the primal inverse.
    pure subroutine inv3_jvp(a, va, ainv, vainv, info)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3), va(3, 3)
        real(dp), intent(out) :: ainv(3, 3), vainv(3, 3)
        integer, intent(out) :: info

        call inv3(a, ainv, info)
        if (info /= LINALG_OK) then
            vainv = 0.0_dp
            return
        end if
        call fortnum_inv3_jvp_kernel( &
            ainv(1, 1), ainv(2, 1), ainv(3, 1), &
            ainv(1, 2), ainv(2, 2), ainv(3, 2), &
            ainv(1, 3), ainv(2, 3), ainv(3, 3), &
            va(1, 1), va(2, 1), va(3, 1), &
            va(1, 2), va(2, 2), va(3, 2), &
            va(1, 3), va(2, 3), va(3, 3), &
            vainv(1, 1), vainv(2, 1), vainv(3, 1), &
            vainv(1, 2), vainv(2, 2), vainv(3, 2), &
            vainv(1, 3), vainv(2, 3), vainv(3, 3))
    end subroutine inv3_jvp

    ! Guarded inverse value and analytical VJP, sharing the primal inverse.
    pure subroutine inv2_vjp(a, u, ainv, abar, info)
        !$acc routine seq
        real(dp), intent(in) :: a(2, 2), u(2, 2)
        real(dp), intent(out) :: ainv(2, 2), abar(2, 2)
        integer, intent(out) :: info

        call inv2(a, ainv, info)
        if (info /= LINALG_OK) then
            abar = 0.0_dp
            return
        end if
        call fortnum_inv2_vjp_kernel( &
            ainv(1, 1), ainv(2, 1), ainv(1, 2), ainv(2, 2), &
            u(1, 1), u(2, 1), u(1, 2), u(2, 2), &
            abar(1, 1), abar(2, 1), abar(1, 2), abar(2, 2))
    end subroutine inv2_vjp

    ! Guarded inverse value and analytical VJP, sharing the primal inverse.
    pure subroutine inv3_vjp(a, u, ainv, abar, info)
        !$acc routine seq
        real(dp), intent(in) :: a(3, 3), u(3, 3)
        real(dp), intent(out) :: ainv(3, 3), abar(3, 3)
        integer, intent(out) :: info

        call inv3(a, ainv, info)
        if (info /= LINALG_OK) then
            abar = 0.0_dp
            return
        end if
        call fortnum_inv3_vjp_kernel( &
            ainv(1, 1), ainv(2, 1), ainv(3, 1), &
            ainv(1, 2), ainv(2, 2), ainv(3, 2), &
            ainv(1, 3), ainv(2, 3), ainv(3, 3), &
            u(1, 1), u(2, 1), u(3, 1), &
            u(1, 2), u(2, 2), u(3, 2), &
            u(1, 3), u(2, 3), u(3, 3), &
            abar(1, 1), abar(2, 1), abar(3, 1), &
            abar(1, 2), abar(2, 2), abar(3, 2), &
            abar(1, 3), abar(2, 3), abar(3, 3))
    end subroutine inv3_vjp

    ! Compact LU factorization with partial pivoting. A is overwritten by L and
    ! U, while ipiv records each row swap for reuse by subsequent solves.
    pure subroutine lu_factor(n, a, ipiv, info)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n, n)
        integer, intent(out) :: ipiv(n)
        integer, intent(out) :: info
        real(dp) :: factor, pivmax, sing_tol, amax, tmp
        integer :: i, j, k, p

        info = LINALG_OK
        do i = 1, n
            ipiv(i) = i
        end do
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
            ipiv(k) = p
            if (p /= k) then
                do j = 1, n
                    tmp = a(k, j)
                    a(k, j) = a(p, j)
                    a(p, j) = tmp
                end do
            end if
            do i = k + 1, n
                factor = a(i, k)/a(k, k)
                a(i, k) = factor
                do j = k + 1, n
                    a(i, j) = a(i, j) - factor*a(k, j)
                end do
            end do
        end do

        if (abs(a(n, n)) <= sing_tol) then
            info = n
            return
        end if
    end subroutine lu_factor

    ! Solve A x = b from a successful lu_factor result without refactorizing A.
    pure subroutine lu_solve_factored(n, a, ipiv, b, info)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n)
        integer, intent(in) :: ipiv(n)
        real(dp), intent(inout) :: b(n)
        integer, intent(out) :: info
        real(dp) :: s
        integer :: i, j, k, p

        info = LINALG_OK
        do k = 1, n - 1
            p = ipiv(k)
            if ((p < k) .or. (p > n)) then
                info = k
                return
            end if
            if (p /= k) then
                s = b(k)
                b(k) = b(p)
                b(p) = s
            end if
        end do

        do i = 2, n
            s = b(i)
            do j = 1, i - 1
                s = s - a(i, j)*b(j)
            end do
            b(i) = s
        end do

        do i = n, 1, -1
            s = b(i)
            do j = i + 1, n
                s = s - a(i, j)*b(j)
            end do
            if (a(i, i) == 0.0_dp) then
                info = i
                return
            end if
            b(i) = s/a(i, i)
        end do
    end subroutine lu_solve_factored

    ! Solve A x = b (n x n, 1 <= n <= LINALG_MAX_N) in place. This compatibility
    ! entry point factors once and immediately consumes the factors.
    pure subroutine lu_solve(n, a, b, info)
        !$acc routine seq
        integer, intent(in) :: n
        real(dp), intent(inout) :: a(n, n)
        real(dp), intent(inout) :: b(n)
        integer, intent(out) :: info
        integer :: ipiv(n)

        call lu_factor(n, a, ipiv, info)
        if (info /= LINALG_OK) return
        call lu_solve_factored(n, a, ipiv, b, info)
    end subroutine lu_solve

    ! Analytical implicit JVP for A x = b: A dx = db - dA x. The caller
    ! supplies the converged primal x, so this neither repeats the primal solve
    ! nor differentiates the elimination algorithm.
    pure subroutine linear_solve_jvp(n, a, x, da, db, dx, info)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n), x(n), da(n, n), db(n)
        real(dp), intent(out) :: dx(n)
        integer, intent(out) :: info
        real(dp) :: work_a(n, n)
        integer :: ipiv(n)

        work_a = a
        call lu_factor(n, work_a, ipiv, info)
        if (info /= LINALG_OK) then
            dx = 0.0_dp
            return
        end if
        call linear_solve_jvp_factored(n, work_a, ipiv, x, da, db, dx, info)
    end subroutine linear_solve_jvp

    ! Analytical JVP using a reusable factorization of the primal matrix.
    pure subroutine linear_solve_jvp_factored(n, a, ipiv, x, da, db, dx, info)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n), x(n), da(n, n), db(n)
        integer, intent(in) :: ipiv(n)
        real(dp), intent(out) :: dx(n)
        integer, intent(out) :: info

        dx = db - matmul(da, x)
        call lu_solve_factored(n, a, ipiv, dx, info)
        if (info /= LINALG_OK) dx = 0.0_dp
    end subroutine linear_solve_jvp_factored

    ! Analytical implicit VJP. Solving A^T lambda = u gives
    ! b_bar=lambda and A_bar=-lambda*x^T.
    pure subroutine linear_solve_vjp(n, a, x, u, abar, bbar, info)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n), x(n), u(n)
        real(dp), intent(out) :: abar(n, n), bbar(n)
        integer, intent(out) :: info
        real(dp) :: work_a(n, n)
        integer :: ipiv(n)

        work_a = transpose(a)
        call lu_factor(n, work_a, ipiv, info)
        if (info /= LINALG_OK) then
            abar = 0.0_dp
            bbar = 0.0_dp
            return
        end if
        call linear_solve_vjp_factored(n, work_a, ipiv, x, u, abar, bbar, info)
    end subroutine linear_solve_vjp

    ! Analytical VJP using a reusable factorization of the transposed matrix.
    pure subroutine linear_solve_vjp_factored(n, at, ipiv, x, u, abar, bbar, info)
        integer, intent(in) :: n
        real(dp), intent(in) :: at(n, n), x(n), u(n)
        integer, intent(in) :: ipiv(n)
        real(dp), intent(out) :: abar(n, n), bbar(n)
        integer, intent(out) :: info
        integer :: i, j

        bbar = u
        call lu_solve_factored(n, at, ipiv, bbar, info)
        if (info /= LINALG_OK) then
            abar = 0.0_dp
            bbar = 0.0_dp
            return
        end if
        do j = 1, n
            do i = 1, n
                abar(i, j) = -bbar(i)*x(j)
            end do
        end do
    end subroutine linear_solve_vjp_factored

end module fortnum_linalg
