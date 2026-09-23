!> Verified linear algebra with norm-wise error bounds.
!>
!> Norms: the spectral norm is bounded above by min(Frobenius,
!> sqrt(norm_1 norm_inf)) with upward rounding, or tightly through a
!> certified largest eigenvalue of A^H A.
!>
!> Products: for real A (m x n) and B, abs(fl(A B) - A B) <= gamma_n
!> abs(A) abs(B) entrywise for any summation order (Higham, Accuracy and
!> Stability of Numerical Algorithms, 2nd ed., sec. 3.5); complex inner
!> products satisfy the same with sqrt(2) gamma_(n+2) (Higham, Problem 3.7).
!> The products are computed by explicit loops so no fast (Strassen-type)
!> BLAS enters.
!>
!> Eigenvalue lower bounds (Cholesky certificate): for symmetric A and a
!> trial shift s, the floating Cholesky factor R of B = fl(A - s I), if it
!> completes, satisfies R^T R = B + dB with abs(dB) <= gamma_(n+1)
!> abs(R^T) abs(R) (Higham, Theorem 10.3), so norm2(dB) <= gamma_(n+1)
!> sum r_ij^2 and lambda_min(A) >= s - gamma_(n+1) sum r_ij^2 - u max abs(b_ii)
!> (the last term is the rounding of the shifted diagonal). A matrix ball
!> (A, rho) lowers the bound by rho (Weyl). Hermitian matrices are certified
!> through the real symmetric embedding [[X, -Y], [Y, X]], which has the same
!> eigenvalues.
!>
!> Verified inverse: an approximate inverse R of Z with theta >= norm2(I - R Z)
!> < 1 gives norm2(Z^-1) <= norm2(R)/(1 - theta) and, for every perturbation
!> norm2(D) <= eta with q = theta + norm2(R) eta < 1,
!> norm2((Z + D)^-1 - R) <= norm2(R) q/(1 - q).
!>
!> This generalises kinetic-compression kc_verified_la and closes its
!> unanalysed roundings (the diagonal of R Z - I, the shifted diagonal of the
!> Cholesky certificate).
module fortnum_verified_linalg
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rounding, only: add_up, mul_up, div_up, &
        sub_down, sqrt_up, gamma_up, unit_roundoff
    use fortnum_interval, only: interval_t, interval, operator(+), operator(*)
    implicit none
    private

    public :: mball_t, fro_up, norm1_up, norminf_up, norm2_up, norm2_tight_up
    public :: matmul_err, verified_inverse, approx_inverse
    public :: sym_lower_bound, eig_lower_bound, eig_upper_bound
    public :: interval_matmul, interval_matvec

    !> Matrix ball {C + E : norm2(E) <= rho}.
    type :: mball_t
        complex(dp), allocatable :: c(:, :)
        real(dp) :: rho = 0.0_dp
    end type mball_t

    interface fro_up
        module procedure fro_up_r, fro_up_c
    end interface fro_up

    interface norm1_up
        module procedure norm1_up_r, norm1_up_c
    end interface norm1_up

    interface norminf_up
        module procedure norminf_up_r, norminf_up_c
    end interface norminf_up

    interface norm2_up
        module procedure norm2_up_r, norm2_up_c
    end interface norm2_up

    interface norm2_tight_up
        module procedure norm2_tight_up_r, norm2_tight_up_c
    end interface norm2_tight_up

    interface matmul_err
        module procedure matmul_err_r, matmul_err_c
    end interface matmul_err

    interface eig_lower_bound
        module procedure eig_lower_bound_r, eig_lower_bound_c
    end interface eig_lower_bound

    interface eig_upper_bound
        module procedure eig_upper_bound_r, eig_upper_bound_c
    end interface eig_upper_bound

    interface
        subroutine dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
            import :: dp
            character, intent(in) :: jobz, uplo
            integer, intent(in) :: n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: w(*), work(*)
            integer, intent(out) :: info
        end subroutine dsyev
    end interface

    real(dp), parameter :: u = unit_roundoff

contains

    !> Upper bound on abs(z) with upward rounding (no hypot).
    pure elemental function cabs_up(z) result(a)
        complex(dp), intent(in) :: z
        real(dp) :: a

        a = sqrt_up(add_up(mul_up(real(z, dp), real(z, dp)), &
            mul_up(aimag(z), aimag(z))))
    end function cabs_up

    pure function fro_up_r(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: s
        integer :: i, j

        s = 0.0_dp
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                s = add_up(s, mul_up(a(i, j), a(i, j)))
            end do
        end do
        s = sqrt_up(s)
    end function fro_up_r

    pure function fro_up_c(a) result(s)
        complex(dp), intent(in) :: a(:, :)
        real(dp) :: s
        integer :: i, j

        s = 0.0_dp
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                s = add_up(s, add_up(mul_up(real(a(i, j), dp), real(a(i, j), dp)), &
                    mul_up(aimag(a(i, j)), aimag(a(i, j)))))
            end do
        end do
        s = sqrt_up(s)
    end function fro_up_c

    pure function norm1_up_r(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: s, t
        integer :: i, j

        s = 0.0_dp
        do j = 1, size(a, 2)
            t = 0.0_dp
            do i = 1, size(a, 1)
                t = add_up(t, abs(a(i, j)))
            end do
            s = max(s, t)
        end do
    end function norm1_up_r

    pure function norm1_up_c(a) result(s)
        complex(dp), intent(in) :: a(:, :)
        real(dp) :: s, t
        integer :: i, j

        s = 0.0_dp
        do j = 1, size(a, 2)
            t = 0.0_dp
            do i = 1, size(a, 1)
                t = add_up(t, cabs_up(a(i, j)))
            end do
            s = max(s, t)
        end do
    end function norm1_up_c

    pure function norminf_up_r(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: s

        s = norm1_up_r(transpose(a))
    end function norminf_up_r

    pure function norminf_up_c(a) result(s)
        complex(dp), intent(in) :: a(:, :)
        real(dp) :: s

        s = norm1_up_c(transpose(a))
    end function norminf_up_c

    !> Upper bound on the spectral norm: min(Frobenius, sqrt(norm1 norminf)).
    pure function norm2_up_r(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: s

        s = min(fro_up_r(a), sqrt_up(mul_up(norm1_up_r(a), norminf_up_r(a))))
    end function norm2_up_r

    pure function norm2_up_c(a) result(s)
        complex(dp), intent(in) :: a(:, :)
        real(dp) :: s

        s = min(fro_up_c(a), sqrt_up(mul_up(norm1_up_c(a), norminf_up_c(a))))
    end function norm2_up_c

    !> P = fl(A B) by explicit loops and err >= norm2(fl(A B) - A B).
    pure subroutine matmul_err_r(a, b, p, err)
        real(dp), intent(in) :: a(:, :), b(:, :)
        real(dp), intent(out) :: p(size(a, 1), size(b, 2))
        real(dp), intent(out) :: err
        integer :: i, j, k, n

        n = size(a, 2)
        p = 0.0_dp
        do j = 1, size(b, 2)
            do k = 1, n
                do i = 1, size(a, 1)
                    p(i, j) = p(i, j) + a(i, k)*b(k, j)
                end do
            end do
        end do
        err = mul_up(gamma_up(n), mul_up(norm2_up_r(abs(a)), norm2_up_r(abs(b))))
        err = add_up(err, underflow_floor(size(a, 1), n, size(b, 2)))
    end subroutine matmul_err_r

    pure subroutine matmul_err_c(a, b, p, err)
        complex(dp), intent(in) :: a(:, :), b(:, :)
        complex(dp), intent(out) :: p(size(a, 1), size(b, 2))
        real(dp), intent(out) :: err
        integer :: i, j, k, n

        n = size(a, 2)
        p = (0.0_dp, 0.0_dp)
        do j = 1, size(b, 2)
            do k = 1, n
                do i = 1, size(a, 1)
                    p(i, j) = p(i, j) + a(i, k)*b(k, j)
                end do
            end do
        end do
        err = mul_up(mul_up(sqrt_up(2.0_dp), gamma_up(n + 2)), &
            mul_up(norm2_up_r(cabs_up(a)), norm2_up_r(cabs_up(b))))
        err = add_up(err, underflow_floor(size(a, 1), n, size(b, 2)))
    end subroutine matmul_err_c

    !> Absolute bound on underflow errors of an m x n by n x q product:
    !> at most 2 n roundings of size 2^-1074 per entry, Frobenius over m q
    !> entries.
    pure function underflow_floor(m, n, q) result(e)
        integer, intent(in) :: m, n, q
        real(dp) :: e

        e = mul_up(mul_up(2.0_dp*real(n, dp), sqrt_up(real(m, dp)*real(q, dp))), &
            tiny(1.0_dp)*epsilon(1.0_dp))
    end function underflow_floor

    !> Cholesky certificate for the symmetric matrix defined by the lower
    !> triangle of a: with trial shift s, ok = .true. and
    !> lambda_min(A) >= lam_lo when the floating Cholesky of fl(A - s I)
    !> completes. An optional rho lowers the bound for the matrix ball
    !> (A, rho).
    pure subroutine sym_lower_bound(a, s, lam_lo, ok, rho)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(in) :: s
        real(dp), intent(out) :: lam_lo
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: rho
        real(dp), allocatable :: l(:, :)
        real(dp) :: d, dmax, rs, delta
        integer :: n, i, j, k

        n = size(a, 1)
        lam_lo = -huge(1.0_dp)
        ok = .false.
        allocate (l(n, n))
        l = 0.0_dp
        dmax = 0.0_dp
        do j = 1, n
            do i = j, n
                l(i, j) = a(i, j)
            end do
            l(j, j) = a(j, j) - s
            dmax = max(dmax, abs(l(j, j)))
        end do
        do j = 1, n
            d = l(j, j)
            do k = 1, j - 1
                d = d - l(j, k)*l(j, k)
            end do
            if (.not. (d > 0.0_dp)) return
            l(j, j) = sqrt(d)
            do i = j + 1, n
                d = l(i, j)
                do k = 1, j - 1
                    d = d - l(i, k)*l(j, k)
                end do
                l(i, j) = d/l(j, j)
            end do
        end do
        rs = 0.0_dp
        do j = 1, n
            do i = j, n
                rs = add_up(rs, mul_up(l(i, j), l(i, j)))
            end do
        end do
        if (.not. (rs <= huge(1.0_dp))) return
        delta = add_up(mul_up(gamma_up(n + 1), rs), mul_up(2.0_dp*u, dmax))
        delta = add_up(delta, underflow_floor(n, n, n))
        if (present(rho)) delta = add_up(delta, rho)
        lam_lo = sub_down(s, delta)
        ok = .true.
    end subroutine sym_lower_bound

    !> Rigorous Gershgorin lower bound for the lower-triangle symmetric matrix.
    pure function gershgorin_lower(a) result(g)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: g, r
        integer :: i, j, n

        n = size(a, 1)
        g = huge(1.0_dp)
        do i = 1, n
            r = 0.0_dp
            do j = 1, n
                if (j == i) cycle
                if (j < i) then
                    r = add_up(r, abs(a(i, j)))
                else
                    r = add_up(r, abs(a(j, i)))
                end if
            end do
            g = min(g, sub_down(a(i, i), r))
        end do
    end function gershgorin_lower

    !> Certified lower bound on the smallest eigenvalue of the symmetric
    !> matrix defined by the lower triangle of a (optionally a matrix ball
    !> (a, rho)). A LAPACK estimate (untrusted) places trial shifts below the
    !> smallest eigenvalue with growing margins; the first completed Cholesky
    !> certificate is compared with the rigorous Gershgorin bound and the
    !> larger one is returned. ok is .false. only when no bound is finite.
    subroutine eig_lower_bound_r(a, lam_lo, ok, rho)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: lam_lo
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: rho
        real(dp), allocatable :: w(:), work(:), c(:, :)
        real(dp) :: est, scale, margin, lo, gersh
        integer :: n, info, it
        logical :: done

        n = size(a, 1)
        if (n == 0) then
            lam_lo = huge(1.0_dp)
            ok = .true.
            return
        end if
        gersh = gershgorin_lower(a)
        if (present(rho)) gersh = sub_down(gersh, rho)
        lam_lo = gersh
        allocate (c(n, n), w(n), work(max(1, 3*n)))
        c = a
        call dsyev("N", "L", n, c, n, w, work, size(work), info)
        if (info == 0) then
            est = w(1)
            scale = max(abs(w(1)), abs(w(n)), tiny(1.0_dp))
            margin = real(n + 1, dp)*epsilon(1.0_dp)*scale
            do it = 1, 40
                call sym_lower_bound(a, est - margin, lo, done, rho)
                if (done) then
                    lam_lo = max(lam_lo, lo)
                    exit
                end if
                margin = 4.0_dp*margin
            end do
        end if
        ok = lam_lo > -huge(1.0_dp)
    end subroutine eig_lower_bound_r

    !> Hermitian version via the real symmetric embedding [[X, -Y], [Y, X]]
    !> of H = X + i Y (lower triangle of h, imaginary diagonal ignored).
    subroutine eig_lower_bound_c(h, lam_lo, ok, rho)
        complex(dp), intent(in) :: h(:, :)
        real(dp), intent(out) :: lam_lo
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: rho
        real(dp), allocatable :: m(:, :)

        call hermitian_embedding(h, m)
        call eig_lower_bound_r(m, lam_lo, ok, rho)
    end subroutine eig_lower_bound_c

    pure subroutine hermitian_embedding(h, m)
        complex(dp), intent(in) :: h(:, :)
        real(dp), allocatable, intent(out) :: m(:, :)
        complex(dp) :: z
        integer :: n, i, j

        n = size(h, 1)
        allocate (m(2*n, 2*n))
        do j = 1, n
            do i = 1, n
                if (i > j) then
                    z = h(i, j)
                else if (i < j) then
                    z = conjg(h(j, i))
                else
                    z = cmplx(real(h(i, i), dp), 0.0_dp, dp)
                end if
                m(i, j) = real(z, dp)
                m(n + i, n + j) = real(z, dp)
                m(n + i, j) = aimag(z)
                m(i, n + j) = -aimag(z)
            end do
        end do
    end subroutine hermitian_embedding

    !> Certified upper bound on the largest eigenvalue: -lower(-A).
    subroutine eig_upper_bound_r(a, lam_hi, ok, rho)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: lam_hi
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: rho

        call eig_lower_bound_r(-a, lam_hi, ok, rho)
        lam_hi = -lam_hi
    end subroutine eig_upper_bound_r

    subroutine eig_upper_bound_c(h, lam_hi, ok, rho)
        complex(dp), intent(in) :: h(:, :)
        real(dp), intent(out) :: lam_hi
        logical, intent(out) :: ok
        real(dp), intent(in), optional :: rho

        call eig_lower_bound_c(-h, lam_hi, ok, rho)
        lam_hi = -lam_hi
    end subroutine eig_upper_bound_c

    !> Tight upper bound on norm2(A): lambda_max(fl(A^T A)) certified plus the
    !> product error, never above norm2_up.
    function norm2_tight_up_r(a) result(s)
        real(dp), intent(in) :: a(:, :)
        real(dp) :: s, err, lam
        real(dp), allocatable :: m(:, :)
        logical :: ok

        s = norm2_up_r(a)
        if (size(a, 2) == 0) return
        allocate (m(size(a, 2), size(a, 2)))
        call matmul_err_r(transpose(a), a, m, err)
        call eig_upper_bound_r(m, lam, ok, err)
        if (ok) s = min(s, sqrt_up(max(lam, 0.0_dp)))
    end function norm2_tight_up_r

    function norm2_tight_up_c(a) result(s)
        complex(dp), intent(in) :: a(:, :)
        real(dp) :: s, err, lam
        complex(dp), allocatable :: m(:, :)
        logical :: ok

        s = norm2_up_c(a)
        if (size(a, 2) == 0) return
        allocate (m(size(a, 2), size(a, 2)))
        call matmul_err_c(conjg(transpose(a)), a, m, err)
        call eig_upper_bound_c(m, lam, ok, err)
        if (ok) s = min(s, sqrt_up(max(lam, 0.0_dp)))
    end function norm2_tight_up_c

    !> Approximate inverse by Gauss-Jordan elimination with partial pivoting
    !> (untrusted; verified_inverse certifies it). ok is .false. for an exactly
    !> singular pivot column.
    pure subroutine approx_inverse(a, r, ok)
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: r(size(a, 1), size(a, 1))
        logical, intent(out) :: ok
        complex(dp), allocatable :: w(:, :), row(:)
        complex(dp) :: f
        integer :: n, k, ip, i

        n = size(a, 1)
        allocate (w(n, 2*n), row(2*n))
        w = (0.0_dp, 0.0_dp)
        w(:, 1:n) = a
        do i = 1, n
            w(i, n + i) = (1.0_dp, 0.0_dp)
        end do
        ok = .true.
        r = (0.0_dp, 0.0_dp)
        do k = 1, n
            ip = k - 1 + maxloc(abs(w(k:n, k)), 1)
            if (.not. (abs(w(ip, k)) > 0.0_dp)) then
                ok = .false.
                return
            end if
            if (ip /= k) then
                row = w(k, :)
                w(k, :) = w(ip, :)
                w(ip, :) = row
            end if
            w(k, :) = w(k, :)/w(k, k)
            do i = 1, n
                if (i == k) cycle
                f = w(i, k)
                w(i, :) = w(i, :) - f*w(k, :)
            end do
        end do
        r = w(:, n + 1:2*n)
    end subroutine approx_inverse

    !> Enclose (Z + D)^-1 for every norm2(D) <= eta + zerr by the matrix ball
    !> g = (R, rho), R an approximate inverse; zerr bounds the error already
    !> present in the computed Z. ok is .false. when the certificate fails
    !> (q >= 1) or Z is numerically singular.
    subroutine verified_inverse(z, zerr, eta, g, ok)
        complex(dp), intent(in) :: z(:, :)
        real(dp), intent(in) :: zerr, eta
        type(mball_t), intent(out) :: g
        logical, intent(out) :: ok
        complex(dp), allocatable :: r(:, :), rz(:, :)
        real(dp) :: theta, e1, e2, nr, q, dmax
        integer :: i, n

        n = size(z, 1)
        allocate (r(n, n), rz(n, n))
        call approx_inverse(z, r, ok)
        if (.not. ok) return
        call matmul_err_c(r, z, rz, e1)
        dmax = 0.0_dp
        do i = 1, n
            rz(i, i) = rz(i, i) - (1.0_dp, 0.0_dp)
            dmax = max(dmax, cabs_up(rz(i, i)))
        end do
        ! fl(p - 1) = (p - 1)(1 + delta) on the real part: the diagonal
        ! rounding is a diagonal matrix of norm <= 2 u max abs(fl(p_ii - 1)).
        e2 = mul_up(2.0_dp*u, dmax)
        theta = add_up(add_up(norm2_tight_up_c(rz), e1), e2)
        nr = norm2_tight_up_c(r)
        q = add_up(theta, mul_up(nr, add_up(eta, zerr)))
        ok = q < 1.0_dp
        if (.not. ok) return
        allocate (g%c(n, n))
        g%c = r
        g%rho = div_up(mul_up(nr, q), sub_down(1.0_dp, q))
    end subroutine verified_inverse

    !> Interval matrix product with outward rounding.
    pure function interval_matmul(a, b) result(c)
        type(interval_t), intent(in) :: a(:, :), b(:, :)
        type(interval_t) :: c(size(a, 1), size(b, 2))
        integer :: i, j, k

        c = interval(0.0_dp)
        do j = 1, size(b, 2)
            do k = 1, size(a, 2)
                do i = 1, size(a, 1)
                    c(i, j) = c(i, j) + a(i, k)*b(k, j)
                end do
            end do
        end do
    end function interval_matmul

    !> Interval matrix-vector product with outward rounding.
    pure function interval_matvec(a, x) result(y)
        type(interval_t), intent(in) :: a(:, :), x(:)
        type(interval_t) :: y(size(a, 1))
        integer :: i, k

        y = interval(0.0_dp)
        do k = 1, size(a, 2)
            do i = 1, size(a, 1)
                y(i) = y(i) + a(i, k)*x(k)
            end do
        end do
    end function interval_matvec
end module fortnum_verified_linalg
