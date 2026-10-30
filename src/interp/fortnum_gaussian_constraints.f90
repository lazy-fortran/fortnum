module fortnum_gaussian_constraints
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: brownian_factor, constrain_discrepancy, constraint_nullspace
    interface
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import dp
            character, intent(in)::jobu, jobvt
            integer, intent(in)::m, n, lda, ldu, ldvt, lwork
            real(dp), intent(inout)::a(lda, *)
            real(dp), intent(out)::s(*), u(ldu, *), vt(ldvt, *), work(*)
            integer, intent(out)::info
        end subroutine
    end interface
contains
    subroutine constraint_nullspace(protected, basis, status)
        !! Euclidean orthonormal complement of a fixed protected tangent space.
        real(dp), intent(in) :: protected(:, :)
        real(dp), allocatable, intent(out) :: basis(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: a(:, :), u(:, :), s(:), work(:)
        real(dp) :: vt(1, 1), query(1)
        integer :: m, n, rank, info
        status = 1
        m = size(protected, 1); n = size(protected, 2)
        if (m < 1 .or. n < 1) return
        if (.not. all(ieee_is_finite(protected))) return
        a = protected
        allocate (u(m, m), s(min(m, n)))
        call dgesvd('A', 'N', m, n, a, m, s, u, m, vt, 1, query, -1, info)
        if (info /= 0) return
        allocate (work(max(1, int(query(1)))))
        call dgesvd('A', 'N', m, n, a, m, s, u, m, vt, 1, work, size(work), info)
        if (info /= 0) return
        rank = 0
        if (s(1) > 0) rank = count(s > max(m, n)*epsilon(1._dp)*s(1))
        basis = u(:, rank + 1:)
        status = 0
    end subroutine

    subroutine brownian_factor(time, initial_sd, drift_sd, factor, status)
        real(dp), intent(in)::time(:), initial_sd, drift_sd
        real(dp), allocatable, intent(out)::factor(:, :)
        integer, intent(out)::status
        integer::i, n
        status = 1
        n = size(time)
        if (n < 1) return
        if (.not. all(ieee_is_finite(time))) return
        if (.not. all(ieee_is_finite([initial_sd, drift_sd]))) return
        if (initial_sd < 0 .or. drift_sd < 0) return
        if (any(time(2:) <= time(:n - 1))) return
        allocate (factor(n, n)); factor = 0
        factor(:, 1) = initial_sd
        do i = 2, n
            factor(i:, i) = drift_sd*sqrt(time(i) - time(i - 1))
        end do
        status = 0
    end subroutine

    subroutine constrain_discrepancy(chol, protected, basis, constrained, status)
        real(dp), intent(in)::chol(:, :), protected(:, :), basis(:, :)
        real(dp), allocatable, intent(out)::constrained(:, :)
        integer, intent(out)::status
        real(dp), allocatable::a(:, :), b(:, :), u(:, :), s(:), work(:)
        real(dp)::vt(1, 1), query(1)
        integer::m, n, k, i, j, info
        status = 1
        m = size(chol, 1); n = size(protected, 2); k = min(m, n)
        if (m < 1 .or. n < 1) return
        if (any(shape(chol) /= [m, m]) .or. size(protected, 1) /= m .or. size(basis, 1) /= m) return
        if (.not. all(ieee_is_finite(chol)) .or. .not. all(ieee_is_finite(protected))) return
        if (.not. all(ieee_is_finite(basis))) return
        if (any([(chol(i, i) <= 0, i=1, m)])) return
        a = protected; b = basis
        do i = 1, m
            do j = 1, n
                a(i, j) = (a(i, j) - dot_product(chol(i, :i - 1), a(:i - 1, j)))/chol(i, i)
            end do
            do j = 1, size(b, 2)
                b(i, j) = (b(i, j) - dot_product(chol(i, :i - 1), b(:i - 1, j)))/chol(i, i)
            end do
        end do
        allocate (u(m, k), s(k))
        call dgesvd('S', 'N', m, n, a, m, s, u, m, vt, 1, query, -1, info)
        if (info /= 0) return
        allocate (work(max(1, int(query(1)))))
        call dgesvd('S', 'N', m, n, a, m, s, u, m, vt, 1, work, size(work), info)
        if (info /= 0) return
        if (s(1) > 0) then
            do i = 1, k
                if (s(i) <= max(m, n)*epsilon(1._dp)*s(1)) cycle
                b = b - matmul(reshape(u(:, i), [m, 1]), &
                    reshape(matmul(u(:, i), b), [1, size(b, 2)]))
            end do
        end if
        constrained = matmul(chol, b)
        status = 0
    end subroutine
end module
