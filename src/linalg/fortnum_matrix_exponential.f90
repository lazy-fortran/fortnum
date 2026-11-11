!> Dense matrix exponential by scaling and squaring with Higham's degree-13
!> diagonal Pade approximant (SIAM J. Matrix Anal. Appl. 26, 1179 (2005)).
module fortnum_matrix_exponential
    use fortnum_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: matrix_exponential

    !> matrix_exponential(a, value, info): value = exp(a) for a square real or
    !> complex matrix; info = 0 on success, nonzero for non-finite input or
    !> result or a singular Pade denominator.
    interface matrix_exponential
        module procedure matrix_exponential_real, matrix_exponential_complex
    end interface

    real(dp), parameter :: b(0:13) = [ &
        64764752532480000.0_dp, 32382376266240000.0_dp, &
        7771770303897600.0_dp, 1187353796428800.0_dp, &
        129060195264000.0_dp, 10559470521600.0_dp, &
        670442572800.0_dp, 33522128640.0_dp, 1323241920.0_dp, &
        40840800.0_dp, 960960.0_dp, 16380.0_dp, 182.0_dp, 1.0_dp]
    real(dp), parameter :: theta13 = 5.371920351148152_dp

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine
        subroutine zgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: dp
            integer, intent(in) :: n, nrhs, lda, ldb
            integer, intent(out) :: ipiv(*)
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    !> Number of squarings that brings the 1-norm below theta_13.
    pure integer function squarings_for(anorm) result(s)
        real(dp), intent(in) :: anorm
        s = 0
        if (anorm > theta13) s = ceiling(log(anorm/theta13)/log(2.0_dp))
    end function

    subroutine matrix_exponential_real(a, value, info)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: value(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: x(:, :), a2(:, :), a4(:, :), a6(:, :)
        real(dp), allocatable :: u(:, :), v(:, :), identity(:, :)
        real(dp) :: anorm
        integer, allocatable :: pivots(:)
        integer :: n, squarings, i

        n = size(a, 1)
        info = 1
        value = 0
        if (size(a, 2) /= n .or. any(shape(value) /= [n, n])) return
        if (.not. all(ieee_is_finite(a))) return
        anorm = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(anorm)) return
        squarings = squarings_for(anorm)
        allocate (x(n, n), a2(n, n), a4(n, n), a6(n, n), u(n, n), v(n, n))
        allocate (identity(n, n), pivots(n))
        x = scale(a, -squarings)
        identity = 0
        do i = 1, n
            identity(i, i) = 1
        end do
        a2 = matmul(x, x)
        a4 = matmul(a2, a2)
        a6 = matmul(a4, a2)
        u = matmul(x, matmul(a6, b(13)*a6 + b(11)*a4 + b(9)*a2) &
                   + b(7)*a6 + b(5)*a4 + b(3)*a2 + b(1)*identity)
        v = matmul(a6, b(12)*a6 + b(10)*a4 + b(8)*a2) &
            + b(6)*a6 + b(4)*a4 + b(2)*a2 + b(0)*identity
        a2 = v - u
        x = v + u
        call dgesv(n, n, a2, n, pivots, x, n, info)
        if (info /= 0) return
        do i = 1, squarings
            x = matmul(x, x)
        end do
        info = 0
        if (.not. all(ieee_is_finite(x))) info = 1
        value = x
    end subroutine

    subroutine matrix_exponential_complex(a, value, info)
        complex(dp), intent(in) :: a(:, :)
        complex(dp), intent(out) :: value(:, :)
        integer, intent(out) :: info
        complex(dp), allocatable :: x(:, :), a2(:, :), a4(:, :), a6(:, :)
        complex(dp), allocatable :: u(:, :), v(:, :), identity(:, :)
        real(dp) :: anorm
        integer, allocatable :: pivots(:)
        integer :: n, squarings, i

        n = size(a, 1)
        info = 1
        value = 0
        if (size(a, 2) /= n .or. any(shape(value) /= [n, n])) return
        if (.not. (all(ieee_is_finite(a%re)) .and. all(ieee_is_finite(a%im)))) return
        anorm = maxval(sum(abs(a), dim=1))
        if (.not. ieee_is_finite(anorm)) return
        squarings = squarings_for(anorm)
        allocate (x(n, n), a2(n, n), a4(n, n), a6(n, n), u(n, n), v(n, n))
        allocate (identity(n, n), pivots(n))
        x = a*2.0_dp**(-squarings)
        identity = 0
        do i = 1, n
            identity(i, i) = 1
        end do
        a2 = matmul(x, x)
        a4 = matmul(a2, a2)
        a6 = matmul(a4, a2)
        u = matmul(x, matmul(a6, b(13)*a6 + b(11)*a4 + b(9)*a2) &
                   + b(7)*a6 + b(5)*a4 + b(3)*a2 + b(1)*identity)
        v = matmul(a6, b(12)*a6 + b(10)*a4 + b(8)*a2) &
            + b(6)*a6 + b(4)*a4 + b(2)*a2 + b(0)*identity
        a2 = v - u
        x = v + u
        call zgesv(n, n, a2, n, pivots, x, n, info)
        if (info /= 0) return
        do i = 1, squarings
            x = matmul(x, x)
        end do
        info = 0
        if (.not. (all(ieee_is_finite(x%re)) .and. all(ieee_is_finite(x%im)))) info = 1
        value = x
    end subroutine
end module fortnum_matrix_exponential
