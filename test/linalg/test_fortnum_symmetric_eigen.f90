program test_fortnum_symmetric_eigen
    !! Independent oracles: matrices built as Q diag(lambda) Q^T from exact
    !! Householder reflectors with prescribed spectra (including degenerate
    !! lowest eigenvalues), plus the closed-form spectrum of the 1D Dirichlet
    !! Laplacian, 2 - 2 cos(j pi/(n+1)).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_symmetric_eigen, only: symmetric_eigen, block_lanczos_lowest
    use fortnum_krylov, only: KRYLOV_OK, KRYLOV_INVALID_ARGUMENT, &
        KRYLOV_MAX_ITERATIONS
    implicit none
    integer, parameter :: n = 40, big = 600
    real(dp) :: a(n, n), values(n), vectors(n, n), lambda(n), u(n), q(n, n)
    real(dp) :: lap(n, n), exact(n), ub(big), lam_big(big)
    real(dp) :: start(big, 4), lv(6), lx(big, 6), lr(6)
    real(dp) :: pi
    integer :: i, j, info, iterations

    pi = acos(-1.0_dp)
    do i = 1, n
        lambda(i) = real(i, dp)**1.5_dp - 3.0_dp
        u(i) = sin(0.37_dp*i) + 0.1_dp*i
    end do
    lambda(2) = lambda(1)
    u = u/norm2(u)
    q = -2.0_dp*spread(u, 2, n)*spread(u, 1, n)
    do i = 1, n
        q(i, i) = q(i, i) + 1.0_dp
    end do
    a = matmul(q, matmul(diag(lambda), transpose(q)))
    call symmetric_eigen(a, values, vectors, info)
    call require(info == KRYLOV_OK, 'dense decomposition succeeds')
    call require(maxval(abs(values - sorted(lambda))) < 2.0e-12_dp*maxval(abs(lambda)), &
                 'prescribed spectrum recovered')
    call require(maxval(abs(matmul(transpose(vectors), vectors) - identity(n))) &
                 < 1.0e-13_dp, 'orthonormal eigenvectors')
    call require(maxval(abs(matmul(a, vectors) - vectors*spread(values, 1, n))) &
                 < 1.0e-11_dp, 'eigen residuals')

    lap = 0.0_dp
    do i = 1, n
        lap(i, i) = 2.0_dp
        if (i < n) then
            lap(i, i + 1) = -1.0_dp
            lap(i + 1, i) = -1.0_dp
        end if
        exact(i) = 2.0_dp - 2.0_dp*cos(i*pi/(n + 1))
    end do
    call symmetric_eigen(lap, values, vectors, info)
    call require(info == KRYLOV_OK .and. maxval(abs(values - exact)) < 1.0e-13_dp, &
                 'Dirichlet Laplacian closed-form spectrum')
    call symmetric_eigen(reshape([1.0_dp], [1, 1]), values(1:1), vectors(1:1, 1:1), info)
    call require(info == KRYLOV_OK .and. abs(values(1) - 1.0_dp) < 1.0e-15_dp, &
                 'one-by-one matrix')

    ! Matrix-free operator Q diag Q^T with a threefold degenerate lowest level.
    do i = 1, big
        lam_big(i) = 1.0_dp + 0.01_dp*i
        ub(i) = cos(0.11_dp*i)
    end do
    lam_big(1:3) = 0.5_dp
    lam_big(4) = 0.75_dp
    lam_big(5) = 0.9_dp
    lam_big(6) = 0.95_dp
    ub = ub/norm2(ub)
    do j = 1, 4
        do i = 1, big
            start(i, j) = sin(0.3_dp*i*j) + 0.01_dp*j
        end do
    end do
    call block_lanczos_lowest(apply_big, start, 6, 240, 1.0e-10_dp, lv, lx, lr, &
                              info, iterations)
    call require(info == KRYLOV_OK, 'block Lanczos converges')
    call require(maxval(abs(lv - [0.5_dp, 0.5_dp, 0.5_dp, 0.75_dp, 0.9_dp, &
                                  0.95_dp])) < 1.0e-9_dp, &
                 'degenerate lowest eigenvalues with multiplicity')
    call require(all(lr <= 1.0e-10_dp*maxval(lam_big)), 'reported residual criterion')
    call require(maxval(abs(matmul(transpose(lx), lx) - identity(6))) < 1.0e-10_dp, &
                 'orthonormal Ritz vectors')
    call require(iterations <= 240, 'operator application count')
    call block_lanczos_lowest(apply_big, start(:, 1:1), 6, 12, 1.0e-12_dp, lv, &
                              lx, lr, info, iterations)
    call require(info == KRYLOV_MAX_ITERATIONS, 'truncated basis reports nonconvergence')
    call block_lanczos_lowest(apply_big, start, 0, 12, 1.0e-12_dp, lv(1:0), &
                              lx(:, 1:0), lr(1:0), info, iterations)
    call require(info == KRYLOV_INVALID_ARGUMENT, 'invalid request rejected')
    print '(a)', 'PASS: dense symmetric eigen and block Lanczos oracles'

contains

    subroutine apply_big(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: z(big)
        z = x - 2.0_dp*dot_product(ub, x)*ub
        z = lam_big*z
        y = z - 2.0_dp*dot_product(ub, z)*ub
    end subroutine apply_big

    function diag(d) result(m)
        real(dp), intent(in) :: d(:)
        real(dp) :: m(size(d), size(d))
        integer :: k
        m = 0.0_dp
        do k = 1, size(d)
            m(k, k) = d(k)
        end do
    end function diag

    function identity(m) result(e)
        integer, intent(in) :: m
        real(dp) :: e(m, m)
        e = diag([(1.0_dp, i=1, m)])
    end function identity

    function sorted(d) result(s)
        real(dp), intent(in) :: d(:)
        real(dp) :: s(size(d)), t
        integer :: k, l
        s = d
        do k = 2, size(s)
            t = s(k)
            l = k - 1
            do while (l >= 1)
                if (s(l) <= t) exit
                s(l + 1) = s(l)
                l = l - 1
            end do
            s(l + 1) = t
        end do
    end function sorted

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            write (error_unit, '(a)') 'FAIL: '//message
            error stop 1
        end if
    end subroutine require
end program test_fortnum_symmetric_eigen
