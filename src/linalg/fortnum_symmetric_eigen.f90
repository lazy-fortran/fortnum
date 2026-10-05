module fortnum_symmetric_eigen
    !! Real symmetric eigenproblems: dense and matrix-free lowest eigenpairs.
    !!
    !! `symmetric_eigen` reduces a dense symmetric matrix to tridiagonal form
    !! by Householder reflections and diagonalizes it with the implicit QL
    !! method (Wilkinson shifts), accumulating the orthogonal eigenvectors
    !! (Golub and Van Loan, Matrix Computations, 4th ed., 8.3; the classical
    !! EISPACK tred2/tql2 schedule). Eigenvalues are returned in ascending
    !! order with orthonormal eigenvector columns.
    !!
    !! `block_lanczos_lowest` computes the lowest eigenpairs of a symmetric
    !! matrix-free operator by block Lanczos with full reorthogonalization and
    !! Rayleigh--Ritz on the whole Krylov basis. A block of b vectors resolves
    !! eigenvalue multiplicities up to b; single-vector Lanczos in exact
    !! arithmetic finds only one vector of a degenerate eigenspace. Residual
    !! norms ||A x - theta x||_2 are recomputed from the stored operator images,
    !! so a reported residual r implies (Weyl/Krylov--Bogoliubov) that an
    !! eigenvalue of A lies within r of theta. It does not prove that the
    !! returned values are the lowest ones; that ordering claim stays a
    !! candidate unless the caller supplies a separate gap or inertia argument.
    !! Memory: two n-by-max_basis real arrays plus O(max_basis^2) workspace.
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_matvec_t, KRYLOV_OK, KRYLOV_MAX_ITERATIONS, &
        KRYLOV_BREAKDOWN, KRYLOV_INVALID_ARGUMENT
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none
    private

    public :: symmetric_eigen, block_lanczos_lowest

contains

    subroutine symmetric_eigen(a, values, vectors, info)
        !! Eigen-decomposition a = vectors diag(values) vectors^T of a real
        !! symmetric matrix (only the lower triangle is read). info = 0 on
        !! success, KRYLOV_INVALID_ARGUMENT for bad shapes or nonfinite data,
        !! KRYLOV_MAX_ITERATIONS if QL fails to converge.
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: values(:), vectors(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: e(:)
        integer :: n, i, j

        info = KRYLOV_INVALID_ARGUMENT
        n = size(a, 1)
        if (n < 1 .or. size(a, 2) /= n .or. size(values) /= n) return
        if (size(vectors, 1) /= n .or. size(vectors, 2) /= n) return
        if (.not. all(ieee_is_finite(a))) return
        allocate (e(n))
        do j = 1, n
            do i = j, n
                vectors(i, j) = a(i, j)
                vectors(j, i) = a(i, j)
            end do
        end do
        call householder_tridiagonal(vectors, values, e)
        call implicit_ql(values, e, vectors, info)
        if (info /= KRYLOV_OK) return
        call sort_ascending(values, vectors)
    end subroutine symmetric_eigen

    subroutine householder_tridiagonal(a, d, e)
        !! Overwrite symmetric a by Q with Q^T a Q = tridiag(e(2:), d, e(2:)).
        real(dp), intent(inout) :: a(:, :)
        real(dp), intent(out) :: d(:), e(:)
        real(dp) :: f, g, h, hh, scale
        integer :: n, i, j, l

        n = size(a, 1)
        do i = n, 2, -1
            l = i - 1
            h = 0.0_dp
            if (l > 1) then
                scale = sum(abs(a(i, 1:l)))
                if (scale == 0.0_dp) then
                    e(i) = a(i, l)
                else
                    a(i, 1:l) = a(i, 1:l)/scale
                    h = sum(a(i, 1:l)**2)
                    f = a(i, l)
                    g = -sign(sqrt(h), f)
                    e(i) = scale*g
                    h = h - f*g
                    a(i, l) = f - g
                    f = 0.0_dp
                    do j = 1, l
                        a(j, i) = a(i, j)/h
                        g = dot_product(a(j, 1:j), a(i, 1:j))
                        if (j < l) g = g + dot_product(a(j + 1:l, j), a(i, j + 1:l))
                        e(j) = g/h
                        f = f + e(j)*a(i, j)
                    end do
                    hh = f/(h + h)
                    do j = 1, l
                        f = a(i, j)
                        g = e(j) - hh*f
                        e(j) = g
                        a(j, 1:j) = a(j, 1:j) - f*e(1:j) - g*a(i, 1:j)
                    end do
                end if
            else
                e(i) = a(i, l)
            end if
            d(i) = h
        end do
        d(1) = 0.0_dp
        e(1) = 0.0_dp
        do i = 1, n
            l = i - 1
            if (d(i) /= 0.0_dp .and. l >= 1) then
                do j = 1, l
                    g = dot_product(a(i, 1:l), a(1:l, j))
                    a(1:l, j) = a(1:l, j) - g*a(1:l, i)
                end do
            end if
            d(i) = a(i, i)
            a(i, i) = 1.0_dp
            if (l >= 1) then
                a(i, 1:l) = 0.0_dp
                a(1:l, i) = 0.0_dp
            end if
        end do
    end subroutine householder_tridiagonal

    subroutine implicit_ql(d, e, z, info)
        !! Diagonalize the tridiagonal (d, e(2:)) and rotate the columns of z.
        real(dp), intent(inout) :: d(:), e(:), z(:, :)
        integer, intent(out) :: info
        real(dp), allocatable :: column(:)
        real(dp) :: b, c, dd, f, g, p, r, s
        integer :: n, l, m, i, iter
        logical :: underflow

        n = size(d)
        info = KRYLOV_OK
        if (n == 1) return
        allocate (column(size(z, 1)))
        e(1:n - 1) = e(2:n)
        e(n) = 0.0_dp
        do l = 1, n
            iter = 0
            do
                do m = l, n - 1
                    dd = abs(d(m)) + abs(d(m + 1))
                    if (abs(e(m)) <= epsilon(1.0_dp)*dd) exit
                end do
                if (m == l) exit
                iter = iter + 1
                if (iter > 60) then
                    info = KRYLOV_MAX_ITERATIONS
                    return
                end if
                g = (d(l + 1) - d(l))/(2.0_dp*e(l))
                r = hypot(g, 1.0_dp)
                g = d(m) - d(l) + e(l)/(g + sign(r, g))
                s = 1.0_dp
                c = 1.0_dp
                p = 0.0_dp
                underflow = .false.
                do i = m - 1, l, -1
                    f = s*e(i)
                    b = c*e(i)
                    r = hypot(f, g)
                    e(i + 1) = r
                    if (r == 0.0_dp) then
                        d(i + 1) = d(i + 1) - p
                        e(m) = 0.0_dp
                        underflow = .true.
                        exit
                    end if
                    s = f/r
                    c = g/r
                    g = d(i + 1) - p
                    r = (d(i) - g)*s + 2.0_dp*c*b
                    p = s*r
                    d(i + 1) = g + p
                    g = c*r - b
                    column = z(:, i + 1)
                    z(:, i + 1) = s*z(:, i) + c*column
                    z(:, i) = c*z(:, i) - s*column
                end do
                if (underflow) cycle
                d(l) = d(l) - p
                e(l) = g
                e(m) = 0.0_dp
            end do
        end do
    end subroutine implicit_ql

    subroutine sort_ascending(values, vectors)
        real(dp), intent(inout) :: values(:), vectors(:, :)
        real(dp), allocatable :: column(:)
        real(dp) :: v
        integer :: i, j, k

        allocate (column(size(vectors, 1)))
        do i = 1, size(values) - 1
            k = i
            v = values(i)
            do j = i + 1, size(values)
                if (values(j) < v) then
                    k = j
                    v = values(j)
                end if
            end do
            if (k /= i) then
                values(k) = values(i)
                values(i) = v
                column = vectors(:, i)
                vectors(:, i) = vectors(:, k)
                vectors(:, k) = column
            end if
        end do
    end subroutine sort_ascending

    subroutine block_lanczos_lowest(matvec, start, nev, max_basis, tolerance, &
                                    values, vectors, residuals, info, iterations)
        !! Lowest nev eigenpairs of the symmetric operator `matvec`.
        !!
        !! start(:, 1:b) is the initial block (b >= 1); rank-deficient or zero
        !! columns are replaced by deterministic pseudo-random vectors. The
        !! basis grows by b columns per step up to max_basis (<= n) columns.
        !! Convergence: every returned residual ||A x - theta x||_2 is at most
        !! tolerance*max(1, max_i |theta_i|). If the Krylov basis becomes
        !! invariant it is extended by pseudo-random vectors; a full basis gives
        !! the exact Rayleigh--Ritz decomposition. info: KRYLOV_OK,
        !! KRYLOV_MAX_ITERATIONS (best available pairs returned),
        !! KRYLOV_INVALID_ARGUMENT, KRYLOV_BREAKDOWN (nonfinite operator data).
        !! iterations counts operator applications.
        procedure(real_matvec_t) :: matvec
        real(dp), intent(in) :: start(:, :)
        integer, intent(in) :: nev, max_basis
        real(dp), intent(in) :: tolerance
        real(dp), intent(out) :: values(:), vectors(:, :), residuals(:)
        integer, intent(out) :: info, iterations
        real(dp), allocatable :: v(:, :), w(:, :), t(:, :), theta(:), y(:, :)
        real(dp), allocatable :: block(:, :), work(:)
        real(dp) :: scale
        integer :: n, b, k, kmax, j, added, seed
        logical :: converged

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        n = size(start, 1)
        b = size(start, 2)
        if (n < 1 .or. b < 1 .or. nev < 1 .or. nev > n) return
        if (size(values) /= nev .or. size(residuals) /= nev) return
        if (size(vectors, 1) /= n .or. size(vectors, 2) /= nev) return
        if (.not. (tolerance > 0.0_dp) .or. .not. all(ieee_is_finite(start))) return
        kmax = min(n, max(max_basis, nev, b))
        allocate (v(n, kmax), w(n, kmax), t(kmax, kmax), block(n, b), work(n))
        seed = 12345
        k = 0
        block = start
        do
            call orthonormal_extension(v, k, block, kmax, added, seed)
            if (added == 0) then
                info = KRYLOV_BREAKDOWN
                if (k < nev) return
                exit
            end if
            do j = k + 1, k + added
                call matvec(v(:, j), w(:, j))
                iterations = iterations + 1
                if (.not. all(ieee_is_finite(w(:, j)))) then
                    info = KRYLOV_BREAKDOWN
                    return
                end if
            end do
            ! Projected matrix rows/columns for the new block.
            t(1:k + added, k + 1:k + added) = &
                matmul(transpose(v(:, 1:k + added)), w(:, k + 1:k + added))
            t(k + 1:k + added, 1:k) = transpose(t(1:k, k + 1:k + added))
            t(k + 1:k + added, k + 1:k + added) = 0.5_dp*( &
                t(k + 1:k + added, k + 1:k + added) + &
                transpose(t(k + 1:k + added, k + 1:k + added)))
            k = k + added
            if (k >= nev) then
                call ritz(v, w, t, k, nev, values, vectors, residuals, &
                          theta, y, work, info)
                if (info /= KRYLOV_OK) return
                scale = max(1.0_dp, maxval(abs(theta)))
                converged = all(residuals <= tolerance*scale)
                if (converged .or. k == n) then
                    info = KRYLOV_OK
                    return
                end if
                if (k >= kmax) then
                    info = KRYLOV_MAX_ITERATIONS
                    return
                end if
            else if (k >= kmax) then
                return
            end if
            ! Next block: operator images of the newest basis block.
            block(:, 1:added) = w(:, k - added + 1:k)
            if (added < b) block(:, added + 1:b) = 0.0_dp
        end do
        info = KRYLOV_MAX_ITERATIONS
    end subroutine block_lanczos_lowest

    subroutine ritz(v, w, t, k, nev, values, vectors, residuals, theta, y, work, &
                    info)
        real(dp), intent(in) :: v(:, :), w(:, :), t(:, :)
        integer, intent(in) :: k, nev
        real(dp), intent(out) :: values(:), vectors(:, :), residuals(:)
        real(dp), allocatable, intent(inout) :: theta(:), y(:, :)
        real(dp), intent(inout) :: work(:)
        integer, intent(out) :: info
        integer :: j

        if (allocated(theta)) deallocate (theta)
        if (allocated(y)) deallocate (y)
        allocate (theta(k), y(k, k))
        call symmetric_eigen(t(1:k, 1:k), theta, y, info)
        if (info /= KRYLOV_OK) return
        do j = 1, nev
            values(j) = theta(j)
            vectors(:, j) = matmul(v(:, 1:k), y(:, j))
            work = matmul(w(:, 1:k), y(:, j)) - theta(j)*vectors(:, j)
            residuals(j) = norm2(work)
        end do
    end subroutine ritz

    subroutine orthonormal_extension(v, k, block, kmax, added, seed)
        !! Append the orthonormalized columns of block to v(:, 1:k) (two
        !! classical Gram--Schmidt passes against v, then modified
        !! Gram--Schmidt within the block). Dependent columns are replaced by
        !! pseudo-random vectors; at most kmax - k columns are added.
        real(dp), intent(inout) :: v(:, :), block(:, :)
        integer, intent(in) :: k, kmax
        integer, intent(out) :: added
        integer, intent(inout) :: seed
        real(dp) :: before, after
        integer :: j, i, pass, attempt, n

        n = size(v, 1)
        added = 0
        do j = 1, size(block, 2)
            if (k + added >= kmax) return
            do attempt = 1, 4
                before = norm2(block(:, j))
                if (before > 0.0_dp) then
                    do pass = 1, 2
                        if (k + added > 0) block(:, j) = block(:, j) - &
                            matmul(v(:, 1:k + added), &
                                   matmul(block(:, j), v(:, 1:k + added)))
                    end do
                    after = norm2(block(:, j))
                    if (after > 1.0e-10_dp*before) exit
                end if
                do i = 1, n
                    seed = int(modulo(48271_int64*int(seed, int64), 2147483647_int64))
                    block(i, j) = real(seed, dp)/2147483647.0_dp - 0.5_dp
                end do
                after = 0.0_dp
            end do
            if (.not. (after > 0.0_dp)) cycle
            added = added + 1
            v(:, k + added) = block(:, j)/after
        end do
    end subroutine orthonormal_extension
end module fortnum_symmetric_eigen
