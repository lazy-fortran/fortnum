program test_fortnum_lapack_svd
    ! Singular value decomposition and minimum-norm least squares by their
    ! defining properties in quad precision: A = U S V^T with orthonormal
    ! factors and descending S for every shape and job, singular values
    ! known by construction, rank-deficient and zero matrices, and least
    ! squares solutions against the pseudo-inverse.
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer :: k
    integer, parameter :: shapes(2, 8) = reshape([1, 1, 5, 3, 3, 5, 40, 25, 25, 40, 60, 60, 100, 7, 7, 100], [2, 8])

    do k = 1, size(shapes, 2)
        call gesvd_jobs(shapes(1, k), shapes(2, k))
        call gesdd_jobs(shapes(1, k), shapes(2, k))
    end do
    call known_values(30, 20)
    call known_values(20, 30)
    call rank_deficient(30, 20, 3)
    call zero_matrix(6, 4)
    do k = 1, 3
        call cyclic_differences(3, 4*2**(k - 1))
    end do
    call least_squares_full_rank(50, 12)
    call least_squares_wide(12, 50)
    call least_squares_rank_deficient(40, 25, 6)
    call least_squares_rank_deficient(25, 40, 6)
    call argument_errors()
    print '(a)', 'SVD and least squares match quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function general(m, n, seed) result(a)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: seed
        real(dp) :: a(m, n)
        integer :: i, j
        do j = 1, n
            do i = 1, m
                a(i, j) = sin(seed*i + 1.91_dp*j + 0.37_dp*i*j)*exp(0.3_dp*cos(1.3_dp*i - 0.7_dp*j))
            end do
        end do
    end function general

    ! max |X^T X - I|
    real(qp) function orth_error(x)
        real(dp), intent(in) :: x(:, :)
        real(qp) :: e(size(x, 2), size(x, 2))
        integer :: i
        e = matmul(transpose(real(x, qp)), real(x, qp))
        do i = 1, size(x, 2)
            e(i, i) = e(i, i) - 1
        end do
        orth_error = maxval(abs(e))
    end function orth_error

    subroutine check_factors(a, s, u, vt, what)
        real(dp), intent(in) :: a(:, :), s(:), u(:, :), vt(:, :)
        character(*), intent(in) :: what
        real(qp) :: r(size(a, 1), size(a, 2))
        integer :: mn, j, m, n
        m = size(a, 1)
        n = size(a, 2)
        mn = min(m, n)
        call require(all(s(1:mn) >= 0), what//': nonnegative')
        do j = 2, mn
            call require(s(j) <= s(j - 1), what//': descending')
        end do
        call require(orth_error(u(:, 1:size(u, 2))) <= 20*max(m, n)*eps, what//': U orthonormal')
        call require(orth_error(transpose(vt)) <= 20*max(m, n)*eps, what//': V orthonormal')
        r = real(a, qp)
        do j = 1, mn
            r = r - s(j)*matmul(real(u(:, j:j), qp), real(vt(j:j, :), qp))
        end do
        call require(maxval(abs(r)) <= 30*max(m, n)*eps*maxval(s(1:mn)), what//': A = U S V^T')
    end subroutine check_factors

    subroutine gesvd_jobs(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), b(m, n), s(min(m, n)), s2(min(m, n)), u(m, m), vt(n, n)
        real(dp) :: us(m, min(m, n)), vts(min(m, n), n), d(1, 1), e(1, 1), q(1)
        real(dp), allocatable :: work(:)
        integer :: info, mn
        mn = min(m, n)
        a = general(m, n, 0.37_dp + 0.01_dp*m)
        b = a
        call fortnum_dgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, q, -1, info)
        call require(info == 0, 'dgesvd query')
        allocate (work(int(q(1))))
        call fortnum_dgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, work, size(work), info)
        call require(info == 0, 'dgesvd A A info')
        call require(orth_error(u) <= 20*max(m, n)*eps, 'dgesvd A: full U orthogonal')
        call require(orth_error(transpose(vt)) <= 20*max(m, n)*eps, 'dgesvd A: full V orthogonal')
        call check_factors(a, s, u(:, 1:mn), vt(1:mn, :), 'dgesvd A A')
        b = a
        call fortnum_dgesvd('S', 'S', m, n, b, m, s2, us, m, vts, mn, work, size(work), info)
        call require(info == 0, 'dgesvd S S info')
        call check_factors(a, s2, us, vts, 'dgesvd S S')
        b = a
        call fortnum_dgesvd('O', 'S', m, n, b, m, s2, d, 1, vts, mn, work, size(work), info)
        call require(info == 0, 'dgesvd O S info')
        call check_factors(a, s2, b(:, 1:mn), vts, 'dgesvd O S')
        b = a
        call fortnum_dgesvd('S', 'O', m, n, b, m, s2, us, m, d, 1, work, size(work), info)
        call require(info == 0, 'dgesvd S O info')
        call check_factors(a, s2, us, b(1:mn, :), 'dgesvd S O')
        b = a
        call fortnum_dgesvd('N', 'N', m, n, b, m, s2, d, 1, e, 1, work, size(work), info)
        call require(info == 0, 'dgesvd N N info')
        call require(maxval(abs(s2 - s)) <= 20*max(m, n)*eps*s(1), 'dgesvd values alone')
    end subroutine gesvd_jobs

    subroutine gesdd_jobs(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), b(m, n), s(min(m, n)), u(m, m), vt(n, n), d(1, 1), e(1, 1), q(1)
        real(dp), allocatable :: work(:)
        integer :: info, mn, iwork(8*min(m, n))
        mn = min(m, n)
        a = general(m, n, 0.83_dp)
        b = a
        call fortnum_dgesdd('O', m, n, b, m, s, u, m, vt, n, q, -1, iwork, info)
        call require(info == 0, 'dgesdd query')
        allocate (work(int(q(1)) + 4*mn*mn + 6*mn + max(m, n)))
        call fortnum_dgesdd('A', m, n, b, m, s, u, m, vt, n, work, size(work), iwork, info)
        call require(info == 0, 'dgesdd A info')
        call check_factors(a, s, u(:, 1:mn), vt(1:mn, :), 'dgesdd A')
        b = a
        call fortnum_dgesdd('S', m, n, b, m, s, u, m, vt, n, work, size(work), iwork, info)
        call require(info == 0, 'dgesdd S info')
        call check_factors(a, s, u(:, 1:mn), vt(1:mn, :), 'dgesdd S')
        b = a
        call fortnum_dgesdd('O', m, n, b, m, s, u, m, vt, n, work, size(work), iwork, info)
        call require(info == 0, 'dgesdd O info')
        if (m >= n) then
            call check_factors(a, s, b(:, 1:n), vt(1:n, :), 'dgesdd O tall')
        else
            call check_factors(a, s, u(:, 1:m), b(1:m, :), 'dgesdd O wide')
        end if
        b = a
        call fortnum_dgesdd('N', m, n, b, m, s, d, 1, e, 1, work, size(work), iwork, info)
        call require(info == 0, 'dgesdd N info')
    end subroutine gesdd_jobs

    ! A = Q1 diag(sigma) Q2^T with reflector products Q1, Q2 in quad
    ! precision; sigma graded over twelve orders of magnitude
    subroutine known_values(m, n)
        integer, intent(in) :: m, n
        real(qp) :: q1(m, m), q2(n, n), sig(min(m, n)), x(m, n)
        real(dp) :: a(m, n), s(min(m, n)), u(m, m), vt(n, n), work(10*(m + n))
        integer :: i, info
        q1 = reflectors(m, 0.3_qp)
        q2 = reflectors(n, 0.7_qp)
        do i = 1, min(m, n)
            sig(i) = 10.0_qp**(-12*real(i - 1, qp)/(min(m, n) - 1))
        end do
        x = 0
        do i = 1, min(m, n)
            x(i, i) = sig(i)
        end do
        a = real(matmul(matmul(q1, x), transpose(q2)), dp)
        call fortnum_dgesvd('A', 'A', m, n, a, m, s, u, m, vt, n, work, size(work), info)
        call require(info == 0, 'known values info')
        call require(maxval(abs(s - sig)) <= 10*max(m, n)*eps, 'singular values known by construction')
    end subroutine known_values

    function reflectors(n, seed) result(q)
        integer, intent(in) :: n
        real(qp), intent(in) :: seed
        real(qp) :: q(n, n), v(n)
        integer :: i, k
        q = 0
        do i = 1, n
            q(i, i) = 1
        end do
        do k = 1, 3
            do i = 1, n
                v(i) = cos(seed*k*i + 0.3_qp*k)
            end do
            q = q - 2*spread(v, 2, n)*spread(matmul(v, q), 1, n)/dot_product(v, v)
        end do
    end function reflectors

    subroutine rank_deficient(m, n, r)
        integer, intent(in) :: m, n, r
        real(dp) :: a(m, n), b(m, n), s(min(m, n)), u(m, m), vt(n, n), work(10*(m + n))
        integer :: info
        a = matmul(general(m, r, 0.21_dp), transpose(general(n, r, 0.57_dp)))
        b = a
        call fortnum_dgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, work, size(work), info)
        call require(info == 0, 'rank deficient info')
        call require(all(s(r + 1:) <= 50*max(m, n)*eps*s(1)), 'rank deficient: trailing values vanish')
        call require(s(r) > 1.0e-3_dp*s(1), 'rank deficient: leading values kept')
        call check_factors(a, s, u(:, 1:min(m, n)), vt(1:min(m, n), :), 'rank deficient')
    end subroutine rank_deficient

    ! Periodic first differences of several interleaved cycles: exact zero
    ! singular values (one per cycle) and pairwise repeated ones, exactly
    ! representable entries apart from the scales
    subroutine cyclic_differences(cycles, planes)
        integer, intent(in) :: cycles, planes
        real(dp) :: a(cycles*planes, cycles*planes), b(cycles*planes, cycles*planes)
        real(dp) :: s(cycles*planes), u(cycles*planes, cycles*planes), vt(cycles*planes, cycles*planes)
        real(dp), allocatable :: work(:)
        real(dp) :: q(1), c
        integer :: n, i, info, iwork(8*cycles*planes)
        n = cycles*planes
        a = 0
        do i = 1, n
            c = 3.5_dp + 0.15_dp*mod(i - 1, cycles)
            a(i, i) = -c
            a(i, mod(i - 1 + cycles, n) + 1) = c
        end do
        b = a
        call fortnum_dgesvd('S', 'A', n, n, b, n, s, u, n, vt, n, q, -1, info)
        allocate (work(int(q(1))))
        call fortnum_dgesvd('S', 'A', n, n, b, n, s, u, n, vt, n, work, size(work), info)
        call require(info == 0, 'dgesvd cyclic info')
        call check_factors(a, s, u, vt, 'dgesvd cyclic')
        call require(count(s <= 8*n*eps*s(1)) == cycles, 'dgesvd cyclic nullity')
        deallocate (work)
        b = a
        call fortnum_dgesdd('A', n, n, b, n, s, u, n, vt, n, q, -1, iwork, info)
        allocate (work(int(q(1))))
        call fortnum_dgesdd('A', n, n, b, n, s, u, n, vt, n, work, size(work), iwork, info)
        call require(info == 0, 'dgesdd cyclic info')
        call check_factors(a, s, u, vt, 'dgesdd cyclic')
        call require(count(s <= 8*n*eps*s(1)) == cycles, 'dgesdd cyclic nullity')
    end subroutine cyclic_differences

    subroutine zero_matrix(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), s(min(m, n)), u(m, m), vt(n, n), work(10*(m + n))
        integer :: info
        a = 0
        call fortnum_dgesvd('A', 'A', m, n, a, m, s, u, m, vt, n, work, size(work), info)
        call require(info == 0 .and. all(s == 0), 'zero matrix values')
        call require(orth_error(u) <= 10*eps .and. orth_error(transpose(vt)) <= 10*eps, &
            'zero matrix: orthonormal factors')
    end subroutine zero_matrix

    ! x of DGELSS and DGELSD against the normal equations, solved in quad
    subroutine least_squares_full_rank(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), b(m, 2), aw(m, n), bw(m, 2), s(n), work(20*(m + n) + 1000)
        real(qp) :: x(n, 2), r(m, 2)
        integer :: info, rank, iwork(20*n + 100)
        a = general(m, n, 0.61_dp)
        b = general(m, 2, 1.3_dp)
        x = solve_qp(matmul(transpose(real(a, qp)), real(a, qp)), matmul(transpose(real(a, qp)), real(b, qp)))
        aw = a
        bw = b
        call fortnum_dgelss(m, n, 2, aw, m, bw, m, s, -1.0_dp, rank, work, size(work), info)
        call require(info == 0 .and. rank == n, 'dgelss full rank')
        call require(maxval(abs(bw(1:n, :) - x)) <= 1.0e4_dp*eps*maxval(abs(x)), 'dgelss solution')
        r = real(b, qp) - matmul(real(a, qp), x)
        call require(abs(sqrt(sum(real(bw(n + 1:m, 1), qp)**2)) - sqrt(sum(r(:, 1)**2))) <= &
            1.0e3_dp*eps*sqrt(sum(real(b(:, 1), qp)**2)), 'dgelss residual rows')
        call require(orth_error(transpose(aw(1:n, 1:n))) <= 20*m*eps, 'dgelss leaves V^T in A')
        aw = a
        bw = b
        call fortnum_dgelsd(m, n, 2, aw, m, bw, m, s, -1.0_dp, rank, work, size(work), iwork, info)
        call require(info == 0 .and. rank == n, 'dgelsd full rank')
        call require(maxval(abs(bw(1:n, :) - x)) <= 1.0e4_dp*eps*maxval(abs(x)), 'dgelsd solution')
    end subroutine least_squares_full_rank

    ! underdetermined: the minimum-norm solution A^T (A A^T)^-1 b
    subroutine least_squares_wide(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), b(n, 1), aw(m, n), s(m), work(20*(m + n) + 1000)
        real(qp) :: x(n, 1), y(m, 1)
        integer :: info, rank
        a = general(m, n, 0.29_dp)
        b = 0
        b(1:m, :) = general(m, 1, 0.9_dp)
        y = solve_qp(matmul(real(a, qp), transpose(real(a, qp))), real(b(1:m, :), qp))
        x = matmul(transpose(real(a, qp)), y)
        aw = a
        call fortnum_dgelss(m, n, 1, aw, m, b, n, s, -1.0_dp, rank, work, size(work), info)
        call require(info == 0 .and. rank == m, 'dgelss wide rank')
        call require(maxval(abs(b - x)) <= 1.0e4_dp*eps*maxval(abs(x)), 'dgelss minimum-norm solution')
        call require(orth_error(transpose(aw(1:m, :))) <= 20*n*eps, 'dgelss wide leaves V^T in A')
    end subroutine least_squares_wide

    ! A = X Y^T of rank r: A^+ = Y (Y^T Y)^-1 (X^T X)^-1 X^T
    subroutine least_squares_rank_deficient(m, n, r)
        integer, intent(in) :: m, n, r
        real(dp) :: xf(m, r), yf(n, r), a(m, n), b(max(m, n), 1), bd(max(m, n), 1), aw(m, n)
        real(dp) :: s(min(m, n)), work(20*(m + n) + 2000)
        real(qp) :: x(n, 1), t(r, 1)
        integer :: info, rank, iwork(20*min(m, n) + 100)
        xf = general(m, r, 0.43_dp)
        yf = general(n, r, 0.77_dp)
        a = matmul(xf, transpose(yf))
        b = 0
        b(1:m, :) = general(m, 1, 0.5_dp)
        t = solve_qp(matmul(transpose(real(xf, qp)), real(xf, qp)), &
            matmul(transpose(real(xf, qp)), real(b(1:m, :), qp)))
        t = solve_qp(matmul(transpose(real(yf, qp)), real(yf, qp)), t)
        x = matmul(real(yf, qp), t)
        bd = b
        aw = a
        call fortnum_dgelss(m, n, 1, aw, m, b, max(m, n), s, 1.0e-10_dp, rank, work, size(work), info)
        call require(info == 0 .and. rank == r, 'dgelss numerical rank')
        call require(maxval(abs(b(1:n, :) - x)) <= 1.0e6_dp*eps*maxval(abs(x)), 'dgelss pseudo-inverse solution')
        aw = a
        call fortnum_dgelsd(m, n, 1, aw, m, bd, max(m, n), s, 1.0e-10_dp, rank, work, size(work), iwork, info)
        call require(info == 0 .and. rank == r, 'dgelsd numerical rank')
        call require(maxval(abs(bd(1:n, :) - x)) <= 1.0e6_dp*eps*maxval(abs(x)), 'dgelsd pseudo-inverse solution')
    end subroutine least_squares_rank_deficient

    ! Gaussian elimination with partial pivoting in quad precision
    function solve_qp(a, b) result(x)
        real(qp), intent(in) :: a(:, :), b(:, :)
        real(qp) :: x(size(b, 1), size(b, 2)), m(size(a, 1), size(a, 2)), f
        real(qp) :: row(size(a, 2)), rrow(size(b, 2))
        integer :: n, i, j, p
        n = size(a, 1)
        m = a
        x = b
        do j = 1, n
            p = j - 1 + maxloc(abs(m(j:, j)), 1)
            row = m(j, :)
            m(j, :) = m(p, :)
            m(p, :) = row
            rrow = x(j, :)
            x(j, :) = x(p, :)
            x(p, :) = rrow
            do i = j + 1, n
                f = m(i, j)/m(j, j)
                m(i, :) = m(i, :) - f*m(j, :)
                x(i, :) = x(i, :) - f*x(j, :)
            end do
        end do
        do j = n, 1, -1
            x(j, :) = (x(j, :) - matmul(m(j, j + 1:), x(j + 1:, :)))/m(j, j)
        end do
    end function solve_qp

    subroutine argument_errors()
        real(dp) :: a(4, 3), s(3), u(4, 4), vt(3, 3), b(4, 1), work(100)
        integer :: info, rank
        a = 1
        call fortnum_dgesvd('X', 'N', 4, 3, a, 4, s, u, 4, vt, 3, work, 100, info)
        call require(info == -1, 'dgesvd jobu')
        call fortnum_dgesvd('O', 'O', 4, 3, a, 4, s, u, 4, vt, 3, work, 100, info)
        call require(info == -2, 'dgesvd jobu and jobvt both O')
        call fortnum_dgesvd('N', 'N', 4, 3, a, 3, s, u, 4, vt, 3, work, 100, info)
        call require(info == -6, 'dgesvd lda')
        call fortnum_dgesvd('A', 'N', 4, 3, a, 4, s, u, 3, vt, 3, work, 100, info)
        call require(info == -9, 'dgesvd ldu')
        call fortnum_dgesvd('N', 'N', 4, 3, a, 4, s, u, 4, vt, 3, work, 12, info)
        call require(info == -13, 'dgesvd lwork')
        call fortnum_dgelss(4, 3, 1, a, 4, b, 3, s, -1.0_dp, rank, work, 100, info)
        call require(info == -7, 'dgelss ldb')
        call fortnum_dgelss(4, 3, 1, a, 4, b, 4, s, -1.0_dp, rank, work, 5, info)
        call require(info == -12, 'dgelss lwork')
    end subroutine argument_errors

end program test_fortnum_lapack_svd
