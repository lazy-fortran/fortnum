program test_fortnum_lapack_gsvd
    ! The generalized SVD by its defining equations in quad precision:
    ! U^T A Q = D1 (0 R), V^T B Q = D2 (0 R) with U, V, Q orthogonal, R upper
    ! triangular and alpha^2 + beta^2 = 1, for tall, wide and rank deficient
    ! pairs; generalized singular values known by construction; the 2 x 2
    ! building blocks by their defining identities; argument errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_gsvd, only: fortnum_dggsvd3, fortnum_dlags2, fortnum_dlasv2
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: eps = epsilon(1.0_dp)

    call svd2x2()
    call triangular_pairs()
    call pair(9, 7, 5, 0.31_dp, 'tall')
    call pair(3, 5, 5, 0.77_dp, 'wide A')
    call pair(8, 2, 6, 1.13_dp, 'low rank B')
    call pair(4, 3, 6, 0.59_dp, 'both short')
    call pair(1, 1, 1, 0.41_dp, 'scalar')
    call deficient()
    call known_values(7)
    call argument_errors()
    print '(a)', 'generalized SVD matches quad precision'

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

    function fro(a) result(s)
        real(qp), intent(in) :: a(:, :)
        real(qp) :: s
        s = sqrt(sum(a**2))
    end function fro

    function eye(n) result(e)
        integer, intent(in) :: n
        real(qp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function eye

    ! the contract of DGGSVD3 for the result of a call on (a0, b0)
    subroutine check_gsvd(m, p, n, a0, b0, k, l, a, b, alpha, beta, u, v, q, iwork, what)
        integer, intent(in) :: m, p, n, k, l, iwork(:)
        real(dp), intent(in) :: a0(m, n), b0(p, n), a(m, n), b(p, n), alpha(n), beta(n)
        real(dp), intent(in) :: u(m, m), v(p, p), q(n, n)
        character(*), intent(in) :: what
        real(qp) :: r(k + l, k + l), c(m, k + l), s(p, k + l), lhs(m, n), rhs(m, n)
        real(qp) :: lhsb(p, n), rhsb(p, n), tol
        real(dp) :: sorted(n), t
        integer :: nr, i, j
        nr = k + l
        call require(nr <= n .and. l <= p, what//': ranks')
        r = 0
        do j = 1, nr
            do i = 1, min(j, m)
                r(i, j) = a(i, n - nr + j)
            end do
            do i = m + 1, j
                r(i, j) = b(i - k, n - nr + j)
            end do
        end do
        c = 0
        s = 0
        do i = 1, nr
            if (i <= m) c(i, i) = alpha(i)
            if (i > k) s(i - k, i) = beta(i)
        end do
        tol = 50*max(m, n, p)*eps
        lhs = matmul(matmul(transpose(real(u, qp)), real(a0, qp)), real(q, qp))
        rhs = 0
        rhs(:, n - nr + 1:n) = matmul(c, r)
        call require(fro(lhs - rhs) <= tol*max(fro(real(a0, qp)), 1.0_qp), what//': U^T A Q')
        lhsb = matmul(matmul(transpose(real(v, qp)), real(b0, qp)), real(q, qp))
        rhsb = 0
        rhsb(:, n - nr + 1:n) = matmul(s, r)
        call require(fro(lhsb - rhsb) <= tol*max(fro(real(b0, qp)), 1.0_qp), what//': V^T B Q')
        call require(fro(matmul(transpose(real(u, qp)), real(u, qp)) - eye(m)) <= tol, &
            what//': U orthogonal')
        call require(fro(matmul(transpose(real(v, qp)), real(v, qp)) - eye(p)) <= tol, &
            what//': V orthogonal')
        call require(fro(matmul(transpose(real(q, qp)), real(q, qp)) - eye(n)) <= tol, &
            what//': Q orthogonal')
        call require(all(abs(alpha(1:nr)**2 + beta(1:nr)**2 - 1) <= 4*eps), &
            what//': alpha^2 + beta^2 = 1')
        call require(all(alpha(1:k) == 1) .and. all(beta(1:k) == 0), what//': first k')
        call require(all(alpha(nr + 1:n) == 0) .and. all(beta(nr + 1:n) == 0), what//': rest')
        call require(all(alpha(1:nr) >= 0) .and. all(beta(1:nr) >= 0), what//': signs')
        do j = 1, nr
            call require(abs(r(j, j)) > 0, what//': R nonsingular')
        end do
        ! iwork sorts alpha(k+1:k+min(l,m-k)) into decreasing order
        sorted = alpha
        do i = k + 1, k + min(l, m - k)
            j = iwork(i)
            t = sorted(i)
            sorted(i) = sorted(j)
            sorted(j) = t
        end do
        do i = k + 1, k + min(l, m - k) - 1
            call require(sorted(i) >= sorted(i + 1), what//': sort order')
        end do
    end subroutine check_gsvd

    subroutine run(m, p, n, a0, b0, k, l, alpha, beta, what)
        integer, intent(in) :: m, p, n
        real(dp), intent(in) :: a0(m, n), b0(p, n)
        integer, intent(out) :: k, l
        real(dp), intent(out) :: alpha(n), beta(n)
        character(*), intent(in) :: what
        real(dp) :: a(m, n), b(p, n), u(m, m), v(p, p), q(n, n), wq(1)
        real(dp), allocatable :: work(:)
        integer :: iwork(n), info
        a = a0
        b = b0
        call fortnum_dggsvd3('U', 'V', 'Q', m, n, p, k, l, a, m, b, p, alpha, beta, u, m, &
            v, p, q, n, wq, -1, iwork, info)
        call require(info == 0, what//': workspace query')
        allocate (work(int(wq(1))))
        call fortnum_dggsvd3('U', 'V', 'Q', m, n, p, k, l, a, m, b, p, alpha, beta, u, m, &
            v, p, q, n, work, size(work), iwork, info)
        call require(info == 0, what//': info')
        call check_gsvd(m, p, n, a0, b0, k, l, a, b, alpha, beta, u, v, q, iwork, what)
    end subroutine run

    subroutine pair(m, p, n, seed, what)
        integer, intent(in) :: m, p, n
        real(dp), intent(in) :: seed
        character(*), intent(in) :: what
        real(dp) :: alpha(n), beta(n)
        integer :: k, l
        call run(m, p, n, general(m, n, seed), general(p, n, seed + 0.83_dp), k, l, &
            alpha, beta, what)
        call require(k + l == min(n, m + p), what//': generic rank')
        call require(l == min(p, n), what//': rank of B')
    end subroutine pair

    ! rank deficient pairs: B = 0, A = 0, a repeated column, and B of rank 2
    subroutine deficient()
        integer, parameter :: m = 6, p = 5, n = 4
        real(dp) :: a(m, n), b(p, n), alpha(n), beta(n)
        integer :: k, l
        a = general(m, n, 0.21_dp)
        b = 0
        call run(m, p, n, a, b, k, l, alpha, beta, 'B = 0')
        call require(k == n .and. l == 0, 'B = 0: k = rank A')
        a = 0
        b = general(p, n, 0.33_dp)
        call run(m, p, n, a, b, k, l, alpha, beta, 'A = 0')
        call require(k == 0 .and. l == n .and. all(alpha(1:n) == 0), 'A = 0: alpha = 0')
        a = general(m, n, 0.45_dp)
        a(:, 3) = a(:, 1)
        b = general(p, n, 0.52_dp)
        b(:, 3) = b(:, 1)
        call run(m, p, n, a, b, k, l, alpha, beta, 'repeated column')
        call require(k + l == n - 1, 'repeated column: rank n - 1')
        b = matmul(general(p, 2, 0.11_dp), general(2, n, 0.71_dp))
        a = general(m, n, 0.63_dp)
        call run(m, p, n, a, b, k, l, alpha, beta, 'rank 2 B')
        call require(l == 2 .and. k == n - 2, 'rank 2 B: k and l')
        ! rank 1 A and rank 2 B: k = 1 < n - l
        a = matmul(general(m, 1, 0.27_dp), general(1, n, 0.93_dp))
        call run(m, p, n, a, b, k, l, alpha, beta, 'rank 1 A, rank 2 B')
        call require(l == 2 .and. k == 1, 'rank 1 A, rank 2 B: k and l')
    end subroutine deficient

    ! A = U0 diag(c) X, B = V0 diag(s) X with X nonsingular: the generalized
    ! singular values alpha/beta are c/s
    subroutine known_values(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), b(n, n), alpha(n), beta(n), x(n, n), cs(n), sn(n), gsv(n), ref(n)
        real(qp) :: u0(n, n), v0(n, n), dc(n, n), ds(n, n)
        integer :: k, l, i
        x = general(n, n, 0.9_dp)
        do i = 1, n
            x(i, i) = x(i, i) + 4
        end do
        u0 = householder(n, 0.3_qp)
        v0 = householder(n, 1.7_qp)
        dc = 0
        ds = 0
        do i = 1, n
            cs(i) = cos(0.2_dp*i)
            sn(i) = sin(0.2_dp*i)
            dc(i, i) = cs(i)
            ds(i, i) = sn(i)
        end do
        a = real(matmul(matmul(u0, dc), real(x, qp)), dp)
        b = real(matmul(matmul(v0, ds), real(x, qp)), dp)
        call run(n, n, n, a, b, k, l, alpha, beta, 'known values')
        call require(k == 0 .and. l == n, 'known values: k and l')
        gsv = alpha/beta
        ref = cs/sn
        call sort(gsv)
        call sort(ref)
        call require(all(abs(gsv - ref) <= 1.0e-12_dp*ref), 'known values: alpha/beta')
    end subroutine known_values

    subroutine sort(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: t
        integer :: i, j
        do i = 2, size(x)
            t = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= t) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = t
        end do
    end subroutine sort

    function householder(n, seed) result(h)
        integer, intent(in) :: n
        real(qp), intent(in) :: seed
        real(qp) :: h(n, n), w(n)
        integer :: i
        do i = 1, n
            w(i) = cos(seed*i + 0.4_qp)
        end do
        h = eye(n) - 2*spread(w, 2, n)*spread(w, 1, n)/sum(w**2)
    end function householder

    ! [csl snl; -snl csl] [f g; 0 h] [csr -snr; snr csr] = diag(ssmax, ssmin)
    subroutine svd2x2()
        real(dp), parameter :: cases(3, 12) = reshape([ &
            1.0_dp, 2.0_dp, 3.0_dp, 3.0_dp, -2.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, -4.0_dp, 1.0e-20_dp, 1.0_dp, 1.0e-20_dp, 2.0_dp, 1.0e20_dp, 3.0_dp, &
            -5.0_dp, 0.5_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 1.0e-300_dp, -1.0_dp, &
            1.0_dp, -10.0_dp, 2.0_dp, 1.0e-20_dp, -1.0_dp, 3.0e-20_dp, -2.0_dp, -7.0_dp, 0.5_dp], [3, 12])
        real(dp) :: f, g, h, smin, smax, snr, csr, snl, csl
        real(qp) :: m(2, 2), lq(2, 2), rq(2, 2), d(2, 2), scale
        integer :: i
        do i = 1, size(cases, 2)
            f = cases(1, i)
            g = cases(2, i)
            h = cases(3, i)
            call fortnum_dlasv2(f, g, h, smin, smax, snr, csr, snl, csl)
            m = reshape([real(f, qp), 0.0_qp, real(g, qp), real(h, qp)], [2, 2])
            lq = reshape([real(csl, qp), -real(snl, qp), real(snl, qp), real(csl, qp)], [2, 2])
            rq = reshape([real(csr, qp), real(snr, qp), -real(snr, qp), real(csr, qp)], [2, 2])
            d = matmul(matmul(lq, m), rq)
            scale = max(abs(f), abs(g), abs(h), tiny(1.0_dp))
            call require(abs(d(1, 1) - smax) <= 8*eps*scale .and. abs(d(2, 2) - smin) <= 8*eps*scale &
                .and. abs(d(1, 2)) <= 8*eps*scale .and. abs(d(2, 1)) <= 8*eps*scale, 'dlasv2: SVD')
            call require(abs(csl**2 + snl**2 - 1) <= 4*eps .and. abs(csr**2 + snr**2 - 1) <= 4*eps, &
                'dlasv2: rotations')
            call require(abs(smax) >= abs(smin), 'dlasv2: order')
        end do
    end subroutine svd2x2

    ! dlags2 zeroes the (1,2) entries of upper and the (2,1) entries of lower
    ! triangular pairs, and the other rows become parallel
    subroutine triangular_pairs()
        real(dp) :: x(6), csu, snu, csv, snv, csq, snq
        real(qp) :: a(2, 2), b(2, 2), uq(2, 2), vq(2, 2), qq(2, 2), ua(2, 2), vb(2, 2), scale
        integer :: i
        logical :: upper
        do i = 1, 24
            x = [sin(1.3_dp*i), cos(0.7_dp*i + 1), sin(2.1_dp*i + 0.3_dp), &
                cos(1.9_dp*i), sin(0.4_dp*i + 2), cos(3.1_dp*i + 0.5_dp)]
            if (mod(i, 5) == 0) x(2) = 0
            if (mod(i, 7) == 0) x(5) = 1.0e-12_dp
            upper = mod(i, 2) == 0
            call fortnum_dlags2(upper, x(1), x(2), x(3), x(4), x(5), x(6), csu, snu, csv, snv, &
                csq, snq)
            if (upper) then
                a = reshape([real(x(1), qp), 0.0_qp, real(x(2), qp), real(x(3), qp)], [2, 2])
                b = reshape([real(x(4), qp), 0.0_qp, real(x(5), qp), real(x(6), qp)], [2, 2])
            else
                a = reshape([real(x(1), qp), real(x(2), qp), 0.0_qp, real(x(3), qp)], [2, 2])
                b = reshape([real(x(4), qp), real(x(5), qp), 0.0_qp, real(x(6), qp)], [2, 2])
            end if
            uq = reshape([real(csu, qp), -real(snu, qp), real(snu, qp), real(csu, qp)], [2, 2])
            vq = reshape([real(csv, qp), -real(snv, qp), real(snv, qp), real(csv, qp)], [2, 2])
            qq = reshape([real(csq, qp), -real(snq, qp), real(snq, qp), real(csq, qp)], [2, 2])
            ua = matmul(matmul(transpose(uq), a), qq)
            vb = matmul(matmul(transpose(vq), b), qq)
            scale = 16*eps*max(maxval(abs(a)), maxval(abs(b)))
            if (upper) then
                call require(abs(ua(1, 2)) <= scale .and. abs(vb(1, 2)) <= scale, 'dlags2 upper: zeros')
                call require(abs(ua(2, 1)*vb(2, 2) - ua(2, 2)*vb(2, 1)) <= 4*scale*maxval(abs(b)) &
                    + 4*scale*maxval(abs(a)), 'dlags2 upper: parallel rows')
            else
                call require(abs(ua(2, 1)) <= scale .and. abs(vb(2, 1)) <= scale, 'dlags2 lower: zeros')
                call require(abs(ua(1, 1)*vb(1, 2) - ua(1, 2)*vb(1, 1)) <= 4*scale*maxval(abs(b)) &
                    + 4*scale*maxval(abs(a)), 'dlags2 lower: parallel rows')
            end if
            call require(abs(csu**2 + snu**2 - 1) <= 4*eps .and. abs(csv**2 + snv**2 - 1) <= 4*eps &
                .and. abs(csq**2 + snq**2 - 1) <= 4*eps, 'dlags2: rotations')
        end do
    end subroutine triangular_pairs

    subroutine argument_errors()
        real(dp) :: a(2, 2), b(2, 2), alpha(2), beta(2), u(2, 2), v(2, 2), q(2, 2), work(64)
        integer :: k, l, iwork(2), info
        a = 1
        b = 1
        call fortnum_dggsvd3('X', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -1, 'dggsvd3: jobu')
        call fortnum_dggsvd3('U', 'X', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -2, 'dggsvd3: jobv')
        call fortnum_dggsvd3('U', 'V', 'X', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -3, 'dggsvd3: jobq')
        call fortnum_dggsvd3('U', 'V', 'Q', -1, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -4, 'dggsvd3: m')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 1, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -10, 'dggsvd3: lda')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 1, alpha, beta, u, 2, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -12, 'dggsvd3: ldb')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 1, v, 2, &
            q, 2, work, 64, iwork, info)
        call require(info == -16, 'dggsvd3: ldu')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 1, &
            q, 2, work, 64, iwork, info)
        call require(info == -18, 'dggsvd3: ldv')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 1, work, 64, iwork, info)
        call require(info == -20, 'dggsvd3: ldq')
        call fortnum_dggsvd3('U', 'V', 'Q', 2, 2, 2, k, l, a, 2, b, 2, alpha, beta, u, 2, v, 2, &
            q, 2, work, 0, iwork, info)
        call require(info == -24, 'dggsvd3: lwork')
    end subroutine argument_errors

end program test_fortnum_lapack_gsvd
