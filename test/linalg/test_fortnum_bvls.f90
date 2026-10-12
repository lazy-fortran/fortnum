program test_fortnum_bvls
    !! Independent oracles for fortnum_bvls and fortnum_gsvd: LAPACK DGELS and
    !! DPOSV for unconstrained and Tikhonov cases, brute-force enumeration of
    !! all active sets (3^n states) for small bounded problems, KKT conditions
    !! evaluated directly from the gradient, analytic elimination of a series
    !! tie, and the GSVD identities plus singular values of A D^-1 from the
    !! independent fortnum_symmetric_eigen.
    use fortnum_kinds, only: dp
    use fortnum_bvls, only: bounded_least_squares, BVLS_OK, BVLS_MAX_ITERATIONS, &
        BVLS_INFEASIBLE, BVLS_INVALID_ARGUMENT
    use fortnum_gsvd, only: generalized_svd
    use fortnum_symmetric_eigen, only: symmetric_eigen
    use, intrinsic :: iso_fortran_env, only: error_unit, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    implicit none

    interface
        subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
            import :: dp
            character, intent(in) :: trans
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: work(*)
            integer, intent(out) :: info
        end subroutine dgels
        subroutine dposv(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            integer, intent(out) :: info
        end subroutine dposv
    end interface

    integer(int64) :: seed = 12345_int64
    real(dp) :: inf

    inf = ieee_value(1.0_dp, ieee_positive_inf)
    call test_unconstrained_vs_dgels()
    call test_known_active_sets_bruteforce()
    call test_kkt_random()
    call test_equality_series_tie()
    call test_equality_general()
    call test_tikhonov()
    call test_status_codes()
    call test_ill_conditioned_flux_gauge()
    call test_gsvd_identities()
    call test_gsvd_rank_deficient_wide()
    write (*, '(a)') 'test_fortnum_bvls: all tests passed'

contains

    real(dp) function urand()
        ! Linear congruential generator; wraparound keeps the high bits usable.
        seed = seed*6364136223846793005_int64 + 1442695040888963407_int64
        urand = real(iand(shiftr(seed, 20), 1048575_int64), dp)/1048576.0_dp
    end function urand

    real(dp) function nrand()
        integer :: i0
        nrand = -6.0_dp
        do i0 = 1, 12
            nrand = nrand + urand()
        end do
    end function nrand

    subroutine fill_random(a)
        real(dp), intent(out) :: a(:, :)
        integer :: i, j
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                a(i, j) = nrand()
            end do
        end do
    end subroutine fill_random

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            write (error_unit, '(a)') 'FAIL: '//message
            error stop 1
        end if
    end subroutine require

    function ls_dgels(a, b) result(x)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp) :: x(size(a, 2))
        real(dp), allocatable :: aw(:, :), bw(:, :), work(:)
        integer :: m, n, info
        m = size(a, 1)
        n = size(a, 2)
        if (n == 0) return
        aw = a
        allocate (bw(max(m, n), 1), work(64*(m + n)))
        bw = 0.0_dp
        bw(1:m, 1) = b
        call dgels('N', m, n, 1, aw, m, bw, max(m, n), work, size(work), info)
        x = bw(1:n, 1)
    end function ls_dgels

    subroutine brute_force(a, b, lo, up, xbest, fbest)
        !! Enumerate free / lower / upper for every variable (finite bounds
        !! only where used); LS on the free set; keep the best feasible one.
        real(dp), intent(in) :: a(:, :), b(:), lo(:), up(:)
        real(dp), intent(out) :: xbest(:), fbest
        integer :: n, code, ncode, j, t, nf
        integer, allocatable :: state(:), idx(:)
        real(dp), allocatable :: x(:), r(:), y(:)
        real(dp) :: f
        logical :: ok

        n = size(a, 2)
        ncode = 3**n
        fbest = huge(1.0_dp)
        xbest = 0.0_dp
        allocate (state(n), x(n))
        do code = 0, ncode - 1
            t = code
            ok = .true.
            do j = 1, n
                state(j) = mod(t, 3)
                t = t/3
                if (state(j) == 1) then
                    if (lo(j) <= -huge(1.0_dp)) ok = .false.
                    x(j) = lo(j)
                else if (state(j) == 2) then
                    if (up(j) >= huge(1.0_dp)) ok = .false.
                    x(j) = up(j)
                else
                    x(j) = 0.0_dp
                end if
            end do
            if (.not. ok) cycle
            idx = pack([(j, j=1, n)], state == 0)
            nf = size(idx)
            r = b - matmul(a, x)
            if (nf > 0) then
                y = ls_dgels(a(:, idx), r)
                x(idx) = y
            end if
            if (any(x < lo - 1.0e-12_dp) .or. any(x > up + 1.0e-12_dp)) cycle
            f = sum((matmul(a, x) - b)**2)
            if (f < fbest) then
                fbest = f
                xbest = x
            end if
        end do
    end subroutine brute_force

    subroutine check_kkt(a, b, lo, up, x, act, tol, label, extra)
        !! Primal feasibility and the sign conditions of the gradient.
        real(dp), intent(in) :: a(:, :), b(:), lo(:), up(:), x(:), tol
        integer, intent(in) :: act(:)
        character(*), intent(in) :: label
        real(dp), intent(in), optional :: extra(:)
        real(dp) :: g(size(x)), gs
        integer :: j
        g = matmul(transpose(a), matmul(a, x) - b)
        if (present(extra)) g = g + extra
        gs = max(1.0_dp, maxval(abs(g)), norm2(a)*norm2(b))
        call require(all(x >= lo) .and. all(x <= up), label//': bounds respected')
        do j = 1, size(x)
            if (act(j) == 0) then
                call require(abs(g(j)) <= tol*gs, label//': free gradient zero')
                call require(x(j) >= lo(j) .and. x(j) <= up(j), &
                    label//': free in range')
            else if (lo(j) == up(j)) then
                call require(x(j) == lo(j), label//': fixed variable at its value')
            else if (act(j) < 0) then
                call require(x(j) == lo(j), label//': lower active at bound')
                call require(g(j) >= -tol*gs, label//': lower multiplier sign')
            else
                call require(x(j) == up(j), label//': upper active at bound')
                call require(g(j) <= tol*gs, label//': upper multiplier sign')
            end if
        end do
    end subroutine check_kkt

    subroutine test_unconstrained_vs_dgels()
        integer, parameter :: m = 20, n = 6
        real(dp) :: a(m, n), b(m), lo(n), up(n), x(n), xr(n), rn
        integer :: info, act(n), its
        call fill_random(a)
        b = [(nrand(), its=1, m)]
        lo = -inf
        up = inf
        call bounded_least_squares(a, b, lo, up, x, info, act, rn, its)
        xr = ls_dgels(a, b)
        call require(info == BVLS_OK, 'unconstrained: status ok')
        call require(maxval(abs(x - xr)) < 1.0e-11_dp, 'unconstrained equals DGELS')
        call require(all(act == 0), 'unconstrained: empty active set')
        call require(abs(rn - norm2(matmul(a, xr) - b)) < 1.0e-11_dp, &
            'unconstrained: residual norm')
    end subroutine test_unconstrained_vs_dgels

    subroutine test_known_active_sets_bruteforce()
        integer, parameter :: m = 9
        integer :: trial, n, j, info, act(6), nactive_seen, ninf
        real(dp), allocatable :: a(:, :)
        real(dp) :: b(m), lo(6), up(6), x(6), xb(6), fb, f
        nactive_seen = 0
        do trial = 1, 60
            n = 2 + mod(trial, 5)
            allocate (a(m, n))
            call fill_random(a)
            b = [(2.0_dp*nrand(), j=1, m)]
            ninf = 0
            do j = 1, n
                lo(j) = -0.5_dp - urand()
                up(j) = 0.3_dp + urand()
                if (mod(trial + j, 7) == 0) lo(j) = -inf
                if (mod(trial + 2*j, 11) == 0) up(j) = inf
            end do
            call bounded_least_squares(a, b, lo(1:n), up(1:n), x(1:n), info, &
                act(1:n))
            call require(info == BVLS_OK, 'bruteforce: status ok')
            call brute_force(a, b, lo(1:n), up(1:n), xb(1:n), fb)
            f = sum((matmul(a, x(1:n)) - b)**2)
            call require(abs(f - fb) <= 1.0e-10_dp*(1.0_dp + fb), &
                'bruteforce: optimal objective')
            call require(maxval(abs(x(1:n) - xb(1:n))) < 1.0e-7_dp, &
                'bruteforce: minimiser (strictly convex)')
            nactive_seen = nactive_seen + count(act(1:n) /= 0)
            deallocate (a)
        end do
        call require(nactive_seen > 40, 'bruteforce: bounds were actually active')
    end subroutine test_known_active_sets_bruteforce

    subroutine test_kkt_random()
        integer, parameter :: m = 50, n = 30
        real(dp) :: a(m, n), b(m), lo(n), up(n), x(n)
        integer :: trial, info, act(n), j, nact
        nact = 0
        do trial = 1, 20
            call fill_random(a)
            b = [(3.0_dp*nrand(), j=1, m)]
            do j = 1, n
                lo(j) = -0.3_dp*urand() - 0.05_dp
                up(j) = 0.3_dp*urand() + 0.05_dp
                if (mod(j + trial, 5) == 0) lo(j) = -inf
                if (mod(j + trial, 6) == 0) up(j) = inf
                if (mod(j + trial, 13) == 0) up(j) = lo(j)
            end do
            call bounded_least_squares(a, b, lo, up, x, info, act)
            call require(info == BVLS_OK, 'kkt: status ok')
            call check_kkt(a, b, lo, up, x, act, 1.0e-8_dp, 'kkt')
            nact = nact + count(act /= 0)
        end do
        call require(nact > 100, 'kkt: many active bounds')
    end subroutine test_kkt_random

    subroutine test_equality_series_tie()
        !! x(1) = x(2) as a tied pair; oracle: eliminate x(2) analytically
        !! (merge columns, intersect bounds) and brute-force the reduced problem.
        integer, parameter :: m = 10, n = 5
        real(dp) :: a(m, n), b(m), lo(n), up(n), x(n), c(1, n), d(1)
        real(dp) :: ar(m, n - 1), lr(n - 1), ur(n - 1), xr(n - 1), fr, f
        integer :: trial, info, act(n), j
        do trial = 1, 25
            call fill_random(a)
            b = [(2.0_dp*nrand(), j=1, m)]
            do j = 1, n
                lo(j) = -0.4_dp - 0.3_dp*urand()
                up(j) = 0.4_dp + 0.3_dp*urand()
            end do
            lo(1) = -0.2_dp
            up(2) = 0.15_dp
            c = 0.0_dp
            c(1, 1) = 1.0_dp
            c(1, 2) = -1.0_dp
            d = 0.0_dp
            call bounded_least_squares(a, b, lo, up, x, info, act, c=c, d=d)
            call require(info == BVLS_OK, 'tie: status ok')
            call require(abs(x(1) - x(2)) <= 1.0e-12_dp, 'tie: equality satisfied')
            ar(:, 1) = a(:, 1) + a(:, 2)
            ar(:, 2:n - 1) = a(:, 3:n)
            lr = [max(lo(1), lo(2)), lo(3:n)]
            ur = [min(up(1), up(2)), up(3:n)]
            call brute_force(ar, b, lr, ur, xr, fr)
            f = sum((matmul(a, x) - b)**2)
            call require(abs(f - fr) <= 1.0e-10_dp*(1.0_dp + fr), &
                'tie: objective equals analytic elimination')
            call require(maxval(abs(x(3:n) - xr(2:n - 1))) < 1.0e-7_dp .and. &
                abs(x(1) - xr(1)) < 1.0e-7_dp, 'tie: minimiser')
        end do
    end subroutine test_equality_series_tie

    subroutine test_equality_general()
        !! Unbounded, two equality rows: oracle is the KKT linear system.
        integer, parameter :: m = 15, n = 6, nc = 2
        real(dp) :: a(m, n), b(m), c(nc, n), d(nc), lo(n), up(n), x(n)
        real(dp) :: kkt(n + nc, n + nc), rhs(n + nc, 1)
        integer :: info, j
        call fill_random(a)
        call fill_random(c)
        b = [(nrand(), j=1, m)]
        d = [1.0_dp, -0.5_dp]
        lo = -inf
        up = inf
        call bounded_least_squares(a, b, lo, up, x, info, c=c, d=d)
        call require(info == BVLS_OK, 'equality: status ok')
        kkt = 0.0_dp
        kkt(1:n, 1:n) = matmul(transpose(a), a)
        kkt(1:n, n + 1:) = transpose(c)
        kkt(n + 1:, 1:n) = c
        rhs(1:n, 1) = matmul(transpose(a), b)
        rhs(n + 1:, 1) = d
        rhs(:, 1) = ls_dgels(kkt, rhs(:, 1))
        info = 0
        call require(info == 0, 'equality: KKT system solvable')
        call require(maxval(abs(x - rhs(1:n, 1))) < 1.0e-10_dp, &
            'equality: equals KKT-system solution')
        call require(maxval(abs(matmul(c, x) - d)) < 1.0e-12_dp, &
            'equality: constraints satisfied')
    end subroutine test_equality_general

    subroutine test_tikhonov()
        integer, parameter :: m = 12, n = 5
        real(dp) :: a(m, n), b(m), lo(n), up(n), x(n), dsc(n), xref(n), lam
        real(dp) :: nm(n, n), rhs(n, 1)
        integer :: info, j
        call fill_random(a)
        b = [(nrand(), j=1, m)]
        dsc = [1.0_dp, 2.0_dp, 0.5_dp, 4.0_dp, 1.5_dp]
        xref = [0.3_dp, -0.2_dp, 0.0_dp, 0.1_dp, 0.4_dp]
        lam = 0.7_dp
        lo = -inf
        up = inf
        call bounded_least_squares(a, b, lo, up, x, info, lambda=lam, dscale=dsc, &
            x_ref=xref)
        call require(info == BVLS_OK, 'tikhonov: status ok')
        nm = matmul(transpose(a), a)
        rhs(:, 1) = matmul(transpose(a), b)
        do j = 1, n
            nm(j, j) = nm(j, j) + lam**2*dsc(j)**2
            rhs(j, 1) = rhs(j, 1) + lam**2*dsc(j)**2*xref(j)
        end do
        call dposv('U', n, 1, nm, n, rhs, n, info)
        call require(info == 0, 'tikhonov: normal equations')
        call require(maxval(abs(x - rhs(:, 1))) < 1.0e-11_dp, &
            'tikhonov: equals normal-equation solution')
        ! Bounded: KKT of the regularized objective.
        lo = -0.1_dp
        up = 0.25_dp
        block
            integer :: act(n)
            real(dp) :: extra(n)
            call bounded_least_squares(a, b, lo, up, x, info, act, lambda=lam, &
                dscale=dsc, x_ref=xref)
            extra = lam**2*dsc**2*(x - xref)
            call require(info == BVLS_OK, 'tikhonov bounded: status ok')
            call check_kkt(a, b, lo, up, x, act, 1.0e-8_dp, 'tikhonov', extra)
        end block
    end subroutine test_tikhonov

    subroutine test_status_codes()
        real(dp) :: a(4, 2), b(4), lo(2), up(2), x(2), c(1, 2), d(1)
        integer :: info, its
        a = reshape([1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, -1.0_dp], &
            [4, 2])
        b = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        lo = [1.0_dp, -1.0_dp]
        up = [0.0_dp, 1.0_dp]
        call bounded_least_squares(a, b, lo, up, x, info)
        call require(info == BVLS_INFEASIBLE, 'status: lower > upper is infeasible')
        up = [2.0_dp, 1.0_dp]
        c(1, :) = [1.0_dp, 1.0_dp]
        d = [10.0_dp]
        call bounded_least_squares(a, b, lo, up, x, info, c=c, d=d)
        call require(info == BVLS_INFEASIBLE, 'status: equality outside box')
        d = [2.0_dp]
        call bounded_least_squares(a, b, lo, up, x, info, c=c, d=d)
        call require(info == BVLS_OK .and. abs(x(1) + x(2) - 2.0_dp) < 1.0e-12_dp, &
            'status: feasible equality')
        call bounded_least_squares(a, b(1:3), lo, up, x, info)
        call require(info == BVLS_INVALID_ARGUMENT, 'status: shape mismatch')
        b(2) = ieee_value(1.0_dp, ieee_positive_inf)
        call bounded_least_squares(a, b, lo, up, x, info)
        call require(info == BVLS_INVALID_ARGUMENT, 'status: nonfinite data')
        ! Iteration limit: a problem that needs several active-set changes.
        block
            integer, parameter :: m = 30, n = 12
            real(dp) :: aa(m, n), bb(m), l2(n), u2(n), x2(n)
            integer :: j
            call fill_random(aa)
            bb = [(5.0_dp*nrand(), j=1, m)]
            l2 = -0.05_dp
            u2 = 0.05_dp
            call bounded_least_squares(aa, bb, l2, u2, x2, info, iterations=its, &
                max_iterations=2)
            call require(info == BVLS_MAX_ITERATIONS .and. its == 2, &
                'status: iteration limit reported')
            call require(all(x2 >= l2) .and. all(x2 <= u2), &
                'status: iterate stays feasible at the limit')
        end block
    end subroutine test_status_codes

    subroutine test_ill_conditioned_flux_gauge()
        !! Five nearly collinear columns (a common flux-like vector plus a
        !! 1e-7 individual part), three ordinary ones. The unbounded optimum
        !! needs huge cancelling coefficients, so bounds bind.
        integer, parameter :: m = 40, n = 8
        real(dp) :: a(m, n), b(m), lo(n), up(n), x(n), xb(n), fb, f, base(m)
        real(dp) :: nm(n, n), rhs(n, 1), dsc(n), xt(n), lam
        integer :: info, act(n), j
        base = [(1.0_dp + 0.1_dp*sin(0.3_dp*j), j=1, m)]
        do j = 1, 5
            a(:, j) = base + 1.0e-7_dp*[(nrand(), info=1, m)]
        end do
        do j = 6, n
            a(:, j) = [(nrand(), info=1, m)]
        end do
        xt = [20.0_dp, -35.0_dp, 10.0_dp, 40.0_dp, -30.0_dp, 0.3_dp, -0.2_dp, 0.5_dp]
        b = matmul(a, xt) + 1.0e-4_dp*[(nrand(), info=1, m)]
        lo = -1.0_dp
        up = 1.0_dp
        call bounded_least_squares(a, b, lo, up, x, info, act)
        call require(info == BVLS_OK, 'flux gauge: status ok')
        call check_kkt(a, b, lo, up, x, act, 1.0e-7_dp, 'flux gauge')
        call brute_force(a, b, lo, up, xb, fb)
        f = sum((matmul(a, x) - b)**2)
        call require(abs(f - fb) <= 1.0e-8_dp*(1.0_dp + fb), &
            'flux gauge: objective equals enumeration')
        call require(count(act /= 0) >= 3, 'flux gauge: bounds are binding')
        ! Tikhonov selects a unique, well-posed solution of the collinear block.
        dsc = 1.0_dp
        lam = 1.0e-3_dp
        lo = -inf
        up = inf
        call bounded_least_squares(a, b, lo, up, x, info, lambda=lam, dscale=dsc)
        nm = matmul(transpose(a), a)
        rhs(:, 1) = matmul(transpose(a), b)
        do j = 1, n
            nm(j, j) = nm(j, j) + lam**2
        end do
        call dposv('U', n, 1, nm, n, rhs, n, info)
        call require(maxval(abs(x - rhs(:, 1))) < 1.0e-6_dp*maxval(abs(rhs)), &
            'flux gauge: Tikhonov equals normal equations')
    end subroutine test_ill_conditioned_flux_gauge

    subroutine test_gsvd_identities()
        integer, parameter :: m = 9, n = 5
        real(dp) :: a(m, n), dm(n, n), dsc(n), ad(m, n), ata(n, n), ev(n)
        real(dp) :: evec(n, n)
        real(dp), allocatable :: alpha(:), beta(:), gsv(:), u(:, :), v(:, :)
        real(dp), allocatable :: xm(:, :), cm(:, :), sm(:, :), sv(:)
        integer :: k, l, info, i, j
        call fill_random(a)
        dsc = [1.0_dp, 0.5_dp, 2.0_dp, 3.0_dp, 0.25_dp]
        dm = 0.0_dp
        do i = 1, n
            dm(i, i) = dsc(i)
        end do
        call generalized_svd(a, dm, alpha, beta, gsv, u, v, xm, cm, sm, k, l, info)
        call require(info == BVLS_OK, 'gsvd: status ok')
        call require(k + l == n, 'gsvd: full rank')
        call require(maxval(abs(a - matmul(u, matmul(cm, transpose(xm))))) &
            < 1.0e-12_dp, &
            'gsvd: A = U C X^T')
        call require(maxval(abs(dm - matmul(v, matmul(sm, transpose(xm))))) &
            < 1.0e-12_dp, &
            'gsvd: D = V S X^T')
        call require(maxval(abs(alpha**2 + beta**2 - 1.0_dp)) < 1.0e-13_dp, &
            'gsvd: C^2 + S^2 = I')
        call require(maxval(abs(matmul(transpose(u), u) - eye(m))) < 1.0e-13_dp, &
            'gsvd: U orthogonal')
        call require(maxval(abs(matmul(transpose(v), v) - eye(n))) < 1.0e-13_dp, &
            'gsvd: V orthogonal')
        ! Independent: gsv are singular values of A D^-1.
        do j = 1, n
            ad(:, j) = a(:, j)/dsc(j)
        end do
        ata = matmul(transpose(ad), ad)
        call symmetric_eigen(ata, ev, evec, info)
        sv = sqrt(max(ev, 0.0_dp))
        call sort_ascending(gsv)
        call require(maxval(abs(gsv - sv)) < 1.0e-11_dp*maxval(sv), &
            'gsvd: generalized values equal singular values of A D^-1')
    end subroutine test_gsvd_identities

    subroutine test_gsvd_rank_deficient_wide()
        !! m < n branch (R split across A and B) with nearly collinear columns.
        integer, parameter :: m = 3, n = 5
        real(dp) :: a(m, n), dm(n, n)
        real(dp), allocatable :: alpha(:), beta(:), gsv(:), u(:, :), v(:, :)
        real(dp), allocatable :: xm(:, :), cm(:, :), sm(:, :)
        integer :: k, l, info, i
        call fill_random(a)
        a(:, 2) = a(:, 1) + 1.0e-8_dp*a(:, 2)
        dm = 0.0_dp
        do i = 1, n
            dm(i, i) = 0.1_dp*real(i, dp)
        end do
        call generalized_svd(a, dm, alpha, beta, gsv, u, v, xm, cm, sm, k, l, info)
        call require(info == BVLS_OK, 'gsvd wide: status ok')
        call require(m - k - l < 0, 'gsvd wide: exercises m < k + l branch')
        call require(maxval(abs(a - matmul(u, matmul(cm, transpose(xm))))) &
            < 1.0e-12_dp, &
            'gsvd wide: A = U C X^T')
        call require(maxval(abs(dm - matmul(v, matmul(sm, transpose(xm))))) &
            < 1.0e-12_dp, &
            'gsvd wide: D = V S X^T')
        call require(maxval(abs(alpha**2 + beta**2 - 1.0_dp)) < 1.0e-13_dp, &
            'gsvd wide: C^2 + S^2 = I')
        call generalized_svd(a(:, 1:4), dm, alpha, beta, gsv, u, v, xm, cm, sm, k, l, &
            info)
        call require(info == BVLS_INVALID_ARGUMENT, 'gsvd: column mismatch rejected')
    end subroutine test_gsvd_rank_deficient_wide

    function eye(n) result(e)
        integer, intent(in) :: n
        real(dp) :: e(n, n)
        integer :: i
        e = 0.0_dp
        do i = 1, n
            e(i, i) = 1.0_dp
        end do
    end function eye

    subroutine sort_ascending(v)
        real(dp), intent(inout) :: v(:)
        real(dp) :: t
        integer :: i, j
        do i = 2, size(v)
            t = v(i)
            j = i - 1
            do while (j >= 1)
                if (v(j) <= t) exit
                v(j + 1) = v(j)
                j = j - 1
            end do
            v(j + 1) = t
        end do
    end subroutine sort_ascending
end program test_fortnum_bvls
