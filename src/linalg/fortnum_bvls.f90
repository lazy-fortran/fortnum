module fortnum_bvls
    !! Dense bounded-variable least squares with optional equality constraints
    !! and a diagonal Tikhonov term:
    !!
    !!     min_x  ||A x - b||_2^2 + lambda^2 ||D (x - x_ref)||_2^2
    !!     s.t.   lower <= x <= upper,   C x = d,
    !!
    !! with D = diag(dscale). Bounds may be infinite (ieee_value infinity).
    !!
    !! Method: primal active-set method for the convex quadratic program
    !! (Nocedal and Wright, Numerical Optimization, 2nd ed., 16.5; for C empty
    !! it is the Lawson-Hanson / Stark-Parker BVLS scheme with two-sided
    !! bounds). The working set is a set of variables held at a bound. Each
    !! iteration solves the equality-constrained least-squares problem on the
    !! free variables by null-space elimination of C (SVD, LAPACK DGESVD) and a
    !! minimum-norm least-squares solve (LAPACK DGELSD, rank decided by
    !! `rcond`), takes the longest feasible step toward it, and adds the
    !! blocking bound to the working set. When no step is blocked the bound
    !! multipliers nu = grad f - C^T mu are computed (mu by least squares on
    !! the free rows); a bound whose multiplier has the wrong sign is released.
    !! With equality constraints a feasible start is found first by solving
    !! min ||C x - d|| subject to the bounds with the same method.
    !!
    !! Status contract (`info`, same convention as fortnum_krylov: 0 is success):
    !!   BVLS_OK                 KKT conditions hold to `tol`.
    !!   BVLS_MAX_ITERATIONS     iteration limit reached; x is the last
    !!                           feasible iterate (objective non-increasing).
    !!   BVLS_INFEASIBLE         lower > upper somewhere, or no x satisfies
    !!                           C x = d within the bounds (x is the point
    !!                           of least constraint violation found).
    !!   BVLS_INVALID_ARGUMENT   bad shapes, nonfinite data, or NaN bounds.
    !! `active(j)` is -1 if x(j) sits at its lower bound (including
    !! lower == upper), +1 at its upper bound, 0 if free.
    !!
    !! Rank-deficient or nearly collinear free subproblems are resolved by the
    !! minimum-norm solution with singular values below rcond*sigma_max
    !! discarded; use lambda > 0 to select a unique solution instead.
    use fortnum_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    implicit none
    private

    integer, parameter, public :: BVLS_OK = 0
    integer, parameter, public :: BVLS_MAX_ITERATIONS = 1
    integer, parameter, public :: BVLS_INFEASIBLE = 2
    integer, parameter, public :: BVLS_INVALID_ARGUMENT = -1

    public :: bounded_least_squares

    interface
        subroutine dgelsd(m, n, nrhs, a, lda, b, ldb, s, rcond, rank, work, lwork, &
                iwork, info)
            import :: dp
            integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
            real(dp), intent(inout) :: a(lda, *), b(ldb, *)
            real(dp), intent(out) :: s(*), work(*)
            real(dp), intent(in) :: rcond
            integer, intent(out) :: rank, iwork(*), info
        end subroutine dgelsd
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, &
                lwork, info)
            import :: dp
            character, intent(in) :: jobu, jobvt
            integer, intent(in) :: m, n, lda, ldu, ldvt, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: s(*), u(ldu, *), vt(ldvt, *), work(*)
            integer, intent(out) :: info
        end subroutine dgesvd
    end interface

contains

    subroutine bounded_least_squares(a, b, lower, upper, x, info, active, &
            residual_norm, iterations, c, d, lambda, &
            dscale, x_ref, tol, max_iterations, rcond)
        !! See the module description. `a` is m-by-n, `c` is nc-by-n. `lambda`
        !! (default 0) with `dscale` (default ones) adds the Tikhonov term;
        !! `x_ref` (default zero) is its centre. `tol` (default 1e-10) scales the
        !! multiplier test; `rcond` (default 1e-12) is the relative singular
        !! value cutoff. `residual_norm` is ||A x - b||_2 without the Tikhonov
        !! rows. `max_iterations` defaults to 50 n + 100 per phase.
        real(dp), intent(in) :: a(:, :), b(:), lower(:), upper(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        integer, intent(out), optional :: active(:), iterations
        real(dp), intent(out), optional :: residual_norm
        real(dp), intent(in), optional :: c(:, :), d(:), lambda, dscale(:), x_ref(:)
        real(dp), intent(in), optional :: tol, rcond
        integer, intent(in), optional :: max_iterations
        real(dp), allocatable :: aa(:, :), bb(:), cc(:, :), dd(:)
        integer, allocatable :: act(:)
        real(dp) :: lam, tl, rc
        integer :: m, n, nc, maxit, iters, i

        info = BVLS_INVALID_ARGUMENT
        m = size(a, 1)
        n = size(a, 2)
        x = 0.0_dp
        if (present(active)) active = 0
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = 0.0_dp
        nc = 0
        if (present(c)) nc = size(c, 1)
        if (n < 1 .or. size(b) /= m .or. size(lower) /= n .or. size(upper) /= n &
            .or. size(x) /= n) return
        if (present(active)) then
            if (size(active) /= n) return
        end if
        if (present(c) .neqv. present(d)) return
        if (present(c)) then
            if (size(c, 2) /= n .or. size(d) /= nc) return
            if (.not. all(ieee_is_finite(c)) .or. .not. all(ieee_is_finite(d))) return
        end if
        lam = 0.0_dp
        if (present(lambda)) lam = lambda
        if (.not. ieee_is_finite(lam)) return
        if (lam < 0.0_dp) return
        if (present(dscale)) then
            if (size(dscale) /= n) return
            if (.not. all(ieee_is_finite(dscale))) return
        end if
        if (present(x_ref)) then
            if (size(x_ref) /= n) return
            if (.not. all(ieee_is_finite(x_ref))) return
        end if
        if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(b))) return
        if (any(ieee_is_nan(lower)) .or. any(ieee_is_nan(upper))) return
        tl = 1.0e-10_dp
        if (present(tol)) tl = tol
        rc = 1.0e-12_dp
        if (present(rcond)) rc = rcond
        maxit = 50*n + 100
        if (present(max_iterations)) maxit = max_iterations
        if (tl <= 0.0_dp .or. rc < 0.0_dp .or. maxit < 1) return

        if (any(lower > upper)) then
            info = BVLS_INFEASIBLE
            return
        end if

        ! Augmented system [A; lambda D] x ~ [b; lambda D x_ref].
        if (lam > 0.0_dp) then
            allocate (aa(m + n, n), bb(m + n))
            aa = 0.0_dp
            aa(1:m, :) = a
            bb(1:m) = b
            do i = 1, n
                aa(m + i, i) = lam
                bb(m + i) = lam
                if (present(dscale)) then
                    aa(m + i, i) = lam*dscale(i)
                    bb(m + i) = lam*dscale(i)
                end if
                if (present(x_ref)) then
                    bb(m + i) = bb(m + i)*x_ref(i)
                else
                    bb(m + i) = 0.0_dp
                end if
            end do
        else
            aa = a
            bb = b
        end if
        if (present(c)) then
            cc = c
            dd = d
        else
            allocate (cc(0, n), dd(0))
        end if
        allocate (act(n))

        call bvls_core(aa, bb, lower, upper, cc, dd, x, act, iters, info, tl, &
            maxit, rc)

        if (present(active)) active = act
        if (present(iterations)) iterations = iters
        if (present(residual_norm)) residual_norm = norm2(matmul(a, x) - b)
    end subroutine bounded_least_squares

    recursive subroutine bvls_core(aa, bb, lo, up, cc, dd, x, act, iters, info, &
            tol, maxit, rcond)
        !! Core active-set loop on the already augmented system.
        real(dp), intent(in) :: aa(:, :), bb(:), lo(:), up(:), cc(:, :), dd(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: act(:), iters, info
        real(dp), intent(in) :: tol, rcond
        integer, intent(in) :: maxit
        real(dp), allocatable :: z(:), p(:), g(:), mu(:), cmu(:), nu(:)
        real(dp), allocatable :: zero_c(:, :), zero_d(:), ctf(:, :)
        logical, allocatable :: free(:)
        real(dp) :: alpha, aj, gscale, tolg, worst, viol, feas
        integer :: n, nc, j, it, jb, jr, nf, rk, it1
        integer, allocatable :: idx(:)
        logical :: blocked

        n = size(aa, 2)
        nc = size(cc, 1)
        allocate (z(n), p(n), g(n), mu(max(nc, 1)), cmu(n), nu(n), free(n))
        iters = 0
        do j = 1, n
            x(j) = min(max(0.0_dp, lo(j)), up(j))
            act(j) = 0
            if (x(j) == lo(j)) then
                act(j) = -1
            else if (x(j) == up(j)) then
                act(j) = 1
            end if
        end do

        if (nc > 0) then
            ! Phase 1: least-violating feasible point of C x = d inside the bounds.
            allocate (zero_c(0, n), zero_d(0))
            call bvls_core(cc, dd, lo, up, zero_c, zero_d, x, act, it1, info, tol, &
                maxit, rcond)
            iters = it1
            if (info /= BVLS_OK) return
            feas = norm2(matmul(cc, x) - dd)
            if (feas > max(1.0e-8_dp, 100.0_dp*tol)*(1.0_dp + norm2(dd))) then
                info = BVLS_INFEASIBLE
                return
            end if
        end if

        gscale = norm2(aa)*norm2(bb)
        if (gscale <= 0.0_dp) gscale = 1.0_dp
        tolg = tol*gscale
        info = BVLS_MAX_ITERATIONS
        do it = 1, maxit
            iters = iters + 1
            free = act == 0
            call free_ls(aa, bb, cc, dd, x, free, z, rcond)
            p = 0.0_dp
            where (free) p = z - x
            alpha = 1.0_dp
            jb = 0
            do j = 1, n
                if (.not. free(j)) cycle
                if (p(j) < 0.0_dp .and. z(j) < lo(j)) then
                    aj = (lo(j) - x(j))/p(j)
                    if (aj < alpha) then
                        alpha = aj
                        jb = -j
                    end if
                else if (p(j) > 0.0_dp .and. z(j) > up(j)) then
                    aj = (up(j) - x(j))/p(j)
                    if (aj < alpha) then
                        alpha = aj
                        jb = j
                    end if
                end if
            end do
            blocked = jb /= 0
            if (blocked) then
                x = x + alpha*p
                j = abs(jb)
                if (jb < 0) then
                    x(j) = lo(j)
                    act(j) = -1
                else
                    x(j) = up(j)
                    act(j) = 1
                end if
                where (free) x = min(max(x, lo), up)
                cycle
            end if
            where (free) x = z
            x = min(max(x, lo), up)

            ! No blocking bound: test the bound multipliers.
            g = matmul(transpose(aa), matmul(aa, x) - bb)
            cmu = 0.0_dp
            nf = count(free)
            if (nc > 0 .and. nf > 0) then
                idx = pack([(j, j=1, n)], free)
                allocate (ctf(nf, nc))
                ctf = transpose(cc(:, idx))
                call min_norm_ls(ctf, g(idx), mu(1:nc), rcond, rk)
                deallocate (ctf)
                cmu = matmul(transpose(cc), mu(1:nc))
            end if
            nu = g - cmu
            jr = 0
            worst = tolg
            do j = 1, n
                if (act(j) == 0 .or. lo(j) == up(j)) cycle
                if (act(j) < 0) then
                    viol = -nu(j)
                else
                    viol = nu(j)
                end if
                if (viol > worst) then
                    worst = viol
                    jr = j
                end if
            end do
            if (jr == 0) then
                info = BVLS_OK
                return
            end if
            act(jr) = 0
        end do
    end subroutine bvls_core

    subroutine free_ls(aa, bb, cc, dd, x, free, z, rcond)
        !! z = argmin ||aa z - bb|| over the free variables (others fixed at x)
        !! subject to cc z = dd; non-free entries of z equal x.
        real(dp), intent(in) :: aa(:, :), bb(:), cc(:, :), dd(:), x(:), rcond
        logical, intent(in) :: free(:)
        real(dp), intent(out) :: z(:)
        real(dp), allocatable :: af(:, :), cf(:, :), r(:), dp_(:), yp(:), y(:)
        real(dp), allocatable :: u(:, :), vt(:, :), s(:), work(:), an(:, :), w(:)
        real(dp) :: wq(1)
        integer, allocatable :: idx(:)
        integer :: n, nf, nc, j, rk, rk2, nn, info

        n = size(x)
        nc = size(cc, 1)
        z = x
        nf = count(free)
        if (nf == 0) return
        idx = pack([(j, j=1, n)], free)
        af = aa(:, idx)
        r = bb - matmul(aa, x) + matmul(af, x(idx))
        if (nc == 0) then
            allocate (y(nf))
            call min_norm_ls(af, r, y, rcond, rk)
            z(idx) = y
            return
        end if

        cf = cc(:, idx)
        dp_ = dd - matmul(cc, x) + matmul(cf, x(idx))
        allocate (u(nc, nc), vt(nf, nf), s(min(nc, nf)))
        call dgesvd('A', 'A', nc, nf, cf, nc, s, u, nc, vt, nf, wq, -1, info)
        allocate (work(max(1, int(wq(1)))))
        call dgesvd('A', 'A', nc, nf, cf, nc, s, u, nc, vt, nf, work, size(work), info)
        rk = 0
        if (info == 0 .and. size(s) > 0) then
            if (s(1) > 0.0_dp) rk = count(s > 1.0e-12_dp*s(1))
        end if
        allocate (yp(nf))
        yp = 0.0_dp
        if (rk > 0) then
            yp = matmul(transpose(vt(1:rk, :)), &
                matmul(transpose(u(:, 1:rk)), dp_)/s(1:rk))
        end if
        nn = nf - rk
        if (nn == 0) then
            z(idx) = yp
            return
        end if
        allocate (an(size(af, 1), nn), w(nn))
        an = matmul(af, transpose(vt(rk + 1:nf, :)))
        call min_norm_ls(an, r - matmul(af, yp), w, rcond, rk2)
        z(idx) = yp + matmul(transpose(vt(rk + 1:nf, :)), w)
    end subroutine free_ls

    subroutine min_norm_ls(a, b, x, rcond, rank)
        !! Minimum-norm least-squares solution of a x ~ b via DGELSD.
        real(dp), intent(in) :: a(:, :), b(:), rcond
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: rank
        real(dp), allocatable :: aw(:, :), bw(:), s(:), work(:)
        integer, allocatable :: iwork(:)
        real(dp) :: wq(1)
        integer :: m, n, ld, info, iq(1)

        m = size(a, 1)
        n = size(a, 2)
        x = 0.0_dp
        rank = 0
        if (n == 0) return
        if (m == 0) return
        ld = max(m, n)
        allocate (aw(m, n), bw(ld), s(min(m, n)))
        aw = a
        bw = 0.0_dp
        bw(1:m) = b
        call dgelsd(m, n, 1, aw, m, bw, ld, s, rcond, rank, wq, -1, iq, info)
        allocate (work(max(1, int(wq(1)))), iwork(max(1, iq(1))))
        call dgelsd(m, n, 1, aw, m, bw, ld, s, rcond, rank, work, size(work), &
            iwork, info)
        if (info == 0) x = bw(1:n)
    end subroutine min_norm_ls

end module fortnum_bvls
