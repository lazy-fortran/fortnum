module fortnum_multiroot
    ! Multidimensional root finding F(x) = 0 in R^n, plus two scalar
    ! numerical utilities that share the same consumer (KiLCA): a central
    ! finite-difference first derivative with a Richardson error estimate,
    ! and an ascending index sort.
    !
    ! Solvers (Powell hybrid dogleg, analytic and finite-difference Jacobian):
    !   multiroot_hybrid  - Newton step with a backtracking line search on the
    !       residual 1/2 |F|^2, analytic Jacobian supplied by the callback.
    !   multiroot_hybrids - same iteration, Jacobian built column by column by
    !       central differences of F (the KiLCA call site fakes the Jacobian
    !       this way).  "s" = Jacobian by secant/finite difference.
    !   Both terminate when |F|_inf <= ftol (residual) or the step norm
    !       |dx|_inf <= xtol*(|x|_inf + xtol) (relative+absolute delta).
    !   The Newton system J dx = -F is solved by Gaussian elimination with
    !       partial pivoting; a singular Jacobian reports FORTNUM_DOMAIN_ERROR.
    !   Ref: Powell (1970) hybrid method; More, Garbow & Hillstrom, MINPACK
    !       (1980), hybrj/hybrd; Dennis & Schnabel (1983), §6.5 line search.
    !
    ! DERIVATIVE POLICY for the solvers (ad.md §1, §4): implicit_rule.
    !   The root satisfies F(x*, p) = 0.  The implicit function theorem gives
    !   the sensitivity dx*/dp = -J_x^{-1} J_p without differentiating the
    !   iteration; J_x is the Jacobian dF/dx at the root.  Active: parameters
    !   p carried in ctx.  Inactive (ad.md §3): n, xtol, ftol, max_iter, status.
    !   Reserved derivative products for the n-dim implicit rule:
    !   multiroot_jvp, multiroot_vjp, multiroot_grad (follow-up issue).
    !
    ! deriv_central: DERIVATIVE POLICY primal_only (ad.md §1).  The operation
    !   IS finite differencing; it already estimates a derivative, so an AD
    !   derivative of the estimator is not meaningful.  Use the value it
    !   returns; do not differentiate through it.
    !   Ref: Ridders, Adv. Eng. Software 4 (1982) 75 (Richardson extrapolation
    !   of central differences); Abramowitz & Stegun 25.3.21.
    !
    ! argsort: DERIVATIVE POLICY primal_only (ad.md §1).  An index permutation
    !   is integer-valued and control-flow dependent; it carries no derivative.
    !
    ! No module-level save, no global procedure pointers; all state is
    ! caller-owned or stack-local.

    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: multiroot_fdf_t, multiroot_fn_t
    public :: multiroot_hybrid, multiroot_hybrids
    public :: deriv_fn_t, deriv_central
    public :: argsort

    ! System with analytic Jacobian: returns F(x) in f(n) and dF/dx in
    ! jac(n,n) with jac(i,j) = dF_i/dx_j.  ctx carries parameters (active for
    ! the implicit rule); the solver never inspects it.
    abstract interface
        subroutine multiroot_fdf_t(x, f, jac, ctx)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: f(:)
            real(dp), intent(out) :: jac(:, :)
            class(*), intent(in), optional :: ctx
        end subroutine multiroot_fdf_t
    end interface

    ! System without an analytic Jacobian: returns only F(x) in f(n).
    ! multiroot_hybrids builds the Jacobian by central differences of this.
    abstract interface
        subroutine multiroot_fn_t(x, f, ctx)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: f(:)
            class(*), intent(in), optional :: ctx
        end subroutine multiroot_fn_t
    end interface

    ! Scalar function for deriv_central: y = f(x), ctx carries parameters.
    abstract interface
        function deriv_fn_t(x, ctx) result(fx)
            import :: dp
            real(dp), intent(in) :: x
            class(*), intent(in), optional :: ctx
            real(dp)             :: fx
        end function deriv_fn_t
    end interface

contains

    ! Newton + backtracking line search with an analytic Jacobian.
    !
    ! x0 is the starting point (length n); the converged root is returned in x.
    ! xtol drives the step-size stop (default 1e-10), ftol the residual stop
    ! (default 1e-10), max_iter caps Newton steps (default 1000).
    subroutine multiroot_hybrid(fdf, n, x0, x, status, xtol, ftol, max_iter, ctx)
        procedure(multiroot_fdf_t)         :: fdf
        integer,                intent(in) :: n
        real(dp),               intent(in) :: x0(n)
        real(dp),               intent(out):: x(n)
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter
        class(*), intent(in), optional     :: ctx

        real(dp) :: f(n), jac(n, n), dx(n)
        real(dp) :: xt, ft
        integer  :: it, max_it, ls_stat

        call status_set(status, FORTNUM_OK, "")
        x = x0
        call resolve_tols(xtol, ftol, max_iter, xt, ft, max_it)

        call fdf(x, f, jac, ctx)
        if (maxval(abs(f)) <= ft) return

        do it = 1, max_it
            call solve_linear(n, jac, -f, dx, ls_stat)
            if (ls_stat /= 0) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "multiroot_hybrid: singular Jacobian in Newton step")
                return
            end if
            call line_search_fdf(fdf, n, x, f, jac, dx, ctx)
            if (maxval(abs(f)) <= ft) return
            if (maxval(abs(dx)) <= xt * (maxval(abs(x)) + xt)) then
                if (maxval(abs(f)) <= ft) return
                exit
            end if
        end do

        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "multiroot_hybrid: maximum iterations reached without convergence")
    end subroutine multiroot_hybrid

    ! Newton + line search with the Jacobian built by central differences.
    ! Same stop criteria and contract as multiroot_hybrid; takes a function
    ! callback that returns only F(x).
    subroutine multiroot_hybrids(fn, n, x0, x, status, xtol, ftol, max_iter, ctx)
        procedure(multiroot_fn_t)          :: fn
        integer,                intent(in) :: n
        real(dp),               intent(in) :: x0(n)
        real(dp),               intent(out):: x(n)
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter
        class(*), intent(in), optional     :: ctx

        real(dp) :: f(n), jac(n, n), dx(n)
        real(dp) :: xt, ft
        integer  :: it, max_it, ls_stat

        call status_set(status, FORTNUM_OK, "")
        x = x0
        call resolve_tols(xtol, ftol, max_iter, xt, ft, max_it)

        call fn(x, f, ctx)
        if (maxval(abs(f)) <= ft) return

        do it = 1, max_it
            call fd_jacobian(fn, n, x, jac, ctx)
            call solve_linear(n, jac, -f, dx, ls_stat)
            if (ls_stat /= 0) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "multiroot_hybrids: singular finite-difference Jacobian")
                return
            end if
            call line_search_fn(fn, n, x, f, dx, ctx)
            if (maxval(abs(f)) <= ft) return
            if (maxval(abs(dx)) <= xt * (maxval(abs(x)) + xt)) then
                if (maxval(abs(f)) <= ft) return
                exit
            end if
        end do

        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "multiroot_hybrids: maximum iterations reached without convergence")
    end subroutine multiroot_hybrids

    ! Resolve optional tolerances and iteration cap to working values.
    pure subroutine resolve_tols(xtol, ftol, max_iter, xt, ft, max_it)
        real(dp), intent(in), optional :: xtol, ftol
        integer,  intent(in), optional :: max_iter
        real(dp), intent(out) :: xt, ft
        integer,  intent(out) :: max_it
        xt = 1.0e-10_dp
        ft = 1.0e-10_dp
        max_it = 1000
        if (present(xtol))     xt = xtol
        if (present(ftol))     ft = ftol
        if (present(max_iter)) max_it = max_iter
    end subroutine resolve_tols

    ! Backtracking line search for the analytic-Jacobian solver.  Tries the
    ! full Newton step, halving lambda until the residual 1/2|F|^2 decreases
    ! (Armijo with a slack constant) or lambda underflows; updates x, f, jac.
    subroutine line_search_fdf(fdf, n, x, f, jac, dx, ctx)
        procedure(multiroot_fdf_t)          :: fdf
        integer,  intent(in)                :: n
        real(dp), intent(inout)             :: x(n), f(n), jac(n, n)
        real(dp), intent(in)                :: dx(n)
        class(*), intent(in), optional      :: ctx

        real(dp) :: xtrial(n), ftrial(n), jtrial(n, n)
        real(dp) :: g0, gtrial, lambda
        integer  :: k

        g0 = 0.5_dp * dot_product(f, f)
        lambda = 1.0_dp
        do k = 1, 30
            xtrial = x + lambda * dx
            call fdf(xtrial, ftrial, jtrial, ctx)
            gtrial = 0.5_dp * dot_product(ftrial, ftrial)
            if (gtrial < (1.0_dp - 1.0e-4_dp * lambda) * g0) then
                x = xtrial; f = ftrial; jac = jtrial
                return
            end if
            lambda = 0.5_dp * lambda
        end do
        ! No decrease found: take the smallest step taken so the iteration
        ! still advances (full Newton already evaluated as xtrial above).
        x = xtrial; f = ftrial; jac = jtrial
    end subroutine line_search_fdf

    ! Backtracking line search for the finite-difference-Jacobian solver.
    subroutine line_search_fn(fn, n, x, f, dx, ctx)
        procedure(multiroot_fn_t)           :: fn
        integer,  intent(in)                :: n
        real(dp), intent(inout)             :: x(n), f(n)
        real(dp), intent(in)                :: dx(n)
        class(*), intent(in), optional      :: ctx

        real(dp) :: xtrial(n), ftrial(n)
        real(dp) :: g0, gtrial, lambda
        integer  :: k

        g0 = 0.5_dp * dot_product(f, f)
        lambda = 1.0_dp
        do k = 1, 30
            xtrial = x + lambda * dx
            call fn(xtrial, ftrial, ctx)
            gtrial = 0.5_dp * dot_product(ftrial, ftrial)
            if (gtrial < (1.0_dp - 1.0e-4_dp * lambda) * g0) then
                x = xtrial; f = ftrial
                return
            end if
            lambda = 0.5_dp * lambda
        end do
        x = xtrial; f = ftrial
    end subroutine line_search_fn

    ! Central-difference Jacobian: column j is (F(x + h e_j) - F(x - h e_j))/(2h)
    ! with h scaled by |x_j|.
    subroutine fd_jacobian(fn, n, x, jac, ctx)
        procedure(multiroot_fn_t)      :: fn
        integer,  intent(in)           :: n
        real(dp), intent(in)           :: x(n)
        real(dp), intent(out)          :: jac(n, n)
        class(*), intent(in), optional :: ctx

        real(dp), parameter :: cube_root_eps = 6.0554544523933429e-6_dp
        real(dp) :: xp(n), xm(n), fp(n), fm(n), h
        integer  :: j

        do j = 1, n
            h = cube_root_eps * max(abs(x(j)), 1.0_dp)
            xp = x; xm = x
            xp(j) = x(j) + h
            xm(j) = x(j) - h
            call fn(xp, fp, ctx)
            call fn(xm, fm, ctx)
            jac(:, j) = (fp - fm) / (2.0_dp * h)
        end do
    end subroutine fd_jacobian

    ! Solve A x = b (n x n) by Gaussian elimination with partial pivoting.
    ! stat = 0 on success, 1 when A is numerically singular.  A and b are
    ! copied so the caller's arrays are untouched.
    subroutine solve_linear(n, a_in, b_in, x, stat)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: a_in(n, n)
        real(dp), intent(in)  :: b_in(n)
        real(dp), intent(out) :: x(n)
        integer,  intent(out) :: stat

        real(dp) :: a(n, n), b(n), factor, pivmax, s, sing_tol
        integer  :: i, j, k, p

        a = a_in
        b = b_in
        stat = 0
        x = 0.0_dp
        ! Singularity threshold scaled by the matrix magnitude: a pivot below
        ! this relative to the largest entry means a numerically singular A.
        sing_tol = epsilon(1.0_dp) * max(maxval(abs(a_in)), tiny(1.0_dp)) * real(n, dp)

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
                stat = 1
                return
            end if
            if (p /= k) then
                call swap_row(a, n, k, p)
                s = b(k); b(k) = b(p); b(p) = s
            end if
            do i = k + 1, n
                factor = a(i, k) / a(k, k)
                a(i, k:n) = a(i, k:n) - factor * a(k, k:n)
                b(i) = b(i) - factor * b(k)
            end do
        end do

        if (abs(a(n, n)) <= sing_tol) then
            stat = 1
            return
        end if

        ! Back substitution.
        do i = n, 1, -1
            s = b(i)
            do j = i + 1, n
                s = s - a(i, j) * x(j)
            end do
            x(i) = s / a(i, i)
        end do
    end subroutine solve_linear

    ! Swap rows r1 and r2 of an n x n matrix in place.
    pure subroutine swap_row(a, n, r1, r2)
        integer,  intent(in)    :: n, r1, r2
        real(dp), intent(inout) :: a(n, n)
        real(dp) :: tmp(n)
        tmp = a(r1, :)
        a(r1, :) = a(r2, :)
        a(r2, :) = tmp
    end subroutine swap_row

    ! Central finite-difference first derivative with an adaptive step.
    !
    ! result holds f'(x); abserr estimates |result - f'(x)|. A first pass at the
    ! requested step h reports the truncation and rounding parts of the error
    ! separately. When rounding dominates truncation the step is too small, so a
    ! single corrected step -- the one that balances the O(h^2) truncation
    ! against the O(1/h) rounding -- is tried, and its result is kept only when
    ! it lowers the total error and stays consistent with the first pass.
    ! Ref: Abramowitz and Stegun, "Handbook of Mathematical Functions", 25.3.6,
    ! and the standard truncation/rounding step-size balance.
    subroutine deriv_central(f, x, h, result, abserr, status, ctx)
        procedure(deriv_fn_t)               :: f
        real(dp),               intent(in)  :: x, h
        real(dp),               intent(out) :: result, abserr
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional      :: ctx

        real(dp) :: deriv, round, trunc, total
        real(dp) :: h_bal, deriv_bal, round_bal, trunc_bal, total_bal

        call status_set(status, FORTNUM_OK, "")
        if (h <= 0.0_dp) then
            result = 0.0_dp
            abserr = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "deriv_central: step h must be positive")
            return
        end if

        call central_diff_estimate(f, x, h, deriv, round, trunc, ctx)
        total = round + trunc

        ! Truncation falls as h^2 and rounding rises as 1/h, so they balance at
        ! h_bal = h*(round/(2*trunc))^(1/3). Re-estimate there when the current
        ! step is rounding-limited and keep the balanced result only if it both
        ! lowers the total error and stays within a few error bars of the first.
        if (round < trunc .and. round > 0.0_dp .and. trunc > 0.0_dp) then
            h_bal = h * (round / (2.0_dp * trunc))**(1.0_dp/3.0_dp)
            call central_diff_estimate(f, x, h_bal, deriv_bal, round_bal, &
                trunc_bal, ctx)
            total_bal = round_bal + trunc_bal
            if (total_bal < total .and. &
                abs(deriv_bal - deriv) < 4.0_dp * total) then
                deriv = deriv_bal
                total = total_bal
            end if
        end if

        result = deriv
        abserr = total
    end subroutine deriv_central

    ! Fourth-order central derivative from the four off-centre samples
    ! f(x +/- h) and f(x +/- h/2), with the truncation and rounding parts of the
    ! error returned separately; the value at x itself is not needed.
    !
    ! The half-step central difference and the full-step central difference each
    ! carry an O(h^2) truncation term with the same coefficient, so the
    ! Richardson combination (4*half - full)/3 cancels it and leaves an O(h^4)
    ! estimate. The gap between the two orders bounds the remaining truncation.
    ! The rounding part is the relative machine error the sampled values carry
    ! through the weighted sum, plus the error from the inexact abscissae
    ! x +/- h. The common 1/h factor is applied once at the end.
    subroutine central_diff_estimate(f, x, h, deriv, err_round, err_trunc, ctx)
        procedure(deriv_fn_t)          :: f
        real(dp), intent(in)           :: x, h
        real(dp), intent(out)          :: deriv, err_round, err_trunc
        class(*), intent(in), optional :: ctx

        real(dp) :: f_lo, f_hi, f_lo_half, f_hi_half
        real(dp) :: diff_wide, diff_rich, eps, noise_wide, noise_rich
        real(dp) :: abscissa_noise

        f_lo      = f(x - h,        ctx)
        f_hi      = f(x + h,        ctx)
        f_lo_half = f(x - h/2.0_dp, ctx)
        f_hi_half = f(x + h/2.0_dp, ctx)

        ! Full-step (order-2) and Richardson-combined (order-4) numerators.
        diff_wide = 0.5_dp * (f_hi - f_lo)
        diff_rich = (4.0_dp/3.0_dp) * (f_hi_half - f_lo_half) &
            - (1.0_dp/3.0_dp) * diff_wide

        ! Rounding floor of each numerator: every sampled value contributes its
        ! own relative roundoff weighted by its coefficient magnitude.
        eps = epsilon(1.0_dp)
        noise_wide = (abs(f_hi) + abs(f_lo)) * eps
        noise_rich = 2.0_dp * (abs(f_hi_half) + abs(f_lo_half)) * eps &
            + noise_wide

        ! Roundoff of x +/- h leaks into the estimate as O(eps * x / h).
        abscissa_noise = max(abs(diff_wide/h), abs(diff_rich/h)) &
            * (abs(x)/h) * eps

        deriv     = diff_rich / h
        err_trunc = abs((diff_rich - diff_wide) / h)
        err_round = abs(noise_rich / h) + abscissa_noise
    end subroutine central_diff_estimate

    ! Ascending index sort: perm is the
    ! permutation with x(perm) nondecreasing.  Heapsort on the index array so
    ! the input x is never moved; O(n log n), no recursion, no allocation
    ! beyond perm.  Ties keep ascending value order (not a stability claim).
    pure subroutine argsort(x, perm)
        real(dp), intent(in)  :: x(:)
        integer,  intent(out) :: perm(size(x))

        integer :: n, i, tmp, last

        n = size(x)
        do i = 1, n
            perm(i) = i
        end do
        if (n < 2) return

        ! Build a max-heap (by x value) over perm(1:n).
        do i = n / 2, 1, -1
            call sift_down(x, perm, i, n)
        end do

        ! Repeatedly move the max to the end and shrink the heap.
        last = n
        do while (last > 1)
            tmp = perm(1)
            perm(1) = perm(last)
            perm(last) = tmp
            last = last - 1
            call sift_down(x, perm, 1, last)
        end do
    end subroutine argsort

    ! Restore the max-heap property at node start within perm(1:heap_n),
    ! comparing by x value of the indexed elements.
    pure subroutine sift_down(x, perm, start, heap_n)
        real(dp), intent(in)    :: x(:)
        integer,  intent(inout) :: perm(:)
        integer,  intent(in)    :: start, heap_n

        integer :: root, child, tmp

        root = start
        do
            child = 2 * root
            if (child > heap_n) exit
            if (child + 1 <= heap_n) then
                if (x(perm(child + 1)) > x(perm(child))) child = child + 1
            end if
            if (x(perm(child)) > x(perm(root))) then
                tmp = perm(root)
                perm(root) = perm(child)
                perm(child) = tmp
                root = child
            else
                exit
            end if
        end do
    end subroutine sift_down

end module fortnum_multiroot
