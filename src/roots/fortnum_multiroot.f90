module fortnum_multiroot
    ! Multidimensional root finding F(x) = 0 in R^n, plus two scalar
    ! numerical utilities that share the same consumer (KiLCA): a fourth-order
    ! central-difference first derivative with an error estimate, and an
    ! ascending index sort.
    !
    ! Solvers:
    !   multiroot_hybrid  - Newton's method with a backtracking line search on
    !       1/2 |F|_2^2, using the analytic Jacobian the callback returns with
    !       every F.  No trust region and no Jacobian updates.
    !   multiroot_hybrids - Powell's hybrid method as in MINPACK hybrd1: dogleg
    !       steps in a trust region, a forward-difference Jacobian at the start
    !       and after repeated poor steps, Broyden rank-one updates of its QR
    !       factors otherwise.  "s": no Jacobian callback.
    !   Both report success only when |F|_inf <= ftol; stagnation below xtol
    !       or the iteration cap is a failure status.
    !   Ref: Powell, in Rabinowitz (ed.), Numerical Methods for Nonlinear
    !       Algebraic Equations (1970); More, Garbow & Hillstrom, User Guide
    !       for MINPACK-1, ANL-80-74 (1980), hybrd; Dennis & Schnabel (1983),
    !       ch. 6 (line search, dogleg), sec. 8.1 (Broyden update).
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
    !   One Richardson step on central differences at h and h/2 (not Ridders'
    !   extrapolation tableau), at most one retry at the error-balancing step.
    !   Ref: Abramowitz & Stegun 25.3.6.
    !
    ! argsort: DERIVATIVE POLICY primal_only (ad.md §1).  An index permutation
    !   is integer-valued and control-flow dependent; it carries no derivative.
    !
    ! No module-level save, no global procedure pointers; all state is
    ! caller-owned or stack-local.

    use fortnum_kinds, only: dp
    use fortnum_linalg, only: LINALG_OK, lu_solve_factored
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: multiroot_fdf_t, multiroot_fn_t, multiroot_residual_jvp_t
    public :: multiroot_state_jacobian_t, multiroot_parameter_vjp_t
    public :: multiroot_preconditioned_solve_t
    public :: multiroot_hybrid, multiroot_hybrids
    public :: multiroot_implicit_jvp, multiroot_implicit_vjp
    public :: multiroot_jvp, multiroot_jvp_factored
    public :: multiroot_vjp, multiroot_vjp_factored, multiroot_grad
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
    ! multiroot_hybrids builds the Jacobian by forward differences of this.
    abstract interface
        subroutine multiroot_fn_t(x, f, ctx)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: f(:)
            class(*), intent(in), optional :: ctx
        end subroutine multiroot_fn_t
    end interface

    ! Contracted residual tangent products at a converged vector root.
    abstract interface
        subroutine multiroot_residual_jvp_t(x, p, tp, jac_x, f_p_tp, context)
            import :: dp
            real(dp), intent(in) :: x(:), p(:), tp(:)
            real(dp), intent(out) :: jac_x(size(x), size(x))
            real(dp), intent(out) :: f_p_tp(size(x))
            class(*), intent(inout), optional :: context
        end subroutine multiroot_residual_jvp_t
    end interface

    abstract interface
        subroutine multiroot_state_jacobian_t(x, p, jac_x, context)
            import :: dp
            real(dp), intent(in) :: x(:), p(:)
            real(dp), intent(out) :: jac_x(size(x), size(x))
            class(*), intent(inout), optional :: context
        end subroutine multiroot_state_jacobian_t
    end interface

    abstract interface
        subroutine multiroot_parameter_vjp_t(x, p, u, f_p_t_u, context)
            import :: dp
            real(dp), intent(in) :: x(:), p(:), u(:)
            real(dp), intent(out) :: f_p_t_u(size(p))
            class(*), intent(inout), optional :: context
        end subroutine multiroot_parameter_vjp_t
    end interface

    ! Optional implicit-product solve hook. Implementations may carry and reuse
    ! a preconditioner through context; the default remains the internal dense
    ! direct solve.
    abstract interface
        subroutine multiroot_preconditioned_solve_t(a, b, x, info, context)
            import :: dp
            real(dp), intent(in) :: a(:, :), b(:)
            real(dp), intent(out) :: x(:)
            integer, intent(out) :: info
            class(*), intent(inout), optional :: context
        end subroutine multiroot_preconditioned_solve_t
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

    ! Newton's method with a backtracking line search, analytic Jacobian.
    !
    ! Unlike multiroot_hybrids this is not Powell's dogleg method: the callback
    ! returns J with every F, so each Newton step uses the exact Jacobian at
    ! the current point and a Broyden update would save no evaluations.
    !
    ! x0 is the starting point (length n); the root is returned in x.  Success
    ! when |F|_inf <= ftol (default 1e-10).  A Newton step with
    ! |dx|_inf <= xtol*(|x|_inf + xtol) (xtol default 1e-10) before that, or
    ! max_iter Newton steps (default 1000), is reported as
    ! FORTNUM_CONVERGENCE_ERROR; a singular Jacobian as FORTNUM_DOMAIN_ERROR.
    subroutine multiroot_hybrid(fdf, n, x0, x, status, xtol, ftol, max_iter, &
            ctx, jacobian)
        procedure(multiroot_fdf_t)         :: fdf
        integer,                intent(in) :: n
        real(dp),               intent(in) :: x0(n)
        real(dp),               intent(out):: x(n)
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter
        class(*), intent(in), optional     :: ctx
        real(dp), intent(out), optional    :: jacobian(n, n)

        real(dp) :: f(n), jac(n, n), dx(n)
        real(dp) :: xt, ft
        integer  :: it, max_it, ls_stat

        call status_set(status, FORTNUM_OK, "")
        x = x0
        call resolve_tols(xtol, ftol, max_iter, xt, ft, max_it)

        call fdf(x, f, jac, ctx)
        if (maxval(abs(f)) <= ft) then
            if (present(jacobian)) jacobian = jac
            return
        end if

        do it = 1, max_it
            call solve_linear(n, jac, -f, dx, ls_stat)
            if (ls_stat /= 0) then
                if (present(jacobian)) jacobian = jac
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "multiroot_hybrid: singular Jacobian in Newton step")
                return
            end if
            call line_search_fdf(fdf, n, x, f, jac, dx, ctx)
            if (maxval(abs(f)) <= ft) then
                if (present(jacobian)) jacobian = jac
                return
            end if
            if (maxval(abs(dx)) <= xt * (maxval(abs(x)) + xt)) then
                exit
            end if
        end do

        if (present(jacobian)) jacobian = jac
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "multiroot_hybrid: maximum iterations reached without convergence")
    end subroutine multiroot_hybrid

    ! Powell hybrid method with a dogleg trust region, as MINPACK hybrd1
    ! (unit scaling, initial bound FACTOR*|x0|_2).  Same defaults and success
    ! criterion as multiroot_hybrid; takes a callback that returns only F(x).
    !
    ! Each iteration minimises the linear model |F + J p|_2 over the region
    ! |p|_2 <= delta along Powell's dogleg path and accepts the trial point
    ! when it achieves at least 1e-4 of the predicted reduction of |F|_2^2.
    ! The region shrinks after a poor ratio of actual to predicted reduction
    ! and grows after good ones.  J is built by forward differences (n
    ! evaluations) at the start and after two consecutive poor ratios; in
    ! between, a Broyden rank-one update is applied directly to J = QR by
    ! Givens rotations, so such a step costs one evaluation of F.
    !
    ! Converged: |F|_inf <= ftol.  Failure (FORTNUM_CONVERGENCE_ERROR, or
    ! FORTNUM_DOMAIN_ERROR when the current Jacobian factor is singular):
    ! the trust radius falls to xtol*(|x|_2 + xtol) after a step that cut
    ! |F|_2 by less than 10% (hybrd stops there regardless), no progress over
    ! five Jacobian evaluations or ten iterations (hybrd's nslow2/nslow1), a
    ! step at rounding level, or max_iter trial points (one F evaluation
    ! each, plus n per Jacobian).
    subroutine multiroot_hybrids(fn, n, x0, x, status, xtol, ftol, max_iter, ctx)
        procedure(multiroot_fn_t)          :: fn
        integer,                intent(in) :: n
        real(dp),               intent(in) :: x0(n)
        real(dp),               intent(out):: x(n)
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter
        class(*), intent(in), optional     :: ctx

        ! MINPACK hybrd constants: initial bound factor and ratio thresholds.
        real(dp), parameter :: FACTOR = 100.0_dp
        real(dp), parameter :: P1 = 0.1_dp, P5 = 0.5_dp, P001 = 1.0e-3_dp, &
            P0001 = 1.0e-4_dp
        real(dp) :: f(n), jac(n, n), q(n, n), r(n, n), qtf(n)
        real(dp) :: p(n), xtrial(n), ftrial(n), pred(n), qtft(n), u(n), v(n)
        real(dp) :: xt, ft, eps, delta, xnorm, fnorm, fnorm1, pnorm
        real(dp) :: actred, prered, ratio
        integer  :: it, max_it, nsucc, ncsuc, ncfail, nslow1, nslow2, j
        logical  :: jeval, sing

        call status_set(status, FORTNUM_OK, "")
        x = x0
        call resolve_tols(xtol, ftol, max_iter, xt, ft, max_it)
        eps = epsilon(1.0_dp)

        call fn(x, f, ctx)
        if (maxval(abs(f)) <= ft) return
        fnorm = norm2(f)
        xnorm = norm2(x)
        it = 0
        nsucc = 0
        ncsuc = 0
        ncfail = 0
        nslow1 = 0
        nslow2 = 0

        do
            jeval = .true.
            call fd_jacobian_forward(fn, n, x, f, jac, ctx)
            ! As in hybrd, the bound restarts until a step has been accepted.
            if (nsucc == 0) then
                delta = FACTOR * xnorm
                if (delta == 0.0_dp) delta = FACTOR
            end if
            call qr_householder(n, jac, q, r)
            sing = any([(r(j, j) == 0.0_dp, j = 1, n)])
            qtf = matmul(f, q)

            do
                it = it + 1
                if (it > max_it) then
                    call fail("maximum iterations reached without convergence")
                    return
                end if
                call dogleg(n, r, qtf, delta, p)
                p = -p
                xtrial = x + p
                pnorm = norm2(p)
                if (nsucc == 0) delta = min(delta, pnorm)
                call fn(xtrial, ftrial, ctx)
                fnorm1 = norm2(ftrial)

                actred = -1.0_dp
                if (fnorm1 < fnorm) actred = 1.0_dp - (fnorm1 / fnorm)**2
                ! Linear-model residual in the Q basis: Q^T F + R p.
                pred = qtf + matmul(r, p)
                prered = 0.0_dp
                if (norm2(pred) < fnorm) prered = 1.0_dp - (norm2(pred) / fnorm)**2
                ratio = 0.0_dp
                if (prered > 0.0_dp) ratio = actred / prered

                if (ratio < P1) then
                    ncsuc = 0
                    ncfail = ncfail + 1
                    delta = P5 * delta
                else
                    ncfail = 0
                    ncsuc = ncsuc + 1
                    if (ratio >= P5 .or. ncsuc > 1) delta = max(delta, pnorm / P5)
                    if (abs(ratio - 1.0_dp) <= P1) delta = pnorm / P5
                end if

                if (ratio >= P0001) then
                    x = xtrial
                    f = ftrial
                    xnorm = norm2(x)
                    fnorm = fnorm1
                    nsucc = nsucc + 1
                    if (maxval(abs(f)) <= ft) return
                end if

                nslow1 = nslow1 + 1
                if (actred >= P001) nslow1 = 0
                if (jeval) nslow2 = nslow2 + 1
                if (actred >= P1) nslow2 = 0

                ! hybrd stops here unconditionally; a step that still cut |F|_2
                ! by 10% is convergence towards ftol, not stagnation.
                if (delta <= xt * (xnorm + xt) .and. actred < P1) then
                    call fail("trust region below xtol without convergence")
                    return
                end if
                if (P1 * max(P1 * delta, pnorm) <= eps * xnorm) then
                    call fail("step at rounding level; xtol too small")
                    return
                end if
                if (nslow2 == 5) then
                    call fail("no progress over five Jacobian evaluations")
                    return
                end if
                if (nslow1 == 10) then
                    call fail("no progress in ten iterations")
                    return
                end if
                if (ncfail == 2) exit

                ! Broyden update J + (dF - J p) p^T / |p|^2 in factored form:
                ! R + u v^T with u = Q^T(dF - J p)/|p|, v = p/|p|.
                qtft = matmul(ftrial, q)
                u = (qtft - pred) / pnorm
                v = p / pnorm
                if (ratio >= P0001) qtf = qtft
                call qr_rank_one_update(n, q, r, qtf, u, v)
                sing = any([(r(j, j) == 0.0_dp, j = 1, n)])
                jeval = .false.
            end do
        end do

    contains

        subroutine fail(reason)
            character(len=*), intent(in) :: reason
            if (sing) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "multiroot_hybrids: singular Jacobian; " // reason)
            else
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "multiroot_hybrids: " // reason)
            end if
        end subroutine fail

    end subroutine multiroot_hybrids

    ! Implicit-rule derivatives for an n-dim root x*(p) of F(x, p) = 0.
    !
    ! Implicit function theorem: differentiating F(x*(p), p) = 0 gives
    !   J_x dx* + J_p dp = 0   =>   dx*/dp = -J_x^{-1} J_p,
    ! where J_x = dF/dx (n x n) and J_p = dF/dp (n x m) at the converged root.
    ! The iteration that found x* is inactive (ad.md §3); only the defining
    ! equation is differentiated.  The module's Gaussian-elimination solve_linear
    ! supplies J_x^{-1} actions; a singular J_x reports FORTNUM_DOMAIN_ERROR.
    !
    ! Active arguments: jac_x (J_x at root), f_p (J_p at root), and the
    ! tangent/adjoint vectors.  Inactive: status, solver internals.
    !
    ! HVP: deferred.  Second derivatives of x*(p) need d^2F second partials that
    ! the current callback interface does not provide.

    ! multiroot_jvp: forward product solves J_x dx = -(J_p tp).
    !   jac_x : J_x at the root (n x n).
    !   f_p   : J_p at the root (n x m), f_p(i, j) = dF_i/dp_j.
    !   tp    : tangent in parameter space (length m).
    !   dx    : output tangent in root space (length n).
    subroutine multiroot_jvp(jac_x, f_p, tp, dx, status, &
            preconditioned_solve, preconditioner_context, reciprocal_condition, &
            minimum_reciprocal_condition)
        real(dp),               intent(in)  :: jac_x(:, :)
        real(dp),               intent(in)  :: f_p(:, :)
        real(dp),               intent(in)  :: tp(:)
        real(dp),               intent(out) :: dx(:)
        type(fortnum_status_t), intent(out) :: status
        procedure(multiroot_preconditioned_solve_t), optional :: preconditioned_solve
        class(*), intent(inout), optional :: preconditioner_context
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition

        real(dp) :: rhs(size(jac_x, 1))
        integer  :: n, ls_stat
        logical :: reliable

        n = size(jac_x, 1)
        call assess_condition(jac_x, reciprocal_condition, &
            minimum_reciprocal_condition, reliable)
        if (.not. reliable) then
            dx = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_jvp: ill-conditioned Jacobian; derivative unreliable")
            return
        end if
        rhs = -matmul(f_p, tp)
        call solve_implicit_product(jac_x, rhs, dx, ls_stat, &
            preconditioned_solve, preconditioner_context)
        if (ls_stat /= 0) then
            dx = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_jvp: singular Jacobian; derivative unreliable")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_jvp

    ! Analytical JVP using a reusable LU factorization of the converged F_x.
    subroutine multiroot_jvp_factored(jacobian_factors, ipiv, f_p, tp, &
            dx, status)
        real(dp), contiguous, intent(in) :: jacobian_factors(:, :)
        real(dp), intent(in) :: f_p(:, :), tp(:)
        integer, contiguous, intent(in) :: ipiv(:)
        real(dp), contiguous, intent(out) :: dx(:)
        type(fortnum_status_t), intent(out) :: status

        integer :: i, info, j, n

        n = size(jacobian_factors, 1)
        dx = 0.0_dp
        do j = 1, size(tp)
            do i = 1, n
                dx(i) = dx(i) - f_p(i, j)*tp(j)
            end do
        end do
        call lu_solve_factored(n, jacobian_factors, ipiv, dx, info)
        if (info /= LINALG_OK) then
            dx = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_jvp_factored: invalid or singular factorization")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_jvp_factored

    ! Generic analytical implicit tangent boundary. The callback evaluates F_x
    ! and the contracted residual tangent F_p*tp at the converged root.
    subroutine multiroot_implicit_jvp(residual_jvp, x, p, tp, dx, status, &
            context, reciprocal_condition, minimum_reciprocal_condition)
        procedure(multiroot_residual_jvp_t) :: residual_jvp
        real(dp), intent(in) :: x(:), p(:), tp(:)
        real(dp), intent(out) :: dx(:)
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(inout), optional :: context
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition

        real(dp) :: jac_x(size(x), size(x)), f_p_tp(size(x))

        call residual_jvp(x, p, tp, jac_x, f_p_tp, context)
        call multiroot_grad(jac_x, f_p_tp, dx, status, &
            reciprocal_condition=reciprocal_condition, &
            minimum_reciprocal_condition=minimum_reciprocal_condition)
    end subroutine multiroot_implicit_jvp

    ! multiroot_vjp: reverse product solves J_x^T lambda = u, then jtu = -J_p^T lambda.
    !   jac_x : J_x at the root (n x n).
    !   f_p   : J_p at the root (n x m).
    !   u     : adjoint on the root output x* (length n).
    !   jtu   : adjoint on the parameter vector (length m).
    subroutine multiroot_vjp(jac_x, f_p, u, jtu, status, &
            preconditioned_solve, preconditioner_context, reciprocal_condition, &
            minimum_reciprocal_condition)
        real(dp),               intent(in)  :: jac_x(:, :)
        real(dp),               intent(in)  :: f_p(:, :)
        real(dp),               intent(in)  :: u(:)
        real(dp),               intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status
        procedure(multiroot_preconditioned_solve_t), optional :: preconditioned_solve
        class(*), intent(inout), optional :: preconditioner_context
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition

        real(dp) :: lambda(size(jac_x, 1))
        integer  :: n, ls_stat
        logical :: reliable

        n = size(jac_x, 1)
        call assess_condition(jac_x, reciprocal_condition, &
            minimum_reciprocal_condition, reliable)
        if (.not. reliable) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_vjp: ill-conditioned Jacobian; derivative unreliable")
            return
        end if
        call solve_implicit_product(transpose(jac_x), u, lambda, ls_stat, &
            preconditioned_solve, preconditioner_context)
        if (ls_stat /= 0) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_vjp: singular Jacobian; derivative unreliable")
            return
        end if
        jtu = -matmul(transpose(f_p), lambda)
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_vjp

    ! Analytical VJP using a reusable LU factorization of the converged F_x^T.
    subroutine multiroot_vjp_factored(transpose_factors, ipiv, f_p, u, &
            jtu, status)
        real(dp), contiguous, intent(in) :: transpose_factors(:, :)
        integer, contiguous, intent(in) :: ipiv(:)
        real(dp), intent(in) :: f_p(:, :), u(:)
        real(dp), contiguous, intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status

        real(dp) :: lambda(size(u))
        integer :: i, info, j, n

        n = size(transpose_factors, 1)
        lambda = u
        call lu_solve_factored(n, transpose_factors, ipiv, lambda, info)
        if (info /= LINALG_OK) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_vjp_factored: invalid or singular factorization")
            return
        end if
        jtu = 0.0_dp
        do j = 1, size(jtu)
            do i = 1, n
                jtu(j) = jtu(j) - f_p(i, j)*lambda(i)
            end do
        end do
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_vjp_factored

    ! Generic analytical implicit adjoint boundary. First solve
    ! F_x^T*lambda=u, then evaluate and negate F_p^T*lambda.
    subroutine multiroot_implicit_vjp(state_jacobian, parameter_vjp, x, p, &
            u, jtu, status, context, reciprocal_condition, &
            minimum_reciprocal_condition)
        procedure(multiroot_state_jacobian_t) :: state_jacobian
        procedure(multiroot_parameter_vjp_t) :: parameter_vjp
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(inout), optional :: context
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition

        real(dp) :: jac_x(size(x), size(x)), jac_x_t(size(x), size(x))
        real(dp) :: lambda(size(x)), f_p_t_lambda(size(p))
        integer :: ls_stat
        logical :: reliable

        call state_jacobian(x, p, jac_x, context)
        call assess_condition(jac_x, reciprocal_condition, &
            minimum_reciprocal_condition, reliable)
        if (.not. reliable) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_implicit_vjp: ill-conditioned Jacobian; " // &
                "derivative unreliable")
            return
        end if
        jac_x_t = transpose(jac_x)
        call solve_implicit_product(jac_x_t, u, lambda, ls_stat)
        if (ls_stat /= 0) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_implicit_vjp: singular Jacobian; derivative unreliable")
            return
        end if
        call parameter_vjp(x, p, lambda, f_p_t_lambda, context)
        jtu = -f_p_t_lambda
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_implicit_vjp

    ! multiroot_grad: scalar-p case dx*/dp = -J_x^{-1} f_p, solved as J_x dxdp = -f_p.
    !   jac_x : J_x at the root (n x n).
    !   f_p   : dF/dp at the root (length n).
    !   dxdp  : sensitivity of the root (length n).
    subroutine multiroot_grad(jac_x, f_p, dxdp, status, &
            preconditioned_solve, preconditioner_context, reciprocal_condition, &
            minimum_reciprocal_condition)
        real(dp),               intent(in)  :: jac_x(:, :)
        real(dp),               intent(in)  :: f_p(:)
        real(dp),               intent(out) :: dxdp(:)
        type(fortnum_status_t), intent(out) :: status
        procedure(multiroot_preconditioned_solve_t), optional :: preconditioned_solve
        class(*), intent(inout), optional :: preconditioner_context
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition

        integer :: n, ls_stat
        logical :: reliable

        n = size(jac_x, 1)
        call assess_condition(jac_x, reciprocal_condition, &
            minimum_reciprocal_condition, reliable)
        if (.not. reliable) then
            dxdp = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_grad: ill-conditioned Jacobian; derivative unreliable")
            return
        end if
        call solve_implicit_product(jac_x, -f_p, dxdp, ls_stat, &
            preconditioned_solve, preconditioner_context)
        if (ls_stat /= 0) then
            dxdp = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "multiroot_grad: singular Jacobian; derivative unreliable")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine multiroot_grad

    subroutine solve_implicit_product(a, b, x, info, &
            preconditioned_solve, preconditioner_context)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        procedure(multiroot_preconditioned_solve_t), optional :: preconditioned_solve
        class(*), intent(inout), optional :: preconditioner_context
        integer :: n

        if (present(preconditioned_solve)) then
            if (present(preconditioner_context)) then
                call preconditioned_solve(a, b, x, info, preconditioner_context)
            else
                call preconditioned_solve(a, b, x, info)
            end if
        else
            n = size(a, 1)
            call solve_linear(n, a, b, x, info)
        end if
    end subroutine solve_implicit_product

    subroutine assess_condition(a, reciprocal_condition, &
            minimum_reciprocal_condition, reliable)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out), optional :: reciprocal_condition
        real(dp), intent(in), optional :: minimum_reciprocal_condition
        logical, intent(out) :: reliable
        real(dp) :: measured

        reliable = .true.
        if (.not. present(reciprocal_condition) .and. &
                .not. present(minimum_reciprocal_condition)) return
        call estimate_rcond_1(a, measured)
        if (present(reciprocal_condition)) reciprocal_condition = measured
        if (present(minimum_reciprocal_condition)) then
            reliable = measured >= minimum_reciprocal_condition
        end if
    end subroutine assess_condition

    subroutine estimate_rcond_1(a, reciprocal_condition)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: reciprocal_condition
        real(dp) :: basis(size(a, 1)), inverse_column(size(a, 1))
        real(dp) :: a_norm, inverse_norm
        integer :: j, info, n

        n = size(a, 1)
        a_norm = 0.0_dp
        do j = 1, n
            a_norm = max(a_norm, sum(abs(a(:, j))))
        end do
        if (a_norm == 0.0_dp) then
            reciprocal_condition = 0.0_dp
            return
        end if

        inverse_norm = 0.0_dp
        do j = 1, n
            basis = 0.0_dp
            basis(j) = 1.0_dp
            call solve_linear(n, a, basis, inverse_column, info)
            if (info /= 0) then
                reciprocal_condition = 0.0_dp
                return
            end if
            inverse_norm = max(inverse_norm, sum(abs(inverse_column)))
        end do
        reciprocal_condition = 1.0_dp/(a_norm*inverse_norm)
    end subroutine estimate_rcond_1

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
        ! No decrease within 30 halvings: accept the last, smallest trial step
        ! so the iteration still advances.
        x = xtrial; f = ftrial; jac = jtrial
    end subroutine line_search_fdf

    ! Powell's dogleg step for min |R p - qtb|_2 subject to |p|_2 <= delta,
    ! with R upper triangular (MINPACK dogleg, unit scaling).  Returns the
    ! Gauss-Newton point R^{-1} qtb when it lies inside the region; otherwise
    ! the point where the path from the Cauchy point of the scaled gradient
    ! R^T qtb to the Gauss-Newton point crosses the boundary, or the scaled
    ! gradient step clipped to the boundary.  A zero diagonal of R is replaced
    ! by eps times its column maximum, so a singular R still yields a step.
    pure subroutine dogleg(n, r, qtb, delta, p)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: r(n, n), qtb(n), delta
        real(dp), intent(out) :: p(n)

        real(dp) :: g(n), rg(n), eps, piv, qnorm, gnorm, sgnorm, bnorm
        real(dp) :: alpha, temp, dq, sd
        integer  :: j

        eps = epsilon(1.0_dp)
        do j = n, 1, -1
            piv = r(j, j)
            if (piv == 0.0_dp) then
                piv = eps * maxval(abs(r(1:j, j)))
                if (piv == 0.0_dp) piv = eps
            end if
            p(j) = (qtb(j) - dot_product(r(j, j+1:n), p(j+1:n))) / piv
        end do
        qnorm = norm2(p)
        if (qnorm <= delta) return

        do j = 1, n
            g(j) = dot_product(r(1:j, j), qtb(1:j))
        end do
        gnorm = norm2(g)
        sgnorm = 0.0_dp
        alpha = delta / qnorm
        if (gnorm /= 0.0_dp) then
            g = g / gnorm
            do j = 1, n
                rg(j) = dot_product(r(j, j:n), g(j:n))
            end do
            temp = norm2(rg)
            ! Distance to the minimiser of the model along the gradient.
            sgnorm = (gnorm / temp) / temp
            alpha = 0.0_dp
            if (sgnorm < delta) then
                bnorm = norm2(qtb)
                dq = delta / qnorm
                sd = sgnorm / delta
                temp = (bnorm / gnorm) * (bnorm / qnorm) * sd
                temp = temp - dq * sd**2 + sqrt((temp - dq)**2 &
                    + (1.0_dp - dq**2) * (1.0_dp - sd**2))
                alpha = dq * (1.0_dp - sd**2) / temp
            end if
        end if
        p = (1.0_dp - alpha) * min(sgnorm, delta) * g + alpha * p
    end subroutine dogleg

    ! Householder QR of a square matrix: a = q r with q orthogonal and r upper
    ! triangular, q formed explicitly for the rank-one updates.
    pure subroutine qr_householder(n, a, q, r)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: a(n, n)
        real(dp), intent(out) :: q(n, n), r(n, n)

        real(dp) :: w(n), alpha, wnorm2
        integer  :: k, j

        r = a
        q = 0.0_dp
        do j = 1, n
            q(j, j) = 1.0_dp
        end do
        do k = 1, n - 1
            alpha = norm2(r(k:n, k))
            if (alpha == 0.0_dp) cycle
            if (r(k, k) > 0.0_dp) alpha = -alpha
            w(k:n) = r(k:n, k)
            w(k) = w(k) - alpha
            wnorm2 = dot_product(w(k:n), w(k:n))
            do j = k, n
                r(k:n, j) = r(k:n, j) &
                    - (2.0_dp * dot_product(w(k:n), r(k:n, j)) / wnorm2) * w(k:n)
            end do
            do j = 1, n
                q(j, k:n) = q(j, k:n) &
                    - (2.0_dp * dot_product(q(j, k:n), w(k:n)) / wnorm2) * w(k:n)
            end do
            r(k+1:n, k) = 0.0_dp
        end do
    end subroutine qr_householder

    ! Refactor q (r + u v^T) as q' r' with r' upper triangular, by Givens
    ! rotations (Golub & Van Loan, Matrix Computations, sec. 6.5.1); qtf is
    ! rotated along so that it stays q'^T F.  Same role as MINPACK
    ! r1updt/r1mpyq, here on dense storage.
    pure subroutine qr_rank_one_update(n, q, r, qtf, u, v)
        integer,  intent(in)    :: n
        real(dp), intent(inout) :: q(n, n), r(n, n), qtf(n)
        real(dp), intent(in)    :: u(n), v(n)

        real(dp) :: w(n), c, s
        integer  :: k

        w = u
        ! Rotate w onto e_1; r becomes upper Hessenberg.
        do k = n - 1, 1, -1
            call givens(w(k), w(k + 1), c, s)
            w(k) = c * w(k) + s * w(k + 1)
            w(k + 1) = 0.0_dp
            call apply_givens(n, k, c, s, q, r, qtf)
        end do
        r(1, :) = r(1, :) + w(1) * v
        ! Restore triangular form.
        do k = 1, n - 1
            call givens(r(k, k), r(k + 1, k), c, s)
            call apply_givens(n, k, c, s, q, r, qtf)
            r(k + 1, k) = 0.0_dp
        end do

    end subroutine qr_rank_one_update

    ! Apply the rotation G = [c -s; s c] in the (k, k+1) plane: r <- G^T r,
    ! qtf <- G^T qtf, q <- q G, so that q r and q qtf are unchanged.
    pure subroutine apply_givens(n, k, c, s, q, r, qtf)
        integer,  intent(in)    :: n, k
        real(dp), intent(in)    :: c, s
        real(dp), intent(inout) :: q(n, n), r(n, n), qtf(n)
        real(dp) :: row(n), col(n), t

        row = r(k, :)
        r(k, :) = c * row + s * r(k + 1, :)
        r(k + 1, :) = -s * row + c * r(k + 1, :)
        col = q(:, k)
        q(:, k) = c * col + s * q(:, k + 1)
        q(:, k + 1) = -s * col + c * q(:, k + 1)
        t = qtf(k)
        qtf(k) = c * t + s * qtf(k + 1)
        qtf(k + 1) = -s * t + c * qtf(k + 1)
    end subroutine apply_givens

    ! Rotation [c s; -s c] that maps (a, b) to (hypot(a, b), 0).
    pure subroutine givens(a, b, c, s)
        real(dp), intent(in)  :: a, b
        real(dp), intent(out) :: c, s
        real(dp) :: h

        h = hypot(a, b)
        if (h == 0.0_dp) then
            c = 1.0_dp
            s = 0.0_dp
        else
            c = a / h
            s = b / h
        end if
    end subroutine givens

    ! Forward-difference Jacobian at x from the known f = F(x), with
    ! MINPACK fdjac1's step sqrt(eps)*|x_j| (sqrt(eps) when x_j = 0).
    subroutine fd_jacobian_forward(fn, n, x, f, jac, ctx)
        procedure(multiroot_fn_t)      :: fn
        integer,  intent(in)           :: n
        real(dp), intent(in)           :: x(n), f(n)
        real(dp), intent(out)          :: jac(n, n)
        class(*), intent(in), optional :: ctx

        real(dp) :: xp(n), fp(n), h, root_eps
        integer  :: j

        root_eps = sqrt(epsilon(1.0_dp))
        do j = 1, n
            h = root_eps * abs(x(j))
            if (h == 0.0_dp) h = root_eps
            xp = x
            xp(j) = x(j) + h
            call fn(xp, fp, ctx)
            jac(:, j) = (fp - f) / h
        end do
    end subroutine fd_jacobian_forward

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
    ! beyond perm.  Not stable: equal values may come in any index order.
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
            ! Check for a leaf before doubling a possibly large index.
            if (root > heap_n / 2) exit
            child = 2 * root
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
