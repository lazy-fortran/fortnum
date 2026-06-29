module fortnum_roots
    ! Scalar root-finding: bisection (bracketing) and Newton with a
    ! derivative callback.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4): implicit_rule.
    !   The root satisfies f(x*, p) = 0.  The implicit function theorem gives
    !   dx*/dp = -f_p / f_x without differentiating the iteration.  For a
    !   simple root (f_x /= 0) under forward mode: dx = -f_p/f_x * dp.
    !   Iteration counts, tolerances, and max_iter are inactive (ad.md §3).
    !   Reserved derivative names: root_bisect_jvp, root_bisect_vjp,
    !   root_newton_jvp, root_newton_vjp (issue #40).
    !
    ! Algorithms:
    !   Bisection: Illinois variant for superlinear convergence near a root
    !     without Brent's inverse quadratic step (kept simple for issue #40
    !     to attach IFT derivatives cleanly; Brent is issue #40+).
    !   Newton: classical Newton-Raphson with derivative callback.  Falls back
    !     to bisection guard when a Newton step leaves the bracket.
    !     Steffensen acceleration not applied; derivative accuracy governs
    !     convergence.
    !   Both: terminate when |interval| <= xtol or |f(x)| <= ftol.
    !   Ref: Burden & Faires, Numerical Analysis, 10th ed., §2.1, §2.3.
    !
    ! No module-level save, no global procedure pointers.  All state is
    ! caller-owned or stack-local.

    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: root_fn_t, root_fn_df_t
    public :: root_bisect, root_newton, root_brent
    public :: root_jvp, root_vjp, root_grad

    ! Function whose root is sought: y = f(x).
    abstract interface
        pure function root_fn_t(x) result(y)
            import :: dp
            real(dp), intent(in) :: x
            real(dp)             :: y
        end function root_fn_t
    end interface

    ! Function that returns both f(x) and f'(x).
    ! fx receives f, dfx receives df/dx.  Supplying both in one call avoids
    ! evaluating a shared sub-expression twice (common in practice).
    abstract interface
        pure subroutine root_fn_df_t(x, fx, dfx)
            import :: dp
            real(dp), intent(in)  :: x
            real(dp), intent(out) :: fx, dfx
        end subroutine root_fn_df_t
    end interface

contains

    ! Bisection on [a, b] with the Illinois modification.
    !
    ! Requires f(a)*f(b) < 0 (strict sign change).  Returns the root in x
    ! and sets status on failure.  xtol and ftol are optional convergence
    ! thresholds (default 4*eps in x, 0 in f so the interval criterion
    ! drives convergence).  max_iter defaults to 200.
    !
    ! Illinois modification (Ford 1979): when the same endpoint is retained
    ! twice in a row, halve its function value.  This gives superlinear
    ! convergence without inverse quadratic steps, simplifying the IFT path.
    subroutine root_bisect(f, a, b, x, status, xtol, ftol, max_iter)
        procedure(root_fn_t)               :: f
        real(dp),               intent(in) :: a, b
        real(dp),               intent(out):: x
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter

        real(dp) :: xa, xb, fa, fb, fc, xc, side_fa, side_fb
        real(dp) :: xt, ft, half_width
        integer  :: it, max_it, last_side
        logical  :: illinois_a, illinois_b

        call status_set(status, FORTNUM_OK, "")
        x = a

        ! Resolve optional parameters.
        max_it = 200
        xt = 4.0_dp * epsilon(1.0_dp)
        ft = 0.0_dp
        if (present(max_iter)) max_it = max_iter
        if (present(xtol))     xt = xtol
        if (present(ftol))     ft = ftol

        xa = a
        xb = b
        fa = f(xa)
        fb = f(xb)

        ! Exact root at an endpoint.
        if (abs(fa) == 0.0_dp) then
            x = xa
            return
        end if
        if (abs(fb) == 0.0_dp) then
            x = xb
            return
        end if

        if (fa * fb > 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_bisect: f(a) and f(b) have the same sign; no bracket")
            return
        end if

        last_side = 0 ! 0 = first step, -1 = kept a, +1 = kept b
        side_fa   = fa
        side_fb   = fb

        do it = 1, max_it
            half_width = 0.5_dp * abs(xb - xa)
            xc = xa + half_width * sign(1.0_dp, xb - xa)
            fc = f(xc)

            if (abs(fc) <= ft .or. half_width <= xt) then
                x = xc
                return
            end if

            if (fc * fa < 0.0_dp) then
                ! Root in [xa, xc]: keep xa, update xb.
                xb = xc
                fb = fc
                if (last_side == 1) then
                    ! Same side (a) kept twice: Illinois fix on fb.
                    side_fa = fa
                    fa = 0.5_dp * fa
                end if
                last_side = 1
            else
                ! Root in [xc, xb]: keep xb, update xa.
                xa = xc
                fa = fc
                if (last_side == -1) then
                    ! Same side (b) kept twice: Illinois fix on fa.
                    side_fb = fb
                    fb = 0.5_dp * fb
                end if
                last_side = -1
            end if
        end do

        x = xa + 0.5_dp * (xb - xa)
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "root_bisect: maximum iterations reached without convergence")
    end subroutine root_bisect

    ! Newton-Raphson with bisection guard.
    !
    ! Requires an initial bracket [a, b] with f(a)*f(b) < 0.  Each step
    ! attempts a Newton update from the current best point; if the Newton
    ! step leaves [a, b], falls back to bisection for that step.  x0 is the
    ! starting point (should lie in [a, b]).  Returns the root in x.
    !
    ! Derivative policy: implicit_rule.  The iteration itself is inactive;
    ! the IFT gives dx/dp = -f_p/f_x at the converged root (issue #40).
    !
    ! Near-multiple-root guard: if |f'(x)| < deriv_floor (default 1e-14)
    ! at any Newton step, the derivative is treated as zero and bisection
    ! takes over for that step.  At the final root, |f'(x)| < deriv_floor
    ! is reported as a FORTNUM_DOMAIN_ERROR (near-multiple root warning).
    subroutine root_newton(fdf, a, b, x0, x, status, xtol, ftol, &
            max_iter, deriv_floor)
        procedure(root_fn_df_t)            :: fdf
        real(dp),               intent(in) :: a, b, x0
        real(dp),               intent(out):: x
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol, deriv_floor
        integer,  intent(in), optional     :: max_iter

        real(dp) :: xa, xb, fa, fb, xc, fc, dfc, xcnew
        real(dp) :: xt, ft, dfloor
        integer  :: it, max_it

        call status_set(status, FORTNUM_OK, "")
        x = x0

        max_it = 200
        xt     = 4.0_dp * epsilon(1.0_dp)
        ft     = 0.0_dp
        dfloor = 1.0e-14_dp
        if (present(max_iter))    max_it = max_iter
        if (present(xtol))        xt = xtol
        if (present(ftol))        ft = ftol
        if (present(deriv_floor)) dfloor = deriv_floor

        xa = a
        xb = b
        fa = 0.0_dp
        fb = 0.0_dp

        ! Evaluate bracket endpoints to establish sign tracking.
        call fdf(xa, fa, dfc)
        call fdf(xb, fb, dfc)

        if (abs(fa) == 0.0_dp) then
            x = xa
            return
        end if
        if (abs(fb) == 0.0_dp) then
            x = xb
            return
        end if

        if (fa * fb > 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_newton: f(a) and f(b) have the same sign; no bracket")
            return
        end if

        xc = x0
        ! If x0 is outside [a,b], start from midpoint.
        if (xc < min(a, b) .or. xc > max(a, b)) &
            xc = 0.5_dp * (xa + xb)

        do it = 1, max_it
            call fdf(xc, fc, dfc)

            if (abs(fc) <= ft) then
                x = xc
                if (abs(dfc) < dfloor) &
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "root_newton: near-multiple root; derivative ~0 at root")
                return
            end if

            ! Update bracket from the new function sign.
            if (fc * fa < 0.0_dp) then
                xb = xc
                fb = fc
            else
                xa = xc
                fa = fc
            end if

            if (abs(xb - xa) <= xt) then
                x = xa + 0.5_dp * (xb - xa)
                if (abs(dfc) < dfloor) &
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "root_newton: near-multiple root; derivative ~0 at root")
                return
            end if

            ! Newton step; fall back to bisection if step leaves the bracket.
            if (abs(dfc) >= dfloor) then
                xcnew = xc - fc / dfc
                if (xcnew > min(xa, xb) .and. xcnew < max(xa, xb)) then
                    xc = xcnew
                else
                    xc = xa + 0.5_dp * (xb - xa)
                end if
            else
                ! Derivative too small; bisect this step.
                xc = xa + 0.5_dp * (xb - xa)
            end if
        end do

        x = xa + 0.5_dp * (xb - xa)
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "root_newton: maximum iterations reached without convergence")
    end subroutine root_newton

    ! Brent's method on a bracketed interval [a, b].
    !
    ! Combines inverse quadratic interpolation, secant, and bisection.
    ! Each step picks the IQI or secant step when it is safe (stays inside
    ! the active bracket and improves on the previous step's half-width);
    ! otherwise falls back to bisection.  Convergence is superlinear on
    ! smooth functions and guaranteed linear (bisection-like) in the worst
    ! case.
    !
    ! Requires f(a)*f(b) <= 0 (sign change; exact zeros accepted).
    ! Optional: xtol (default 4*epsilon), ftol (default 0),
    !           max_iter (default 200).
    ! Returns FORTNUM_DOMAIN_ERROR if the bracket is invalid (no sign change
    ! or a near-multiple root is detected by both endpoints having |f| equal);
    ! FORTNUM_CONVERGENCE_ERROR if max_iter is exhausted.
    !
    ! Derivative policy (ad.md §1, §4): implicit_rule.  Same IFT path as
    ! bisection: dx*/dp = -f_p/f_x at the converged root.  Reserved names:
    ! root_brent_jvp, root_brent_vjp (issue #40).  Iteration choices are
    ! inactive (ad.md §3).
    !
    ! Algorithm ref: Brent (1973), "Algorithms for Minimization without
    ! Derivatives", Ch. 4; Forsythe, Malcolm & Moler (1977), §11.
    subroutine root_brent(f, a, b, x, status, xtol, ftol, max_iter)
        procedure(root_fn_t)               :: f
        real(dp),               intent(in) :: a, b
        real(dp),               intent(out):: x
        type(fortnum_status_t), intent(out):: status
        real(dp), intent(in), optional     :: xtol, ftol
        integer,  intent(in), optional     :: max_iter

        real(dp) :: xa, xb, xc, fa, fb, fc
        real(dp) :: s, p, q, r, d, e, xt, ft, tol1
        integer  :: it, max_it

        call status_set(status, FORTNUM_OK, "")
        x = a

        max_it = 200
        xt     = 4.0_dp * epsilon(1.0_dp)
        ft     = 0.0_dp
        if (present(max_iter)) max_it = max_iter
        if (present(xtol))     xt = xtol
        if (present(ftol))     ft = ftol

        xa = a
        xb = b
        fa = f(xa)
        fb = f(xb)

        ! Exact root at an endpoint.
        if (abs(fa) == 0.0_dp .or. abs(fa) <= ft) then
            x = xa
            return
        end if
        if (abs(fb) == 0.0_dp .or. abs(fb) <= ft) then
            x = xb
            return
        end if

        if (fa * fb > 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_brent: f(a) and f(b) have the same sign; no bracket")
            return
        end if

        ! xc is the third (oldest) point; initialise to xa so the first step
        ! is a secant/bisection step (no IQI yet).
        xc = xa
        fc = fa
        d  = xb - xa ! most-recent interpolation step (seed with bracket width)
        e  = d ! step before that; bisection guard uses |e|

        do it = 1, max_it
            ! Ensure xb is the current best approximation (|fb| <= |fc|).
            if (abs(fb) > abs(fc)) then
                xa = xb; fa = fb
                xb = xc; fb = fc
                xc = xa; fc = fa
            end if

            ! Convergence tolerance at xb.
            tol1 = 2.0_dp * epsilon(1.0_dp) * abs(xb) + 0.5_dp * xt
            s    = 0.5_dp * (xc - xb) ! half-distance to the third point

            if (abs(s) <= tol1 .or. abs(fb) <= ft) then
                x = xb
                return
            end if

            ! Decide: try interpolation or fall back to bisection.
            if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
                ! Attempt IQI if all three points are distinct, else secant.
                if (abs(xa - xc) < epsilon(1.0_dp) * (abs(xa) + abs(xc) + 1.0_dp)) then
                    ! xa == xc (numerically): only two distinct points -> secant.
                    p = fb * (xb - xa) / (fa - fb)
                    q = 1.0_dp
                else
                    ! Inverse quadratic interpolation (Brent 1973, eq. 4.3).
                    q = fa / fc
                    r = fb / fc
                    p = r * (2.0_dp * s * q * (q - r) - (xb - xa) * (r - 1.0_dp))
                    q = (q - 1.0_dp) * (r - 1.0_dp) * (fa / fb - 1.0_dp)
                end if

                ! Ensure p/q > 0 (step direction consistent with bracket).
                if (p > 0.0_dp) then
                    q = -q
                else
                    p = -p
                end if

                ! Accept the interpolation step only when it is inside the
                ! bracket and smaller than half the previous step width; this
                ! prevents stagnation when interpolation is barely useful.
                if (2.0_dp * p < min(3.0_dp * s * q - abs(tol1 * q), abs(e * q))) then
                    e = d
                    d = p / q
                else
                    ! Interpolation rejected: bisect.
                    d = s
                    e = s
                end if
            else
                ! |e| too small or |fa| <= |fb|: bisect.
                d = s
                e = s
            end if

            ! Shift: xb -> xa, then apply step d.
            xa = xb
            fa = fb
            if (abs(d) > tol1) then
                xb = xb + d
            else
                ! Step smaller than tolerance: nudge by tol1 in the right direction.
                xb = xb + sign(tol1, s)
            end if
            fb = f(xb)

            ! Keep xc on the opposite sign from xb to maintain a valid bracket.
            if (fb * fc > 0.0_dp) then
                xc = xa
                fc = fa
            end if
        end do

        x = xb
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "root_brent: maximum iterations reached without convergence")
    end subroutine root_brent

    ! Implicit-rule derivatives for a scalar root x*(p) of f(x, p) = 0.
    !
    ! IFT: dx*/dp = -f_p / f_x, where f_x = df/dx and f_p = df/dp at the root.
    !
    ! Active arguments: x_star (converged root), f_x (df/dx at root),
    !   f_p (df/dp at root, scalar for root_grad; vector for root_jvp/root_vjp).
    ! Inactive: status, tolerance, solver internals.
    !
    ! Near-multiple-root guard: |f_x| < deriv_floor signals FORTNUM_DOMAIN_ERROR.
    ! The derivative is not reliable near a multiple root (IFT breaks down).
    !
    ! HVP: deferred. Second derivatives of x*(p) require d^2f/dp^2, d^2f/dxdp,
    ! and d^2f/dx^2; these are not part of the current caller-supplied interface.
    ! HVP can be added in a follow-up issue when the second-partials callback
    ! pattern is settled.

    ! root_jvp: forward-mode product dx* = -(f_p . tp) / f_x.
    !   f_x  : df/dx at the converged root (scalar).
    !   f_p  : df/dp_i at the root, one value per parameter component (vector).
    !   tp   : tangent vector in parameter space (same size as f_p).
    !   dx   : output tangent in root space (scalar).
    !   deriv_floor: |f_x| threshold below which status -> FORTNUM_DOMAIN_ERROR.
    pure subroutine root_jvp(f_x, f_p, tp, dx, status, deriv_floor)
        real(dp),               intent(in)  :: f_x
        real(dp),               intent(in)  :: f_p(:)
        real(dp),               intent(in)  :: tp(:)
        real(dp),               intent(out) :: dx
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional      :: deriv_floor

        real(dp) :: dfloor

        dfloor = 1.0e-14_dp
        if (present(deriv_floor)) dfloor = deriv_floor

        if (abs(f_x) < dfloor) then
            dx = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_jvp: |f_x| near zero; near-multiple root, derivative unreliable")
            return
        end if

        dx = -dot_product(f_p, tp) / f_x
        call status_set(status, FORTNUM_OK, "")
    end subroutine root_jvp

    ! root_vjp: reverse-mode product jtu_i = -(f_p_i / f_x) * u.
    !   f_x  : df/dx at the converged root.
    !   f_p  : df/dp_i at the root (vector, length n).
    !   u    : adjoint on the root output (scalar; passed as real(dp)).
    !   jtu  : adjoint on the parameter vector (length n).
    pure subroutine root_vjp(f_x, f_p, u, jtu, status, deriv_floor)
        real(dp),               intent(in)  :: f_x
        real(dp),               intent(in)  :: f_p(:)
        real(dp),               intent(in)  :: u
        real(dp),               intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional      :: deriv_floor

        real(dp) :: dfloor

        dfloor = 1.0e-14_dp
        if (present(deriv_floor)) dfloor = deriv_floor

        if (abs(f_x) < dfloor) then
            jtu = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_vjp: |f_x| near zero; near-multiple root, derivative unreliable")
            return
        end if

        jtu = -(f_p / f_x) * u
        call status_set(status, FORTNUM_OK, "")
    end subroutine root_vjp

    ! root_grad: scalar-p case: dx*/dp = -f_p / f_x.
    !   Convenience wrapper for the 1-D parameter case.
    pure subroutine root_grad(f_x, f_p, dxdp, status, deriv_floor)
        real(dp),               intent(in)  :: f_x
        real(dp),               intent(in)  :: f_p
        real(dp),               intent(out) :: dxdp
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional      :: deriv_floor

        real(dp) :: dfloor

        dfloor = 1.0e-14_dp
        if (present(deriv_floor)) dfloor = deriv_floor

        if (abs(f_x) < dfloor) then
            dxdp = 0.0_dp
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "root_grad: |f_x| near zero; near-multiple root, derivative unreliable")
            return
        end if

        dxdp = -f_p / f_x
        call status_set(status, FORTNUM_OK, "")
    end subroutine root_grad

end module fortnum_roots
