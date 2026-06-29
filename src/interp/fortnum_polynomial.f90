module fortnum_polynomial
    ! Lagrange interpolation weights on an arbitrary support node set.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4):
    !   Default class: analytic_rule.
    !   lagrange_weights: the barycentric-form weights are linear in the nodal
    !     values; value weights are transparent (linear map). The derivative
    !     weights are a closed-form analytic rule: see lagrange_deriv_weights.
    !   Active arguments (ad.md §3): nodal values f(1:n) are active when the
    !     caller differentiates through the interpolated result. The evaluation
    !     point x and support nodes xp are active when the caller differentiates
    !     with respect to position; that path is deferred to issue #40.
    !   Inactive: n (size, selects problem dimension).
    !   Derivative entry points (ad.md §2): lagrange_weights_jvp -- not
    !     implemented here; will be added in the AD milestone without touching
    !     the primal signatures below.
    !
    ! The two exported routines compute weights coef(1:n) and dcoef(1:n) such
    ! that for any function f sampled at the support nodes xp(1:n):
    !
    !   p(x)  = sum_{i=1}^{n} f(xp(i)) * coef(i)     [interpolant at x]
    !   p'(x) = sum_{i=1}^{n} f(xp(i)) * dcoef(i)    [derivative of interpolant]
    !
    ! Weights depend only on x and xp, not on f, so they can be precomputed
    ! once and reused over many function-value vectors (the pattern issue #40
    ! exploits for JVPs).
    !
    ! Algorithm: the direct Lagrange formula without barycentric precomputation
    ! is used because n is assumed small (typically 2--8). The derivative weights
    ! follow the product-rule expansion of the Lagrange basis, matching the
    ! algorithm in libneo plag_coeff.f90 in behavior (clean-room reimplementation).
    !
    ! Nodes xp(1:n) must be distinct; behavior is undefined if any two coincide.

    use fortnum_kinds, only: dp
    implicit none
    private

    public :: lagrange_weights
    public :: lagrange_deriv_weights
    public :: lagrange_weights_jvp ! d p(x)/d x  . v_x   (x active)
    public :: lagrange_weights_vjp ! (d p(x)/d x)^T . u   (x active)
    public :: lagrange_fval_jvp ! d p(x)/d f  . v_f   (f values active)
    public :: lagrange_fval_vjp ! (d p(x)/d f)^T . u  (f values active)

contains

    ! Compute Lagrange value weights coef(1:n) at point x over nodes xp(1:n).
    ! coef(i) = prod_{k /= i} (x - xp(k)) / (xp(i) - xp(k)).
    pure subroutine lagrange_weights(n, x, xp, coef)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(out) :: coef(n)

        integer  :: i, k
        real(dp) :: w

        do i = 1, n
            w = 1.0_dp
            do k = 1, n
                if (k == i) cycle
                w = w*(x - xp(k))/(xp(i) - xp(k))
            end do
            coef(i) = w
        end do
    end subroutine lagrange_weights


    ! Compute Lagrange derivative weights dcoef(1:n) at point x over nodes xp.
    ! dcoef(i) = d/dx L_i(x) where L_i is the i-th Lagrange basis polynomial.
    !
    ! Each basis derivative is the sum over terms that drop one factor from the
    ! product defining L_i:
    !   d/dx L_i(x) = sum_{k /= i} [ prod_{j /= i, j /= k} (x-xp(j))/(xp(i)-xp(j)) ]
    !                               * 1/(xp(i)-xp(k))
    !
    ! Implementation keeps an auxiliary array to accumulate the contribution from
    ! dropping each factor k, matching the loop structure of plag_coeff.f90 in
    ! behavior without copying code.
    pure subroutine lagrange_deriv_weights(n, x, xp, dcoef)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(out) :: dcoef(n)

        integer  :: i, k, j
        real(dp) :: fac
        real(dp) :: tmp(n) ! running partial products for each j

        do i = 1, n
            ! tmp(j) accumulates the product over k /= i, updated factor by factor.
            ! Dropping factor k from L_i contributes 1/(xp(i)-xp(k)) to tmp(k)
            ! and (x-xp(k))/(xp(i)-xp(k)) to all tmp(j), j /= k.
            tmp    = 1.0_dp
            tmp(i) = 0.0_dp ! L_i never contributes to its own derivative slot
            do k = 1, n
                if (k == i) cycle
                fac = (x - xp(k))/(xp(i) - xp(k))
                do j = 1, n
                    if (j == k) then
                        ! dropping factor k: replace (x-xp(k))/(xp(i)-xp(k)) with
                        ! its derivative 1/(xp(i)-xp(k))
                        tmp(j) = tmp(j)/(xp(i) - xp(k))
                    else
                        tmp(j) = tmp(j)*fac
                    end if
                end do
            end do
            dcoef(i) = sum(tmp)
        end do
    end subroutine lagrange_deriv_weights


    ! JVP of p(x) = sum_i f(i)*coef_i(x) with respect to the evaluation
    ! point x.  Policy: analytic_rule (ad.md §4).
    !
    ! Active: x (scalar). Inactive: f (nodal values), xp (support), n.
    ! Valid only inside a fixed cell; crossing a cell boundary is non-smooth
    ! (ad.md §4 interp note) and the caller must hold the cell index fixed.
    !
    ! jv = (d p / d x) * vx = [sum_i f(i) * dcoef(i)] * vx
    !
    ! Deferred: HVP (d^2 p / dx^2 would need second-order weights; out of
    ! scope for the current sweep).
    pure subroutine lagrange_weights_jvp(n, x, xp, f, vx, jv)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(in)  :: f(n) ! nodal values (inactive in this map)
        real(dp), intent(in)  :: vx ! tangent for x
        real(dp), intent(out) :: jv ! d p(x)/d x * vx
        real(dp) :: dcoef(n)
        call lagrange_deriv_weights(n, x, xp, dcoef)
        jv = dot_product(f, dcoef) * vx
    end subroutine lagrange_weights_jvp


    ! VJP of p(x) w.r.t. x.  Scalar output -> VJP collapses to the gradient
    ! times the output cotangent u.
    !
    ! jtu = u * (d p / d x) = u * [sum_i f(i) * dcoef(i)]
    pure subroutine lagrange_weights_vjp(n, x, xp, f, u, jtu)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(in)  :: f(n)
        real(dp), intent(in)  :: u ! output cotangent
        real(dp), intent(out) :: jtu ! input cotangent for x
        real(dp) :: dcoef(n)
        call lagrange_deriv_weights(n, x, xp, dcoef)
        jtu = u * dot_product(f, dcoef)
    end subroutine lagrange_weights_vjp


    ! JVP of p(x) w.r.t. the nodal values f(1:n).  Policy: transparent
    ! (p is linear in f; value weights coef are the Jacobian row).
    !
    ! Active: f(1:n). Inactive: x, xp, n.
    ! Smooth everywhere (the linear dependence on f has no branch).
    !
    ! jv = J_f p . vf = sum_i coef(i) * vf(i)
    pure subroutine lagrange_fval_jvp(n, x, xp, vf, jv)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(in)  :: vf(n) ! tangent for f values
        real(dp), intent(out) :: jv ! scalar tangent for p(x)
        real(dp) :: coef(n)
        call lagrange_weights(n, x, xp, coef)
        jv = dot_product(coef, vf)
    end subroutine lagrange_fval_jvp


    ! VJP of p(x) w.r.t. f(1:n).  J_f p = coef^T (1 x n), so
    ! (J_f p)^T u = coef * u  (scalar u times the weight vector).
    !
    ! jtu(i) = coef(i) * u
    pure subroutine lagrange_fval_vjp(n, x, xp, u, jtu)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x
        real(dp), intent(in)  :: xp(n)
        real(dp), intent(in)  :: u ! output cotangent (scalar)
        real(dp), intent(out) :: jtu(n) ! input cotangents for f(1:n)
        real(dp) :: coef(n)
        call lagrange_weights(n, x, xp, coef)
        jtu = coef * u
    end subroutine lagrange_fval_vjp

end module fortnum_polynomial
