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
        real(dp) :: tmp(n)  ! running partial products for each j

        do i = 1, n
            ! tmp(j) accumulates the product over k /= i, updated factor by factor.
            ! Dropping factor k from L_i contributes 1/(xp(i)-xp(k)) to tmp(k)
            ! and (x-xp(k))/(xp(i)-xp(k)) to all tmp(j), j /= k.
            tmp    = 1.0_dp
            tmp(i) = 0.0_dp   ! L_i never contributes to its own derivative slot
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

end module fortnum_polynomial
