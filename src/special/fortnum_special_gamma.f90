module fortnum_special_gamma
    ! Lower incomplete gamma function and the regularized variant P(a,x).
    !
    ! Convention (DLMF 8.2.1):
    !   gamma_lower(a, x) = integral_0^x  t^{a-1} exp(-t) dt
    ! This is the UNNORMALIZED lower incomplete gamma; equivalently,
    !   gamma_lower(a, x) = P(a, x) * Gamma(a)
    ! where P is the regularized form returned by gamma_reg_p.
    !
    ! Derivative policy (docs/design/ad.md sec. 1, 4):
    !   Default class: analytic_rule.
    !   Reason: d/dx gamma_lower(a,x) = x^{a-1} exp(-x) (DLMF 8.8.1) and
    !           d/da involves the derivative of P w.r.t. a which requires
    !           the derivative of log_gamma; closed-form recurrences exist
    !           and are cheaper/more stable than differentiating the
    !           iteration branches. Derivative entry points (gamma_lower_grad,
    !           gamma_reg_p_grad) are reserved; not yet implemented (issue #40).
    !   Active arguments: a (shape parameter), x (integration limit).
    !   Inactive arguments: none beyond standard control flow.
    !
    ! Algorithms (DLMF 8.11.4, 8.9.2):
    !   series expansion converges for x < a+1; continued fraction for x >= a+1.
    !
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: gamma_lower
    public :: gamma_reg_p
    ! analytic_rule derivatives (ad.md §2):
    !   gamma_lower_jvp: d/dx gamma_lower(a,x) * v = x^{a-1} e^{-x} * v
    !                    (DLMF 8.8.1; a fixed, x active).
    !   d/da: deferred. Requires d/da log_gamma (digamma) and a series for
    !         d/da P(a,x); no trivial closed form. Documented in the module
    !         header and reserved as gamma_lower_jvp_da for a future issue.
    public :: gamma_lower_jvp

contains

    ! Unnormalized lower incomplete gamma: gamma_lower(a,x) = P(a,x)*Gamma(a).
    ! a > 0, x >= 0.  Returns 0 for x == 0.
    ! Stops with an error message on bad domain (mirrors the libneo convention;
    ! a graceful status path can be added when fortnum_status integration is
    ! needed without changing the primal signature).
    pure function gamma_lower(a, x) result(g)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: x
        real(dp) :: g

        if (a <= 0.0_dp) error stop "gamma_lower: a must be positive"
        if (x < 0.0_dp)  error stop "gamma_lower: x must be non-negative"

        if (x == 0.0_dp) then
            g = 0.0_dp
        else if (x < a + 1.0_dp) then
            g = gamma_reg_p_series(a, x) * gamma(a)
        else
            g = (1.0_dp - gamma_reg_q_contfrac(a, x)) * gamma(a)
        end if
    end function gamma_lower

    ! Regularized lower incomplete gamma P(a,x) = gamma_lower(a,x)/Gamma(a).
    ! a > 0, x >= 0.
    pure function gamma_reg_p(a, x) result(p)
        real(dp), intent(in) :: a
        real(dp), intent(in) :: x
        real(dp) :: p

        if (a <= 0.0_dp) error stop "gamma_reg_p: a must be positive"
        if (x < 0.0_dp)  error stop "gamma_reg_p: x must be non-negative"

        if (x == 0.0_dp) then
            p = 0.0_dp
        else if (x < a + 1.0_dp) then
            p = gamma_reg_p_series(a, x)
        else
            p = 1.0_dp - gamma_reg_q_contfrac(a, x)
        end if
    end function gamma_reg_p

    ! Forward product w.r.t. x (a inactive): jv = x^{a-1} e^{-x} * v.
    ! x(1) = x (integration limit), v(1) = direction, jv(1) = product.
    ! a is passed as x(2) (inactive, read-only; harness interface uses x(:)).
    ! Convention: x array = [x_val, a_val] so the harness can vary x(1).
    subroutine gamma_lower_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)   ! x(1)=limit, x(2)=shape a
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: xv, av
        xv = x(1)
        av = x(2)
        ! d/dx gamma_lower(a,x) = x^{a-1} exp(-x)  (DLMF 8.8.1)
        jv(1) = xv**(av - 1.0_dp) * exp(-xv) * v(1)
    end subroutine gamma_lower_jvp

    ! Series expansion of P(a,x) (DLMF 8.11.4).
    ! Converges for all x but is most efficient for x < a+1.
    pure function gamma_reg_p_series(a, x) result(p)
        real(dp), intent(in) :: a, x
        real(dp) :: p

        integer,  parameter :: max_iter = 500
        real(dp) :: term, total
        integer  :: n

        term  = 1.0_dp / a
        total = term
        do n = 1, max_iter
            term  = term * x / (a + real(n, dp))
            total = total + term
            if (abs(term) < abs(total) * epsilon(1.0_dp)) then
                p = total * exp(-x + a * log(x) - log_gamma(a))
                return
            end if
        end do
        error stop "gamma_reg_p_series: series did not converge"
    end function gamma_reg_p_series

    ! Regularized Q(a,x) = 1 - P(a,x) via modified Lentz continued fraction
    ! (DLMF 8.9.2).  Efficient for x >= a+1.
    pure function gamma_reg_q_contfrac(a, x) result(q)
        real(dp), intent(in) :: a, x
        real(dp) :: q

        integer,  parameter :: max_iter = 500
        real(dp), parameter :: tiny_val = 1.0e-300_dp
        real(dp) :: b, c, d, h, an, del
        integer  :: n

        b = x + 1.0_dp - a
        c = 1.0_dp / tiny_val
        d = 1.0_dp / max(b, tiny_val)
        h = d
        do n = 1, max_iter
            an = -real(n, dp) * (real(n, dp) - a)
            b  = b + 2.0_dp
            d  = an * d + b
            if (abs(d) < tiny_val) d = tiny_val
            c  = b + an / c
            if (abs(c) < tiny_val) c = tiny_val
            d   = 1.0_dp / d
            del = d * c
            h   = h * del
            if (abs(del - 1.0_dp) < epsilon(1.0_dp)) then
                q = h * exp(-x + a * log(x) - log_gamma(a))
                return
            end if
        end do
        error stop "gamma_reg_q_contfrac: continued fraction did not converge"
    end function gamma_reg_q_contfrac

end module fortnum_special_gamma
