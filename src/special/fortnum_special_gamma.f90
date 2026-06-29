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
    !           iteration branches. Forward products gamma_lower_jvp_da and
    !           gamma_reg_p_jvp and the reverse product gamma_reg_p_grad
    !           implement the d/da pieces (issue #40).
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
    ! analytic_rule derivatives (ad.md §2). Packing x = [x_val, a_val],
    ! v = [dx, da] so the square test harness can vary each input.
    !   gamma_lower_jvp:    d/dx gamma_lower(a,x) (DLMF 8.8.1; a inactive).
    !   gamma_lower_jvp_da: full forward product in [x, a]. The a-direction
    !                       uses d/da gamma_lower = x^a e^{-x} (ln x . S - T)
    !                       with P-series sums S, T; the Gamma(a) factor
    !                       cancels psi(a), so no digamma is needed here.
    !   gamma_reg_p_jvp:    full forward product for P(a,x) = gamma_lower/Gamma(a).
    !   gamma_reg_p_grad:   reverse product (scalar output) for P.
    ! d/da P keeps psi(a): dP/da = A[(ln x - psi(a)) S - T], A = x^a e^{-x}/Gamma(a).
    public :: gamma_lower_jvp
    public :: gamma_lower_jvp_da
    public :: gamma_reg_p_jvp
    public :: gamma_reg_p_grad

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
        real(dp), intent(in)  :: x(:) ! x(1)=limit, x(2)=shape a
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: xv, av
        xv = x(1)
        av = x(2)
        ! d/dx gamma_lower(a,x) = x^{a-1} exp(-x)  (DLMF 8.8.1)
        jv(1) = xv**(av - 1.0_dp) * exp(-xv) * v(1)
    end subroutine gamma_lower_jvp

    ! Forward product for gamma_lower in both inputs. x = [x_val, a_val],
    ! v = [dx, da]; jv(1) = d/dx g * dx + d/da g * da.
    !   d/dx g = x^{a-1} e^{-x}                 (DLMF 8.8.1)
    !   d/da g = x^a e^{-x} (ln x . S - T)
    ! with S = sum_n term_n, T = sum_n term_n H_n, term_n = x^n/prod_{k=0}^n(a+k),
    ! H_n = sum_{k=0}^n 1/(a+k). The Gamma(a) factor of P cancels psi(a) here.
    subroutine gamma_lower_jvp_da(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: xv, av, s, t, dgdx, dgda
        xv = x(1)
        av = x(2)
        call gamma_p_series_sums(av, xv, s, t)
        dgdx = xv**(av - 1.0_dp) * exp(-xv)
        dgda = xv**av * exp(-xv) * (log(xv) * s - t)
        jv = 0.0_dp
        jv(1) = dgdx * v(1) + dgda * v(2)
    end subroutine gamma_lower_jvp_da

    ! Forward product for P(a,x) = gamma_lower(a,x)/Gamma(a). Packing as above.
    !   d/dx P = x^{a-1} e^{-x} / Gamma(a)
    !   d/da P = A[(ln x - psi(a)) S - T],  A = x^a e^{-x} / Gamma(a)
    subroutine gamma_reg_p_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: dpdx, dpda
        call gamma_reg_p_partials(x(1), x(2), dpdx, dpda)
        jv = 0.0_dp
        jv(1) = dpdx * v(1) + dpda * v(2)
    end subroutine gamma_reg_p_jvp

    ! Reverse product for the scalar map P(a,x). u in R^1, jtu in R^2:
    !   jtu(1) = d/dx P * u(1),  jtu(2) = d/da P * u(1).
    subroutine gamma_reg_p_grad(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        real(dp) :: dpdx, dpda
        call gamma_reg_p_partials(x(1), x(2), dpdx, dpda)
        jtu(1) = dpdx * u(1)
        jtu(2) = dpda * u(1)
    end subroutine gamma_reg_p_grad

    ! Partials of P(a,x) at (xv, av) shared by the P forward and reverse
    ! products. dP/da carries the digamma psi(a) (DLMF 5.2.2) that the
    ! Gamma(a) normalization introduces.
    pure subroutine gamma_reg_p_partials(xv, av, dpdx, dpda)
        real(dp), intent(in)  :: xv, av
        real(dp), intent(out) :: dpdx, dpda
        real(dp) :: s, t, ga, psi, a_pref
        call gamma_p_series_sums(av, xv, s, t)
        ga = gamma(av)
        psi = digamma(av)
        a_pref = xv**av * exp(-xv) / ga
        dpdx = xv**(av - 1.0_dp) * exp(-xv) / ga
        dpda = a_pref * ((log(xv) - psi) * s - t)
    end subroutine gamma_reg_p_partials

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

    ! Sums driving d/da of the lower incomplete gamma and P (DLMF 8.11.4):
    !   term_n = x^n / (a(a+1)...(a+n)),  H_n = sum_{k=0}^n 1/(a+k)
    !   S = sum_n term_n,  T = sum_n term_n H_n.
    ! Both series converge for all x > 0, a > 0; the terms are positive.
    pure subroutine gamma_p_series_sums(a, x, s, t)
        real(dp), intent(in)  :: a, x
        real(dp), intent(out) :: s, t

        integer,  parameter :: max_iter = 500
        real(dp) :: term, hsum
        integer  :: n

        term = 1.0_dp / a
        hsum = 1.0_dp / a
        s = term
        t = term * hsum
        do n = 1, max_iter
            term = term * x / (a + real(n, dp))
            hsum = hsum + 1.0_dp / (a + real(n, dp))
            s = s + term
            t = t + term * hsum
            if (abs(term) < abs(s) * epsilon(1.0_dp)) then
                if (abs(term * hsum) < abs(t) * epsilon(1.0_dp)) return
            end if
        end do
        error stop "gamma_p_series_sums: series did not converge"
    end subroutine gamma_p_series_sums

    ! Digamma psi(a) = d/da log Gamma(a) for a > 0 (DLMF 5.2.2): recurrence
    ! psi(a) = psi(a+1) - 1/a lifts the argument past 10, then the asymptotic
    ! series psi(x) ~ ln x - 1/(2x) - 1/(12 x^2) + 1/(120 x^4) - 1/(252 x^6)
    ! + 1/(240 x^8) - 1/(132 x^10), accurate to ~1e-13 once x >= 10.
    pure function digamma(a) result(psi)
        real(dp), intent(in) :: a
        real(dp) :: psi
        real(dp) :: x, f, f2
        psi = 0.0_dp
        x = a
        do while (x < 10.0_dp)
            psi = psi - 1.0_dp / x
            x = x + 1.0_dp
        end do
        f  = 1.0_dp / x
        f2 = f * f
        psi = psi + log(x) - 0.5_dp * f - f2 * (1.0_dp/12.0_dp &
            - f2 * (1.0_dp/120.0_dp - f2 * (1.0_dp/252.0_dp &
            - f2 * (1.0_dp/240.0_dp - f2/132.0_dp))))
    end function digamma

end module fortnum_special_gamma
