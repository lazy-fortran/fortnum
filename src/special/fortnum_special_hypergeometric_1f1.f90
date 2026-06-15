module fortnum_special_hypergeometric_1f1
    ! Confluent hypergeometric function 1F1(a;b;z) (Kummer's function M),
    ! complex a, b, z.  Clean-room implementation from DLMF chapter 13.
    !
    ! Definition (DLMF 13.2.2):
    !   M(a,b,z) = sum_{k>=0} (a)_k / (b)_k * z^k / k!
    ! where (x)_k is the Pochhammer symbol.  Holomorphic in a and z; has poles
    ! in b at the non-positive integers b = 0,-1,-2,...
    !
    ! Replaces the >6 self-tuned algorithm variants of MEPHIT src/hyper1F1.c and
    ! KAMEL KiLCA math/hyper/hyper1F1.cpp (Kummer series, inverse continued
    ! fractions, Levin-u acceleration, quadrature) with one robust selector that
    ! exposes a single converged primal.  Those consumers fix a = 1 and pass
    ! complex b = 1 + t2, z = t1 from the FLR / plasma-dispersion sums, with
    ! |z| moderate and b possibly near integers; hyperg_1f1_a1 is the matching
    ! specialization.
    !
    ! Method selection (internal trace; only the converged value is exposed):
    !   Re(z) < 0          Kummer transformation M(a,b,z)=e^z M(b-a,b,-z)
    !                      (DLMF 13.2.39) so the series argument has Re >= 0 and
    !                      the terms do not alternate with growing magnitude.
    !   |z| <= z_series    Taylor series (DLMF 13.2.2) via the term ratio
    !                      t_{k+1}/t_k = (a+k) z / ((b+k)(k+1)).
    !   |z| >  z_series    Asymptotic expansion for large z (DLMF 13.7.2):
    !                      M(a,b,z) ~ Gamma(b)/Gamma(a) e^z z^{a-b} 2F0(b-a,1-a;;1/z).
    !
    ! Derivative policy (docs/design/ad.md sec. 1, 4): analytic_rule.
    !   d/dz M(a,b,z) = (a/b) M(a+1,b+1,z)  (DLMF 13.3.15).
    !   Active argument: z.  Inactive arguments: a, b (order/shape selectors).
    !   The harness exercises the real and imaginary parts of z as a 2-vector
    !   map R^2 -> R^2 (the Cauchy-Riemann Jacobian of the analytic map), with
    !   hyperg_1f1_a1_jvp / hyperg_1f1_a1_vjp the forward and reverse products.
    !
    ! References: DLMF 13.2 (series), 13.2.39 (Kummer transformation),
    !   13.3.15 (derivative), 13.7.2 (large-z asymptotics);
    !   Abramowitz & Stegun chapter 13.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: hyperg_1f1
    public :: hyperg_1f1_a1
    public :: hyperg_1f1m_a1
    ! analytic_rule derivatives (ad.md sec. 2): d/dz M = (a/b) M(a+1,b+1,z).
    !   hyperg_1f1_a1_jvp: forward product J(z) v, z = (Re z, Im z).
    !   hyperg_1f1_a1_vjp: reverse product J(z)^T u.
    public :: hyperg_1f1_a1_jvp
    public :: hyperg_1f1_a1_vjp

    ! Series/asymptotic crossover on |z|.  Below this the Taylor series
    ! converges in a bounded number of terms at float64 precision; above it
    ! the large-z asymptotic series is the accurate and cheaper branch.
    real(dp), parameter :: z_series_max = 60.0_dp

    ! Distance from a non-positive integer within which b is treated as a pole.
    real(dp), parameter :: pole_guard = 1.0e-12_dp

    integer,  parameter :: series_max_terms = 5000
    integer,  parameter :: asy_max_terms    = 60

contains

    ! 1F1(a;b;z) for complex a, b, z.  Returns the converged value in result
    ! and a status: FORTNUM_DOMAIN_ERROR when b is at (or within pole_guard of)
    ! a non-positive integer, FORTNUM_CONVERGENCE_ERROR when no branch reaches
    ! float64 agreement, FORTNUM_OK otherwise.
    subroutine hyperg_1f1(a, b, z, result, status)
        complex(dp), intent(in)  :: a, b, z
        complex(dp), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: m
        logical     :: ok

        if (b_near_nonpositive_integer(b)) then
            result = cmplx(0.0_dp, 0.0_dp, dp)
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "hyperg_1f1: b at or near a non-positive integer (pole)")
            return
        end if

        if (real(z, dp) < 0.0_dp) then
            ! Kummer transformation: M(a,b,z) = e^z M(b-a,b,-z) (DLMF 13.2.39).
            call eval_positive_re(b - a, b, -z, m, ok)
            m = exp(z) * m
        else
            call eval_positive_re(a, b, z, m, ok)
        end if

        result = m
        if (ok) then
            call status_set(status, FORTNUM_OK, "")
        else
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "hyperg_1f1: no branch reached float64 tolerance")
        end if
    end subroutine hyperg_1f1

    ! Specialization a = 1 used by the MEPHIT/KiLCA consumers.
    subroutine hyperg_1f1_a1(b, z, result, status)
        complex(dp), intent(in)  :: b, z
        complex(dp), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        call hyperg_1f1(cmplx(1.0_dp, 0.0_dp, dp), b, z, result, status)
    end subroutine hyperg_1f1_a1

    ! Modified form for a = 1 used by the KiLCA FLR reduction, which writes
    !   M(1,b,z) = 1 + z/b + z^2/(b(b+1)) * (1 + F11m)
    ! and consumes F11m, not M.  Reconstructing F11m from M cancels two ~1
    ! quantities and divides by z^2 (tiny at small z), amplifying the float64
    ! error in M by |M| |b(b+1)/z^2| / |F11m|.  The closed form
    !   F11m = M(1, b+2, z) - 1
    ! (a = 1 telescoping of (b)_k against the modified-0 reduction) computes
    ! the modified value directly with no cancellation.
    subroutine hyperg_1f1m_a1(b, z, result, status)
        complex(dp), intent(in)  :: b, z
        complex(dp), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        complex(dp) :: m
        call hyperg_1f1(cmplx(1.0_dp, 0.0_dp, dp), b + 2.0_dp, z, m, status)
        result = m - 1.0_dp
    end subroutine hyperg_1f1m_a1

    ! Evaluate M(a,b,z) assuming Re(z) >= 0: pick the Taylor series for
    ! moderate |z| and the large-z asymptotic series otherwise.
    !
    ! The large-z asymptotic (DLMF 13.7.2) is valid only when |z| dominates
    ! the order: its 2F0 tail carries (b-a)_s (1-a)_s / z^s, which diverges
    ! when |b| is comparable to |z|.  The Taylor series, by contrast, has
    ! ratio (a+k) z / ((b+k)(k+1)), so |b| >= |z| keeps the terms contracting
    ! and machine-accurate without large intermediate magnitudes.  The KiLCA
    ! modified-form consumer evaluates M(1, b+2, z) with b ~ 1 + z, i.e.
    ! |b| > |z| at large z; routing it through the asymptotic branch (the old
    ! |z|-only gate) returned a wholly wrong value.  Take the asymptotic
    ! branch only when |z| is large AND |b| < |z| (the small-order regime it
    ! is built for).
    subroutine eval_positive_re(a, b, z, m, ok)
        complex(dp), intent(in)  :: a, b, z
        complex(dp), intent(out) :: m
        logical,     intent(out) :: ok

        if (abs(z) <= z_series_max .or. abs(b) >= abs(z)) then
            call series_1f1(a, b, z, m, ok)
        else
            call asymptotic_1f1(a, b, z, m, ok)
        end if
    end subroutine eval_positive_re

    ! Taylor series M(a,b,z) = sum_k (a)_k/(b)_k z^k/k! via the term ratio
    !   term_{k+1} = term_k * (a+k) z / ((b+k)(k+1))   (DLMF 13.2.2).
    ! Converges absolutely; the loop stops when the running term is below the
    ! float64 epsilon relative to the partial sum and at least two further
    ! terms are negligible, guarding the alternating-sign transient.
    subroutine series_1f1(a, b, z, m, ok)
        complex(dp), intent(in)  :: a, b, z
        complex(dp), intent(out) :: m
        logical,     intent(out) :: ok

        complex(dp) :: term, total
        real(dp)    :: tol
        integer     :: k, small_run

        term      = cmplx(1.0_dp, 0.0_dp, dp)
        total     = term
        ok        = .false.
        small_run = 0
        tol       = epsilon(1.0_dp)

        do k = 0, series_max_terms - 1
            term  = term * (a + real(k, dp)) * z &
                / ((b + real(k, dp)) * real(k + 1, dp))
            total = total + term
            if (abs(term) <= tol * abs(total)) then
                small_run = small_run + 1
                if (small_run >= 2) then
                    ok = .true.
                    exit
                end if
            else
                small_run = 0
            end if
        end do

        m = total
    end subroutine series_1f1

    ! Large-z asymptotic expansion for Re(z) > 0 (DLMF 13.7.2):
    !   M(a,b,z) ~ Gamma(b)/Gamma(a) e^z z^{a-b} sum_s (b-a)_s (1-a)_s / s! z^{-s}
    ! The 2F0 tail is divergent; it is summed until the terms stop shrinking
    ! (optimal truncation).  Convergence is reported only when the smallest
    ! term reached float64 epsilon relative to the running sum.
    subroutine asymptotic_1f1(a, b, z, m, ok)
        complex(dp), intent(in)  :: a, b, z
        complex(dp), intent(out) :: m
        logical,     intent(out) :: ok

        complex(dp) :: term, tail, prefac
        real(dp)    :: prev_mag, mag, tol
        integer     :: s

        term     = cmplx(1.0_dp, 0.0_dp, dp)
        tail     = term
        prev_mag = abs(term)
        ok       = .false.
        tol      = epsilon(1.0_dp)

        do s = 0, asy_max_terms - 1
            term = term * (b - a + real(s, dp)) * (1.0_dp - a + real(s, dp)) &
                / (real(s + 1, dp) * z)
            mag = abs(term)
            if (mag > prev_mag) exit          ! optimal-truncation reached
            tail = tail + term
            if (mag <= tol * abs(tail)) then
                ok = .true.
                exit
            end if
            prev_mag = mag
        end do

        ! Gamma(b)/Gamma(a) e^z z^{a-b}.  log_gamma keeps the ratio finite for
        ! the moderate complex orders the consumers use.
        prefac = exp(clog_gamma(b) - clog_gamma(a) + z + (a - b) * log(z))
        m = prefac * tail
    end subroutine asymptotic_1f1

    ! .true. when b is at or within pole_guard of a non-positive integer, where
    ! M(a,b,z) has a pole in b.
    pure function b_near_nonpositive_integer(b) result(near)
        complex(dp), intent(in) :: b
        logical :: near
        real(dp) :: br, nearest

        near = .false.
        if (abs(aimag(b)) > pole_guard) return
        br = real(b, dp)
        if (br > 0.5_dp) return
        nearest = anint(br)
        if (nearest <= 0.0_dp .and. abs(br - nearest) <= pole_guard) near = .true.
    end function b_near_nonpositive_integer

    ! Complex log-Gamma via the Lanczos approximation (g = 7, n = 9), reflecting
    ! across Re < 0.5 with the reflection formula.  Sufficient accuracy for the
    ! Gamma(b)/Gamma(a) prefactor of the large-z branch at float64.
    pure recursive function clog_gamma(zin) result(lg)
        complex(dp), intent(in) :: zin
        complex(dp) :: lg

        real(dp), parameter :: g = 7.0_dp
        real(dp), parameter :: pi = 3.14159265358979324_dp
        real(dp), parameter :: c(0:8) = [ &
            0.99999999999980993_dp, &
            676.5203681218851_dp, &
            -1259.1392167224028_dp, &
            771.32342877765313_dp, &
            -176.61502916214059_dp, &
            12.507343278686905_dp, &
            -0.13857109526572012_dp, &
            9.9843695780195716e-6_dp, &
            1.5056327351493116e-7_dp]
        complex(dp) :: z, x, t, sum_c
        integer :: i

        if (real(zin, dp) < 0.5_dp) then
            ! Reflection: log Gamma(z) = log(pi/sin(pi z)) - log Gamma(1-z).
            lg = log(pi / sin(pi * zin)) - clog_gamma(1.0_dp - zin)
            return
        end if

        z = zin - 1.0_dp
        x = cmplx(c(0), 0.0_dp, dp)
        do i = 1, 8
            x = x + c(i) / (z + real(i, dp))
        end do
        t = z + g + 0.5_dp
        sum_c = 0.5_dp * log(2.0_dp * pi) + (z + 0.5_dp) * log(t) - t + log(x)
        lg = sum_c
    end function clog_gamma

    ! Forward product for a = 1: J(z) v with z packed as (Re z, Im z), v its
    ! tangent, jv the tangent of (Re M, Im M).  The analytic derivative is
    ! d/dz M(1,b,z) = (1/b) M(2,b+1,z) (DLMF 13.3.15); the real Jacobian is the
    ! Cauchy-Riemann matrix [[u_x,-v_x],[v_x,u_x]] with u_x+i v_x = dM/dz.
    subroutine hyperg_1f1_a1_jvp(z, b, v, jv)
        real(dp), intent(in)  :: z(:)   ! z(1)=Re z, z(2)=Im z
        complex(dp), intent(in) :: b
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: deriv, zc
        zc = cmplx(z(1), z(2), dp)
        deriv = hyperg_1f1_a1_deriv(b, zc)
        jv(1) = real(deriv, dp) * v(1) - aimag(deriv) * v(2)
        jv(2) = aimag(deriv) * v(1) + real(deriv, dp) * v(2)
    end subroutine hyperg_1f1_a1_jvp

    ! Reverse product J(z)^T u for the same map.  J^T transposes the
    ! Cauchy-Riemann matrix, i.e. flips the sign of the off-diagonal term.
    subroutine hyperg_1f1_a1_vjp(z, b, u, jtu)
        real(dp), intent(in)  :: z(:)
        complex(dp), intent(in) :: b
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        complex(dp) :: deriv, zc
        zc = cmplx(z(1), z(2), dp)
        deriv = hyperg_1f1_a1_deriv(b, zc)
        jtu(1) = real(deriv, dp) * u(1) + aimag(deriv) * u(2)
        jtu(2) = -aimag(deriv) * u(1) + real(deriv, dp) * u(2)
    end subroutine hyperg_1f1_a1_vjp

    ! d/dz M(1,b,z) = (1/b) M(2,b+1,z)  (DLMF 13.3.15 with a = 1).
    function hyperg_1f1_a1_deriv(b, z) result(d)
        complex(dp), intent(in) :: b, z
        complex(dp) :: d, m
        type(fortnum_status_t) :: st
        call hyperg_1f1(cmplx(2.0_dp, 0.0_dp, dp), b + 1.0_dp, z, m, st)
        d = m / b
    end function hyperg_1f1_a1_deriv

end module fortnum_special_hypergeometric_1f1
