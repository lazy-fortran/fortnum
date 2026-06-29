module fortnum_special_complex_bessel
    ! Complex-argument Bessel functions for integer order: J_n(z), I_n(z),
    ! K_n(z).  Clean-room from DLMF chapter 10; replaces the AMOS zbesj/zbesi/
    ! zbesk that KiLCA (KAMEL) bundles.
    !
    ! Methods (DLMF 10):
    !   J_n(z)  power series 10.2.2 off the real axis and for small |z|;
    !           Hankel asymptotic 10.17.3 for |z| above the crossover; downward
    !           (Miller) recurrence near the real axis at moderate |z|, where
    !           the oscillatory series cancels.  Backward recurrence is unstable
    !           for J off the real axis, so it is restricted to the near-real
    !           strip where it is well conditioned.
    !   I_n(z)  power series 10.25.2 for small |z|; downward (Miller) recurrence
    !           10.29.1 with the e^z normalization 10.41.4 for moderate |z|;
    !           asymptotic 10.40.1 for large |z|.  KODE scaling factors out
    !           e^{Re z}, mandatory for the conductivity grid where Re z is large.
    !   K_n(z)  K_0 and K_1 from the integral representation 10.32.18,
    !             K_m(z) = int_0^inf e^{-z cosh t} cosh(m t) dt   (Re z > 0),
    !           by an adaptive trapezoid below |z| = k_asym_min and the DLMF
    !           10.40.2 asymptotic with optimal truncation at or above it; the
    !           higher orders follow by the stable upward recurrence 10.29.4.
    !           The integrand of 10.32.18 oscillates as e^{-i Im(z) cosh t}, so
    !           the trapezoid panel count tracks Im(z) sinh(T) and the truncation
    !           T tracks the e^{-Re(z) cosh t} decay; near the imaginary axis the
    !           oscillation count grows with |z|, where the asymptotic takes over.
    !           The scaled form integrates e^{-z(cosh t - 1)} to avoid the
    !           e^{Re z} overflow.
    !
    ! Negative integer order: J_{-n} = (-1)^n J_n, I_{-n} = I_n, K_{-n} = K_n;
    ! handled internally so callers may pass any integer order.
    !
    ! DERIVATIVE POLICY (ad.md sec 1, sec 4): analytic_rule.
    !   Active argument: z (complex).  Inactive: order, scaling flag.
    !   J_n'(z) =  (J_{n-1}(z) - J_{n+1}(z)) / 2   (DLMF 10.6.1)
    !   I_n'(z) =  (I_{n-1}(z) + I_{n+1}(z)) / 2   (DLMF 10.29.2)
    !   K_n'(z) = -(K_{n-1}(z) + K_{n+1}(z)) / 2   (DLMF 10.29.4)
    !   Forward complex products are bessel_*_complex_jvp.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR

    implicit none
    private

    public :: bessel_j_complex
    public :: bessel_i_complex, bessel_i_complex_array
    public :: bessel_k_complex, bessel_k_complex_array
    public :: bessel_j_complex_jvp
    public :: bessel_i_complex_jvp
    public :: bessel_k_complex_jvp

    real(dp), parameter :: pi      = 3.14159265358979323846_dp
    real(dp), parameter :: two_pi  = 6.28318530717958647692_dp

    ! J: Hankel asymptotic above a phase-dependent onset (asym_j_min on the
    ! positive real axis, rising to ~20.6 near arg z = pi where 10.17.3
    ! converges slower); below it the power series (off the real axis it stays
    ! machine-accurate well past this |z|) or, on the near-real strip beyond
    ! miller_j_min where the series cancels, Miller recurrence.
    real(dp), parameter :: asym_j_min   = 13.0_dp
    real(dp), parameter :: miller_j_min = 9.0_dp
    real(dp), parameter :: near_real_ratio = 0.25_dp
    ! I: power series below this |z|, recurrence/asymptotic above.
    real(dp), parameter :: series_i_max = 10.0_dp
    real(dp), parameter :: asym_i_min   = 25.0_dp

    real(dp), parameter :: rescale_big   = 1.0e250_dp
    real(dp), parameter :: rescale_small = 1.0e-250_dp
    complex(dp), parameter :: miller_seed = (1.0e-30_dp, 0.0_dp)

    ! K: DLMF 10.40.2 asymptotic at and above this |z| (machine-accurate with
    ! optimal truncation), adaptive trapezoid of 10.32.18 below it.
    real(dp), parameter :: k_asym_min = 14.0_dp
    ! K trapezoid: floor panel count, points resolving one oscillation
    ! wavelength, panel cap, and the decay margin e^{-k_decay} for truncation.
    integer,  parameter :: k_panels   = 600
    real(dp), parameter :: k_ppw      = 20.0_dp
    integer,  parameter :: k_panel_cap = 1000000
    real(dp), parameter :: k_decay    = 40.0_dp

contains

    ! ------------------------------------------------------------------ J_n(z)

    ! J_n(z), integer order, complex z.  Single value, unscaled (KiLCA KODE=1).
    subroutine bessel_j_complex(order, z, result, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        complex(dp),            intent(out) :: result
        type(fortnum_status_t), intent(out) :: status

        integer  :: n
        real(dp) :: az

        n  = abs(order)
        az = abs(z)
        ! The Hankel asymptotic 10.17.3 converges slower as |arg z| -> pi, so
        ! the |z| onset rises with phase; the power series covers the off-axis
        ! interior (accurate there to |z| >= 30).
        if (az > asym_j_min + 11.0_dp* &
            (min(abs(atan2(aimag(z), real(z, dp))), pi)/pi)**2) then
            result = j_asymptotic(n, z)
        else if (az > miller_j_min .and. &
                abs(aimag(z)) <= near_real_ratio*abs(real(z, dp))) then
            result = j_miller(n, z)
        else
            result = j_series(n, z)
        end if
        if (order < 0 .and. mod(n, 2) == 1) result = -result
        call status_set(status, FORTNUM_OK, "")
    end subroutine bessel_j_complex

    ! DLMF 10.2.2: J_n(z) = (z/2)^n sum_k (-z^2/4)^k / (k! (n+k)!).
    pure function j_series(n, z) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        complex(dp) :: value

        complex(dp) :: pref, q, term, total
        integer     :: j, k

        q    = -0.25_dp*z*z
        pref = (1.0_dp, 0.0_dp)
        do j = 1, n
            pref = pref*(0.5_dp*z)/real(j, dp)
        end do
        term  = (1.0_dp, 0.0_dp)
        total = term
        do k = 1, 400
            term  = term*q/(real(k, dp)*real(n + k, dp))
            total = total + term
            if (abs(term) <= epsilon(1.0_dp)*abs(total) .and. k > 2) exit
        end do
        value = pref*total
    end function j_series

    ! DLMF 10.17.3 Hankel asymptotic for -pi < arg z < pi.
    pure function j_asymptotic(n, z) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        complex(dp) :: value

        complex(dp) :: omega, s_even, s_odd, aval
        real(dp)    :: mu, prev
        integer     :: k

        mu     = 4.0_dp*real(n, dp)**2
        omega  = z - (0.5_dp*real(n, dp) + 0.25_dp)*pi
        s_even = (1.0_dp, 0.0_dp)
        s_odd  = (0.0_dp, 0.0_dp)
        aval   = (1.0_dp, 0.0_dp)
        prev   = huge(1.0_dp)
        do k = 1, 80
            aval = aval*(mu - real(2*k - 1, dp)**2)/(real(k, dp)*8.0_dp*z)
            if (abs(aval) >= prev) exit
            if (mod(k, 2) == 0) then
                s_even = s_even + alt_sign(k/2)*aval
            else
                s_odd  = s_odd  + alt_sign((k - 1)/2)*aval
            end if
            prev = abs(aval)
        end do
        value = sqrt(2.0_dp/(pi*z))*(cos(omega)*s_even - sin(omega)*s_odd)
    end function j_asymptotic

    pure function alt_sign(m) result(s)
        integer, intent(in) :: m
        real(dp) :: s
        s = real(1 - 2*mod(m, 2), dp)
    end function alt_sign

    ! Downward recurrence J_{k-1} = (2k/z) J_k - J_{k+1} (DLMF 10.6.1) with the
    ! normalization 1 = J_0 + 2 sum_{k>=1} J_{2k} (DLMF 10.12.4).  Stable on the
    ! near-real strip; the caller restricts it there.
    pure function j_miller(n, z) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        complex(dp) :: value

        complex(dp) :: p_hi, p_cur, p_lo, ssum, j0
        integer     :: m, k

        m = n + int(40.0_dp + 8.0_dp*sqrt(abs(z)))
        if (mod(m, 2) == 1) m = m + 1
        p_hi  = (0.0_dp, 0.0_dp)
        p_cur = miller_seed
        ssum  = (0.0_dp, 0.0_dp)
        value = (0.0_dp, 0.0_dp)
        do k = m, 1, -1
            p_lo = (2.0_dp*real(k, dp)/z)*p_cur - p_hi
            if (mod(k, 2) == 0) ssum = ssum + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_big) then
                p_lo  = p_lo*rescale_small
                p_cur = p_cur*rescale_small
                ssum  = ssum*rescale_small
                value = value*rescale_small
            end if
            p_hi  = p_cur
            p_cur = p_lo
            if (k - 1 == n) value = p_cur
        end do
        j0 = p_cur
        if (n == 0) value = j0
        value = value/(j0 + ssum)
    end function j_miller

    ! ------------------------------------------------------------------ I_n(z)

    ! I_n(z), integer order, complex z.  scaled = .true. returns e^{-Re z} I.
    subroutine bessel_i_complex(order, z, scaled, result, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(out) :: result
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: seq(0:abs(order))

        call bessel_i_complex_array(0, abs(order) + 1, z, scaled, seq, status)
        result = seq(abs(order))
    end subroutine bessel_i_complex

    ! Contiguous order sequence I_{order0}(z) .. I_{order0+nseq-1}(z).  Matches
    ! the AMOS (FNU, N) convention KiLCA uses; order0 >= 0.
    subroutine bessel_i_complex_array(order0, nseq, z, scaled, result, status)
        integer,                intent(in)  :: order0
        integer,                intent(in)  :: nseq
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(out) :: result(nseq)
        type(fortnum_status_t), intent(out) :: status

        integer     :: ntop, k
        complex(dp) :: full(0:order0 + nseq - 1)

        if (nseq < 1 .or. order0 < 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "i array bad range")
            result = (0.0_dp, 0.0_dp)
            return
        end if
        ntop = order0 + nseq - 1
        if (abs(z) <= series_i_max) then
            do k = 0, ntop
                full(k) = i_series(k, z)
            end do
            if (scaled) full = full*exp(-real(z, dp))
        else
            call i_miller_array(ntop, z, scaled, full)
        end if
        result = full(order0:ntop)
        call status_set(status, FORTNUM_OK, "")
    end subroutine bessel_i_complex_array

    ! DLMF 10.25.2: I_n(z) = (z/2)^n sum_k (z^2/4)^k / (k! (n+k)!).
    pure function i_series(n, z) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        complex(dp) :: value

        complex(dp) :: pref, q, term, total
        integer     :: j, k

        q    = 0.25_dp*z*z
        pref = (1.0_dp, 0.0_dp)
        do j = 1, n
            pref = pref*(0.5_dp*z)/real(j, dp)
        end do
        term  = (1.0_dp, 0.0_dp)
        total = term
        do k = 1, 400
            term  = term*q/(real(k, dp)*real(n + k, dp))
            total = total + term
            if (abs(term) <= epsilon(1.0_dp)*abs(total) .and. k > 2) exit
        end do
        value = pref*total
    end function i_series

    ! Downward recurrence I_{k-1} = I_{k+1} + (2k/z) I_k (DLMF 10.29.1) with the
    ! normalization e^z = I_0 + 2 sum_{k>=1} I_k (DLMF 10.41.4).  Stable for I.
    subroutine i_miller_array(ntop, z, scaled, full)
        integer,     intent(in)  :: ntop
        complex(dp), intent(in)  :: z
        logical,     intent(in)  :: scaled
        complex(dp), intent(out) :: full(0:ntop)

        complex(dp) :: p_hi, p_cur, p_lo, ssum, norm
        integer     :: m, k

        m     = ntop + int(40.0_dp + 12.0_dp*sqrt(abs(z)))
        p_hi  = (0.0_dp, 0.0_dp)
        p_cur = miller_seed
        ssum  = (0.0_dp, 0.0_dp)
        full  = (0.0_dp, 0.0_dp)
        do k = m, 1, -1
            p_lo = p_hi + (2.0_dp*real(k, dp)/z)*p_cur
            ssum = ssum + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_big) then
                p_lo  = p_lo*rescale_small
                p_cur = p_cur*rescale_small
                ssum  = ssum*rescale_small
                if (k <= ntop) full(k:ntop) = full(k:ntop)*rescale_small
            end if
            p_hi  = p_cur
            p_cur = p_lo
            if (k - 1 <= ntop) full(k - 1) = p_cur
        end do
        ssum = ssum + p_cur
        if (scaled) then
            norm = exp(z - real(z, dp))/ssum
        else
            norm = exp(z)/ssum
        end if
        full = full*norm
    end subroutine i_miller_array

    ! ------------------------------------------------------------------ K_n(z)

    ! K_n(z), integer order, complex z with Re z > 0.  scaled returns e^{z} K.
    subroutine bessel_k_complex(order, z, scaled, result, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(out) :: result
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: seq(0:abs(order))

        call bessel_k_complex_array(0, abs(order) + 1, z, scaled, seq, status)
        result = seq(abs(order))
    end subroutine bessel_k_complex

    ! Contiguous order sequence K_{order0}(z) .. K_{order0+nseq-1}(z).
    subroutine bessel_k_complex_array(order0, nseq, z, scaled, result, status)
        integer,                intent(in)  :: order0
        integer,                intent(in)  :: nseq
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(out) :: result(nseq)
        type(fortnum_status_t), intent(out) :: status

        integer     :: ntop, k
        complex(dp) :: full(0:order0 + nseq - 1)

        if (nseq < 1 .or. order0 < 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "k array bad range")
            result = (0.0_dp, 0.0_dp)
            return
        end if
        if (real(z, dp) <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "k needs Re z > 0")
            result = (0.0_dp, 0.0_dp)
            return
        end if
        ntop = order0 + nseq - 1
        full(0) = k_value(0, z, scaled)
        if (ntop >= 1) full(1) = k_value(1, z, scaled)
        ! DLMF 10.29.4: K_{m+1}(z) = K_{m-1}(z) + (2m/z) K_m(z); stable upward.
        do k = 1, ntop - 1
            full(k + 1) = full(k - 1) + (2.0_dp*real(k, dp)/z)*full(k)
        end do
        result = full(order0:ntop)
        call status_set(status, FORTNUM_OK, "")
    end subroutine bessel_k_complex_array

    ! K_n(z) for the base orders n = 0, 1: DLMF 10.40.2 asymptotic at large |z|
    ! (machine-accurate with optimal truncation), adaptive 10.32.18 trapezoid
    ! below.  Higher orders come from the upward recurrence in the caller.
    pure function k_value(n, z, scaled) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        logical,     intent(in) :: scaled
        complex(dp) :: value

        if (abs(z) >= k_asym_min) then
            value = k_asymptotic(n, z, scaled)
        else
            value = k_trapezoid(n, z, scaled)
        end if
    end function k_value

    ! DLMF 10.40.2: K_n(z) ~ sqrt(pi/2z) e^{-z} sum_k a_k(n)/z^k for large |z|,
    ! a_k(n) = prod_{j=1}^{k} (4n^2 - (2j-1)^2) / (k! 8^k).  Optimal truncation
    ! (stop before the smallest term grows) reaches machine precision across the
    ! right half-plane once |z| >= k_asym_min.  scaled drops the e^{-z} factor.
    pure function k_asymptotic(n, z, scaled) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        logical,     intent(in) :: scaled
        complex(dp) :: value

        complex(dp) :: pref, term, total
        real(dp)    :: mu, prev, at
        integer     :: k

        mu   = 4.0_dp*real(n, dp)**2
        pref = sqrt(pi/(2.0_dp*z))
        if (.not. scaled) pref = pref*exp(-z)
        term  = (1.0_dp, 0.0_dp)
        total = term
        prev  = huge(1.0_dp)
        do k = 1, 300
            term = term*(mu - real(2*k - 1, dp)**2)/(real(k, dp)*8.0_dp*z)
            at   = abs(term)
            if (at > prev) exit
            total = total + term
            prev  = at
        end do
        value = pref*total
    end function k_asymptotic

    ! DLMF 10.32.18: K_n(z) = int_0^inf e^{-z cosh t} cosh(n t) dt, Re z > 0.
    ! Trapezoid on [0, T].  The integrand decays like e^{-Re(z) cosh t}, fixing
    ! T from the k_decay margin, and oscillates like e^{-i Im(z) cosh t}, so the
    ! panel count resolves Im(z) sinh(T) at k_ppw points per wavelength.  scaled
    ! removes the e^{-z} factor (cosh t - 1 in the exponent) to avoid overflow.
    pure function k_trapezoid(n, z, scaled) result(value)
        integer,     intent(in) :: n
        complex(dp), intent(in) :: z
        logical,     intent(in) :: scaled
        complex(dp) :: value

        complex(dp) :: acc, integrand
        real(dp)    :: rz, tmax, h, t, weight, shift, phase
        integer     :: k, panels

        rz   = max(real(z, dp), 1.0e-3_dp)
        tmax = k_trunc(rz, n)
        phase = abs(aimag(z))*sinh(tmax)
        panels = max(k_panels, &
            min(int(k_ppw*phase/two_pi) + 1, k_panel_cap))
        h    = tmax/real(panels, dp)
        shift = 0.0_dp
        if (scaled) shift = 1.0_dp
        acc = (0.0_dp, 0.0_dp)
        do k = 0, panels
            t = real(k, dp)*h
            integrand = exp(-z*(cosh(t) - shift))*cosh(real(n, dp)*t)
            weight = 1.0_dp
            if (k == 0 .or. k == panels) weight = 0.5_dp
            acc = acc + weight*integrand
        end do
        value = acc*h
    end function k_trapezoid

    ! Upper limit T for the 10.32.18 trapezoid: the decay exponent
    ! g(t) = Re(z) cosh t - n t must drop k_decay below its peak at the saddle
    ! t* = asinh(n/Re z).  Newton on g(T) - g(t*) = k_decay.
    pure function k_trunc(rz, n) result(tmax)
        real(dp), intent(in) :: rz
        integer,  intent(in) :: n
        real(dp) :: tmax

        real(dp) :: tpk, gpk, target, f, fp
        integer  :: it

        tpk    = asinh(real(n, dp)/rz)
        gpk    = rz*cosh(tpk) - real(n, dp)*tpk
        target = gpk + k_decay
        tmax   = log(2.0_dp*target/rz + 2.0_dp) + 1.0_dp
        do it = 1, 80
            f  = rz*cosh(tmax) - real(n, dp)*tmax - target
            fp = rz*sinh(tmax) - real(n, dp)
            if (abs(fp) < 1.0e-12_dp) fp = sign(1.0e-12_dp, fp)
            tmax = tmax - f/fp
            if (tmax < 0.05_dp) tmax = 0.05_dp
        end do
        tmax = max(tmax, 2.0_dp)
    end function k_trunc

    ! ----------------------------------------------------- analytic derivatives

    ! J_n'(z) v = (J_{n-1}(z) - J_{n+1}(z))/2 * v   (DLMF 10.6.1).
    subroutine bessel_j_complex_jvp(order, z, v, jv, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        complex(dp),            intent(in)  :: v
        complex(dp),            intent(out) :: jv
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: jm, jp
        type(fortnum_status_t) :: s1, s2

        call bessel_j_complex(order - 1, z, jm, s1)
        call bessel_j_complex(order + 1, z, jp, s2)
        jv = 0.5_dp*(jm - jp)*v
        if (s1%code /= FORTNUM_OK) then
            status = s1
        else
            status = s2
        end if
    end subroutine bessel_j_complex_jvp

    ! I_n'(z) v = (I_{n-1}(z) + I_{n+1}(z))/2 * v   (DLMF 10.29.2).
    subroutine bessel_i_complex_jvp(order, z, scaled, v, jv, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(in)  :: v
        complex(dp),            intent(out) :: jv
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: im, ip
        type(fortnum_status_t) :: s1, s2

        call bessel_i_complex(abs(order - 1), z, scaled, im, s1)
        call bessel_i_complex(order + 1, z, scaled, ip, s2)
        jv = 0.5_dp*(im + ip)*v
        if (s1%code /= FORTNUM_OK) then
            status = s1
        else
            status = s2
        end if
    end subroutine bessel_i_complex_jvp

    ! K_n'(z) v = -(K_{n-1}(z) + K_{n+1}(z))/2 * v   (DLMF 10.29.4).
    subroutine bessel_k_complex_jvp(order, z, scaled, v, jv, status)
        integer,                intent(in)  :: order
        complex(dp),            intent(in)  :: z
        logical,                intent(in)  :: scaled
        complex(dp),            intent(in)  :: v
        complex(dp),            intent(out) :: jv
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: km, kp
        type(fortnum_status_t) :: s1, s2

        call bessel_k_complex(abs(order - 1), z, scaled, km, s1)
        call bessel_k_complex(order + 1, z, scaled, kp, s2)
        jv = -0.5_dp*(km + kp)*v
        if (s1%code /= FORTNUM_OK) then
            status = s1
        else
            status = s2
        end if
    end subroutine bessel_k_complex_jvp

end module fortnum_special_complex_bessel
