!> Radix-2 complex FFT with certified twiddle factors and a rigorous,
!> non-asymptotic a priori error bound for FFT-based convolution.
!>
!> Twiddles. Each w_k = exp(-2 pi i k/n) is enclosed in binary128 interval
!> arithmetic: k is reduced exactly by quarter- and eighth-turn symmetries to
!> an angle phi in [0, pi/4], cos phi and sin phi are Taylor polynomials of
!> degree 33 and 32 with Lagrange remainders (all derivatives bounded by 1),
!> and pi is enclosed by the binary128 value nearest pi widened by one ulp.
!> The stored binary64 twiddle is the enclosure midpoint rounded to binary64,
!> and mu = max_k abs(w_stored - w_exact) is bounded from the enclosure
!> radii. No libm cos/sin enters; mu is certified (about 0.7 u) instead of
!> assumed.
!>
!> FFT bound (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed.,
!> Theorem 24.2). For y = F_n x by radix-2 FFT with per-twiddle error mu:
!>   norm2(y_hat - y) <= kappa norm2(y),  kappa = p eta/(1 - p eta),
!>   eta = mu + gamma_4 (sqrt 2 + mu),  p = log2 n,
!> with no O(u^2) term dropped.
!>
!> Convolution. c = IFFT(FFT(a) .* FFT(b)) with n larger than the linear
!> convolution length (no wraparound). With A = F a, B = F b exact,
!> norm2(A) = sqrt(n) norm2(a), and the pointwise product rounded with
!> rho = sqrt 5 u/(1 - sqrt 5 u) (Brent, Percival, Zimmermann 2007):
!>   norm2(p_hat - P) <= n norm2(a) norm2(b) sigma,
!>   sigma = 2 kappa + kappa^2 + rho (1 + kappa)^2,
!> and the inverse transform (Higham bound plus u for the 1/n scaling) gives
!>   norm2(c_hat - c) <= sqrt(n) [(kappa + u)(1 + sigma) + sigma]
!>                       norm2(a) norm2(b) =: Gamma(n) norm2(a) norm2(b),
!> which bounds every entry. l1 norms are accepted since norm2 <= l1. The 2D
!> transform (rows of length pn, then columns of length pm) composes as
!> kappa2 = kappa_m + kappa_n + kappa_m kappa_n with length pm pn. Every
!> constant is rounded upward (denominators downward) with fortnum_rounding.
!> This is the kinetic-compression kc_cfft derivation with its l1 sums and
!> rho quotient rounded upward.
module fortnum_cfft_rigorous
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    use fortnum_rounding, only: round_up, mul_up, add_up, div_up, &
        sub_down, sqrt_up, sum_up, gamma_up, unit_roundoff
    use fortnum_interval_qp, only: qp, qinterval_t, qinterval, qinterval_pi, &
        qsqr, qmid, qrad, operator(+), operator(-), operator(*), operator(/)
    implicit none
    private

    public :: rigorous_fft_plan_t, rigorous_fft_plan_init, rigorous_fft_apply
    public :: rigorous_fft2_apply
    public :: certified_twiddles, fft_kappa, fft_kappa2, conv_error_bound
    public :: conv_error_bound_2d, conv_nonneg, conv_nonneg_2d, next_pow2

    !> Plan for length n = 2^p: certified twiddles w(0:n/2-1) and their
    !> error bound mu, plus the Higham factor kappa.
    type :: rigorous_fft_plan_t
        integer :: n = 0
        integer :: p = 0
        complex(dp), allocatable :: w(:)
        real(dp) :: mu = 0.0_dp
        real(dp) :: kappa = 0.0_dp
    end type rigorous_fft_plan_t

    real(dp), parameter :: u = unit_roundoff
    integer, parameter :: taylor_m = 16

contains

    !> Smallest power of two >= max(m, 1).
    pure function next_pow2(m) result(n)
        integer, intent(in) :: m
        integer :: n

        n = 1
        do while (n < m)
            n = 2*n
        end do
    end function next_pow2

    pure function log2_pow2(n) result(p)
        integer, intent(in) :: n
        integer :: p, m

        p = 0
        m = n
        do while (m > 1)
            m = m/2
            p = p + 1
        end do
    end function log2_pow2

    !> Binary128 enclosures of cos(phi), sin(phi) for phi = 2 pi m/n, with
    !> 0 <= m <= n/8.
    pure subroutine qsincos_octant(m, n, c, s)
        integer, intent(in) :: m, n
        type(qinterval_t), intent(out) :: c, s
        type(qinterval_t) :: phi, phi2
        real(qp) :: rm, ec, es
        integer :: j

        phi = qinterval_pi()*qinterval(2*m)/qinterval(n)
        phi2 = qsqr(phi)
        c = qinterval(1)
        s = qinterval(1)
        do j = taylor_m, 1, -1
            s = qinterval(1) - phi2*s/qinterval((2*j)*(2*j + 1))
            c = qinterval(1) - phi2*c/qinterval((2*j - 1)*(2*j))
        end do
        s = phi*s
        rm = phi%hi
        ec = 1.0_qp
        do j = 1, 2*taylor_m + 2
            ec = nearest(nearest(ec*rm, 1.0_qp)/real(j, qp), 1.0_qp)
        end do
        es = nearest(nearest(ec*rm, 1.0_qp)/real(2*taylor_m + 3, qp), 1.0_qp)
        c%lo = nearest(c%lo - ec, -1.0_qp)
        c%hi = nearest(c%hi + ec, 1.0_qp)
        s%lo = nearest(s%lo - es, -1.0_qp)
        s%hi = nearest(s%hi + es, 1.0_qp)
    end subroutine qsincos_octant

    !> Stored twiddles w(k) ~ exp(-2 pi i k/n), k = 0..n/2-1, n = 2^p >= 2,
    !> and mu >= max_k abs(w(k) - exp(-2 pi i k/n)).
    pure subroutine certified_twiddles(n, w, mu)
        integer, intent(in) :: n
        complex(dp), intent(out) :: w(0:)
        real(dp), intent(out) :: mu
        type(qinterval_t) :: c, s, cr, sr
        type(qinterval_t), allocatable :: co(:), so(:)
        real(qp) :: wr, wi, er, ei, e2, e2max
        integer :: k, m, q, quarter

        e2max = 0.0_qp
        quarter = max(n/4, 1)
        allocate (co(0:n/8), so(0:n/8))
        do m = 0, n/8
            call qsincos_octant(m, n, co(m), so(m))
        end do
        do k = 0, n/2 - 1
            q = k/quarter
            m = k - q*quarter
            if (n < 4) then
                c = qinterval(1)
                s = qinterval(0)
            else if (8*m <= n) then
                c = co(m)
                s = so(m)
            else
                c = so(quarter - m)
                s = co(quarter - m)
            end if
            ! Rotate by q quarter turns: (cos, sin)(q pi/2 + phi).
            if (q == 1) then
                cr = qinterval(0) - s
                sr = c
            else
                cr = c
                sr = s
            end if
            ! exp(-i theta) = cos theta - i sin theta.
            wr = qmid(cr)
            wi = -qmid(sr)
            w(k) = cmplx(real(wr, dp), real(wi, dp), dp)
            er = nearest(abs(real(real(w(k), dp), qp) - wr), 1.0_qp)
            er = nearest(er + qrad(cr), 1.0_qp)
            ei = nearest(abs(real(aimag(w(k)), qp) - wi), 1.0_qp)
            ei = nearest(ei + qrad(sr), 1.0_qp)
            e2 = nearest(nearest(er*er, 1.0_qp) + nearest(ei*ei, 1.0_qp), 1.0_qp)
            e2max = max(e2max, e2)
        end do
        mu = round_up(real(e2max, dp))
        mu = sqrt_up(mu)
    end subroutine certified_twiddles

    !> Higham's kappa(n) for n = 2^p and per-twiddle error mu.
    pure function fft_kappa(n, mu) result(kappa)
        integer, intent(in) :: n
        real(dp), intent(in) :: mu
        real(dp) :: kappa, eta, num, den

        eta = add_up(mu, mul_up(gamma_up(4), add_up(sqrt_up(2.0_dp), mu)))
        num = mul_up(real(log2_pow2(max(n, 2)), dp), eta)
        den = sub_down(1.0_dp, num)
        if (den <= 0.0_dp) then
            kappa = huge(1.0_dp)
        else
            kappa = div_up(num, den)
        end if
    end function fft_kappa

    !> Composed 2D factor kappa_m + kappa_n + kappa_m kappa_n.
    pure function fft_kappa2(km, kn) result(kappa2)
        real(dp), intent(in) :: km, kn
        real(dp) :: kappa2

        if (km >= huge(1.0_dp)/4.0_dp .or. kn >= huge(1.0_dp)/4.0_dp) then
            kappa2 = huge(1.0_dp)
        else
            kappa2 = add_up(add_up(km, kn), mul_up(km, kn))
        end if
    end function fft_kappa2

    !> Initialise a plan for length n (a power of two, n >= 1).
    pure subroutine rigorous_fft_plan_init(plan, n, status)
        type(rigorous_fft_plan_t), intent(out) :: plan
        integer, intent(in) :: n
        type(fortnum_status_t), intent(out), optional :: status

        if (present(status)) call status_set(status, FORTNUM_OK, "")
        if (n < 1 .or. next_pow2(n) /= n) then
            if (present(status)) call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rigorous_fft_plan_init: n must be a power of two")
            return
        end if
        plan%n = n
        plan%p = log2_pow2(n)
        allocate (plan%w(0:max(n/2 - 1, 0)))
        plan%mu = 0.0_dp
        plan%w = (1.0_dp, 0.0_dp)
        if (n >= 2) call certified_twiddles(n, plan%w, plan%mu)
        plan%kappa = fft_kappa(n, plan%mu)
    end subroutine rigorous_fft_plan_init

    !> In-place radix-2 FFT of x(0:n-1): forward y_j = sum_k x_k
    !> exp(-2 pi i jk/n), or inverse (conjugate twiddles, scaled by 1/n).
    pure subroutine rigorous_fft_apply(plan, x, inverse)
        type(rigorous_fft_plan_t), intent(in) :: plan
        complex(dp), intent(inout) :: x(0:)
        logical, intent(in) :: inverse
        complex(dp) :: t, wk
        integer :: n, i, j, m, mh, k, istep, stride

        n = plan%n
        if (n <= 1) return
        j = 0
        do i = 0, n - 2
            if (i < j) then
                t = x(i)
                x(i) = x(j)
                x(j) = t
            end if
            m = n/2
            do while (m >= 1)
                if (j < m) exit
                j = j - m
                m = m/2
            end do
            j = j + m
        end do
        mh = 1
        do while (mh < n)
            istep = 2*mh
            stride = n/istep
            do i = 0, mh - 1
                wk = plan%w(i*stride)
                if (inverse) wk = conjg(wk)
                do k = 0, n - 1, istep
                    t = wk*x(k + i + mh)
                    x(k + i + mh) = x(k + i) - t
                    x(k + i) = x(k + i) + t
                end do
            end do
            mh = istep
        end do
        if (inverse) x = x/real(n, dp)
    end subroutine rigorous_fft_apply

    !> 2D transform of x(0:pm-1, 0:pn-1): rows (length pn) then columns
    !> (length pm), the order the kappa2 composition assumes.
    pure subroutine rigorous_fft2_apply(plan_m, plan_n, x, inverse)
        type(rigorous_fft_plan_t), intent(in) :: plan_m, plan_n
        complex(dp), intent(inout) :: x(0:, 0:)
        logical, intent(in) :: inverse
        complex(dp) :: row(0:plan_n%n - 1)
        integer :: i, j

        do i = 0, plan_m%n - 1
            row = x(i, :)
            call rigorous_fft_apply(plan_n, row, inverse)
            x(i, :) = row
        end do
        do j = 0, plan_n%n - 1
            call rigorous_fft_apply(plan_m, x(:, j), inverse)
        end do
    end subroutine rigorous_fft2_apply

    !> Shared tail of the convolution bound for a transform with factor
    !> kappa and length xlen (n in 1D, pm pn in 2D).
    pure function compose_conv_err(kappa, xlen, l1a, l1b) result(e)
        real(dp), intent(in) :: kappa, xlen, l1a, l1b
        real(dp) :: e, sigma, rho, s5u, onepk

        if (kappa >= huge(1.0_dp)/4.0_dp) then
            e = huge(1.0_dp)
            return
        end if
        s5u = mul_up(sqrt_up(5.0_dp), u)
        rho = div_up(s5u, sub_down(1.0_dp, s5u))
        onepk = add_up(1.0_dp, kappa)
        sigma = add_up(add_up(mul_up(2.0_dp, kappa), mul_up(kappa, kappa)), &
            mul_up(rho, mul_up(onepk, onepk)))
        e = mul_up(sqrt_up(xlen), add_up(mul_up(add_up(kappa, u), &
            add_up(1.0_dp, sigma)), sigma))
        e = add_up(mul_up(e, mul_up(l1a, l1b)), 2.0_dp*tiny(1.0_dp))
    end function compose_conv_err

    !> Entrywise bound on abs(c_hat - c) for the FFT convolution (forward,
    !> forward, pointwise product, inverse) of two length-n sequences with
    !> l1 (or 2-) norms at most l1a and l1b, when n exceeds the linear
    !> convolution length.
    pure function conv_error_bound(plan, l1a, l1b) result(e)
        type(rigorous_fft_plan_t), intent(in) :: plan
        real(dp), intent(in) :: l1a, l1b
        real(dp) :: e

        e = compose_conv_err(plan%kappa, real(plan%n, dp), l1a, l1b)
    end function conv_error_bound

    !> 2D analogue of conv_error_bound for a pm x pn transform.
    pure function conv_error_bound_2d(plan_m, plan_n, l1a, l1b) result(e)
        type(rigorous_fft_plan_t), intent(in) :: plan_m, plan_n
        real(dp), intent(in) :: l1a, l1b
        real(dp) :: e

        e = compose_conv_err(fft_kappa2(plan_m%kappa, plan_n%kappa), &
            real(plan_m%n, dp)*real(plan_n%n, dp), l1a, l1b)
    end function conv_error_bound_2d

    !> Rigorous upper bound on the linear convolution of nonnegative
    !> sequences aa(-na:na), bb(-nb:nb): cc(modulo(k, n)) >= sum_i aa(i)
    !> bb(k - i) for abs(k) <= na + nb, where n = plan%n > 2 (na + nb).
    pure subroutine conv_nonneg(plan, aa, na, bb, nb, cc)
        type(rigorous_fft_plan_t), intent(in) :: plan
        integer, intent(in) :: na, nb
        real(dp), intent(in) :: aa(-na:na), bb(-nb:nb)
        real(dp), intent(out) :: cc(0:plan%n - 1)
        complex(dp) :: xa(0:plan%n - 1), xb(0:plan%n - 1)
        real(dp) :: margin
        integer :: i, n

        n = plan%n
        xa = (0.0_dp, 0.0_dp)
        xb = (0.0_dp, 0.0_dp)
        do i = -na, na
            xa(modulo(i, n)) = cmplx(aa(i), 0.0_dp, dp)
        end do
        do i = -nb, nb
            xb(modulo(i, n)) = cmplx(bb(i), 0.0_dp, dp)
        end do
        call rigorous_fft_apply(plan, xa, .false.)
        call rigorous_fft_apply(plan, xb, .false.)
        xa = xa*xb
        call rigorous_fft_apply(plan, xa, .true.)
        margin = conv_error_bound(plan, sum_up(aa), sum_up(bb))
        do i = 0, n - 1
            cc(i) = add_up(max(real(xa(i), dp), 0.0_dp), margin)
        end do
    end subroutine conv_nonneg

    !> 2D analogue of conv_nonneg for aa(-ma:ma, -na:na), bb(-mb:mb, -nb:nb)
    !> on a pm x pn grid with pm > 2 (ma + mb), pn > 2 (na + nb).
    pure subroutine conv_nonneg_2d(plan_m, plan_n, aa, ma, na, bb, mb, nb, cc)
        type(rigorous_fft_plan_t), intent(in) :: plan_m, plan_n
        integer, intent(in) :: ma, na, mb, nb
        real(dp), intent(in) :: aa(-ma:ma, -na:na), bb(-mb:mb, -nb:nb)
        real(dp), intent(out) :: cc(0:plan_m%n - 1, 0:plan_n%n - 1)
        complex(dp) :: xa(0:plan_m%n - 1, 0:plan_n%n - 1)
        complex(dp) :: xb(0:plan_m%n - 1, 0:plan_n%n - 1)
        real(dp) :: margin, l1a, l1b
        integer :: i, j, pm, pn

        pm = plan_m%n
        pn = plan_n%n
        xa = (0.0_dp, 0.0_dp)
        xb = (0.0_dp, 0.0_dp)
        l1a = 0.0_dp
        l1b = 0.0_dp
        do j = -na, na
            do i = -ma, ma
                xa(modulo(i, pm), modulo(j, pn)) = cmplx(aa(i, j), 0.0_dp, dp)
            end do
            l1a = add_up(l1a, sum_up(aa(:, j)))
        end do
        do j = -nb, nb
            do i = -mb, mb
                xb(modulo(i, pm), modulo(j, pn)) = cmplx(bb(i, j), 0.0_dp, dp)
            end do
            l1b = add_up(l1b, sum_up(bb(:, j)))
        end do
        call rigorous_fft2_apply(plan_m, plan_n, xa, .false.)
        call rigorous_fft2_apply(plan_m, plan_n, xb, .false.)
        xa = xa*xb
        call rigorous_fft2_apply(plan_m, plan_n, xa, .true.)
        margin = conv_error_bound_2d(plan_m, plan_n, l1a, l1b)
        do j = 0, pn - 1
            do i = 0, pm - 1
                cc(i, j) = add_up(max(real(xa(i, j), dp), 0.0_dp), margin)
            end do
        end do
    end subroutine conv_nonneg_2d
end module fortnum_cfft_rigorous
