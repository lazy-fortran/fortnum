!> Rigorous Fourier series with ball coefficients and l1 tails.
!>
!> 1D: F(theta) = sum_{abs(n)<=k} a_n exp(i n theta) + T with ball
!> coefficients a_n and tau >= the l1 norm of the coefficients of T, so
!> abs(T(theta)) <= tau everywhere.
!>
!> 2D: F(theta, phi) = sum_{abs(m)<=mk, abs(n)<=nk} a_mn
!> exp(i (m theta - n nper phi)) + T, with the tail measured in the
!> exponentially weighted l1 norm norm_sigma(c) = sum abs(c_mn)
!> exp(sigma (abs(m) + abs(n))), sigma >= 0. The weighted norm is a Banach
!> algebra norm (norm_sigma(F G) <= norm_sigma(F) norm_sigma(G)), which makes
!> truncated products and the Wiener-algebra reciprocal rigorous; sigma = 0
!> is the plain l1 norm. Weights are upper bounds from the libm-free interval
!> exp.
!>
!> Every operation returns a series whose exact function set contains the
!> exact result. Operations that cannot bound a result (a derivative of a
!> series with a tail, incompatible periods or weights, a reciprocal whose
!> Neumann series is not certified) return tau = huge, never a wrong bound.
!> The fast products replace the O(K^2) ball double sums by FFT convolutions
!> of the midpoints plus the rigorous bounds of fortnum_cfft_rigorous; the
!> input radii are propagated per coefficient by rigorous nonnegative
!> convolutions. This module generalises kinetic-compression kc_fseries and
!> kc_fseries2d.
module fortnum_fseries
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rounding, only: round_up, add_up, mul_up, div_up, &
        sub_down
    use fortnum_interval, only: interval_t, interval, exp
    use fortnum_ball, only: ball_t, ball, bpoint, badd, bsub, bmul, bneg, &
        babs_hi, cabs_hi, bconjg
    use fortnum_cfft_rigorous, only: rigorous_fft_plan_t, rigorous_fft_plan_init, &
        rigorous_fft_apply, rigorous_fft2_apply, conv_error_bound, &
        conv_error_bound_2d, &
        conv_nonneg, conv_nonneg_2d, next_pow2
    implicit none
    private

    public :: fseries_t, fs_zero, fs_trig, fs_add, fs_sub, fs_scale, fs_axpy
    public :: fs_l1_hi, fs_mul, fs_mul_fast, fs_dtheta, fs_inner_w
    public :: fs_inner_w_fast, fs_inverse, fs_fft_size
    public :: fseries2d_t, fs2_zero, fs2_trig, fs2_trig_tail, fs2_add, fs2_scale
    public :: fs2_axpy, fs2_sup_hi, fs2_wnorm_hi, fs2_mul, fs2_mul_fast
    public :: fs2_dtheta, fs2_dphi, fs2_dpar, fs2_inner_w, fs2_inner_w_fast
    public :: fs2_inverse

    type :: fseries_t
        integer :: k = 0
        type(ball_t), allocatable :: a(:)
        real(dp) :: tau = 0.0_dp
    end type fseries_t

    type :: fseries2d_t
        integer :: mk = 0
        integer :: nk = 0
        integer :: nper = 1
        real(dp) :: sigma = 0.0_dp
        type(ball_t), allocatable :: a(:, :)
        real(dp) :: tau = 0.0_dp
    end type fseries2d_t

    real(dp), parameter :: eta = 2.0_dp*tiny(1.0_dp)

contains

    !> FFT length that makes the product of series with cutoffs kf and kg
    !> free of wraparound.
    pure function fs_fft_size(kf, kg) result(n)
        integer, intent(in) :: kf, kg
        integer :: n

        n = next_pow2(2*(kf + kg) + 2)
    end function fs_fft_size

    pure function fs_zero(k) result(f)
        integer, intent(in) :: k
        type(fseries_t) :: f

        f%k = k
        allocate (f%a(-k:k))
        f%a = bpoint((0.0_dp, 0.0_dp))
        f%tau = 0.0_dp
    end function fs_zero

    !> Exact trigonometric polynomial from point coefficients c(-n:n).
    pure function fs_trig(c, n) result(f)
        integer, intent(in) :: n
        complex(dp), intent(in) :: c(-n:n)
        type(fseries_t) :: f

        f = fs_zero(n)
        f%a = bpoint(c)
    end function fs_trig

    pure function fs_add(f, g) result(h)
        type(fseries_t), intent(in) :: f, g
        type(fseries_t) :: h
        integer :: n

        h = fs_zero(max(f%k, g%k))
        do n = -f%k, f%k
            h%a(n) = badd(h%a(n), f%a(n))
        end do
        do n = -g%k, g%k
            h%a(n) = badd(h%a(n), g%a(n))
        end do
        h%tau = add_up(f%tau, g%tau)
    end function fs_add

    pure function fs_sub(f, g) result(h)
        type(fseries_t), intent(in) :: f, g
        type(fseries_t) :: h
        integer :: n

        h = fs_zero(max(f%k, g%k))
        do n = -f%k, f%k
            h%a(n) = badd(h%a(n), f%a(n))
        end do
        do n = -g%k, g%k
            h%a(n) = bsub(h%a(n), g%a(n))
        end do
        h%tau = add_up(f%tau, g%tau)
    end function fs_sub

    pure function fs_scale(f, s) result(h)
        type(fseries_t), intent(in) :: f
        type(ball_t), intent(in) :: s
        type(fseries_t) :: h
        integer :: n

        h = fs_zero(f%k)
        do n = -f%k, f%k
            h%a(n) = bmul(f%a(n), s)
        end do
        h%tau = mul_up(f%tau, babs_hi(s))
    end function fs_scale

    !> h = h + s f.
    pure subroutine fs_axpy(h, s, f)
        type(fseries_t), intent(inout) :: h
        type(ball_t), intent(in) :: s
        type(fseries_t), intent(in) :: f

        h = fs_add(h, fs_scale(f, s))
    end subroutine fs_axpy

    !> Upper bound on the l1 norm of all coefficients (hence on sup abs(F)).
    pure function fs_l1_hi(f) result(s)
        type(fseries_t), intent(in) :: f
        real(dp) :: s
        integer :: n

        s = f%tau
        do n = -f%k, f%k
            s = add_up(s, babs_hi(f%a(n)))
        end do
    end function fs_l1_hi

    !> Product truncated to abs(n) <= kmax: kept coefficients are exact ball
    !> double sums; truncated terms and tail cross terms go into tau:
    !> tau_h <= tau_f l1(g) + tau_g l1(f) + dropped (l1 includes the tails).
    pure function fs_mul(f, g, kmax) result(h)
        type(fseries_t), intent(in) :: f, g
        integer, intent(in) :: kmax
        type(fseries_t) :: h
        type(ball_t) :: p
        integer :: i, j, n
        real(dp) :: dropped

        h = fs_zero(kmax)
        dropped = 0.0_dp
        do j = -g%k, g%k
            do i = -f%k, f%k
                n = i + j
                p = bmul(f%a(i), g%a(j))
                if (abs(n) <= kmax) then
                    h%a(n) = badd(h%a(n), p)
                else
                    dropped = add_up(dropped, babs_hi(p))
                end if
            end do
        end do
        h%tau = add_up(add_up(mul_up(f%tau, fs_l1_hi(g)), &
            mul_up(g%tau, fs_l1_hi(f))), dropped)
    end function fs_mul

    !> Plan of length >= nmin: the caller's plan when it is large enough,
    !> otherwise a new one.
    pure subroutine pick_plan(nmin, plan_in, plan)
        integer, intent(in) :: nmin
        type(rigorous_fft_plan_t), intent(in), optional :: plan_in
        type(rigorous_fft_plan_t), intent(out) :: plan
        logical :: reuse

        reuse = present(plan_in)
        if (reuse) reuse = plan_in%n >= nmin
        if (reuse) then
            plan = plan_in
        else
            call rigorous_fft_plan_init(plan, next_pow2(nmin))
        end if
    end subroutine pick_plan

    !> Fast product truncated to abs(n) <= kmax. Midpoints by FFT
    !> convolution; per output index the exact coefficient differs from the
    !> FFT midpoint by (i) the FFT convolution error (conv_error_bound) and
    !> (ii) the radius propagation sum_i (abs(f_i) rad(g) + rad(f) abs(g_i) +
    !> rad(f) rad(g)), a convolution of nonnegative upper bounds computed by
    !> conv_nonneg, so each radius stays local to its coefficient. Tails and
    !> truncation are handled as in fs_mul. An optional plan of length at
    !> least 2 (f%k + g%k) + 1 is reused; otherwise one is built.
    pure function fs_mul_fast(f, g, kmax, plan) result(h)
        type(fseries_t), intent(in) :: f, g
        integer, intent(in) :: kmax
        type(rigorous_fft_plan_t), intent(in), optional :: plan
        type(fseries_t) :: h
        type(rigorous_fft_plan_t) :: pl
        integer :: nf, ng, nlin, npad, i, n
        real(dp) :: dropped, l1am, l1bm, merr, r
        complex(dp), allocatable :: am(:), bm(:)
        real(dp), allocatable :: aa(:), ra(:), bb(:), rb(:), c1(:), c2(:), c3(:)

        nf = f%k
        ng = g%k
        nlin = nf + ng
        call pick_plan(2*nlin + 1, plan, pl)
        npad = pl%n
        allocate (am(0:npad - 1), bm(0:npad - 1))
        allocate (aa(-nf:nf), ra(-nf:nf), bb(-ng:ng), rb(-ng:ng))
        allocate (c1(0:npad - 1), c2(0:npad - 1), c3(0:npad - 1))
        am = (0.0_dp, 0.0_dp)
        bm = (0.0_dp, 0.0_dp)
        l1am = 0.0_dp
        l1bm = 0.0_dp
        do i = -nf, nf
            am(modulo(i, npad)) = f%a(i)%c
            aa(i) = cabs_hi(f%a(i)%c)
            ra(i) = f%a(i)%r
            l1am = add_up(l1am, aa(i))
        end do
        do i = -ng, ng
            bm(modulo(i, npad)) = g%a(i)%c
            bb(i) = cabs_hi(g%a(i)%c)
            rb(i) = g%a(i)%r
            l1bm = add_up(l1bm, bb(i))
        end do
        call rigorous_fft_apply(pl, am, .false.)
        call rigorous_fft_apply(pl, bm, .false.)
        am = am*bm
        call rigorous_fft_apply(pl, am, .true.)
        merr = conv_error_bound(pl, l1am, l1bm)
        if (any(ra > 0.0_dp) .or. any(rb > 0.0_dp)) then
            call conv_nonneg(pl, aa, nf, rb, ng, c1)
            call conv_nonneg(pl, ra, nf, bb, ng, c2)
            call conv_nonneg(pl, ra, nf, rb, ng, c3)
        else
            c1 = 0.0_dp
            c2 = 0.0_dp
            c3 = 0.0_dp
        end if
        h = fs_zero(kmax)
        dropped = 0.0_dp
        do n = -nlin, nlin
            i = modulo(n, npad)
            r = add_up(add_up(add_up(add_up(merr, c1(i)), c2(i)), c3(i)), eta)
            if (abs(n) <= kmax) then
                h%a(n) = ball(am(i), r)
            else
                dropped = add_up(dropped, add_up(cabs_hi(am(i)), r))
            end if
        end do
        h%tau = add_up(add_up(mul_up(f%tau, fs_l1_hi(g)), &
            mul_up(g%tau, fs_l1_hi(f))), dropped)
    end function fs_mul_fast

    !> d/dtheta of an exact trigonometric polynomial; a series with a tail
    !> gives tau = huge.
    pure function fs_dtheta(f) result(h)
        type(fseries_t), intent(in) :: f
        type(fseries_t) :: h
        integer :: n

        h = fs_zero(f%k)
        do n = -f%k, f%k
            h%a(n) = bmul(f%a(n), bpoint(cmplx(0.0_dp, real(n, dp), dp)))
        end do
        if (f%tau > 0.0_dp) h%tau = huge(1.0_dp)
    end function fs_dtheta

    !> Enclosure of (1/2 pi) int conj(F) G w dtheta for a weight series w
    !> with sup abs(w) <= wmax: exact ball sums on the finite parts, tails
    !> bounded by sup norms.
    pure function fs_inner_w(f, g, w, wmax) result(s)
        type(fseries_t), intent(in) :: f, g, w
        real(dp), intent(in) :: wmax
        type(ball_t) :: s
        integer :: i, j
        real(dp) :: rad, sf, sg

        s = bpoint((0.0_dp, 0.0_dp))
        do i = -f%k, f%k
            do j = max(-g%k, i - w%k), min(g%k, i + w%k)
                s = badd(s, bmul(bmul(bconjg(f%a(i)), g%a(j)), w%a(i - j)))
            end do
        end do
        sf = round_up(fs_l1_hi(f) - f%tau)
        sg = round_up(fs_l1_hi(g) - g%tau)
        rad = mul_up(mul_up(sf, sg), w%tau)
        rad = add_up(rad, mul_up(add_up(mul_up(f%tau, fs_l1_hi(g)), &
            mul_up(g%tau, sf)), wmax))
        s%r = add_up(s%r, rad)
    end function fs_inner_w

    !> Fast version of fs_inner_w: P = conj(F) G by fs_mul_fast truncated to
    !> abs(n) <= kmax, then sum_n P_n w_{-n} on the overlap, the remainder
    !> bounded by l1(P) tau_w + tau_P wmax.
    pure function fs_inner_w_fast(f, g, w, wmax, kmax, plan) result(s)
        type(fseries_t), intent(in) :: f, g, w
        real(dp), intent(in) :: wmax
        integer, intent(in) :: kmax
        type(rigorous_fft_plan_t), intent(in), optional :: plan
        type(ball_t) :: s
        type(fseries_t) :: fbar, p
        integer :: n, kk
        real(dp) :: pl1, rad

        fbar = fs_zero(f%k)
        do n = -f%k, f%k
            fbar%a(n) = bconjg(f%a(-n))
        end do
        fbar%tau = f%tau
        p = fs_mul_fast(fbar, g, kmax, plan)
        s = bpoint((0.0_dp, 0.0_dp))
        kk = min(p%k, w%k)
        do n = -kk, kk
            s = badd(s, bmul(p%a(n), w%a(-n)))
        end do
        pl1 = round_up(fs_l1_hi(p) - p%tau)
        rad = add_up(mul_up(pl1, w%tau), mul_up(p%tau, wmax))
        s%r = add_up(s%r, rad)
    end function fs_inner_w_fast

    !> Wiener-algebra reciprocal. An approximate inverse y (trigonometric
    !> polynomial of degree kinv from samples of 1/F by FFT, not trusted) is
    !> certified by rho = 1 - F y computed as a ball series: if
    !> q = l1(rho) < 1, then 1/F = y + E with l1(E) <= l1(y) q/(1 - q)
    !> (Neumann series in the Banach algebra). The result is y with that tail;
    !> q >= 1 gives tau = huge.
    pure function fs_inverse(f, kinv, plan) result(h)
        type(fseries_t), intent(in) :: f
        integer, intent(in) :: kinv
        type(rigorous_fft_plan_t), intent(in), optional :: plan
        type(fseries_t) :: h, rho
        type(rigorous_fft_plan_t) :: ps
        complex(dp), allocatable :: x(:)
        integer :: ns, n
        real(dp) :: q, ly

        ns = next_pow2(4*max(f%k, kinv, 1))
        call rigorous_fft_plan_init(ps, ns)
        allocate (x(0:ns - 1))
        x = (0.0_dp, 0.0_dp)
        do n = -f%k, f%k
            x(modulo(n, ns)) = f%a(n)%c
        end do
        call rigorous_fft_apply(ps, x, .true.)
        x = 1.0_dp/(x*real(ns, dp))
        call rigorous_fft_apply(ps, x, .false.)
        x = x/real(ns, dp)
        h = fs_zero(kinv)
        do n = -kinv, kinv
            h%a(n) = bpoint(x(modulo(n, ns)))
        end do
        rho = fs_mul_fast(f, h, f%k + kinv, plan)
        do n = -rho%k, rho%k
            rho%a(n) = bneg(rho%a(n))
        end do
        rho%a(0) = badd(rho%a(0), bpoint((1.0_dp, 0.0_dp)))
        q = fs_l1_hi(rho)
        ly = fs_l1_hi(h)
        if (q < 1.0_dp) then
            h%tau = div_up(mul_up(ly, q), sub_down(1.0_dp, q))
        else
            h%tau = huge(1.0_dp)
        end if
    end function fs_inverse

    !> Upper bound on exp(sigma).
    pure function ebase_hi(sigma) result(e)
        real(dp), intent(in) :: sigma
        real(dp) :: e
        type(interval_t) :: t

        if (sigma > 0.0_dp) then
            t = exp(interval(sigma))
            e = t%hi
        else
            e = 1.0_dp
        end if
    end function ebase_hi

    !> Upper bounds w(j) on exp(sigma j), j = 0..jmax.
    pure function weight_table(sigma, jmax) result(w)
        real(dp), intent(in) :: sigma
        integer, intent(in) :: jmax
        real(dp) :: w(0:max(jmax, 0))
        real(dp) :: eb
        integer :: j

        eb = ebase_hi(sigma)
        w(0) = 1.0_dp
        do j = 1, jmax
            w(j) = w(j - 1)
            if (eb > 1.0_dp) w(j) = mul_up(w(j - 1), eb)
        end do
    end function weight_table

    pure function compatible(f, g) result(ok)
        type(fseries2d_t), intent(in) :: f, g
        logical :: ok

        ok = f%nper == g%nper .and. .not. (f%sigma < g%sigma .or. f%sigma > g%sigma)
    end function compatible

    pure function fs2_zero(mk, nk, nper, sigma) result(f)
        integer, intent(in) :: mk, nk, nper
        real(dp), intent(in) :: sigma
        type(fseries2d_t) :: f

        f%mk = mk
        f%nk = nk
        f%nper = nper
        f%sigma = sigma
        allocate (f%a(-mk:mk, -nk:nk))
        f%a = bpoint((0.0_dp, 0.0_dp))
        f%tau = 0.0_dp
    end function fs2_zero

    !> Exact trigonometric polynomial from point coefficients.
    pure function fs2_trig(c, mk, nk, nper, sigma) result(f)
        integer, intent(in) :: mk, nk, nper
        complex(dp), intent(in) :: c(-mk:mk, -nk:nk)
        real(dp), intent(in) :: sigma
        type(fseries2d_t) :: f

        f = fs2_zero(mk, nk, nper, sigma)
        f%a = bpoint(c)
    end function fs2_trig

    !> Exact spectrum c(-mfull:mfull, -nfull:nfull) truncated to the kept
    !> window abs(m) <= mk, abs(n) <= nk; every coefficient outside the window
    !> goes into the weighted tail, so later certificates hold for the full
    !> spectrum.
    pure function fs2_trig_tail(c, mfull, nfull, mk, nk, nper, sigma) result(f)
        integer, intent(in) :: mfull, nfull, mk, nk, nper
        complex(dp), intent(in) :: c(-mfull:mfull, -nfull:nfull)
        real(dp), intent(in) :: sigma
        type(fseries2d_t) :: f
        integer :: m, n
        real(dp) :: w(0:mfull + nfull)

        f = fs2_zero(mk, nk, nper, sigma)
        w = weight_table(sigma, mfull + nfull)
        do n = -nfull, nfull
            do m = -mfull, mfull
                if (abs(m) <= mk .and. abs(n) <= nk) then
                    f%a(m, n) = bpoint(c(m, n))
                else
                    f%tau = add_up(f%tau, mul_up(cabs_hi(c(m, n)), w(abs(m) + abs(n))))
                end if
            end do
        end do
    end function fs2_trig_tail

    pure function fs2_add(f, g) result(h)
        type(fseries2d_t), intent(in) :: f, g
        type(fseries2d_t) :: h
        integer :: m, n

        h = fs2_zero(max(f%mk, g%mk), max(f%nk, g%nk), f%nper, f%sigma)
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                h%a(m, n) = badd(h%a(m, n), f%a(m, n))
            end do
        end do
        do n = -g%nk, g%nk
            do m = -g%mk, g%mk
                h%a(m, n) = badd(h%a(m, n), g%a(m, n))
            end do
        end do
        h%tau = add_up(f%tau, g%tau)
        if (.not. compatible(f, g)) h%tau = huge(1.0_dp)
    end function fs2_add

    pure function fs2_scale(f, s) result(h)
        type(fseries2d_t), intent(in) :: f
        type(ball_t), intent(in) :: s
        type(fseries2d_t) :: h
        integer :: m, n

        h = fs2_zero(f%mk, f%nk, f%nper, f%sigma)
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                h%a(m, n) = bmul(f%a(m, n), s)
            end do
        end do
        h%tau = mul_up(f%tau, babs_hi(s))
    end function fs2_scale

    pure subroutine fs2_axpy(h, s, f)
        type(fseries2d_t), intent(inout) :: h
        type(ball_t), intent(in) :: s
        type(fseries2d_t), intent(in) :: f

        h = fs2_add(h, fs2_scale(f, s))
    end subroutine fs2_axpy

    !> Upper bound on the unweighted l1 norm, hence on sup abs(F) (the
    !> weighted tail bound dominates the unweighted one since sigma >= 0).
    pure function fs2_sup_hi(f) result(s)
        type(fseries2d_t), intent(in) :: f
        real(dp) :: s
        integer :: m, n

        s = f%tau
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                s = add_up(s, babs_hi(f%a(m, n)))
            end do
        end do
    end function fs2_sup_hi

    !> Upper bound on the full sigma-weighted l1 norm.
    pure function fs2_wnorm_hi(f) result(s)
        type(fseries2d_t), intent(in) :: f
        real(dp) :: s
        real(dp) :: w(0:f%mk + f%nk)
        integer :: m, n

        w = weight_table(f%sigma, f%mk + f%nk)
        s = f%tau
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                s = add_up(s, mul_up(babs_hi(f%a(m, n)), w(abs(m) + abs(n))))
            end do
        end do
    end function fs2_wnorm_hi

    !> Product truncated to abs(m) <= mmax, abs(n) <= nmax. Dropped terms are
    !> weighted at their output index; tails contribute
    !> tau_f norm_sigma(g) + tau_g norm_sigma(f).
    pure function fs2_mul(f, g, mmax, nmax) result(h)
        type(fseries2d_t), intent(in) :: f, g
        integer, intent(in) :: mmax, nmax
        type(fseries2d_t) :: h
        type(ball_t) :: p
        integer :: m1, n1, m2, n2, m, n
        real(dp) :: dropped
        real(dp) :: w(0:f%mk + g%mk + f%nk + g%nk)

        h = fs2_zero(mmax, nmax, f%nper, f%sigma)
        w = weight_table(f%sigma, f%mk + g%mk + f%nk + g%nk)
        dropped = 0.0_dp
        do n2 = -g%nk, g%nk
            do m2 = -g%mk, g%mk
                do n1 = -f%nk, f%nk
                    do m1 = -f%mk, f%mk
                        m = m1 + m2
                        n = n1 + n2
                        p = bmul(f%a(m1, n1), g%a(m2, n2))
                        if (abs(m) <= mmax .and. abs(n) <= nmax) then
                            h%a(m, n) = badd(h%a(m, n), p)
                        else
                            dropped = add_up(dropped, mul_up(babs_hi(p), &
                                w(abs(m) + abs(n))))
                        end if
                    end do
                end do
            end do
        end do
        h%tau = add_up(add_up(mul_up(f%tau, fs2_wnorm_hi(g)), &
            mul_up(g%tau, fs2_wnorm_hi(f))), dropped)
        if (.not. compatible(f, g)) h%tau = huge(1.0_dp)
    end function fs2_mul

    !> d/dtheta (factor i m); a series with a tail gives tau = huge.
    pure function fs2_dtheta(f) result(h)
        type(fseries2d_t), intent(in) :: f
        type(fseries2d_t) :: h
        integer :: m, n

        h = fs2_zero(f%mk, f%nk, f%nper, f%sigma)
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                h%a(m, n) = bmul(f%a(m, n), bpoint(cmplx(0.0_dp, real(m, dp), dp)))
            end do
        end do
        if (f%tau > 0.0_dp) h%tau = huge(1.0_dp)
    end function fs2_dtheta

    !> d/dphi (factor -i n nper); a series with a tail gives tau = huge.
    pure function fs2_dphi(f) result(h)
        type(fseries2d_t), intent(in) :: f
        type(fseries2d_t) :: h
        integer :: m, n
        real(dp) :: fac

        h = fs2_zero(f%mk, f%nk, f%nper, f%sigma)
        do n = -f%nk, f%nk
            fac = -real(n, dp)*real(f%nper, dp)
            do m = -f%mk, f%mk
                h%a(m, n) = bmul(f%a(m, n), bpoint(cmplx(0.0_dp, fac, dp)))
            end do
        end do
        if (f%tau > 0.0_dp) h%tau = huge(1.0_dp)
    end function fs2_dphi

    !> Parallel derivative iota d/dtheta + d/dphi with a ball iota.
    pure function fs2_dpar(f, iota) result(h)
        type(fseries2d_t), intent(in) :: f
        type(ball_t), intent(in) :: iota
        type(fseries2d_t) :: h

        h = fs2_add(fs2_scale(fs2_dtheta(f), iota), fs2_dphi(f))
    end function fs2_dpar

    !> Enclosure of the mean of conj(F) G w over the torus for a weight
    !> series w with sup abs(w) <= wmax; tails bounded by sup norms.
    pure function fs2_inner_w(f, g, w, wmax) result(s)
        type(fseries2d_t), intent(in) :: f, g, w
        real(dp), intent(in) :: wmax
        type(ball_t) :: s
        integer :: mi, ni, mj, nj
        real(dp) :: rad, sf, sg

        s = bpoint((0.0_dp, 0.0_dp))
        do ni = -f%nk, f%nk
            do mi = -f%mk, f%mk
                do nj = max(-g%nk, ni - w%nk), min(g%nk, ni + w%nk)
                    do mj = max(-g%mk, mi - w%mk), min(g%mk, mi + w%mk)
                        s = badd(s, bmul(bmul(bconjg(f%a(mi, ni)), g%a(mj, nj)), &
                            w%a(mi - mj, ni - nj)))
                    end do
                end do
            end do
        end do
        sf = round_up(fs2_sup_hi(f) - f%tau)
        sg = round_up(fs2_sup_hi(g) - g%tau)
        rad = mul_up(mul_up(sf, sg), w%tau)
        rad = add_up(rad, mul_up(add_up(mul_up(f%tau, fs2_sup_hi(g)), &
            mul_up(g%tau, sf)), wmax))
        s%r = add_up(s%r, rad)
        if (.not. (compatible(f, g) .and. compatible(f, w))) s%r = huge(1.0_dp)
    end function fs2_inner_w

    !> 2D plans of sizes at least (mmin, nmin).
    pure subroutine pick_plans(mmin, nmin, plan_m, plan_n, pm, pn)
        integer, intent(in) :: mmin, nmin
        type(rigorous_fft_plan_t), intent(in), optional :: plan_m, plan_n
        type(rigorous_fft_plan_t), intent(out) :: pm, pn

        call pick_plan(mmin, plan_m, pm)
        call pick_plan(nmin, plan_n, pn)
    end subroutine pick_plans

    !> Fast 2D product truncated to abs(m) <= mmax, abs(n) <= nmax: 2D FFT
    !> convolution of the midpoints with conv_error_bound_2d, radii by
    !> conv_nonneg_2d, dropped coefficients weighted at their output index
    !> (weight(m1, n1) weight(m2, n2) >= weight(m1 + m2, n1 + n2)).
    pure function fs2_mul_fast(f, g, mmax, nmax, plan_m, plan_n) result(h)
        type(fseries2d_t), intent(in) :: f, g
        integer, intent(in) :: mmax, nmax
        type(rigorous_fft_plan_t), intent(in), optional :: plan_m, plan_n
        type(fseries2d_t) :: h
        type(rigorous_fft_plan_t) :: pm, pn
        integer :: ml, nl, m, n, i, j, mm, nn
        real(dp) :: dropped, l1a, l1b, merr, r
        complex(dp), allocatable :: am(:, :), bm(:, :)
        real(dp), allocatable :: aa(:, :), ra(:, :), bb(:, :), rb(:, :)
        real(dp), allocatable :: c1(:, :), c2(:, :), c3(:, :), w(:)

        ml = f%mk + g%mk
        nl = f%nk + g%nk
        call pick_plans(2*ml + 1, 2*nl + 1, plan_m, plan_n, pm, pn)
        mm = pm%n
        nn = pn%n
        allocate (am(0:mm - 1, 0:nn - 1), bm(0:mm - 1, 0:nn - 1))
        allocate (aa(-f%mk:f%mk, -f%nk:f%nk), ra(-f%mk:f%mk, -f%nk:f%nk))
        allocate (bb(-g%mk:g%mk, -g%nk:g%nk), rb(-g%mk:g%mk, -g%nk:g%nk))
        allocate (c1(0:mm - 1, 0:nn - 1), c2(0:mm - 1, 0:nn - 1))
        allocate (c3(0:mm - 1, 0:nn - 1), w(0:ml + nl))
        am = (0.0_dp, 0.0_dp)
        bm = (0.0_dp, 0.0_dp)
        l1a = 0.0_dp
        l1b = 0.0_dp
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                am(modulo(m, mm), modulo(n, nn)) = f%a(m, n)%c
                aa(m, n) = cabs_hi(f%a(m, n)%c)
                ra(m, n) = f%a(m, n)%r
                l1a = add_up(l1a, aa(m, n))
            end do
        end do
        do n = -g%nk, g%nk
            do m = -g%mk, g%mk
                bm(modulo(m, mm), modulo(n, nn)) = g%a(m, n)%c
                bb(m, n) = cabs_hi(g%a(m, n)%c)
                rb(m, n) = g%a(m, n)%r
                l1b = add_up(l1b, bb(m, n))
            end do
        end do
        call rigorous_fft2_apply(pm, pn, am, .false.)
        call rigorous_fft2_apply(pm, pn, bm, .false.)
        am = am*bm
        call rigorous_fft2_apply(pm, pn, am, .true.)
        merr = conv_error_bound_2d(pm, pn, l1a, l1b)
        if (any(ra > 0.0_dp) .or. any(rb > 0.0_dp)) then
            call conv_nonneg_2d(pm, pn, aa, f%mk, f%nk, rb, g%mk, g%nk, c1)
            call conv_nonneg_2d(pm, pn, ra, f%mk, f%nk, bb, g%mk, g%nk, c2)
            call conv_nonneg_2d(pm, pn, ra, f%mk, f%nk, rb, g%mk, g%nk, c3)
        else
            c1 = 0.0_dp
            c2 = 0.0_dp
            c3 = 0.0_dp
        end if
        w = weight_table(f%sigma, ml + nl)
        h = fs2_zero(mmax, nmax, f%nper, f%sigma)
        dropped = 0.0_dp
        do n = -nl, nl
            j = modulo(n, nn)
            do m = -ml, ml
                i = modulo(m, mm)
                r = add_up(add_up(add_up(add_up(merr, c1(i, j)), c2(i, j)), &
                    c3(i, j)), eta)
                if (abs(m) <= mmax .and. abs(n) <= nmax) then
                    h%a(m, n) = ball(am(i, j), r)
                else
                    dropped = add_up(dropped, mul_up(add_up(cabs_hi(am(i, j)), r), &
                        w(abs(m) + abs(n))))
                end if
            end do
        end do
        h%tau = add_up(add_up(mul_up(f%tau, fs2_wnorm_hi(g)), &
            mul_up(g%tau, fs2_wnorm_hi(f))), dropped)
        if (.not. compatible(f, g)) h%tau = huge(1.0_dp)
    end function fs2_mul_fast

    !> Fast version of fs2_inner_w via P = conj(F) G from fs2_mul_fast.
    pure function fs2_inner_w_fast(f, g, w, wmax, mmax, nmax, plan_m, plan_n) &
            result(s)
        type(fseries2d_t), intent(in) :: f, g, w
        real(dp), intent(in) :: wmax
        integer, intent(in) :: mmax, nmax
        type(rigorous_fft_plan_t), intent(in), optional :: plan_m, plan_n
        type(ball_t) :: s
        type(fseries2d_t) :: fbar, p
        integer :: m, n, km, kn
        real(dp) :: pl1, rad

        fbar = fs2_zero(f%mk, f%nk, f%nper, f%sigma)
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                fbar%a(m, n) = bconjg(f%a(-m, -n))
            end do
        end do
        fbar%tau = f%tau
        p = fs2_mul_fast(fbar, g, mmax, nmax, plan_m, plan_n)
        s = bpoint((0.0_dp, 0.0_dp))
        km = min(p%mk, w%mk)
        kn = min(p%nk, w%nk)
        do n = -kn, kn
            do m = -km, km
                s = badd(s, bmul(p%a(m, n), w%a(-m, -n)))
            end do
        end do
        pl1 = round_up(fs2_sup_hi(p) - p%tau)
        rad = add_up(mul_up(pl1, w%tau), mul_up(p%tau, wmax))
        s%r = add_up(s%r, rad)
        if (.not. (compatible(f, g) .and. compatible(f, w))) s%r = huge(1.0_dp)
    end function fs2_inner_w_fast

    !> 2D Wiener-algebra reciprocal in the sigma-weighted norm, as
    !> fs_inverse: 1/F = y + E, norm_sigma(E) <= norm_sigma(y) q/(1 - q),
    !> q = norm_sigma(1 - F y) < 1; otherwise tau = huge.
    pure function fs2_inverse(f, minv, ninv, plan_m, plan_n) result(h)
        type(fseries2d_t), intent(in) :: f
        integer, intent(in) :: minv, ninv
        type(rigorous_fft_plan_t), intent(in), optional :: plan_m, plan_n
        type(fseries2d_t) :: h, rho
        type(rigorous_fft_plan_t) :: sm, sn
        complex(dp), allocatable :: x(:, :)
        integer :: ms, ns, m, n
        real(dp) :: q, ly

        ms = next_pow2(4*max(f%mk, minv, 1))
        ns = next_pow2(4*max(f%nk, ninv, 1))
        call rigorous_fft_plan_init(sm, ms)
        call rigorous_fft_plan_init(sn, ns)
        allocate (x(0:ms - 1, 0:ns - 1))
        x = (0.0_dp, 0.0_dp)
        do n = -f%nk, f%nk
            do m = -f%mk, f%mk
                x(modulo(m, ms), modulo(n, ns)) = f%a(m, n)%c
            end do
        end do
        call rigorous_fft2_apply(sm, sn, x, .true.)
        x = 1.0_dp/(x*(real(ms, dp)*real(ns, dp)))
        call rigorous_fft2_apply(sm, sn, x, .false.)
        x = x/(real(ms, dp)*real(ns, dp))
        h = fs2_zero(minv, ninv, f%nper, f%sigma)
        do n = -ninv, ninv
            do m = -minv, minv
                h%a(m, n) = bpoint(x(modulo(m, ms), modulo(n, ns)))
            end do
        end do
        rho = fs2_mul_fast(f, h, f%mk + minv, f%nk + ninv, plan_m, plan_n)
        do n = -rho%nk, rho%nk
            do m = -rho%mk, rho%mk
                rho%a(m, n) = bneg(rho%a(m, n))
            end do
        end do
        rho%a(0, 0) = badd(rho%a(0, 0), bpoint((1.0_dp, 0.0_dp)))
        q = fs2_wnorm_hi(rho)
        ly = fs2_wnorm_hi(h)
        if (q < 1.0_dp) then
            h%tau = div_up(mul_up(ly, q), sub_down(1.0_dp, q))
        else
            h%tau = huge(1.0_dp)
        end if
    end function fs2_inverse
end module fortnum_fseries
