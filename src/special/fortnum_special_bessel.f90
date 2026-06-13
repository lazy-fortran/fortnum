module fortnum_special_bessel
    ! Modified Bessel functions of integer order.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4):
    !   Default class: analytic_rule.
    !   Derivatives follow known recurrences:
    !     d/dx I_n(x) = (I_{n-1}(x) + I_{n+1}(x)) / 2
    !     d/dx K_n(x) = -(K_{n-1}(x) + K_{n+1}(x)) / 2
    !   Active arguments (ad.md §3): x only.
    !   Inactive: n (integer order, selects branch).
    !   Derivative entry points (ad.md §2): bessel_in_grad, bessel_kn_grad —
    !   not implemented here; signatures will be added in the AD milestone
    !   without touching the primal signatures below.
    !
    ! I_n(x) — bessel_in, bessel_in_array
    !   Power series DLMF 10.25.2 for small |x|.
    !   Downward (Miller) recurrence with normalization e^x = I_0(x) + 2 sum_{k>=1}
    !   I_k(x) (DLMF 10.29.1, A&S 9.6.36) for moderate |x|.
    !   Asymptotic expansion DLMF 10.40.1 for large |x|.
    !
    ! K_n(x) — bessel_kn
    !   K_0, K_1 for x <= 3: power series DLMF 10.31.1/10.31.2.
    !   K_0, K_1 for x > 3: Chebyshev fits of sqrt(2x/pi) e^x K_nu(x).
    !   K_n: upward recurrence DLMF 10.29.1 (forward-stable for K).

    use fortnum_kinds, only: dp

    implicit none
    private

    public :: bessel_in, bessel_in_array, bessel_kn
    ! analytic_rule derivatives (ad.md §2):
    !   bessel_in_jvp: d/dx I_n(x) * v = (I_{n-1}(x)+I_{n+1}(x))/2 * v  (DLMF 10.29.2)
    !   bessel_kn_jvp: d/dx K_n(x) * v = -(K_{n-1}(x)+K_{n+1}(x))/2 * v (DLMF 10.29.4)
    ! Active argument: x. Inactive: n (integer order selector, ad.md §3).
    ! HVP not provided (second-order recurrences involve more orders; deferred).
    public :: bessel_in_jvp
    public :: bessel_kn_jvp

    ! --- I_n parameters ---
    real(dp), parameter :: series_x_max    = 2.0_dp
    real(dp), parameter :: rescale_limit   = 1.0e250_dp
    real(dp), parameter :: rescale_factor  = 1.0e-250_dp
    real(dp), parameter :: start_decay_tol = 1.0e-18_dp
    real(dp), parameter :: seed            = 1.0e-30_dp
    real(dp), parameter :: two_pi          = 6.2831853071795864769_dp
    real(dp), parameter :: asym_x_min      = 20.0_dp
    real(dp), parameter :: asym_n_factor   = 0.7_dp

    ! --- K_n parameters ---
    real(dp), parameter :: euler_gamma     = 0.57721566490153286060651209008240243_dp
    real(dp), parameter :: half_pi         = 1.5707963267948966192313216916397514_dp
    real(dp), parameter :: k_series_x_max  = 3.0_dp
    real(dp), parameter :: cheb_hi_x_min   = 16.0_dp
    real(dp), parameter :: cheb_mid_s_lo   = 0.0625_dp
    real(dp), parameter :: cheb_mid_s_hi   = 0.33333333333333333333333333333333333_dp

    ! I_0 series coefficients: 1/(k!)^2, DLMF 10.25.2
    real(dp), parameter :: poly_i0(18) = [ &
        1.0_dp, 1.0_dp, &
        0.25_dp, 0.027777777777777777778_dp, &
        0.0017361111111111111111_dp, 0.000069444444444444444444_dp, &
        1.9290123456790123457e-6_dp, 3.9367598891408415218e-8_dp, &
        6.1511873267825648778e-10_dp, 7.5940584281266233059e-12_dp, &
        7.5940584281266233059e-14_dp, 6.2760813455591928148e-16_dp, &
        4.3583898233049950103e-18_dp, 2.5789288895295828463e-20_dp, &
        1.3157800456783585951e-22_dp, 5.8479113141260382003e-25_dp, &
        2.284340357080483672e-27_dp, 7.9042918930120542283e-30_dp]

    ! K_0 series: h_k/(k!)^2 with h_k = sum_{j=1..k} 1/j
    real(dp), parameter :: poly_s0(18) = [ &
        0.0_dp, 1.0_dp, &
        0.375_dp, 0.050925925925925925926_dp, &
        0.0036168981481481481481_dp, 0.00015856481481481481481_dp, &
        4.7260802469135802469e-6_dp, 1.0207455998272324803e-7_dp, &
        1.6718048413148328114e-9_dp, 2.1483350211950276805e-11_dp, &
        2.2242756054762939135e-13_dp, 1.8952995870061529211e-15_dp, &
        1.3525001839484811536e-17_dp, 8.2013388136826374592e-20_dp, &
        4.2783408265702079911e-22_dp, 1.9404708872364882507e-24_dp, &
        7.7227356755850624588e-27_dp, 2.7187227120298503047e-29_dp]

    ! I_1 series coefficients: 1/(k!(k+1)!)
    real(dp), parameter :: poly_i1(18) = [ &
        1.0_dp, 0.5_dp, &
        0.083333333333333333333_dp, 0.0069444444444444444444_dp, &
        0.00034722222222222222222_dp, 0.000011574074074074074074_dp, &
        2.7557319223985890653e-7_dp, 4.9209498614260519022e-9_dp, &
        6.8346525853139609753e-11_dp, 7.5940584281266233059e-13_dp, &
        6.9036894801151120963e-15_dp, 5.2300677879659940123e-17_dp, &
        3.3526075563884577002e-19_dp, 1.8420920639497020331e-21_dp, &
        8.7718669711890573004e-24_dp, 3.6549445713287738752e-26_dp, &
        1.3437296218120492188e-28_dp, 4.3912732738955856824e-31_dp]

    ! K_1 series: (h_k + h_{k+1} - 2*gamma)/(k!(k+1)!)
    real(dp), parameter :: poly_s1(18) = [ &
        -0.15443132980306572121_dp, 0.67278433509846713939_dp, &
        0.18157516696085563434_dp, 0.019182189839330562121_dp, &
        0.0011153594919665281061_dp, 0.000041422476892711430696_dp, &
        1.0715459140911808686e-6_dp, 2.0452860035938779413e-8_dp, &
        3.0020487465891878859e-10_dp, 3.4959287296928819208e-12_dp, &
        3.309914735250272068e-14_dp, 2.5986411321011287351e-16_dp, &
        1.7195232826992565241e-18_dp, 9.7212075188236180165e-21_dp, &
        4.7502817433276669896e-23_dp, 2.0264937604328579082e-25_dp, &
        7.6133707277671159281e-28_dp, 2.5382566314880592899e-30_dp]

    ! Chebyshev fits in s = 1/x on (1/16, 1/3) for K_0, K_1
    real(dp), parameter :: cheb_mid_k0(18) = [ &
        1.9559376047895837408_dp, -0.014059053153789941523_dp, &
        0.00039802254448432713646_dp, -0.000018657916479681936657_dp, &
        1.1670921100603393265e-6_dp, -8.8637049995049946834e-8_dp, &
        7.7618435419763736638e-9_dp, -7.5903717486808993347e-10_dp, &
        8.1130444797671898315e-11_dp, -9.3346884671836175433e-12_dp, &
        1.1431343028648823935e-12_dp, -1.4771049122256338483e-13_dp, &
        2.0002339096669601218e-14_dp, -2.8230963387131460904e-15_dp, &
        4.1342702911603966579e-16_dp, -6.2586795525484043848e-17_dp, &
        9.7637329920200614695e-18_dp, -1.5654516590192795234e-18_dp]
    real(dp), parameter :: cheb_mid_k1(18) = [ &
        2.1390906750838713041_dp, 0.045796323598856923076_dp, &
        -0.00071639260360370820452_dp, 0.000028086590879043600368_dp, &
        -1.6076666325349922967e-6_dp, 1.1571456080375915070e-7_dp, &
        -9.7722547348374030052e-9_dp, 9.3086831666491458465e-10_dp, &
        -9.7528215763417439059e-11_dp, 1.1045807616711904736e-11_dp, &
        -1.3354794780368619141e-12_dp, 1.7074202453340075633e-13_dp, &
        -2.2914753368530012268e-14_dp, 3.2093797833233986927e-15_dp, &
        -4.6686970748489713955e-16_dp, 7.0264558327001338014e-17_dp, &
        -1.0904794030447221775e-17_dp, 1.7403276745873914004e-18_dp]

    ! Chebyshev fits in u = 32/x - 1 for x >= 16 for K_0, K_1
    real(dp), parameter :: cheb_hi_k0(12) = [ &
        1.9923831602004083696_dp, -0.0037766317031612104235_dp, &
        0.000031310052233667130168_dp, -4.6784426011890146e-7_dp, &
        1.0013945577406131153e-8_dp, -2.7653023813480403757e-10_dp, &
        9.2831416580801534697e-12_dp, -3.6467062746741891765e-13_dp, &
        1.6325975548816327494e-14_dp, -8.1701412478988243465e-16_dp, &
        4.5035599694971361505e-17_dp, -2.7029444926470642288e-18_dp]
    real(dp), parameter :: cheb_hi_k1(12) = [ &
        2.0231087348982674842_dp, 0.011500735467006774991_dp, &
        -0.000052954102844829022231_dp, 6.6446695239227473635e-7_dp, &
        -1.3058062241211819979e-8_dp, 3.4269978500370871047e-10_dp, &
        -1.1121538753031974548e-11_dp, 4.2645352569241138993e-13_dp, &
        -1.8748591483665248044e-14_dp, 9.2508042350153445531e-16_dp, &
        -5.0417867748703588125e-17_dp, 2.9979881728938231118e-18_dp]

contains

    ! I_n(x) for integer n. Handles negative n (I_{-n} = I_n) and negative x.
    elemental function bessel_in(n, x) result(value)
        integer,  intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: value

        integer  :: na
        real(dp) :: ax

        na = abs(n)
        ax = abs(x)
        if (ax <= series_x_max) then
            value = series_in(na, ax)
        else if (ax >= max(asym_x_min, asym_n_factor*real(na, dp)**2)) then
            value = asym_in(na, ax)
        else
            value = miller_in(na, ax)
        end if
        if (x < 0.0_dp .and. mod(na, 2) == 1) value = -value
    end function bessel_in


    ! Fill values(0:nmax) with I_0(x)...I_nmax(x) in one pass.
    pure subroutine bessel_in_array(nmax, x, values)
        integer,  intent(in)  :: nmax
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: ax

        ax = abs(x)
        if (ax <= series_x_max) then
            call series_in_array(nmax, ax, values)
        else
            call miller_in_array(nmax, ax, values)
        end if
        if (x < 0.0_dp .and. nmax >= 1) then
            values(1:nmax:2) = -values(1:nmax:2)
        end if
    end subroutine bessel_in_array


    ! K_n(x) for integer n >= 0 and x > 0.  K_{-n} = K_n by symmetry.
    pure function bessel_kn(n, x) result(kn)
        integer,  intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: kn

        real(dp) :: k_prev, k_curr, k_next
        integer  :: n_abs, m

        n_abs = abs(n)
        call bessel_k01(x, k_prev, k_curr)
        if (n_abs == 0) then
            kn = k_prev
        else if (n_abs == 1) then
            kn = k_curr
        else
            do m = 1, n_abs - 1
                k_next = k_prev + (2.0_dp*real(m, dp)/x)*k_curr
                k_prev = k_curr
                k_curr = k_next
            end do
            kn = k_curr
        end if
    end function bessel_kn

    ! Forward product for I_n w.r.t. x (n inactive).
    ! dI_n/dx = (I_{n-1}(x) + I_{n+1}(x)) / 2  (DLMF 10.29.2).
    ! Uses I_{-1} = I_1 (symmetry) so n=0 is handled correctly.
    pure subroutine bessel_in_jvp(n, x, v, jv)
        integer,  intent(in)  :: n    ! inactive order
        real(dp), intent(in)  :: x    ! active argument
        real(dp), intent(in)  :: v    ! tangent
        real(dp), intent(out) :: jv   ! directional derivative
        jv = 0.5_dp*(bessel_in(n - 1, x) + bessel_in(n + 1, x)) * v
    end subroutine bessel_in_jvp

    ! Forward product for K_n w.r.t. x (n inactive).
    ! dK_n/dx = -(K_{n-1}(x) + K_{n+1}(x)) / 2  (DLMF 10.29.4).
    ! Uses K_{-1} = K_1 (symmetry) so n=0 is handled correctly.
    pure subroutine bessel_kn_jvp(n, x, v, jv)
        integer,  intent(in)  :: n    ! inactive order
        real(dp), intent(in)  :: x    ! active argument (must be > 0)
        real(dp), intent(in)  :: v    ! tangent
        real(dp), intent(out) :: jv   ! directional derivative
        jv = -0.5_dp*(bessel_kn(n - 1, x) + bessel_kn(n + 1, x)) * v
    end subroutine bessel_kn_jvp

    ! ------------------------------------------------------------------ I_n internals

    elemental function series_in(n, ax) result(value)
        integer,  intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: pref, halfx
        integer  :: j

        halfx = 0.5_dp*ax
        pref  = 1.0_dp
        do j = 1, n
            pref = pref*halfx/real(j, dp)
            if (pref == 0.0_dp) exit
        end do
        if (pref == 0.0_dp) then
            value = 0.0_dp
        else
            value = pref*series_sum(n, 0.25_dp*ax*ax)
        end if
    end function series_in


    elemental function series_sum(n, q) result(total)
        ! DLMF 10.25.2: sum over k of q^k / (k! (n+k)!/n!) with (x/2)^n/n! split off.
        integer,  intent(in) :: n
        real(dp), intent(in) :: q
        real(dp) :: total

        real(dp) :: term
        integer  :: k

        term  = 1.0_dp
        total = 1.0_dp
        k = 0
        do
            k = k + 1
            term = term*q/(real(k, dp)*real(n + k, dp))
            total = total + term
            if (term <= epsilon(1.0_dp)*total) exit
        end do
    end function series_sum


    elemental function miller_in(n, ax) result(value)
        integer,  intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: p_hi, p_cur, p_lo, s, res
        integer  :: m, k

        m     = start_order(max(n, int(ax) + 1), ax)
        p_hi  = 0.0_dp
        p_cur = seed
        s     = 0.0_dp
        res   = 0.0_dp
        do k = m, 1, -1
            p_lo = p_hi + (2.0_dp*real(k, dp)/ax)*p_cur
            s    = s + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_limit) then
                p_lo  = p_lo*rescale_factor
                p_cur = p_cur*rescale_factor
                s     = s*rescale_factor
                res   = res*rescale_factor
            end if
            p_hi  = p_cur
            p_cur = p_lo
            if (k - 1 == n) res = p_cur
        end do
        s     = s + p_cur
        value = (res/s)*exp(ax)
    end function miller_in


    elemental function asym_in(n, ax) result(value)
        ! DLMF 10.40.1: I_n(x) ~ e^x/sqrt(2 pi x) sum_k (-1)^k a_k(n)/x^k.
        integer,  intent(in) :: n
        real(dp), intent(in) :: ax
        real(dp) :: value

        real(dp) :: mu, term, total, prev, r8x
        integer  :: k

        mu    = 4.0_dp*real(n, dp)**2
        r8x   = 1.0_dp/(8.0_dp*ax)
        term  = 1.0_dp
        total = 1.0_dp
        prev  = 1.0_dp
        do k = 1, 40
            term = -term*(mu - real(2*k - 1, dp)**2)*r8x/real(k, dp)
            if (abs(term) >= prev) exit
            total = total + term
            prev  = abs(term)
            if (prev <= epsilon(1.0_dp)*abs(total)) exit
        end do
        value = exp(ax)/sqrt(two_pi*ax)*total
    end function asym_in


    pure subroutine series_in_array(nmax, ax, values)
        integer,  intent(in)  :: nmax
        real(dp), intent(in)  :: ax
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: pref, halfx, q
        integer  :: j

        halfx = 0.5_dp*ax
        q     = 0.25_dp*ax*ax
        pref  = 1.0_dp
        do j = 0, nmax
            if (j > 0) pref = pref*halfx/real(j, dp)
            if (pref == 0.0_dp) then
                values(j:nmax) = 0.0_dp
                return
            end if
            values(j) = pref*series_sum(j, q)
        end do
    end subroutine series_in_array


    pure subroutine miller_in_array(nmax, ax, values)
        integer,  intent(in)  :: nmax
        real(dp), intent(in)  :: ax
        real(dp), intent(out) :: values(0:nmax)

        real(dp) :: p_hi, p_cur, p_lo, s
        integer  :: m, k

        m     = start_order(max(nmax, int(ax) + 1), ax)
        p_hi  = 0.0_dp
        p_cur = seed
        s     = 0.0_dp
        do k = m, 1, -1
            p_lo = p_hi + (2.0_dp*real(k, dp)/ax)*p_cur
            s    = s + 2.0_dp*p_cur
            if (abs(p_lo) > rescale_limit) then
                p_lo  = p_lo*rescale_factor
                p_cur = p_cur*rescale_factor
                s     = s*rescale_factor
                if (k <= nmax) values(k:nmax) = values(k:nmax)*rescale_factor
            end if
            p_hi  = p_cur
            p_cur = p_lo
            if (k - 1 <= nmax) values(k - 1) = p_cur
        end do
        s      = s + p_cur
        values = values*(exp(ax)/s)
    end subroutine miller_in_array


    elemental function start_order(k0, ax) result(m)
        ! Estimate starting order for Miller recurrence so the seed is negligible
        ! relative to the true value at k0.
        integer,  intent(in) :: k0
        real(dp), intent(in) :: ax
        integer :: m

        real(dp) :: decay

        m     = k0
        decay = 1.0_dp
        do while (decay > start_decay_tol)
            m     = m + 1
            decay = decay*ax/(real(m, dp) + sqrt(real(m, dp)**2 + ax*ax))
        end do
        m = m + 2
    end function start_order

    ! ------------------------------------------------------------------ K_n internals

    pure subroutine bessel_k01(x, k0, k1)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: k0, k1

        real(dp) :: u, pref

        if (x <= k_series_x_max) then
            call k01_series(x, k0, k1)
        else
            pref = sqrt(half_pi/x)*exp(-x)
            if (x < cheb_hi_x_min) then
                u  = (2.0_dp/x - cheb_mid_s_lo - cheb_mid_s_hi) &
                     /(cheb_mid_s_hi - cheb_mid_s_lo)
                k0 = pref*cheb_eval(cheb_mid_k0, u)
                k1 = pref*cheb_eval(cheb_mid_k1, u)
            else
                u  = 32.0_dp/x - 1.0_dp
                k0 = pref*cheb_eval(cheb_hi_k0, u)
                k1 = pref*cheb_eval(cheb_hi_k1, u)
            end if
        end if
    end subroutine bessel_k01


    ! DLMF 10.31.2 for K_0; DLMF 10.31.1 for K_1.
    pure subroutine k01_series(x, k0, k1)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: k0, k1

        real(dp) :: q, log_half_x, i0, s0, i1s, s1

        q          = 0.25_dp*x*x
        log_half_x = log(0.5_dp*x)
        i0         = poly_eval(poly_i0, q)
        s0         = poly_eval(poly_s0, q)
        i1s        = poly_eval(poly_i1, q)
        s1         = poly_eval(poly_s1, q)
        k0         = -(log_half_x + euler_gamma)*i0 + s0
        k1         = 1.0_dp/x + 0.5_dp*x*(log_half_x*i1s - 0.5_dp*s1)
    end subroutine k01_series


    pure function poly_eval(c, q) result(p)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: q
        real(dp) :: p

        integer :: k

        p = c(size(c))
        do k = size(c) - 1, 1, -1
            p = p*q + c(k)
        end do
    end function poly_eval


    ! Clenshaw recurrence; c(1) is the halved T_0 coefficient.
    pure function cheb_eval(c, u) result(f)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: u
        real(dp) :: f

        real(dp) :: b0, b1, b2, two_u
        integer  :: j

        two_u = 2.0_dp*u
        b1    = 0.0_dp
        b2    = 0.0_dp
        do j = size(c), 2, -1
            b0 = two_u*b1 - b2 + c(j)
            b2 = b1
            b1 = b0
        end do
        f = u*b1 - b2 + 0.5_dp*c(1)
    end function cheb_eval

end module fortnum_special_bessel
