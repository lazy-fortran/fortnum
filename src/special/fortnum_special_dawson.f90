module fortnum_special_dawson
    ! Dawson integral F(x) = exp(-x^2) * integral_0^x exp(t^2) dt.
    !
    ! Derivative policy (ad.md sec 1, sec 4): analytic_rule.
    !   F'(x) = 1 - 2*x*F(x).
    ! Active argument: x (real scalar). Result F(x) is the active output.
    ! Derivative entry point (not yet implemented): dawson_grad(x, dF_dout).
    !
    ! Algorithm:
    !   |x| < 1   Maclaurin series (DLMF 7.6.4).
    !   1 <= |x| < 10   Rybicki sampling-theorem method (Rybicki 1989,
    !                    Computers in Physics 3, 85; DLMF 7.10).
    !   |x| >= 10   Asymptotic expansion (DLMF 7.12.2).
    !
    ! Numerics are verbatim from the libneo math-kit clean-room port
    ! (src/math/dawson.f90, git branch math-kit) with literal style
    ! conformance: 0.0d0 -> _dp.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: dawson
    ! analytic_rule derivatives (ad.md §2):
    !   dawson_jvp: forward product F'(x) * v = (1 - 2*x*F(x)) * v
    !   dawson_grad: scalar gradient (same coefficient, d/dx applied to scalar)
    ! Active argument: x. No HVP here (would require F'' = -2F - 2xF').
    public :: dawson_jvp
    public :: dawson_grad

    real(dp), parameter :: inv_sqrt_pi = 0.564189583547756287_dp

    ! Maclaurin coefficients (-2)^n / (2n+1)!!, n = 0..20 (DLMF 7.6.4)
    integer, parameter :: nmac = 20
    real(dp), parameter :: mac(0:nmac) = [ &
        1.0_dp, &
        -0.666666666666666667_dp, &
        0.266666666666666667_dp, &
        -0.0761904761904761905_dp, &
        0.0169312169312169312_dp, &
        -0.0030784030784030784_dp, &
        0.0004736004736004736_dp, &
        -0.0000631467298133964801_dp, &
        7.42902703687017413e-6_dp, &
        -7.82002845986334118e-7_dp, &
        7.44764615225080113e-8_dp, &
        -6.47621404543547924e-9_dp, &
        5.18097123634838339e-10_dp, &
        -3.83775647136917288e-11_dp, &
        2.64672860094425716e-12_dp, &
        -1.70756683931887559e-13_dp, &
        1.03488899352659127e-14_dp, &
        -5.91365139158052152e-16_dp, &
        3.19656831977325487e-17_dp, &
        -1.63926580501192558e-18_dp, &
        7.9964185610337833e-20_dp]

    ! Rybicki sampling step h and precomputed weights exp(-((2k-1)*h)^2),
    ! k = 1..14; sampling error ~ exp(-(pi/(2h))^2) ~ 7e-18 for h = 1/4.
    real(dp), parameter :: h_ryb = 0.25_dp
    integer,  parameter :: nryb = 14
    real(dp), parameter :: cryb(nryb) = [ &
        0.939413062813475786_dp, &
        0.56978282473092301_dp, &
        0.209611387151097823_dp, &
        0.0467706223839589837_dp, &
        0.00632971542748574658_dp, &
        0.000519574682154838482_dp, &
        0.0000258681002226541213_dp, &
        7.8114894083044908e-7_dp, &
        1.43072419185676883e-8_dp, &
        1.58939100945163665e-10_dp, &
        1.07092323825080765e-12_dp, &
        4.37661850287084989e-15_dp, &
        1.0848552640429378e-17_dp, &
        1.63101392267018568e-20_dp]

    ! Asymptotic coefficients (2m-1)!!, m = 0..14 (DLMF 7.12.2)
    integer, parameter :: nasy = 14
    real(dp), parameter :: asy(0:nasy) = [ &
        1.0_dp, 1.0_dp, 3.0_dp, 15.0_dp, 105.0_dp, 945.0_dp, 10395.0_dp, &
        135135.0_dp, 2027025.0_dp, 34459425.0_dp, 654729075.0_dp, &
        13749310575.0_dp, 316234143225.0_dp, 7905853580625.0_dp, &
        213458046676875.0_dp]

contains

    elemental function dawson(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f

        real(dp) :: a

        a = abs(x)
        if (a < 1.0_dp) then
            f = dawson_series(x)
        else if (a < 10.0_dp) then
            f = sign(dawson_rybicki(a), x)
        else
            f = sign(dawson_asymptotic(a), x)
        end if
    end function dawson

    pure function dawson_series(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f

        real(dp) :: t, s
        integer  :: n

        t = x*x
        s = mac(nmac)
        do n = nmac - 1, 0, -1
            s = mac(n) + t*s
        end do
        f = x*s
    end function dawson_series

    pure function dawson_rybicki(a) result(f)
        ! Rybicki (1989): F(a) via Gaussian sampling on a shifted grid so the
        ! denominator factors never hit zero.
        real(dp), intent(in) :: a
        real(dp) :: f

        integer  :: n0, k
        real(dp) :: xp, gauss, e1, e2, einv2, ep, em, s

        n0 = 2*nint(0.5_dp*a/h_ryb)
        xp = a - real(n0, dp)*h_ryb
        gauss = exp(-xp*xp)
        e1 = exp(2.0_dp*h_ryb*xp)
        e2 = e1*e1
        einv2 = 1.0_dp/e2
        ep = e1
        em = 1.0_dp/e1
        s = 0.0_dp
        do k = 1, nryb
            s = s + cryb(k)*(ep/real(n0 + 2*k - 1, dp) + em/real(n0 - 2*k + 1, dp))
            ep = ep*e2
            em = em*einv2
        end do
        f = s*gauss*inv_sqrt_pi
    end function dawson_rybicki

    pure function dawson_asymptotic(a) result(f)
        ! Asymptotic: F(a) ~ (1/2a) * sum_{m=0}^{nasy} (2m-1)!!/(2a^2)^m.
        ! t underflows to 0 for a > ~1e154, leaving the exact 1/(2a) tail.
        real(dp), intent(in) :: a
        real(dp) :: f

        real(dp) :: t, s
        integer  :: m

        t = 0.5_dp/(a*a)
        s = asy(nasy)
        do m = nasy - 1, 0, -1
            s = asy(m) + t*s
        end do
        f = 0.5_dp/a*s
    end function dawson_asymptotic

    ! Forward product: jv = F'(x) v, where F'(x) = 1 - 2*x*F(x).
    ! x and v are length-1 arrays conforming to the harness interface.
    subroutine dawson_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = (1.0_dp - 2.0_dp*x(1)*dawson(x(1))) * v(1)
    end subroutine dawson_jvp

    ! Scalar gradient: grad = dF/dx = 1 - 2*x*F(x).
    ! Returns the gradient as a length-1 array so callers may use it uniformly.
    subroutine dawson_grad(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        ! For scalar output the VJP and JVP coincide; u is the output cotangent.
        jtu(1) = u(1) * (1.0_dp - 2.0_dp*x(1)*dawson(x(1)))
    end subroutine dawson_grad

end module fortnum_special_dawson
