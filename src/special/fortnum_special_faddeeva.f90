module fortnum_special_faddeeva
    ! Faddeeva function w(z) = exp(-z^2) erfc(-i z) (DLMF 7.2.3) and the
    ! Fried-Conte plasma dispersion function Z(zeta) = i sqrt(pi) w(zeta)
    ! (DLMF 7.19.2 in Z form).
    !
    ! Z(zeta) = pi^(-1/2) int exp(-t^2)/(t - zeta) dt for Im zeta > 0; w and
    ! Z are entire, so the same formula is the Landau analytic continuation
    ! to Im zeta <= 0. Z(0) = i sqrt(pi).
    !
    ! Algorithm (Gautschi, SIAM J. Numer. Anal. 7 (1970) 187-198, with the
    ! region and term-count constants of Poppe and Wijers, ACM TOMS 16 (1990)
    ! 38-46, algorithm 680): first-quadrant evaluation by a Taylor series
    ! of erf near the origin and by Laplace continued-fraction/Taylor
    ! combinations elsewhere, then w(-conj z) = conj w(z) and
    ! w(-z) = 2 exp(-z^2) - w(z) (DLMF 7.4.3, 7.4.4). Relative accuracy is
    ! about 1e-14; for |z| > 1e7 the first-quadrant value is the two-term
    ! asymptotic DLMF 7.12.1. In the lower half plane 2 exp(-z^2) overflows
    ! to IEEE infinity once Im(z)^2 - Re(z)^2 exceeds about 709.
    !
    ! Derivative candidate (ad.md sec. 1): analytical,
    !   w'(z) = 2 i/sqrt(pi) - 2 z w(z),  Z'(zeta) = -2 (1 + zeta Z(zeta)).
    ! plasma_dispersion_z_derivative returns Z'; higher derivatives follow
    ! from Z'' = -2 (Z + zeta Z').
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: faddeeva_w
    public :: plasma_dispersion_z, plasma_dispersion_z_derivative

    real(dp), parameter :: two_over_sqrt_pi = 1.12837916709551257390_dp
    real(dp), parameter :: sqrt_pi = 1.77245385090551602730_dp

contains

    elemental function faddeeva_w(z) result(w)
        complex(dp), intent(in) :: z
        complex(dp) :: w
        real(dp) :: xabs, yabs, xs, ys, qrho, xquad, yquad, u, v, u2, v2
        real(dp) :: scale
        logical :: series

        xabs = abs(real(z))
        yabs = abs(aimag(z))
        xquad = (xabs - yabs)*(xabs + yabs)
        yquad = 2.0_dp*xabs*yabs
        series = .false.
        u2 = 0.0_dp
        v2 = 0.0_dp
        if (max(xabs, yabs) > 1.0e7_dp) then
            w = first_quadrant_asymptotic(cmplx(xabs, yabs, dp))
            u = real(w)
            v = aimag(w)
        else
            xs = xabs/6.3_dp
            ys = yabs/4.4_dp
            qrho = xs*xs + ys*ys
            series = qrho < 0.085264_dp
            if (series) then
                call first_quadrant_series(xabs, yabs, xquad, yquad, &
                    (1.0_dp - 0.85_dp*ys)*sqrt(qrho), u, v, u2, v2)
            else
                call first_quadrant_fraction(xabs, yabs, ys, qrho, u, v)
            end if
        end if

        if (aimag(z) < 0.0_dp) then
            if (series) then
                u2 = 2.0_dp*u2
                v2 = 2.0_dp*v2
            else
                scale = 2.0_dp*exp(-xquad)
                u2 = scale*cos(yquad)
                v2 = -scale*sin(yquad)
            end if
            u = u2 - u
            v = v2 - v
            if (real(z) > 0.0_dp) v = -v
        else if (real(z) < 0.0_dp) then
            v = -v
        end if
        w = cmplx(u, v, dp)
    end function faddeeva_w

    elemental function plasma_dispersion_z(zeta) result(z)
        complex(dp), intent(in) :: zeta
        complex(dp) :: z

        z = cmplx(0.0_dp, sqrt_pi, dp)*faddeeva_w(zeta)
    end function plasma_dispersion_z

    elemental function plasma_dispersion_z_derivative(zeta) result(dz)
        complex(dp), intent(in) :: zeta
        complex(dp) :: dz

        dz = -2.0_dp*(1.0_dp + zeta*plasma_dispersion_z(zeta))
    end function plasma_dispersion_z_derivative

    pure subroutine first_quadrant_series(xabs, yabs, xquad, yquad, rho, &
            u, v, u2, v2)
        ! w = exp(-z^2) (1 + (2i/sqrt(pi)) z sum_n (z^2)^n/(n!(2n+1)))
        ! for small |z| in the first quadrant; u2 + i v2 = exp(-z^2).
        real(dp), intent(in) :: xabs, yabs, xquad, yquad, rho
        real(dp), intent(out) :: u, v, u2, v2
        real(dp) :: xsum, ysum, xaux, u1, v1, decay
        integer :: n, i, j

        n = nint(6.0_dp + 72.0_dp*rho)
        j = 2*n + 1
        xsum = 1.0_dp/real(j, dp)
        ysum = 0.0_dp
        do i = n, 1, -1
            j = j - 2
            xaux = (xsum*xquad - ysum*yquad)/real(i, dp)
            ysum = (xsum*yquad + ysum*xquad)/real(i, dp)
            xsum = xaux + 1.0_dp/real(j, dp)
        end do
        u1 = 1.0_dp - two_over_sqrt_pi*(xsum*yabs + ysum*xabs)
        v1 = two_over_sqrt_pi*(xsum*xabs - ysum*yabs)
        decay = exp(-xquad)
        u2 = decay*cos(yquad)
        v2 = -decay*sin(yquad)
        u = u1*u2 - v1*v2
        v = u1*v2 + v1*u2
    end subroutine first_quadrant_series

    pure subroutine first_quadrant_fraction(xabs, yabs, ys, qrho, u, v)
        ! Laplace continued fraction (qrho > 1) or Gautschi's truncated
        ! Taylor/continued-fraction combination with step h (qrho <= 1).
        real(dp), intent(in) :: xabs, yabs, ys, qrho
        real(dp), intent(out) :: u, v
        real(dp) :: rho, h, h2, qlambda, rx, ry, sx, sy, tx, ty, c
        integer :: kapn, nu, n

        if (qrho > 1.0_dp) then
            h = 0.0_dp
            kapn = 0
            rho = sqrt(qrho)
            nu = int(3.0_dp + 1442.0_dp/(26.0_dp*rho + 77.0_dp))
        else
            rho = (1.0_dp - ys)*sqrt(1.0_dp - qrho)
            h = 1.88_dp*rho
            kapn = nint(7.0_dp + 34.0_dp*rho)
            nu = nint(16.0_dp + 26.0_dp*rho)
        end if
        h2 = 2.0_dp*h
        qlambda = 0.0_dp
        if (h > 0.0_dp) qlambda = h2**kapn
        rx = 0.0_dp
        ry = 0.0_dp
        sx = 0.0_dp
        sy = 0.0_dp
        do n = nu, 0, -1
            tx = yabs + h + real(n + 1, dp)*rx
            ty = xabs - real(n + 1, dp)*ry
            c = 0.5_dp/(tx*tx + ty*ty)
            rx = c*tx
            ry = c*ty
            if (h > 0.0_dp .and. n <= kapn) then
                tx = qlambda + sx
                sx = rx*tx - ry*sy
                sy = ry*tx + rx*sy
                qlambda = qlambda/h2
            end if
        end do
        if (h > 0.0_dp) then
            u = two_over_sqrt_pi*sx
            v = two_over_sqrt_pi*sy
        else
            u = two_over_sqrt_pi*rx
            v = two_over_sqrt_pi*ry
        end if
        if (yabs == 0.0_dp) u = exp(-xabs*xabs)
    end subroutine first_quadrant_fraction

    pure function first_quadrant_asymptotic(z) result(w)
        ! w(z) ~ (i/(sqrt(pi) z)) (1 + 1/(2 z^2)), remainder 3/(4|z|^4).
        complex(dp), intent(in) :: z
        complex(dp) :: w
        complex(dp) :: r

        r = 1.0_dp/z
        w = cmplx(0.0_dp, 1.0_dp/sqrt_pi, dp)*r*(1.0_dp + 0.5_dp*r*r)
        if (aimag(z) == 0.0_dp) w = cmplx(exp(-real(z)**2), aimag(w), dp)
    end function first_quadrant_asymptotic

end module fortnum_special_faddeeva
