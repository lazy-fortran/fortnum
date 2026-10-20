! Bounded, real l=0 Coulomb functions, DLMF 33.2 and 33.6.
module fortnum_special_coulomb
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_special_gamma, only: digamma_complex
    use fortnum_coulomb_series_kernel, only: coulomb_series_step
    implicit none
    private
    public :: coulomb_fg0, coulomb_c0_squared, coulomb_h
contains
    pure function coulomb_h(eta) result(h)
        real(dp), intent(in) :: eta
        real(dp) :: h
        ! Only meaningful for nonzero eta; the caller handles eta=0.
        h = real(digamma_complex(cmplx(1.0_dp, eta, dp)), dp) - log(abs(eta))
    end function

    pure function coulomb_c0_squared(eta) result(c2)
        real(dp), intent(in) :: eta
        real(dp) :: c2, x
        x = 2*acos(-1.0_dp)*eta
        if (abs(x) < 1.0e-5_dp) then
            c2 = 1 - x/2 + x*x/12 - x**4/720
        else if (x > 700) then
            c2 = 0
        else
            c2 = x/(exp(x) - 1)
        end if
    end function

    ! Derivatives fp,gp are with respect to rho. No unreported extrapolation:
    ! supported domain |eta|<=30, 0<rho<=8, |eta*rho|<=2.
    ! Internal binary128 summation suppresses cancellation of logarithmic terms.
    ! Returned success is convergence/finite status, NOT a rigorous enclosure.
    pure subroutine coulomb_fg0(eta, rho, f, g, fp, gp, ok)
        real(dp), intent(in) :: eta, rho
        real(dp), intent(out) :: f, g, fp, gp
        logical, intent(out) :: ok
        real(qp), parameter :: euler = &
                               0.577215664901532860606512090082402431_qp
        real(qp) :: z, t, f1, f2, v1, v2, fn, vn, fs, vs, df, dv, c, h, lg
        integer :: n
        f = 0; g = 0; fp = 0; gp = 0; ok = .false.
        if (.not. ieee_is_finite(eta)) return
        if (.not. ieee_is_finite(rho)) return
        if (abs(eta) > 30 .or. rho <= 0 .or. rho > 8) return
        if (abs(eta*rho) > 2) return
        z = real(rho, qp); t = real(eta, qp)
        c = sqrt(real(coulomb_c0_squared(eta), qp))
        h = real(real(digamma_complex(cmplx(1.0_dp, eta, dp)), dp), qp)
        f2 = 1; f1 = t*z
        v2 = 1; v1 = 2*t*z*(h + 2*euler - 1)
        fs = f2 + f1; vs = v2 + v1; df = f1/z; dv = v1/z
        do n = 2, 256
            call coulomb_series_step(real(n, qp), t, z, f1, f2, v1, v2, fn, vn)
            fs = fs + fn; vs = vs + vn
            df = df + n*fn/z; dv = dv + n*vn/z
            f2 = f1; f1 = fn; v2 = v1; v1 = vn
            if (n < 32) cycle
            if (max(abs(fn), abs(vn), abs(f2), abs(v2)) < 1.0e-32_qp) exit
        end do
        lg = log(2*z)
        f = real(c*z*fs, dp); fp = real(c*(fs + z*df), dp)
        g = real((vs + 2*t*z*fs*lg)/c, dp)
        gp = real((dv + 2*t*((fs + z*df)*lg + fs))/c, dp)
        ok = n <= 256 .and. all(ieee_is_finite([f, g, fp, gp]))
    end subroutine
end module
