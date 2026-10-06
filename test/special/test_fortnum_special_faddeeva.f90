program test_fortnum_special_faddeeva
    ! Faddeeva w(z) and plasma dispersion Z(zeta) against independent
    ! oracles:
    !   * quad-precision Maclaurin series w(z) = sum (iz)^n/Gamma(n/2+1)
    !     (DLMF 7.6.3) for |z| <= 5 in all quadrants, used only where the
    !     series cancellation sum|t_n|/|w| stays below 1e12;
    !   * quad-precision asymptotic series DLMF 7.12.1 at optimal truncation
    !     for |z| >= 8, Im z >= 0, and in the lower half plane through
    !     w(z) = 2 exp(-z^2) - w(-z) (DLMF 7.4.3) with quad exp(-z^2);
    !   * the real axis: Re w(x) = exp(-x^2), Im w(x) = (2/sqrt(pi)) F(x)
    !     with FortNum's Dawson integral F;
    !   * Z(0) = i sqrt(pi), and Z'(zeta) = -2(1 + zeta Z) against a
    !     Cauchy-integral derivative of plasma_dispersion_z itself;
    !   * the Landau-damped Langmuir root of 1 - Z'(zeta)/(2 k^2) = 0 at
    !     k lambda_D = 0.5, omega/omega_p = 1.4156 - 0.1533 i (Fried and
    !     Conte, The Plasma Dispersion Function, 1961).
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_special, only: faddeeva_w, plasma_dispersion_z, &
        plasma_dispersion_z_derivative, dawson
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: tol = 2.0e-14_dp
    real(dp), parameter :: radii(5) = [8.0_dp, 12.5_dp, 40.0_dp, 3.0e3_dp, &
        2.0e7_dp]
    complex(dp), parameter :: lower(6) = [(8.0_dp, -2.0_dp), &
        (-9.0_dp, -1.5_dp), (3.0_dp, -7.0_dp), (-5.0_dp, -9.0_dp), &
        (20.0_dp, -19.0_dp), (0.0_dp, -10.0_dp)]
    complex(dp), parameter :: zetas(6) = [(0.0_dp, 0.0_dp), (0.7_dp, 0.2_dp), &
        (2.0_dp, -0.2_dp), (-1.3_dp, -1.1_dp), (4.5_dp, 3.0_dp), &
        (6.2_dp, -0.05_dp)]
    integer :: nfail, i, j, k
    real(dp) :: x, radius, angle
    complex(dp) :: z, ref, omega, zeta, f, df
    complex(qp) :: zq
    logical :: ok

    nfail = 0
    do i = -10, 10
        do j = -10, 10
            z = cmplx(0.5_dp*i, 0.5_dp*j, dp)
            if (abs(z) > 5.0_dp) cycle
            call series_oracle(cmplx(z, kind=qp), ref, ok)
            if (ok) call check("series", z, faddeeva_w(z), ref)
        end do
    end do

    do k = 0, 23
        angle = pi*k/23.0_dp
        do i = 1, 5
            radius = radii(i)
            z = radius*cmplx(cos(angle), sin(angle), dp)
            call check("asymptotic", z, faddeeva_w(z), &
                cmplx(asymptotic_oracle(cmplx(z, kind=qp)), kind=dp))
        end do
    end do
    do i = 1, 6
        z = lower(i)
        zq = cmplx(z, kind=qp)
        ref = cmplx(2.0_qp*exp(-zq*zq) - asymptotic_oracle(-zq), kind=dp)
        call check("lower half plane", z, faddeeva_w(z), ref)
    end do

    do i = -60, 60
        x = 0.5_dp*i
        z = cmplx(x, 0.0_dp, dp)
        call check("real axis", z, faddeeva_w(z), &
            cmplx(exp(-x*x), 2.0_dp/sqrt(pi)*dawson(x), dp))
    end do

    call check("Z(0)", (0.0_dp, 0.0_dp), plasma_dispersion_z((0.0_dp, 0.0_dp)), &
        cmplx(0.0_dp, sqrt(pi), dp))
    do i = 1, 6
        zeta = zetas(i)
        ! Z' = -2(1 + zeta Z) cancels like 1/(2 zeta^2) at large zeta.
        call check("Z' Cauchy", zeta, plasma_dispersion_z_derivative(zeta), &
            cauchy_derivative(zeta), 1.0e-13_dp)
    end do

    omega = (1.4_dp, -0.15_dp)
    do k = 1, 50
        zeta = omega/(sqrt(2.0_dp)*0.5_dp)
        f = 1.0_dp - plasma_dispersion_z_derivative(zeta)/(2.0_dp*0.25_dp)
        df = 2.0_dp*(plasma_dispersion_z(zeta) + zeta* &
            plasma_dispersion_z_derivative(zeta))/(2.0_dp*0.25_dp)/ &
            (sqrt(2.0_dp)*0.5_dp)
        omega = omega - f/df
        if (abs(f/df) < 1.0e-14_dp) exit
    end do
    if (.not. (abs(real(omega) - 1.4156_dp) < 1.0e-4_dp .and. &
            abs(aimag(omega) + 0.1533_dp) < 1.0e-4_dp)) then
        nfail = nfail + 1
        write (error_unit, "(a,2es24.16)") "Langmuir root", omega
    end if

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine series_oracle(zq, w, ok)
        complex(qp), intent(in) :: zq
        complex(dp), intent(out) :: w
        logical, intent(out) :: ok
        complex(qp) :: term_even, term_odd, total
        real(qp) :: magnitude
        integer :: n

        ! (iz)^n/Gamma(n/2+1), even and odd chains with Gamma(n/2+1)
        ! advanced by n/2 + 1 every two steps.
        term_even = (1.0_qp, 0.0_qp)
        term_odd = (0.0_qp, 1.0_qp)*zq/(0.5_qp*sqrt(acos(-1.0_qp)))
        total = term_even + term_odd
        magnitude = abs(term_even) + abs(term_odd)
        do n = 2, 600, 2
            term_even = term_even*(0.0_qp, 1.0_qp)**2*zq*zq/real(n/2, qp)
            term_odd = term_odd*(0.0_qp, 1.0_qp)**2*zq*zq/(real(n + 1, qp)/2.0_qp)
            total = total + term_even + term_odd
            magnitude = magnitude + abs(term_even) + abs(term_odd)
        end do
        w = cmplx(total, kind=dp)
        ok = magnitude < 1.0e12_qp*abs(total)
    end subroutine series_oracle

    function asymptotic_oracle(zq) result(w)
        ! i/(sqrt(pi) z) sum_k (2k-1)!!/(2 z^2)^k, stopped at the smallest
        ! term (optimal truncation, error about that term).
        complex(qp), intent(in) :: zq
        complex(qp) :: w
        complex(qp) :: term, total
        real(qp) :: previous
        integer :: kk

        term = (1.0_qp, 0.0_qp)
        total = term
        previous = huge(1.0_qp)
        do kk = 1, 400
            term = term*real(2*kk - 1, qp)/(2.0_qp*zq*zq)
            if (abs(term) >= previous .or. abs(term) < 1.0e-40_qp*abs(total)) &
                exit
            previous = abs(term)
            total = total + term
        end do
        w = (0.0_qp, 1.0_qp)/(sqrt(acos(-1.0_qp))*zq)*total
    end function asymptotic_oracle

    function cauchy_derivative(zeta0) result(d)
        ! Trapezoidal Cauchy integral on |t - zeta0| = 0.1: spectrally
        ! accurate for the entire function Z.
        complex(dp), intent(in) :: zeta0
        complex(dp) :: d
        complex(dp) :: dz
        integer :: m

        d = (0.0_dp, 0.0_dp)
        do m = 0, 63
            dz = 0.1_dp*exp(cmplx(0.0_dp, 2.0_dp*pi*m/64.0_dp, dp))
            d = d + plasma_dispersion_z(zeta0 + dz)/dz
        end do
        d = d/64.0_dp
    end function cauchy_derivative

    subroutine check(label, zin, got, expected, custom_tol)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: zin, got, expected
        real(dp), intent(in), optional :: custom_tol
        real(dp) :: t
        t = tol
        if (present(custom_tol)) t = custom_tol
        if (.not. abs(got - expected) <= t*abs(expected) + tiny(1.0_dp)) then
            nfail = nfail + 1
            write (error_unit, "(a,' z=',2es12.4,' got=',2es24.16,' rel=',es9.2)") &
                label, zin, got, abs(got - expected)/abs(expected)
        end if
    end subroutine check

end program test_fortnum_special_faddeeva
