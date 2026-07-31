module fortnum_special_spherical
    ! Standard orthonormal complex spherical harmonics Y_l^m(theta, phi).
    !
    ! The Condon-Shortley phase and the standard orthonormal normalization
    ! follow DLMF 14.30. The angular derivatives are analytic derivatives
    ! with respect to theta and phi; they are defined away from the poles.
    use, intrinsic :: ieee_arithmetic, only: &
        ieee_is_finite, ieee_quiet_nan, ieee_value
    use fortnum_kinds, only: dp
    use fortnum_special_legendre, only: &
        legendre_p, legendre_p_derivative
    implicit none
    private

    public :: spherical_harmonic
    public :: spherical_harmonic_theta_derivative
    public :: spherical_harmonic_phi_derivative

contains

    elemental function spherical_harmonic( &
            degree, order, theta, phi) result(value)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: theta, phi
        complex(dp) :: value
        complex(dp) :: theta_derivative

        call evaluate_spherical_harmonic( &
            degree, order, theta, phi, value, theta_derivative)
    end function spherical_harmonic

    elemental function spherical_harmonic_theta_derivative( &
            degree, order, theta, phi) result(value)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: theta, phi
        complex(dp) :: value
        complex(dp) :: function_value

        call evaluate_spherical_harmonic( &
            degree, order, theta, phi, function_value, value)
    end function spherical_harmonic_theta_derivative

    elemental function spherical_harmonic_phi_derivative( &
            degree, order, theta, phi) result(value)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: theta, phi
        complex(dp) :: value
        complex(dp) :: function_value

        function_value = spherical_harmonic(degree, order, theta, phi)
        value = cmplx(0.0_dp, real(order, dp), dp)*function_value
    end function spherical_harmonic_phi_derivative

    pure elemental subroutine evaluate_spherical_harmonic( &
            degree, order, theta, phi, value, theta_derivative)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: theta, phi
        complex(dp), intent(out) :: value, theta_derivative
        integer :: absolute_order
        real(dp) :: cosine, sine, normalization, p_value, p_derivative
        complex(dp) :: phase, positive_value, positive_theta

        if (.not. valid_arguments(degree, order, theta, phi)) then
            value = complex_nan(theta)
            theta_derivative = value
            return
        end if

        absolute_order = abs(order)
        cosine = min(1.0_dp, max(-1.0_dp, cos(theta)))
        sine = sin(theta)
        normalization = spherical_normalization(degree, absolute_order)
        p_value = legendre_p(degree, absolute_order, cosine)
        phase = cmplx( &
            cos(real(absolute_order, dp)*phi), &
            sin(real(absolute_order, dp)*phi), dp)
        positive_value = normalization*p_value*phase
        if (sine == 0.0_dp) then
            positive_theta = complex_nan(theta)
        else
            p_derivative = legendre_p_derivative( &
                degree, absolute_order, cosine)
            positive_theta = normalization*(-sine*p_derivative)*phase
        end if

        if (order >= 0) then
            value = positive_value
            theta_derivative = positive_theta
        else
            if (mod(absolute_order, 2) == 1) then
                value = -conjg(positive_value)
                theta_derivative = -conjg(positive_theta)
            else
                value = conjg(positive_value)
                theta_derivative = conjg(positive_theta)
            end if
        end if
    end subroutine evaluate_spherical_harmonic

    pure elemental function spherical_normalization(degree, order) result(value)
        integer, intent(in) :: degree, order
        real(dp) :: value
        integer :: k
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

        value = (2.0_dp*real(degree, dp) + 1.0_dp)/(4.0_dp*pi)
        do k = degree - order + 1, degree + order
            value = value/real(k, dp)
        end do
        value = sqrt(value)
    end function spherical_normalization

    pure elemental function valid_arguments( &
            degree, order, theta, phi) result(valid)
        integer, intent(in) :: degree, order
        real(dp), intent(in) :: theta, phi
        logical :: valid
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

        valid = degree >= 0 .and. abs(order) <= degree .and. &
            ieee_is_finite(theta) .and. ieee_is_finite(phi) .and. &
            theta >= 0.0_dp .and. theta <= pi
    end function valid_arguments

    pure elemental function complex_nan(reference) result(value)
        real(dp), intent(in) :: reference
        complex(dp) :: value
        real(dp) :: nan

        nan = ieee_value(reference, ieee_quiet_nan)
        value = cmplx(nan, nan, dp)
    end function complex_nan

end module fortnum_special_spherical
