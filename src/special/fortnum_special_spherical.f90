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
    public :: spherical_harmonic_product_coefficient

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

    pure function spherical_harmonic_product_coefficient( &
            degree_1, order_1, degree_2, order_2, degree_out, order_out) &
            result(value)
        ! Return the coefficient of Y_degree_out^order_out in
        ! Y_degree_1^order_1 * Y_degree_2^order_2.
        !
        ! The Gaunt coefficient is evaluated from the two Wigner-3j symbols
        ! in DLMF 34.3. The finite Racah sum is independent of the pointwise
        ! Legendre implementation and is intended for moderate degrees.
        integer, intent(in) :: degree_1, order_1, degree_2, order_2
        integer, intent(in) :: degree_out, order_out
        real(dp) :: value
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
        real(dp) :: first_symbol, second_symbol, prefactor

        value = ieee_value(0.0_dp, ieee_quiet_nan)
        if (.not. valid_product_indices( &
            degree_1, order_1, degree_2, order_2, degree_out, order_out)) then
            return
        end if
        if (order_out /= order_1 + order_2) then
            value = 0.0_dp
            return
        end if
        if (degree_out < abs(degree_1 - degree_2) .or. &
            degree_out > degree_1 + degree_2 .or. &
            mod(degree_1 + degree_2 + degree_out, 2) /= 0) then
            value = 0.0_dp
            return
        end if
        first_symbol = wigner_3j(degree_1, degree_2, degree_out, &
            0, 0, 0)
        second_symbol = wigner_3j(degree_1, degree_2, degree_out, &
            order_1, order_2, -order_out)
        prefactor = sqrt((2.0_dp*real(degree_1, dp) + 1.0_dp)* &
            (2.0_dp*real(degree_2, dp) + 1.0_dp)* &
            (2.0_dp*real(degree_out, dp) + 1.0_dp)/(4.0_dp*pi))
        if (mod(abs(order_out), 2) == 0) then
            value = prefactor*first_symbol*second_symbol
        else
            value = -prefactor*first_symbol*second_symbol
        end if
    end function spherical_harmonic_product_coefficient

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

    pure function wigner_3j(j_one, j_two, j_three, m_one, m_two, m_three) &
            result(value)
        integer, intent(in) :: j_one, j_two, j_three
        integer, intent(in) :: m_one, m_two, m_three
        real(dp) :: value
        real(dp) :: log_scale, maximum_log_term, term_sum, log_term
        integer :: lower, upper, k

        value = 0.0_dp
        if (m_one + m_two + m_three /= 0) return
        if (j_three < abs(j_one - j_two) .or. &
            j_three > j_one + j_two) return
        if (abs(m_one) > j_one .or. abs(m_two) > j_two .or. &
            abs(m_three) > j_three) return
        lower = max(0, j_two - j_three - m_one, j_one - j_three + m_two)
        upper = min(j_one + j_two - j_three, j_one - m_one, j_two + m_two)
        if (lower > upper) return

        maximum_log_term = -huge(1.0_dp)
        do k = lower, upper
            log_term = -log_factorial(k) - &
                log_factorial(j_one + j_two - j_three - k) - &
                log_factorial(j_one - m_one - k) - &
                log_factorial(j_two + m_two - k) - &
                log_factorial(j_three - j_two + m_one + k) - &
                log_factorial(j_three - j_one - m_two + k)
            maximum_log_term = max(maximum_log_term, log_term)
        end do
        term_sum = 0.0_dp
        do k = lower, upper
            log_term = -log_factorial(k) - &
                log_factorial(j_one + j_two - j_three - k) - &
                log_factorial(j_one - m_one - k) - &
                log_factorial(j_two + m_two - k) - &
                log_factorial(j_three - j_two + m_one + k) - &
                log_factorial(j_three - j_one - m_two + k)
            if (mod(k, 2) == 0) then
                term_sum = term_sum + exp(log_term - maximum_log_term)
            else
                term_sum = term_sum - exp(log_term - maximum_log_term)
            end if
        end do
        log_scale = 0.5_dp*(log_factorial(j_one + j_two - j_three) + &
            log_factorial(j_one - j_two + j_three) + &
            log_factorial(-j_one + j_two + j_three) - &
            log_factorial(j_one + j_two + j_three + 1) + &
            log_factorial(j_one + m_one) + log_factorial(j_one - m_one) + &
            log_factorial(j_two + m_two) + log_factorial(j_two - m_two) + &
            log_factorial(j_three + m_three) + log_factorial(j_three - m_three))
        value = exp(log_scale + maximum_log_term)*term_sum
        if (mod(j_one - j_two - m_three, 2) /= 0) value = -value
    end function wigner_3j

    pure function log_factorial(integer_value) result(value)
        integer, intent(in) :: integer_value
        real(dp) :: value

        value = log_gamma(real(integer_value + 1, dp))
    end function log_factorial

    pure elemental function valid_product_indices( &
            degree_1, order_1, degree_2, order_2, degree_out, order_out) &
            result(valid)
        integer, intent(in) :: degree_1, order_1, degree_2, order_2
        integer, intent(in) :: degree_out, order_out
        logical :: valid

        valid = degree_1 >= 0 .and. abs(order_1) <= degree_1 .and. &
            degree_2 >= 0 .and. abs(order_2) <= degree_2 .and. &
            degree_out >= 0 .and. abs(order_out) <= degree_out
    end function valid_product_indices

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
