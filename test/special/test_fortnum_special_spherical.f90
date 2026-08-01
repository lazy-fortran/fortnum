program test_fortnum_special_spherical
    ! Independent pointwise checks for normalized complex spherical harmonics.
    ! The degree-one values and their angular derivatives are closed forms,
    ! rather than a second evaluation of the Legendre recurrence.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special, only: &
        spherical_harmonic, spherical_harmonic_phi_derivative, &
        spherical_harmonic_product_coefficient, &
        spherical_harmonic_theta_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: tol = 3.0e-13_dp
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: theta = 0.9_dp
    real(dp), parameter :: phi = 0.4_dp
    real(dp) :: amplitude
    complex(dp) :: expected
    real(dp) :: coefficient

    nfail = 0
    call check_complex("Y_0^0", spherical_harmonic(0, 0, theta, phi), &
        cmplx(1.0_dp/sqrt(4.0_dp*pi), 0.0_dp, dp))
    call check_complex("Y_1^0", spherical_harmonic(1, 0, theta, phi), &
        cmplx(sqrt(3.0_dp/(4.0_dp*pi))*cos(theta), 0.0_dp, dp))

    amplitude = sqrt(3.0_dp/(8.0_dp*pi))*sin(theta)
    expected = cmplx(-amplitude*cos(phi), -amplitude*sin(phi), dp)
    call check_complex("Y_1^1", spherical_harmonic(1, 1, theta, phi), expected)
    expected = cmplx(amplitude*cos(phi), -amplitude*sin(phi), dp)
    call check_complex("Y_1^-1", spherical_harmonic(1, -1, theta, phi), expected)

    expected = cmplx(-sqrt(3.0_dp/(8.0_dp*pi))*cos(theta)*cos(phi), &
        -sqrt(3.0_dp/(8.0_dp*pi))*cos(theta)*sin(phi), dp)
    call check_complex("dY_1^1/dtheta", &
        spherical_harmonic_theta_derivative(1, 1, theta, phi), expected)
    call check_complex("dY_1^1/dphi", &
        spherical_harmonic_phi_derivative(1, 1, theta, phi), &
        cmplx(0.0_dp, 1.0_dp, dp)*spherical_harmonic(1, 1, theta, phi))

    call check_true("negative-order conjugacy", &
        abs(spherical_harmonic(1, -1, theta, phi) - &
        (-conjg(spherical_harmonic(1, 1, theta, phi)))) <= tol)
    call check_complex("north-pole Y_2^2", spherical_harmonic(2, 2, 0.0_dp, phi), &
        cmplx(0.0_dp, 0.0_dp, dp))

    call check_true("negative degree", &
        complex_is_nan(spherical_harmonic(-1, 0, theta, phi)))
    call check_true("order out of range", &
        complex_is_nan(spherical_harmonic(1, 2, theta, phi)))
    call check_true("theta below domain", &
        complex_is_nan(spherical_harmonic(1, 0, -1.0e-6_dp, phi)))
    call check_true("theta above domain", &
        complex_is_nan(spherical_harmonic(1, 0, pi + 1.0e-6_dp, phi)))

    coefficient = spherical_harmonic_product_coefficient(0, 0, 0, 0, 0, 0)
    call check_real("Y00*Y00 -> Y00", coefficient, 1.0_dp/sqrt(4.0_dp*pi))
    coefficient = spherical_harmonic_product_coefficient(1, 0, 1, 0, 0, 0)
    call check_real("Y10*Y10 -> Y00", coefficient, 1.0_dp/sqrt(4.0_dp*pi))
    call check_true("product m selection rule", &
        spherical_harmonic_product_coefficient(2, 1, 1, -1, 3, 1) == 0.0_dp)
    call check_true("product parity selection rule", &
        spherical_harmonic_product_coefficient(1, 0, 1, 0, 1, 0) == 0.0_dp)
    call check_true("product rejects invalid indices", &
        ieee_is_nan(spherical_harmonic_product_coefficient(1, 2, 1, 0, 1, 0)))
    call check_quadrature(2, 1, 1, -1, 3, 0)
    call check_product_reconstruction(1, 1, 1, -1, 0.73_dp, 0.41_dp)

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_complex(label, got, expected, custom_tol)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: got, expected
        real(dp), intent(in), optional :: custom_tol
        real(dp) :: local_tol

        local_tol = tol
        if (present(custom_tol)) local_tol = custom_tol
        if (.not. abs(got - expected) <= local_tol*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", abs(got), " expected=", abs(expected)
        end if
    end subroutine check_complex

    logical function complex_is_nan(value)
        complex(dp), intent(in) :: value

        complex_is_nan = ieee_is_nan(real(value, dp)) .or. &
            ieee_is_nan(aimag(value))
    end function complex_is_nan

    subroutine check_true(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: "//label
        end if
    end subroutine check_true

    subroutine check_real(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected

        if (abs(got - expected) > 3.0e-12_dp*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check_real

    subroutine check_quadrature(degree_1, order_1, degree_2, order_2, &
            degree_out, order_out)
        integer, intent(in) :: degree_1, order_1, degree_2, order_2
        integer, intent(in) :: degree_out, order_out
        real(dp), parameter :: nodes(8) = [ &
            -0.9602898564975363_dp, -0.7966664774136267_dp, &
            -0.5255324099163290_dp, -0.1834346424956498_dp, &
            0.1834346424956498_dp,  0.5255324099163290_dp, &
            0.7966664774136267_dp,  0.9602898564975363_dp]
        real(dp), parameter :: weights(8) = [ &
            0.1012285362903763_dp, 0.2223810344533745_dp, &
            0.3137066458778873_dp, 0.3626837833783620_dp, &
            0.3626837833783620_dp, 0.3137066458778873_dp, &
            0.2223810344533745_dp, 0.1012285362903763_dp]
        integer, parameter :: phi_count = 32
        real(dp), parameter :: two_pi = 2.0_dp*pi
        real(dp) :: quadrature, theta, azimuth, dphi
        complex(dp) :: integrand
        integer :: node, angle

        quadrature = 0.0_dp
        dphi = two_pi/real(phi_count, dp)
        do node = 1, size(nodes)
            theta = acos(nodes(node))
            do angle = 0, phi_count - 1
                azimuth = dphi*real(angle, dp)
                integrand = spherical_harmonic(degree_1, order_1, theta, azimuth)* &
                    spherical_harmonic(degree_2, order_2, theta, azimuth)* &
                    conjg(spherical_harmonic(degree_out, order_out, theta, azimuth))
                quadrature = quadrature + weights(node)*dphi*real(integrand, dp)
            end do
        end do
        call check_real("Gaunt tensor-product quadrature", quadrature, &
            spherical_harmonic_product_coefficient( &
            degree_1, order_1, degree_2, order_2, degree_out, order_out))
    end subroutine check_quadrature

    subroutine check_product_reconstruction( &
            degree_1, order_1, degree_2, order_2, theta, azimuth)
        integer, intent(in) :: degree_1, order_1, degree_2, order_2
        real(dp), intent(in) :: theta, azimuth
        complex(dp) :: reconstructed, direct
        integer :: degree_out, order_out

        order_out = order_1 + order_2
        reconstructed = cmplx(0.0_dp, 0.0_dp, dp)
        do degree_out = abs(degree_1 - degree_2), degree_1 + degree_2
            reconstructed = reconstructed + &
                spherical_harmonic_product_coefficient( &
                degree_1, order_1, degree_2, order_2, degree_out, order_out)* &
                spherical_harmonic(degree_out, order_out, theta, azimuth)
        end do
        direct = spherical_harmonic(degree_1, order_1, theta, azimuth)* &
            spherical_harmonic(degree_2, order_2, theta, azimuth)
        call check_true("spherical product reconstruction", &
            abs(reconstructed - direct) < 2.0e-12_dp)
    end subroutine check_product_reconstruction

end program test_fortnum_special_spherical
