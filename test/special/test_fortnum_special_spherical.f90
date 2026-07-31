program test_fortnum_special_spherical
    ! Independent pointwise checks for normalized complex spherical harmonics.
    ! The degree-one values and their angular derivatives are closed forms,
    ! rather than a second evaluation of the Legendre recurrence.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special, only: &
        spherical_harmonic, spherical_harmonic_phi_derivative, &
        spherical_harmonic_theta_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: tol = 3.0e-13_dp
    real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
    real(dp), parameter :: theta = 0.9_dp
    real(dp), parameter :: phi = 0.4_dp
    real(dp) :: amplitude
    complex(dp) :: expected

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
        -conjg(spherical_harmonic(1, 1, theta, phi))) <= tol)
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

end program test_fortnum_special_spherical
