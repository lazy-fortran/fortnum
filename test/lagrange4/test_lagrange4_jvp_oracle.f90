program test_lagrange4_jvp_oracle
    !! Independent-oracle check of the committed Fortran Lagrange-4 JVP leaf.
    !!
    !! fortnum_lagrange4_jvp_kernel interpolates the cubic that passes
    !! through the four samples at nodes {-1, 0, 1, 2} and returns its JVP
    !! along the tangent direction (tx, ty1..ty4). This test evaluates the
    !! interpolant and its directional derivative directly from the closed
    !! form (cubic coefficients solved from the samples) and requires the
    !! generated leaf to agree to the 1e-13 relative tolerance the kernel's
    !! own equivalence test uses.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_generated_lagrange4_jvp, only: fortnum_lagrange4_jvp_kernel
    implicit none

    integer, parameter :: n = 4096
    real(dp), parameter :: tolerance = 1.0e-13_dp
    real(dp), parameter :: nodes(4) = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    real(dp) :: x, samples(4), tx, ty(4)
    real(dp) :: leaf_value, leaf_jvp, expected_value, expected_jvp
    real(dp) :: value_scale, jvp_scale
    real(dp) :: c0, c1, c2, c3
    integer :: i, j, failures

    failures = 0
    do i = 1, n
        x = -0.8_dp + 2.6_dp*real(mod(17*i, 4093), dp)/4092.0_dp
        tx = -0.7_dp + 1.4_dp*real(mod(19*i, 4091), dp)/4090.0_dp
        do j = 1, 4
            samples(j) = primal_polynomial(nodes(j), i)
            ty(j) = tangent_polynomial(nodes(j), i)
        end do

        call fortnum_lagrange4_jvp_kernel( &
            x, samples(1), samples(2), samples(3), samples(4), &
            tx, ty(1), ty(2), ty(3), ty(4), leaf_value, leaf_jvp)

        ! Independent oracle: solve for the cubic interpolant coefficients
        ! from the four samples, then differentiate.
        call cubic_coefficients(samples, c0, c1, c2, c3)
        expected_value = c0 + x*(c1 + x*(c2 + x*c3))
        expected_jvp = tx*(c1 + x*(2.0_dp*c2 + 3.0_dp*x*c3)) + &
            tangent_polynomial(x, i)

        value_scale = max(1.0_dp, abs(expected_value))
        jvp_scale = max(1.0_dp, abs(expected_jvp))
        if (abs(leaf_value - expected_value) > tolerance*value_scale) then
            write (error_unit, "(a,i0,3es24.16)") &
                "FAIL value ", i, leaf_value, expected_value, &
                abs(leaf_value - expected_value)
            failures = failures + 1
        end if
        if (abs(leaf_jvp - expected_jvp) > tolerance*jvp_scale) then
            write (error_unit, "(a,i0,3es24.16)") &
                "FAIL jvp ", i, leaf_jvp, expected_jvp, &
                abs(leaf_jvp - expected_jvp)
            failures = failures + 1
        end if
    end do

    if (failures /= 0) then
        write (error_unit, "(a,i0,a)") &
            "test_lagrange4_jvp_oracle: ", failures, " mismatches"
        error stop 1
    end if
    print *, "test_lagrange4_jvp_oracle: Fortran leaf agrees with cubic oracle"

contains

    pure function primal_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 17), dp)/31.0_dp
        value = (0.3_dp + shift) - 0.4_dp*z + 0.2_dp*z*z - 0.05_dp*z*z*z
    end function primal_polynomial

    pure function tangent_polynomial(z, index) result(value)
        real(dp), intent(in) :: z
        integer, intent(in) :: index
        real(dp) :: value, shift

        shift = real(mod(index, 13), dp)/29.0_dp
        value = (-0.2_dp + shift) + 0.3_dp*z - 0.1_dp*z*z + 0.04_dp*z*z*z
    end function tangent_polynomial

    pure subroutine cubic_coefficients(samples, c0, c1, c2, c3)
        !! Solve the 4x4 Vandermonde system for the cubic interpolant through
        !! nodes {-1, 0, 1, 2}. Closed form derived from the Lagrange basis.
        real(dp), intent(in) :: samples(4)
        real(dp), intent(out) :: c0, c1, c2, c3

        ! value(z) = y1*L1 + y2*L2 + y3*L3 + y4*L4 with
        !   L1 = z(z-1)(z-2)/(-6) = -z^3/6 + z^2/2 - z/3
        !   L2 = (z+1)(z-1)(z-2)/2  = z^3/2 - z^2 - z/2 + 1
        !   L3 = -(z+1)z(z-2)/2     = -z^3/2 + z^2/2 + z
        !   L4 = (z+1)z(z-1)/6      = z^3/6 - z/6
        c0 = samples(2)
        c1 = -samples(1)/3.0_dp - samples(2)/2.0_dp + samples(3) &
            - samples(4)/6.0_dp
        c2 = samples(1)/2.0_dp - samples(2) + samples(3)/2.0_dp
        c3 = -samples(1)/6.0_dp + samples(2)/2.0_dp - samples(3)/2.0_dp &
            + samples(4)/6.0_dp
    end subroutine cubic_coefficients

end program test_lagrange4_jvp_oracle
