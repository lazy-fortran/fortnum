program test_fortnum_special_toroidal
    ! High-precision oracle values use mpmath legenp/legenq(type=3), 50 digits.
    ! The ODE residual is an independent centered-difference behavioral check.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_toroidal, only: &
        toroidal_p, toroidal_q, &
        toroidal_p_derivative, toroidal_q_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: value_tol = 3.0e-12_dp

    nfail = 0
    call check("P_-1/2^0(2)", toroidal_p(0, 0, 2.0_dp), &
        0.90128629936044729874_dp)
    call check("Q_-1/2^0(2)", toroidal_q(0, 0, 2.0_dp), &
        1.6566381702365941664_dp)
    call check("P_3/2^0(2)", toroidal_p(2, 0, 2.0_dp), &
        3.2439396660408049155_dp)
    call check("Q_3/2^0(2)", toroidal_q(2, 0, 2.0_dp), &
        0.045158724151576976637_dp)
    call check("P_-1/2^1(2)", toroidal_p(0, 1, 2.0_dp), &
        -0.13666874968871549533_dp)
    call check("Q_-1/2^1(2)", toroidal_q(0, 1, 2.0_dp), &
        -0.89179313740019260390_dp)
    call check("P_5/2^2(3)", toroidal_p(3, 2, 3.0_dp), &
        100.13371998025088164_dp)
    call check("Q_5/2^2(3)", toroidal_q(3, 2, 3.0_dp), &
        0.033773724240355700700_dp)
    call check("dP_5/2^2(2)", toroidal_p_derivative(3, 2, 2.0_dp), &
        48.525686151003148510_dp)
    call check("dQ_3/2^1(2.5)", &
        toroidal_q_derivative(2, 1, 2.5_dp), &
        0.067949860642696846331_dp)

    call check_finite("near-cut Q_-1/2^0", toroidal_q(0, 0, 1.000001_dp))
    call check("near-cut Q_-1/2^0", toroidal_q(0, 0, 1.000001_dp), &
        8.6406222753454926764_dp)
    call check_finite("near-cut Q_-1/2^1", toroidal_q(0, 1, 1.000001_dp))
    call check_scaled("near-cut Q_-1/2^1", &
        toroidal_q(0, 1, 1.000001_dp), -707.10822028768059246_dp, 5.0e-11_dp)
    call check_finite("near-cut Q_5/2^2", toroidal_q(3, 2, 1.0001_dp))
    call check_scaled("near-cut Q_5/2^2", &
        toroidal_q(3, 2, 1.0001_dp), 9996.1311322726379956_dp, 5.0e-11_dp)

    call check_ode("P toroidal ODE", .true., 3, 2, 2.0_dp)
    call check_ode("Q toroidal ODE", .false., 2, 1, 2.5_dp)

    call check_finite("large-degree P", toroidal_p(80, 4, 1.2_dp))
    call check_finite("large-degree Q", toroidal_q(80, 4, 1.2_dp))
    call check_scaled("large-degree P reference", &
        toroidal_p(80, 4, 1.2_dp), 7.8574629022740074243e27_dp, 5.0e-12_dp)
    call check_scaled("large-degree Q reference", &
        toroidal_q(80, 4, 1.2_dp), 1.9996464751570121823e-15_dp, 5.0e-10_dp)
    call check_degree_recurrence("large-degree P recurrence", .true., &
        80, 4, 1.2_dp)
    call check_degree_recurrence("large-degree Q recurrence", .false., &
        80, 4, 1.2_dp)

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected

        if (.not. abs(got - expected) <= &
            value_tol*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check

    subroutine check_scaled(label, got, expected, relative_tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, relative_tolerance

        if (abs(got - expected) > relative_tolerance*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check_scaled

    subroutine check_finite(label, got)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got

        if (.not. ieee_is_finite(got)) then
            nfail = nfail + 1
            write (error_unit, "(a,a)") "FAIL: ", label//" is not finite"
        end if
    end subroutine check_finite

    subroutine check_ode(label, first_kind, degree_index, order, x)
        character(*), intent(in) :: label
        logical, intent(in) :: first_kind
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: x
        real(dp), parameter :: h = 2.0e-4_dp
        real(dp) :: center, left, right, first, second, degree, residual

        if (first_kind) then
            left = toroidal_p(degree_index, order, x - h)
            center = toroidal_p(degree_index, order, x)
            right = toroidal_p(degree_index, order, x + h)
        else
            left = toroidal_q(degree_index, order, x - h)
            center = toroidal_q(degree_index, order, x)
            right = toroidal_q(degree_index, order, x + h)
        end if
        first = (right - left)/(2.0_dp*h)
        second = (right - 2.0_dp*center + left)/(h*h)
        degree = real(degree_index, dp) - 0.5_dp
        residual = (1.0_dp - x*x)*second - 2.0_dp*x*first + &
            (degree*(degree + 1.0_dp) - &
            real(order*order, dp)/(1.0_dp - x*x))*center
        if (abs(residual) > 2.0e-5_dp*(1.0_dp + abs(center))) then
            nfail = nfail + 1
            write (error_unit, "(a,a,es22.14)") &
                "FAIL: ", label//" residual=", residual
        end if
    end subroutine check_ode

    subroutine check_degree_recurrence(label, first_kind, degree_index, &
            order, x)
        character(*), intent(in) :: label
        logical, intent(in) :: first_kind
        integer, intent(in) :: degree_index, order
        real(dp), intent(in) :: x
        real(dp) :: previous, current, following, degree, residual

        if (first_kind) then
            previous = toroidal_p(degree_index - 1, order, x)
            current = toroidal_p(degree_index, order, x)
            following = toroidal_p(degree_index + 1, order, x)
        else
            previous = toroidal_q(degree_index - 1, order, x)
            current = toroidal_q(degree_index, order, x)
            following = toroidal_q(degree_index + 1, order, x)
        end if
        degree = real(degree_index, dp) - 0.5_dp
        residual = (degree - real(order, dp) + 1.0_dp)*following - &
            (2.0_dp*degree + 1.0_dp)*x*current + &
            (degree + real(order, dp))*previous
        if (abs(residual) > 2.0e-12_dp*(1.0_dp + abs(current))) then
            nfail = nfail + 1
            write (error_unit, "(a,a,es22.14)") &
                "FAIL: ", label//" residual=", residual
        end if
    end subroutine check_degree_recurrence

end program test_fortnum_special_toroidal
