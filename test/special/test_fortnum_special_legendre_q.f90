program test_fortnum_special_legendre_q
    ! Behavioral checks for the real Legendre functions of the second kind.
    ! The values through degree three are closed forms, independent of the
    ! recurrence used by the implementation.  Derivatives are also checked
    ! against a centered difference of the primal routine.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special, only: legendre_q, legendre_q_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: tol = 3.0e-13_dp
    real(dp), parameter :: x = 2.3_dp
    real(dp) :: q0, q1, q2, q3, h, finite_difference

    nfail = 0
    q0 = 0.5_dp*log((x + 1.0_dp)/(x - 1.0_dp))
    q1 = x*q0 - 1.0_dp
    q2 = 0.5_dp*(3.0_dp*x*x - 1.0_dp)*q0 - 1.5_dp*x
    q3 = 0.5_dp*(5.0_dp*x*x*x - 3.0_dp*x)*q0 - &
        2.5_dp*x*x + 2.0_dp/3.0_dp

    call check("Q_0", legendre_q(0, x), q0)
    call check("Q_1", legendre_q(1, x), q1)
    call check("Q_2", legendre_q(2, x), q2)
    call check("Q_3", legendre_q(3, x), q3)

    call check("dQ_0", legendre_q_derivative(0, x), -1.0_dp/(x*x - 1.0_dp))
    call check("dQ_1", legendre_q_derivative(1, x), q0 + &
        x*legendre_q_derivative(0, x))

    h = 2.0e-5_dp
    finite_difference = (legendre_q(4, x + h) - legendre_q(4, x - h))/(2.0_dp*h)
    call check("dQ_4 centered difference", legendre_q_derivative(4, x), &
        finite_difference, 2.0e-9_dp)

    call check_true("degree below zero", ieee_is_nan(legendre_q(-1, x)))
    call check_true("singular endpoint", ieee_is_nan(legendre_q(0, 1.0_dp)))
    call check_true("wrong side of cut", ieee_is_nan(legendre_q(2, 0.5_dp)))
    call check_true("invalid derivative degree", &
        ieee_is_nan(legendre_q_derivative(-1, x)))

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check(label, got, expected, custom_tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        real(dp), intent(in), optional :: custom_tol
        real(dp) :: local_tol

        local_tol = tol
        if (present(custom_tol)) local_tol = custom_tol
        if (.not. abs(got - expected) <= local_tol*(1.0_dp + abs(expected))) then
            nfail = nfail + 1
            write (error_unit, "(a,2(a,es22.14))") &
                "FAIL: "//label, " got=", got, " expected=", expected
        end if
    end subroutine check

    subroutine check_true(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition

        if (.not. condition) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: "//label
        end if
    end subroutine check_true

end program test_fortnum_special_legendre_q
