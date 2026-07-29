program test_fortnum_special_legendre
    ! Independent behavioral checks from Rodrigues' formula and parity.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use fortnum_special_legendre, only: legendre_p, legendre_p_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: tol = 2.0e-13_dp
    real(dp), parameter :: x = 0.37_dp

    nfail = 0

    call check("P_0", legendre_p(0, 0, x), 1.0_dp)
    call check("P_1", legendre_p(1, 0, x), x)
    call check("P_2", legendre_p(2, 0, x), 0.5_dp*(3.0_dp*x*x - 1.0_dp))
    call check("P_3", legendre_p(3, 0, x), &
        0.5_dp*(5.0_dp*x*x*x - 3.0_dp*x))
    call check("P_2^1", legendre_p(2, 1, x), &
        -3.0_dp*x*sqrt(1.0_dp - x*x))
    call check("P_3^2", legendre_p(3, 2, x), &
        15.0_dp*x*(1.0_dp - x*x))
    call check("P_4^4", legendre_p(4, 4, x), &
        105.0_dp*(1.0_dp - x*x)**2)

    call check("negative-order normalization", legendre_p(3, -2, x), &
        legendre_p(3, 2, x)/120.0_dp)
    call check("parity", legendre_p(5, 2, -x), -legendre_p(5, 2, x))
    call check("degree below order", legendre_p(2, 3, x), 0.0_dp)
    call check("P_l(1)", legendre_p(8, 0, 1.0_dp), 1.0_dp)
    call check("P_l(-1)", legendre_p(7, 0, -1.0_dp), -1.0_dp)
    call check("P_l^m endpoint", legendre_p(7, 3, 1.0_dp), 0.0_dp)
    call check_true("outside Ferrers domain", ieee_is_nan(legendre_p(2, 0, 1.1_dp)))

    call check("dP_1", legendre_p_derivative(1, 0, x), 1.0_dp)
    call check("dP_2", legendre_p_derivative(2, 0, x), 3.0_dp*x)
    call check("dP_3^2", legendre_p_derivative(3, 2, x), &
        15.0_dp*(1.0_dp - 3.0_dp*x*x))

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        real(dp) :: error

        error = abs(got - expected)
        if (.not. error <= tol*(1.0_dp + abs(expected))) then
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

end program test_fortnum_special_legendre
