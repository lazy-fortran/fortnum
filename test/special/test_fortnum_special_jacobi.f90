program test_fortnum_special_jacobi
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_jacobi, only: jacobi_p, jacobi_p_derivative
    implicit none

    integer :: nfail
    real(dp), parameter :: tolerance = 3.0e-13_dp
    real(dp), parameter :: x = 0.37_dp

    nfail = 0
    call check("P0", jacobi_p(0, 2.0_dp, 1.0_dp, x), 1.0_dp)
    call check("P1", jacobi_p(1, 2.0_dp, 1.0_dp, x), &
        0.5_dp*(1.0_dp + 5.0_dp*x))
    call check("P4(1)", jacobi_p(4, 2.0_dp, 1.0_dp, 1.0_dp), 15.0_dp)
    call check("P5(-1)", jacobi_p(5, 2.0_dp, 1.0_dp, -1.0_dp), -6.0_dp)
    call check("Legendre specialization", &
        jacobi_p(6, 0.0_dp, 0.0_dp, x), &
        (231.0_dp*x**6 - 315.0_dp*x**4 + 105.0_dp*x*x - 5.0_dp)/16.0_dp)
    call check("derivative identity", &
        jacobi_p_derivative(6, 2.0_dp, 1.0_dp, x), &
        5.0_dp*jacobi_p(5, 3.0_dp, 2.0_dp, x))
    call check_true("negative degree rejected", &
        jacobi_p(-1, 0.0_dp, 0.0_dp, x) == 0.0_dp)
    call check_weighted_orthogonality()

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine check_weighted_orthogonality()
        integer, parameter :: polynomial_count = 9
        integer, parameter :: quadrature_order = 24
        real(dp) :: nodes(quadrature_order), weights(quadrature_order)
        real(dp) :: computed, expected
        integer :: first, left, right

        call gauss_legendre(quadrature_order, nodes, weights)
        do left = 0, polynomial_count - 1
            do right = 0, polynomial_count - 1
                computed = 0.0_dp
                do first = 1, quadrature_order
                    computed = computed + weights(first)* &
                        (1.0_dp - nodes(first))**2* &
                        jacobi_p(left, 2.0_dp, 0.0_dp, nodes(first))* &
                        jacobi_p(right, 2.0_dp, 0.0_dp, nodes(first))
                end do
                expected = 0.0_dp
                if (left == right) expected = 8.0_dp/real(2*left + 3, dp)
                call check("weighted orthogonality", computed, expected)
            end do
        end do
    end subroutine check_weighted_orthogonality

    subroutine check(label, got, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        real(dp) :: error

        error = abs(got - expected)
        if (.not. error <= tolerance*(1.0_dp + abs(expected))) then
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

end program test_fortnum_special_jacobi
