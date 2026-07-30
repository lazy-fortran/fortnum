program test_fortnum_special_jacobi
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_special_jacobi, only: jacobi_p, jacobi_p_derivative, &
        scaled_jacobi_p, tetrahedron_koornwinder, triangle_dubiner
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
    call check("scaled definition", &
        scaled_jacobi_p(7, 2.0_dp, 1.0_dp, 0.31_dp, 0.73_dp), &
        0.73_dp**7*jacobi_p(7, 2.0_dp, 1.0_dp, 0.31_dp/0.73_dp))
    call check("scaled Legendre polynomial", &
        scaled_jacobi_p(3, 0.0_dp, 0.0_dp, 0.31_dp, 0.73_dp), &
        (5.0_dp*0.31_dp**3 - 3.0_dp*0.31_dp*0.73_dp**2)/2.0_dp)
    call check("scaled removable limit", &
        scaled_jacobi_p(3, 0.0_dp, 0.0_dp, 0.31_dp, 0.0_dp), &
        2.5_dp*0.31_dp**3)
    call check("scaled zero degree", &
        scaled_jacobi_p(0, 2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp), 1.0_dp)
    call check("triangle constant", &
        triangle_dubiner(0, 0, 0.2_dp, 0.3_dp), 1.0_dp)
    call check("triangle first x mode", &
        triangle_dubiner(1, 0, 0.2_dp, 0.3_dp), -0.3_dp)
    call check("triangle first y mode", &
        triangle_dubiner(0, 1, 0.2_dp, 0.3_dp), -0.1_dp)
    call check("tetrahedron constant", &
        tetrahedron_koornwinder(0, 0, 0, 0.1_dp, 0.2_dp, 0.3_dp), &
        1.0_dp)
    call check("tetrahedron first x mode", &
        tetrahedron_koornwinder(1, 0, 0, 0.1_dp, 0.2_dp, 0.3_dp), &
        -0.3_dp)
    call check("tetrahedron first y mode", &
        tetrahedron_koornwinder(0, 1, 0, 0.1_dp, 0.2_dp, 0.3_dp), &
        -0.1_dp)
    call check("tetrahedron first z mode", &
        tetrahedron_koornwinder(0, 0, 1, 0.1_dp, 0.2_dp, 0.3_dp), &
        0.2_dp)
    call check_true("negative degree rejected", &
        jacobi_p(-1, 0.0_dp, 0.0_dp, x) == 0.0_dp)
    call check_weighted_orthogonality()
    call check_simplex_orthogonality()

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

    subroutine check_simplex_orthogonality()
        integer, parameter :: order = 10
        integer, parameter :: triangle_count = 10, tetrahedron_count = 20
        real(dp) :: nodes(order), weights(order)
        real(dp) :: triangle_gram(triangle_count, triangle_count)
        real(dp) :: tetrahedron_gram(tetrahedron_count, tetrahedron_count)
        real(dp) :: triangle_values(triangle_count)
        real(dp) :: tetrahedron_values(tetrahedron_count)
        real(dp) :: r, s, t, weight, xq, yq, zq
        integer :: first, i, j, k, left, mode, right, second, third

        call gauss_legendre(order, nodes, weights)
        triangle_gram = 0.0_dp
        tetrahedron_gram = 0.0_dp
        do first = 1, order
            r = 0.5_dp*(nodes(first) + 1.0_dp)
            do second = 1, order
                s = 0.5_dp*(nodes(second) + 1.0_dp)
                xq = r
                yq = (1.0_dp - r)*s
                weight = 0.25_dp*weights(first)*weights(second)* &
                    (1.0_dp - r)
                mode = 0
                do i = 0, 3
                    do j = 0, 3 - i
                        mode = mode + 1
                        triangle_values(mode) = triangle_dubiner( &
                            i, j, xq, yq)
                    end do
                end do
                do left = 1, triangle_count
                    do right = 1, triangle_count
                        triangle_gram(left, right) = &
                            triangle_gram(left, right) + weight* &
                            triangle_values(left)*triangle_values(right)
                    end do
                end do
                do third = 1, order
                    t = 0.5_dp*(nodes(third) + 1.0_dp)
                    zq = (1.0_dp - r)*(1.0_dp - s)*t
                    weight = 0.125_dp*weights(first)*weights(second)* &
                        weights(third)*(1.0_dp - r)**2*(1.0_dp - s)
                    mode = 0
                    do i = 0, 3
                        do j = 0, 3 - i
                            do k = 0, 3 - i - j
                                mode = mode + 1
                                tetrahedron_values(mode) = &
                                    tetrahedron_koornwinder( &
                                    i, j, k, xq, yq, zq)
                            end do
                        end do
                    end do
                    do left = 1, tetrahedron_count
                        do right = 1, tetrahedron_count
                            tetrahedron_gram(left, right) = &
                                tetrahedron_gram(left, right) + weight* &
                                tetrahedron_values(left)* &
                                tetrahedron_values(right)
                        end do
                    end do
                end do
            end do
        end do
        call check_orthogonal_gram("triangle modal orthogonality", &
            triangle_gram)
        call check_orthogonal_gram("tetrahedron modal orthogonality", &
            tetrahedron_gram)
    end subroutine check_simplex_orthogonality

    subroutine check_orthogonal_gram(label, gram)
        character(*), intent(in) :: label
        real(dp), intent(in) :: gram(:, :)
        real(dp) :: correlation
        integer :: left, right

        do left = 1, size(gram, 1)
            call check_true(label//" positive norm", gram(left, left) > 0.0_dp)
            do right = 1, left - 1
                correlation = abs(gram(left, right))/ &
                    sqrt(gram(left, left)*gram(right, right))
                call check_true(label, correlation < 2.0e-12_dp)
            end do
        end do
    end subroutine check_orthogonal_gram

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
