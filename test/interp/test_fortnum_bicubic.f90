program test_fortnum_bicubic
    use fortnum_kinds, only: dp
    use fortnum_bicubic, only: bicubic_jet, bicubic_jet_jvp, bicubic_jet_vjp
    use fortnum_interval, only: interval_t, interval, mag, &
                                operator(+), operator(-), operator(*), operator(**)
    use fortnum_rounding, only: round_up, round_down
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), parameter :: points(2, 5) = reshape( &
                           [0.0_dp, 0.0_dp, 0.125_dp, -0.3_dp, -0.8_dp, &
                            0.4_dp, 1.0_dp, 1.0_dp, 2.0_dp, -1.25_dp], [2, 5])
    real(dp) :: c(4, 4), dc(4, 4), bc(4, 4), expected_bc(4, 4)
    real(dp) :: jet(3), expected(3), hess(3), djet(3), expected_djet(3)
    real(dp) :: cj(3), ch(3), weights(3), bx, by, expected_bx, expected_by
    real(dp) :: x, y, dx, dy, plus(3), minus(3), step, scale
    integer :: i, j, k, m

    do k = 1, size(points, 2)
        do j = 1, 4
            do i = 1, 4
                c = 0.0_dp
                c(i, j) = 1.0_dp
                call check_case(c, points(1, k), points(2, k))
            end do
        end do
    end do
    do j = 1, 4
        do i = 1, 4
            c(i, j) = real((-1)**(i + j)*(i + 3*j), dp)/17.0_dp
            dc(i, j) = real(2*i - j, dp)/13.0_dp
        end do
    end do
    do k = 1, size(points, 2)
        call check_case(c, points(1, k), points(2, k))
        call check_case(1.0e100_dp*c, points(1, k), points(2, k))
        call check_case(1.0e-100_dp*c, points(1, k), points(2, k))
    end do
    ! Exact cancellation is judged by an absolute coefficient-weight budget.
    c = 0.0_dp
    c(:, 1) = [-1.0_dp, 3.0_dp, -3.0_dp, 1.0_dp]
    call check_case(c, 1.0_dp, 0.25_dp)

    do j = 1, 4
        do i = 1, 4
            c(i, j) = real((-1)**(i + j)*(i + 3*j), dp)/17.0_dp
        end do
    end do
    x = -0.37_dp
    y = 0.62_dp
    weights = [0.7_dp, -0.4_dp, 1.3_dp]
    dx = 0.3_dp
    dy = -0.2_dp
    call reference(c, x, y, expected, hess)
    call reference(dc, x, y, cj, ch)
    expected_djet = cj + [expected(2)*dx + expected(3)*dy, &
                          hess(1)*dx + hess(2)*dy, hess(2)*dx + hess(3)*dy]
    call bicubic_jet_jvp(c, x, y, dc, dx, dy, jet, djet)
    call close_vector(jet, expected, 2.0e-13_dp, 'JVP primal')
    call close_vector(djet, expected_djet, 2.0e-13_dp, 'mixed active JVP')
    step = 1.0e-5_dp
    call bicubic_jet(c + step*dc, x + step*dx, y + step*dy, plus)
    call bicubic_jet(c - step*dc, x - step*dx, y - step*dy, minus)
    call close_vector(djet, (plus - minus)/(2*step), 2.0e-9_dp, 'central JVP')
    call bicubic_jet_vjp(c, x, y, weights, jet, bc, bx, by)
    call close_vector(jet, expected, 2.0e-13_dp, 'VJP primal')
    expected_bx = dot_product(weights, [expected(2), hess(1), hess(2)])
    expected_by = dot_product(weights, [expected(3), hess(2), hess(3)])
    scale = sum(dc*bc) + dx*bx + dy*by
    call close_vector([dot_product(weights, djet)], [scale], &
                      2.0e-13_dp, 'mixed JVP/VJP pairing')
    do j = 1, 4
        do i = 1, 4
            dc = 0.0_dp
            dc(i, j) = 1.0_dp
            call reference(dc, x, y, cj, ch)
            expected_bc(i, j) = dot_product(weights, cj)
            call bicubic_jet_jvp(c, x, y, dc, 0.0_dp, 0.0_dp, jet, djet)
            call close_vector(djet, cj, 2.0e-13_dp, 'coefficient JVP')
        end do
    end do
    call close_vector(reshape(bc, [16]), reshape(expected_bc, [16]), &
                      2.0e-13_dp, 'coefficient VJP')
    call close_vector([bx, by], [expected_bx, expected_by], &
                      2.0e-13_dp, 'coordinate VJP')
    ! Separate pure-coordinate directions catch omitted second partials.
    dc = 0.0_dp
    do m = 1, 2
        dx = real(2 - m, dp)
        dy = real(m - 1, dp)
        call bicubic_jet_jvp(c, x, y, dc, dx, dy, jet, djet)
        expected_djet = [expected(2)*dx + expected(3)*dy, &
                         hess(1)*dx + hess(2)*dy, hess(2)*dx + hess(3)*dy]
        call close_vector(djet, expected_djet, 2.0e-13_dp, 'coordinate JVP')
        scale = sum(dc*bc) + dx*bx + dy*by
        call close_vector([dot_product(weights, djet)], [scale], &
                          2.0e-13_dp, 'JVP/VJP Euclidean pairing')
    end do
    print '(a)', 'PASS bicubic basis, independent partials, JVP/VJP and enclosures'
contains
    subroutine reference(coeff, a, b, value, second)
        ! Independent direct monomial sums, never the production Horner schedule.
        real(dp), intent(in) :: coeff(4, 4), a, b
        real(dp), intent(out) :: value(3), second(3)
        integer :: r, s

        value = 0.0_dp
        second = 0.0_dp
        do s = 0, 3
            do r = 0, 3
                value(1) = value(1) + coeff(r + 1, s + 1)*a**r*b**s
                if (r > 0) then
                    value(2) = value(2) &
                               + real(r, dp)*coeff(r + 1, s + 1)*a**(r - 1)*b**s
                end if
                if (s > 0) then
                    value(3) = value(3) &
                               + real(s, dp)*coeff(r + 1, s + 1)*a**r*b**(s - 1)
                end if
                if (r > 1) then
                    second(1) = second(1) &
                               + real(r*(r - 1), dp)*coeff(r + 1, s + 1)*a**(r - 2)*b**s
                end if
                if (r > 0 .and. s > 0) then
                    second(2) = second(2) &
                               + real(r*s, dp)*coeff(r + 1, s + 1)*a**(r - 1)*b**(s - 1)
                end if
                if (s > 1) then
                    second(3) = second(3) &
                               + real(s*(s - 1), dp)*coeff(r + 1, s + 1)*a**r*b**(s - 2)
                end if
            end do
        end do
    end subroutine reference

    subroutine check_case(coeff, a, b)
        real(dp), intent(in) :: coeff(4, 4), a, b
        type(interval_t) :: exact(3), size_sum(3), term, error_bound
        real(dp) :: value(3), ref(3), second(3), gamma, radius
        integer :: r, s, q

        call bicubic_jet(coeff, a, b, value)
        call reference(coeff, a, b, ref, second)
        exact = interval(0.0_dp)
        size_sum = interval(0.0_dp)
        ! Strict outward arithmetic on the independent monomial representation.
        do s = 0, 3
            do r = 0, 3
                term = interval(coeff(r + 1, s + 1))*interval(a)**r*interval(b)**s
                exact(1) = exact(1) + term
                size_sum(1) = size_sum(1) + interval(mag(term))
                if (r > 0) then
                    term = interval(real(r, dp))*interval(coeff(r + 1, s + 1)) &
                           *interval(a)**(r - 1)*interval(b)**s
                    exact(2) = exact(2) + term
                    size_sum(2) = size_sum(2) + interval(mag(term))
                end if
                if (s > 0) then
                    term = interval(real(s, dp))*interval(coeff(r + 1, s + 1)) &
                           *interval(a)**r*interval(b)**(s - 1)
                    exact(3) = exact(3) + term
                    size_sum(3) = size_sum(3) + interval(mag(term))
                end if
            end do
        end do
        ! 128 bounds all 71 emitted primal operations; see docs/bicubic.md.
        gamma = round_up((128*epsilon(1.0_dp)/2)/(1 - 128*epsilon(1.0_dp)/2))
        do q = 1, 3
            if (.not. ieee_is_finite(value(q))) error stop 'nonfinite finite case'
            error_bound = interval(gamma)*size_sum(q)
            radius = round_up(error_bound%hi + 512*nearest(0.0_dp, 1.0_dp))
            if (abs(value(q) - ref(q)) > 2*radius) &
                error stop 'direct monomial disagreement'
            ! The rigorous oracle's entire enclosure fits the center budget.
            if (round_down(value(q) - radius) > exact(q)%lo) &
                error stop 'lower enclosure failed'
            if (round_up(value(q) + radius) < exact(q)%hi) &
                error stop 'upper enclosure failed'
        end do
    end subroutine check_case

    subroutine close_vector(actual, wanted, tolerance, label)
        real(dp), intent(in) :: actual(:), wanted(:), tolerance
        character(*), intent(in) :: label
        if (any(.not. ieee_is_finite(actual))) error stop 'nonfinite product'
        if (any(abs(actual - wanted) > tolerance*max(1.0_dp, abs(wanted)))) then
            print '(a)', label
            print *, actual
            print *, wanted
            error stop 'bicubic derivative oracle'
        end if
    end subroutine close_vector
end program test_fortnum_bicubic
