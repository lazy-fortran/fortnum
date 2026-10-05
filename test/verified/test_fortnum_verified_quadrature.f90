module verified_midpoint_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use fortnum_interval, only: interval_t, interval, sqr, exp, mag, &
        operator(+), operator(-), operator(*)
    implicit none
    integer :: fixture = 1, evaluations = 0
contains
    subroutine evaluate(cell, midpoint, value, curvature, ok)
        type(interval_t), intent(in) :: cell, midpoint
        type(interval_t), intent(out) :: value, curvature
        logical, intent(out) :: ok
        type(interval_t) :: derivative_range

        evaluations = evaluations + 1
        ok = .true.
        select case (fixture)
        case (1)
            value = sqr(midpoint) + 2.0_dp*midpoint + 1.0_dp
            curvature = interval(2.0_dp)
        case (2)
            value = sqr(sqr(midpoint))
            curvature = interval(0.0_dp, mag(12.0_dp*sqr(cell)))
        case (3)
            value = exp(midpoint)
            derivative_range = exp(cell)
            curvature = interval(0.0_dp, derivative_range%hi)
        case (4)
            value = interval(1.0_dp)
            curvature = interval(0.0_dp)
        case (5)
            value = 1.0e12_dp*(midpoint - 0.5_dp)
            curvature = interval(0.0_dp)
        case (6)
            value%lo = ieee_value(0.0_dp, ieee_quiet_nan)
            value%hi = value%lo
            curvature = interval(1.0_dp)
        case (7)
            value = interval(1.0_dp)
            curvature = interval(-1.0_dp, 1.0_dp)
        case (8)
            value = interval(1.0_dp)
            curvature = interval(1.0_dp)
            ok = .false.
        case (9)
            value = interval(huge(1.0_dp))
            curvature = interval(0.0_dp)
        case (10)
            value = interval(1.0_dp)
            curvature = interval(ieee_value(0.0_dp, ieee_positive_inf))
        end select
    end subroutine evaluate

    subroutine evaluate_batch(cell, midpoint, value, curvature, ok)
        type(interval_t), intent(in) :: cell, midpoint
        type(interval_t), intent(out) :: value(:), curvature(:)
        logical, intent(out) :: ok
        type(interval_t) :: common
        integer :: component

        evaluations = evaluations + 1
        common = sqr(midpoint) + 2.0_dp*midpoint + 1.0_dp
        do component = 1, size(value)
            value(component) = real(component, dp)*common
            curvature(component) = interval(2.0_dp*real(component, dp))
        end do
        ok = .true.
        if (fixture == 7) curvature(size(curvature)) = interval(-1.0_dp)
        if (fixture == 8) ok = .false.
    end subroutine evaluate_batch
end module verified_midpoint_fixture

program test_fortnum_verified_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use fortnum_interval, only: interval_t, rad, is_empty
    use fortnum_verified_quadrature, only: composite_midpoint, composite_midpoint_batch
    use verified_midpoint_fixture, only: evaluate, evaluate_batch, fixture, evaluations
    implicit none
    type(interval_t) :: result, coarse
    type(interval_t) :: batch(8)
    real(dp) :: edges(6), grid(65), invalid(2)
    real(qp) :: expected
    logical :: ok
    integer :: j

    edges = [-1.0_dp, -0.98_dp, -0.25_dp, 0.0_dp, 0.2_dp, 2.0_dp]
    fixture = 1
    call composite_midpoint(evaluate, edges, result, ok)
    call check(ok, 'nonuniform quadratic')
    ! Integral of (x+1)^2 from -1 to2 is9. Positive curvature makes
    ! the midpoint sum an underestimate; a full integral enclosure must include9.
    call enclosed(result, 9.0_qp)
    if (result%hi > 9.0_dp + 1.0e-12_dp) error stop 'quadratic remainder too broad'
    evaluations = 0
    call composite_midpoint_batch(evaluate_batch, edges, batch, ok)
    call check(ok, 'batched nonuniform quadratic')
    if (evaluations /= 5) error stop 'batch repeated shared callback work'
    do j = 1, 8
        call enclosed(batch(j), 9.0_qp*real(j, qp))
    end do
    fixture = 2
    call composite_midpoint(evaluate, edges, result, ok)
    call check(ok, 'cellwise quartic curvature')
    call enclosed(result, 33.0_qp/5.0_qp)
    fixture = 3
    edges = [-0.5_dp, -0.49_dp, -0.125_dp, 0.0_dp, 0.25_dp, 1.125_dp]
    call composite_midpoint(evaluate, edges, result, ok)
    call check(ok, 'exponential cell bounds')
    expected = exp(1.125_qp) - exp(-0.5_qp)
    call enclosed(result, expected)

    ! Exact stored edges have mathematical midpoints that may not be binary64.
    ! Zero-curvature integral checks expose midpoint/weight/sum roundoff loss.
    fixture = 4
    edges = [1.0_dp, nearest(1.0_dp, 1.0_dp), &
             nearest(nearest(1.0_dp, 1.0_dp), 1.0_dp), 1.1_dp, 1.25_dp, 1.5_dp]
    call composite_midpoint(evaluate, edges, result, ok)
    call check(ok, 'rounded midpoint geometry')
    call enclosed(result, 0.5_qp)
    if (rad(result) > 1.0e-13_dp) error stop 'rounding width needlessly broad'
    fixture = 5
    edges = [0.0_dp, 0.1_dp, 0.25_dp, 0.5_dp, 0.9_dp, 1.0_dp]
    call composite_midpoint(evaluate, edges, result, ok)
    call check(ok, 'signed cancellation')
    call enclosed(result, 0.0_qp)

    fixture = 1
    do j = 1, 65
        grid(j) = real(j - 1, dp)/64.0_dp
    end do
    call composite_midpoint(evaluate, grid(1:65:2), coarse, ok)
    call check(ok, 'coarse quadratic')
    call enclosed(coarse, 7.0_qp/3.0_qp)
    call composite_midpoint(evaluate, grid, result, ok)
    call check(ok, 'fine quadratic')
    call enclosed(result, 7.0_qp/3.0_qp)
    if (rad(result) >= 0.26_dp*rad(coarse)) error stop 'midpoint order lost'

    ! Preconditions must fail before invoking an integrand on an invalid domain.
    invalid = [1.0_dp, 0.0_dp]
    evaluations = 0
    call composite_midpoint(evaluate, invalid, result, ok)
    call rejected(ok, result)
    if (evaluations /= 0) error stop 'invalid domain reached callback'
    invalid = [0.0_dp, 0.0_dp]
    call composite_midpoint(evaluate, invalid, result, ok)
    call rejected(ok, result)
    call composite_midpoint(evaluate, invalid(1:1), result, ok)
    call rejected(ok, result)
    invalid(2) = ieee_value(0.0_dp, ieee_quiet_nan)
    call composite_midpoint(evaluate, invalid, result, ok)
    call rejected(ok, result)
    invalid = [0.0_dp, 2.0_dp]
    do fixture = 6, 10
        call composite_midpoint(evaluate, invalid, result, ok)
        call rejected(ok, result)
    end do
    fixture = 7
    call composite_midpoint_batch(evaluate_batch, invalid, batch, ok)
    do j = 1, 8
        call rejected(ok, batch(j))
    end do
    fixture = 8
    call composite_midpoint_batch(evaluate_batch, invalid, batch, ok)
    do j = 1, 8
        call rejected(ok, batch(j))
    end do
    call composite_midpoint_batch(evaluate_batch, invalid, batch(1:0), ok)
    if (ok) error stop 'zero-component batch accepted'
    print '(a)', 'Verified midpoint independent polynomial/rounding oracles PASS'

contains
    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) error stop label
    end subroutine check

    subroutine enclosed(packet, exact)
        type(interval_t), intent(in) :: packet
        real(qp), intent(in) :: exact
        if (is_empty(packet)) error stop 'empty oracle packet'
        if (.not. ieee_is_finite(packet%lo)) error stop 'nonfinite lower oracle'
        if (.not. ieee_is_finite(packet%hi)) error stop 'nonfinite upper oracle'
        if (real(packet%lo, qp) > exact .or. real(packet%hi, qp) < exact) &
            error stop 'binary128 integral excluded'
    end subroutine enclosed

    subroutine rejected(accepted, packet)
        logical, intent(in) :: accepted
        type(interval_t), intent(in) :: packet
        if (accepted) error stop 'invalid quadrature accepted'
        if (.not. is_empty(packet)) error stop 'failure returned a usable packet'
    end subroutine rejected
end program test_fortnum_verified_quadrature
