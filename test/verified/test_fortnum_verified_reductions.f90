program test_fortnum_verified_reductions
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
                                                       ieee_positive_inf, ieee_is_finite
    use fortnum_verified_linalg, only: dot_enclosure, weighted_norm2_enclosure, norm2_up
    use fortnum_interval, only: interval_t
    implicit none
    real(dp) :: x(8), y(8), w(8), bound
    real(dp) :: many(16384), weights(16384), other(16384)
    real(qp) :: expected, energy
    type(interval_t) :: actual
    integer :: i, k
    x = [1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 0.5_dp, 0.25_dp, -0.125_dp, 0.0_dp]
    y = [4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, -0.5_dp, 0.125_dp, 0.25_dp, 0.0_dp]
    w = [0.0_dp, 0.125_dp, 0.25_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    do k = 1, 4
        expected = sum(real(x, qp)*real(y, qp))
        actual = dot_enclosure(x, y)
        call contains('dot dyadic/scale/cancellation', actual, expected)
        energy = sum(real(w, qp)*real(x, qp)**2)
        actual = weighted_norm2_enclosure(w, x)
        call contains('weighted energy', actual, energy)
        if (actual%lo < 0) error stop 'negative valid energy lower endpoint'
        bound = norm2_up(x)
        expected = sqrt(sum(real(x, qp)**2))
        if (real(bound, qp) < expected) error stop 'norm upper bound'
        select case (k)
        case (1)
            x = [2.0_dp**500, 1.0_dp, -2.0_dp**500, 2.0_dp**(-500), &
                 -2.0_dp**(-500), 0.0_dp, 0.0_dp, 0.0_dp]
            y = 1
            actual = dot_enclosure(x, y)
            call contains('exact cancellation leaves one', actual, 1.0_qp)
        case (2)
            x = nearest(0.0_dp, 1.0_dp)
            y = 0.5_dp
        case (3)
            x = huge(1.0_dp)/2
            y = 4
        end select
    end do
    actual = dot_enclosure(x(:0), y(:0))
    if (actual%lo /= 0 .or. actual%hi /= 0) error stop 'empty dot'
    actual = weighted_norm2_enclosure(w(:0), x(:0))
    if (actual%lo /= 0 .or. actual%hi /= 0) error stop 'empty energy'
    if (norm2_up(x(:0)) /= 0) error stop 'empty norm'
    actual = dot_enclosure(x(:2), y(:1))
    call uninformative(actual)
    w(1) = -1
    actual = weighted_norm2_enclosure(w, x)
    call uninformative(actual)
    w = 1
    x = 1
    x(2) = ieee_value(0.0_dp, ieee_quiet_nan)
    actual = dot_enclosure(x, y)
    call uninformative(actual)
    actual = weighted_norm2_enclosure(w, x)
    call uninformative(actual)
    if (ieee_is_finite(norm2_up(x))) error stop 'nonfinite norm input'
    x(2) = ieee_value(0.0_dp, ieee_positive_inf)
    actual = dot_enclosure(x, y)
    call uninformative(actual)
    ! Exercise the OpenMP reduction tree; all products and the total are exact
    ! dyadic values in the independent binary128 oracle.
    do i = 1, size(many)
        many(i) = real(mod(i, 17) - 8, dp)/16
        other(i) = real(mod(i, 13) - 6, dp)/8
        weights(i) = real(mod(i, 5), dp)/4
    end do
    expected = sum(real(many, qp)*real(other, qp))
    actual = dot_enclosure(many, other)
    call contains('parallel dot', actual, expected)
    expected = sum(real(weights, qp)*real(many, qp)**2)
    actual = weighted_norm2_enclosure(weights, many)
    call contains('parallel energy', actual, expected)
    expected = sqrt(sum(real(many, qp)**2))
    if (real(norm2_up(many), qp) < expected) error stop 'parallel norm'
    print *, 'verified allocation-free vector reductions PASS'
contains
    subroutine contains(label, iv, value)
        character(*), intent(in) :: label
        type(interval_t), intent(in) :: iv
        real(qp), intent(in) :: value
        if (real(iv%lo, qp) <= value .and. value <= real(iv%hi, qp)) return
        print *, label, iv, value
        error stop 'enclosure excludes independent oracle'
    end subroutine
    subroutine uninformative(iv)
        type(interval_t), intent(in) :: iv
        if (ieee_is_finite(iv%lo) .or. ieee_is_finite(iv%hi)) &
            error stop 'invalid input returned plausible finite bound'
        if (.not. (iv%lo < 0 .and. iv%hi > 0)) error stop 'invalid interval'
    end subroutine
end program
