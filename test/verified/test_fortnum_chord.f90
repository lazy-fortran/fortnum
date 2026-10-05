program test_fortnum_chord
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
        ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, is_empty, contains, &
        linear_interpolation_remainder
    implicit none
    type(interval_t) :: left, right, curvature, r, bad, endpoints(2), bounds(2)
    real(dp) :: eta, lo, hi, m
    real(qp) :: midpoint, chord, actual, exact
    integer :: j

    ! g(x)=x^2 on [-1,1], M=2: the midpoint error is exactly one.
    left = interval(-1.0_dp)
    right = interval(1.0_dp)
    curvature = interval(2.0_dp)
    r = linear_interpolation_remainder(left, right, curvature)
    midpoint = (real(left%lo, qp) + real(right%hi, qp))/2.0_qp
    chord = (real(left%lo, qp)**2 + real(right%hi, qp)**2)/2.0_qp
    actual = abs(midpoint**2 - chord)
    if (actual /= 1.0_qp) error stop 'independent quadratic fixture'
    call encloses(r, actual)
    if (r%hi > 1.00000000000001_dp) error stop 'quadratic bound too broad'
    ! The same sharp error occurs for -x^2; an affine function has zero error.
    actual = abs(-midpoint**2 + chord)
    call encloses(r, actual)
    r = linear_interpolation_remainder(left, right, interval(0.0_dp))
    call exact_zero(r)
    r = linear_interpolation_remainder(left, left, curvature)
    call exact_zero(r)

    ! Small-significand dyadic fixtures: the independent binary128 oracle is
    ! exact for every subtraction/square/product (well under 113 bits).
    do j = 1, 127
        lo = real(j, dp)/16.0_dp
        hi = lo + real(mod(j, 7) + 1, dp)/32.0_dp
        m = real(j + 1, dp)/8.0_dp
        exact = real(m, qp)*(real(hi, qp) - real(lo, qp))**2/8.0_qp
        r = linear_interpolation_remainder(interval(lo), interval(hi), interval(m))
        call encloses(r, exact)
    end do
    ! Promotion preserves the exact binary64 input; division by 32 is exact.
    r = linear_interpolation_remainder(interval(0.0_dp), interval(0.5_dp), &
                                       interval(0.1_dp))
    call encloses(r, real(0.1_dp, qp)/32.0_qp)
    ! A 53-bit squared width and a four-bit curvature require at most 110
    ! significand bits; this oracle is exact despite inexact binary64 products.
    r = linear_interpolation_remainder(interval(0.0_dp), interval(0.1_dp), &
                                       interval(1.375_dp))
    call encloses(r, 1.375_qp*real(0.1_dp, qp)**2/8.0_qp)
    ! Endpoint/curvature uncertainty includes the exact product extrema.
    r = linear_interpolation_remainder(interval(0.0_dp, 0.25_dp), &
                                       interval(0.5_dp, 1.0_dp), &
                                       interval(1.0_dp, 2.0_dp))
    call encloses(r, 1.0_qp/128.0_qp)
    call encloses(r, 1.0_qp/4.0_qp)
    ! Gradual underflow: eta^2/4 is positive and exact in binary128.
    eta = tiny(1.0_dp)*epsilon(1.0_dp)
    r = linear_interpolation_remainder(interval(0.0_dp), interval(eta), curvature)
    exact = real(eta, qp)**2/4.0_qp
    if (exact <= 0.0_qp) error stop 'underflow oracle unavailable'
    call encloses(r, exact)
    if (r%lo < 0.0_dp) error stop 'negative remainder radius'
    ! Elemental use with independent exact answers.
    endpoints(1) = interval(1.0_dp)
    endpoints(2) = interval(2.0_dp)
    bounds = linear_interpolation_remainder(interval(0.0_dp), endpoints, curvature)
    if (.not. contains(bounds(1), 0.25_dp)) error stop 'elemental first'
    if (.not. contains(bounds(2), 1.0_dp)) error stop 'elemental second'

    bad%lo = 1.0_dp
    bad%hi = 0.0_dp
    call rejects(bad, right, curvature)
    call rejects(left, bad, curvature)
    call rejects(left, right, bad)
    call rejects(right, left, curvature)
    call rejects(interval(0.0_dp, 0.75_dp), interval(0.5_dp, 1.0_dp), curvature)
    call rejects(left, right, interval(-1.0_dp))
    call rejects(left, right, interval(-1.0_dp, 1.0_dp))
    bad = interval(ieee_value(1.0_dp, ieee_quiet_nan))
    call rejects(bad, right, curvature)
    call rejects(left, bad, curvature)
    call rejects(left, right, bad)
    bad = interval(ieee_value(1.0_dp, ieee_positive_inf))
    call rejects(bad, right, curvature)
    call rejects(left, bad, curvature)
    call rejects(left, right, bad)
    call rejects(interval(-huge(1.0_dp)), interval(huge(1.0_dp)), curvature)
    call rejects(left, right, interval(huge(1.0_dp)))
    print '(a)', 'C2 chord remainder independent fixtures PASS'

contains

    subroutine encloses(value, reference)
        type(interval_t), intent(in) :: value
        real(qp), intent(in) :: reference

        if (is_empty(value)) error stop 'valid remainder rejected'
        if (real(value%lo, qp) > reference .or. &
            real(value%hi, qp) < reference) error stop 'exact remainder excluded'
    end subroutine encloses

    subroutine exact_zero(value)
        type(interval_t), intent(in) :: value

        if (is_empty(value)) error stop 'zero remainder rejected'
        if (value%lo /= 0.0_dp .or. value%hi /= 0.0_dp) &
            error stop 'exact zero remainder'
    end subroutine exact_zero

    subroutine rejects(a, b, bound)
        type(interval_t), intent(in) :: a, b, bound
        type(interval_t) :: value

        value = linear_interpolation_remainder(a, b, bound)
        if (.not. is_empty(value)) error stop 'invalid remainder accepted'
    end subroutine rejects

end program test_fortnum_chord
