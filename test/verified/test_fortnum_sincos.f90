program test_fortnum_sincos
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
        ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, interval_pi, sin, cos, sincos, &
        contains, operator(+), operator(-), operator(*), operator(/)
    implicit none
    real(dp), parameter :: exact_s(0:3) = [0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp]
    real(dp), parameter :: exact_c(0:3) = [1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp]
    type(interval_t) :: a, s, c, separate_s, separate_c
    type(interval_t) :: arguments(9), sines(9), cosines(9)
    real(dp) :: lo, hi, eta
    real(qp) :: x
    integer :: i, j, quadrant

    ! Exact quadrant values do not depend on a transcendental reference library.
    do i = -16, 16
        a = interval_pi()*real(i, dp)/2.0_dp
        call sincos(a, s, c)
        quadrant = modulo(i, 4)
        if (.not. contains(s, exact_s(quadrant))) error stop 'quadrant sine excluded'
        if (.not. contains(c, exact_c(quadrant))) error stop 'quadrant cosine excluded'
        a = a + interval(-1.0e-12_dp, 1.0e-12_dp)
        call sincos(a, s, c)
        if (.not. contains(s, exact_s(quadrant))) error stop 'interior sine extremum'
        if (.not. contains(c, exact_c(quadrant))) error stop 'interior cosine extremum'
    end do

    ! Deterministic broad/narrow intervals in positive and negative quadrants.
    ! Binary128 samples are independent accuracy oracles, not the range proof.
    do i = -128, 128
        lo = real(i, dp)/8.0_dp
        hi = lo + real(modulo(i, 11) + 1, dp)/16.0_dp
        a = interval(lo, hi)
        call sincos(a, s, c)
        do j = 0, 16
            x = real(lo, qp) + (real(hi, qp) - real(lo, qp))*real(j, qp)/16.0_qp
            call encloses(s, sin(x))
            call encloses(c, cos(x))
        end do
        ! Output equality checks only the optimization's preservation property;
        ! the independent identities/samples above supply behavioral oracles.
        separate_s = sin(a)
        separate_c = cos(a)
        call same_interval(s, separate_s)
        call same_interval(c, separate_c)
    end do

    arguments(1) = interval(-3.0_dp, 4.0_dp)
    arguments(2) = interval(-1.0e16_dp, 1.0e16_dp)
    arguments(3) = interval(-1.0e12_dp)
    arguments(4) = interval(1.0e12_dp)
    arguments(5) = interval(-0.01_dp, 0.01_dp)
    arguments(6) = interval(0.0_dp)
    arguments(7) = interval(1.0e-300_dp)
    arguments(8) = interval(-1.0e-300_dp)
    arguments(9) = interval(-1000.0_dp, 1000.0_dp)
    call sincos(arguments, sines, cosines)
    do i = 1, size(arguments)
        call same_interval(sines(i), sin(arguments(i)))
        call same_interval(cosines(i), cos(arguments(i)))
        x = (real(arguments(i)%lo, qp) + real(arguments(i)%hi, qp))/2.0_qp
        call encloses(sines(i), sin(x))
        call encloses(cosines(i), cos(x))
    end do
    call whole_unit(sines(1), cosines(1))
    call whole_unit(sines(2), cosines(2))
    call whole_unit(sines(9), cosines(9))
    eta = tiny(1.0_dp)*epsilon(1.0_dp)
    call sincos(interval(eta), s, c)
    call encloses(s, sin(real(eta, qp)))
    call encloses(c, cos(real(eta, qp)))

    ! Nonfinite and malformed arguments use the explicit conservative fallback.
    a%lo = -ieee_value(1.0_dp, ieee_positive_inf)
    a%hi = ieee_value(1.0_dp, ieee_positive_inf)
    call sincos(a, s, c)
    call whole_unit(s, c)
    a%lo = ieee_value(1.0_dp, ieee_quiet_nan)
    a%hi = 1.0_dp
    call sincos(a, s, c)
    call whole_unit(s, c)
    a%lo = 1.0_dp
    a%hi = ieee_value(1.0_dp, ieee_quiet_nan)
    call sincos(a, s, c)
    call whole_unit(s, c)
    a%lo = 2.0_dp
    a%hi = 1.0_dp
    call sincos(a, s, c)
    call whole_unit(s, c)
    print '(a)', 'paired interval sincos independent fixtures PASS'

contains

    subroutine encloses(value, exact)
        type(interval_t), intent(in) :: value
        real(qp), intent(in) :: exact

        if (real(value%lo, qp) > exact .or. real(value%hi, qp) < exact) &
            error stop 'paired range misses independent binary128 sample'
    end subroutine encloses

    subroutine same_interval(left, right)
        type(interval_t), intent(in) :: left, right

        if (left%lo < right%lo .or. left%lo > right%lo .or. &
            left%hi < right%hi .or. left%hi > right%hi) &
            error stop 'paired result differs from separate range'
    end subroutine same_interval

    subroutine whole_unit(sine, cosine)
        type(interval_t), intent(in) :: sine, cosine

        if (sine%lo > -1.0_dp .or. sine%hi < 1.0_dp .or. &
            cosine%lo > -1.0_dp .or. cosine%hi < 1.0_dp) &
            error stop 'conservative unit range fallback failed'
    end subroutine whole_unit

end program test_fortnum_sincos
