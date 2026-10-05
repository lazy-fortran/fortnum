module fortnum_reduction_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_round_type, ieee_get_rounding_mode, &
        ieee_set_rounding_mode, ieee_down, ieee_support_rounding, ieee_set_flag, ieee_all
    use fortnum_certified_reductions, only: reduction_t, reduction_ok, &
        reduction_input, reduction_range, reduction_environment_supported, &
        aggregate_sum, aggregate_dot, interval_sum, interval_dot
    implicit none
    private
    public :: validate_reductions, make_reduction_workload, require_encloses

contains

    subroutine require_encloses(result, exact)
        type(reduction_t), intent(in) :: result
        real(qp), intent(in) :: exact

        if (result%status /= reduction_ok) error stop "unexpected reduction failure"
        if (result%radius < 0.0_dp) error stop "negative reduction radius"
        if (real(result%center, qp) - real(result%radius, qp) > exact) &
            error stop "reduction lower bound misses exact dyadic oracle"
        if (real(result%center, qp) + real(result%radius, qp) < exact) &
            error stop "reduction upper bound misses exact dyadic oracle"
    end subroutine require_encloses

    subroutine check_sum(x, exact)
        real(dp), intent(in) :: x(:)
        real(qp), intent(in) :: exact

        call require_encloses(aggregate_sum(x), exact)
        call require_encloses(interval_sum(x), exact)
    end subroutine check_sum

    subroutine check_dot(x, y, exact)
        real(dp), intent(in) :: x(:), y(:)
        real(qp), intent(in) :: exact

        call require_encloses(aggregate_dot(x, y), exact)
        call require_encloses(interval_dot(x, y), exact)
    end subroutine check_dot

    ! The integer accumulator is exact for the bounded benchmark sizes.
    ! No tested floating reduction is used to construct its own oracle.
    subroutine make_reduction_workload(x, y, cancel, exact_sum, exact_dot)
        real(dp), intent(out) :: x(:), y(:)
        logical, intent(in) :: cancel
        real(qp), intent(out) :: exact_sum, exact_dot
        integer(int64) :: sum_integer, dot_integer, a, b
        integer :: i

        if (size(x) /= size(y)) error stop "workload shape mismatch"
        if (size(x) > 1048576) error stop "workload exceeds integer oracle domain"
        sum_integer = 0_int64
        dot_integer = 0_int64
        do i = 1, size(x)
            a = int(1 + mod(i, 17), int64)
            if (cancel .and. mod(i, 2) == 0) a = -a
            b = int(1 + mod(i, 7), int64)
            x(i) = real(a, dp)/1024.0_dp
            y(i) = real(b, dp)/128.0_dp
            sum_integer = sum_integer + a
            dot_integer = dot_integer + a*b
        end do
        exact_sum = real(sum_integer, qp)/1024.0_qp
        exact_dot = real(dot_integer, qp)/131072.0_qp
    end subroutine make_reduction_workload

    subroutine validate_reductions()
        real(dp) :: x(4096), y(4096), eta, bad(2)
        real(qp) :: exact_sum, exact_dot
        type(reduction_t) :: result
        type(ieee_round_type) :: mode
        integer :: i

        if (.not. reduction_environment_supported()) &
            error stop "unsupported FP environment for certified reductions"
        call check_sum([real(dp) ::], 0.0_qp)
        call check_dot([real(dp) ::], [real(dp) ::], 0.0_qp)
        call check_sum([0.5_dp, -0.25_dp, 0.125_dp], 0.375_qp)
        call check_dot([0.5_dp, -0.25_dp], [0.5_dp, 0.5_dp], 0.125_qp)
        call check_sum([2.0_dp**53, 1.0_dp, -2.0_dp**53], 1.0_qp)
        call check_dot([2.0_dp**53, 1.0_dp, -2.0_dp**53], &
            [1.0_dp, 1.0_dp, 1.0_dp], 1.0_qp)
        call check_sum([2.0_dp**53, -2.0_dp**53, 1.0_dp], 1.0_qp)
        call check_dot([1.0_dp + 2.0_dp**(-52)], &
            [1.0_dp - 2.0_dp**(-52)], 1.0_qp - 2.0_qp**(-104))
        do i = 1, 2
            call make_reduction_workload(x, y, i == 2, exact_sum, exact_dot)
            call check_sum(x, exact_sum)
            call check_dot(x, y, exact_dot)
        end do
        eta = tiny(1.0_dp)*epsilon(1.0_dp)
        call check_sum([eta, eta, -eta], real(eta, qp))
        call check_dot([eta, eta], [0.5_dp, 0.5_dp], real(eta, qp))
        call check_dot([eta], [0.5_dp], real(eta, qp)/2.0_qp)
        call check_dot([tiny(1.0_dp)], [2.0_dp**(-78)], 2.0_qp**(-1100))

        bad = [huge(1.0_dp), huge(1.0_dp)]
        result = aggregate_sum(bad)
        if (result%status /= reduction_range) error stop "sum overflow not rejected"
        result = aggregate_sum([huge(1.0_dp), -huge(1.0_dp)])
        if (result%status /= reduction_range) error stop "magnitude overflow not rejected"
        result = aggregate_dot(bad, [2.0_dp, 2.0_dp])
        if (result%status /= reduction_range) error stop "product overflow not rejected"
        result = interval_sum(bad)
        if (result%status /= reduction_range) error stop "interval overflow not reported"
        result = aggregate_dot([1.0_dp], [1.0_dp, 2.0_dp])
        if (result%status /= reduction_input) error stop "shape mismatch not rejected"
        bad = [ieee_value(1.0_dp, ieee_quiet_nan), 1.0_dp]
        result = aggregate_sum(bad)
        if (result%status /= reduction_input) error stop "NaN not rejected"
        bad(1) = ieee_value(1.0_dp, ieee_positive_inf)
        result = aggregate_dot(bad, [1.0_dp, 1.0_dp])
        if (result%status /= reduction_input) error stop "infinity not rejected"
        if (ieee_support_rounding(ieee_down)) then
            call ieee_get_rounding_mode(mode)
            call ieee_set_rounding_mode(ieee_down)
            if (reduction_environment_supported()) error stop "wrong rounding accepted"
            call ieee_set_rounding_mode(mode)
        end if
        ! Overflow/underflow are expected fixture events, not benchmark failures.
        call ieee_set_flag(ieee_all, .false.)
    end subroutine validate_reductions

end module fortnum_reduction_fixtures
