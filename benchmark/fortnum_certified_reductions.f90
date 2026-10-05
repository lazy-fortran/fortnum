!> Experimental benchmark candidates; no production selection is implied.
module fortnum_certified_reductions
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_positive_inf, ieee_support_datatype, ieee_support_denormal, &
        ieee_get_rounding_mode, ieee_round_type, ieee_nearest, operator(==)
    use fortnum_rounding, only: gamma_up, add_up, sub_down, sub_up, mul_up, div_up
    use fortnum_interval, only: interval_t, interval, operator(+), operator(*)
    implicit none
    private
    public :: reduction_t, reduction_ok, reduction_input, reduction_range
    public :: reduction_environment_supported, aggregate_sum, aggregate_dot
    public :: interval_sum, interval_dot

    integer, parameter :: reduction_ok = 0, reduction_input = 1, reduction_range = 2
    real(dp), parameter :: eta = tiny(1.0_dp)*epsilon(1.0_dp)

    type :: reduction_t
        real(dp) :: center = 0.0_dp
        real(dp) :: radius = 0.0_dp
        integer :: status = reduction_ok
    end type reduction_t

contains

    ! IEEE metadata checks cannot detect all unsafe compiler transformations.
    ! Caller must also preserve operation order and gradual underflow.
    function reduction_environment_supported() result(ok)
        logical :: ok
        type(ieee_round_type) :: mode
        real(dp), volatile :: normal, smallest

        call ieee_get_rounding_mode(mode)
        normal = tiny(1.0_dp)
        smallest = eta
        ok = ieee_support_datatype(1.0_dp) .and. ieee_support_denormal(1.0_dp)
        ok = ok .and. radix(1.0_dp) == 2 .and. digits(1.0_dp) == 53
        ok = ok .and. mode == ieee_nearest
        ok = ok .and. normal/2.0_dp > 0.0_dp
        ok = ok .and. smallest + smallest > smallest
    end function reduction_environment_supported

    function aggregate_sum(x) result(result)
        real(dp), intent(in) :: x(:)
        type(reduction_t) :: result
        real(dp) :: magnitude
        integer :: i

        if (.not. all(ieee_is_finite(x))) then
            result = failed(reduction_input)
            return
        end if
        result%center = 0.0_dp
        magnitude = 0.0_dp
        do i = 1, size(x)
            result%center = result%center + x(i)
            magnitude = magnitude + abs(x(i))
        end do
        call finish(result, magnitude, size(x))
    end function aggregate_sum

    function aggregate_dot(x, y) result(result)
        real(dp), intent(in) :: x(:), y(:)
        type(reduction_t) :: result
        real(dp) :: magnitude, product
        integer :: i

        if (size(x) /= size(y)) then
            result = failed(reduction_input)
            return
        end if
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
            result = failed(reduction_input)
            return
        end if
        if (size(x) > shiftr(huge(i), 1)) then
            result = failed(reduction_range)
            return
        end if
        result%center = 0.0_dp
        magnitude = 0.0_dp
        do i = 1, size(x)
            product = x(i)*y(i)
            result%center = result%center + product
            magnitude = magnitude + abs(product)
        end do
        call finish(result, magnitude, 2*size(x))
    end function aggregate_dot

    ! For k rounding operations and exact absolute sum A, both computed
    ! center and magnitude have error <= g*A+b, with g=gamma_k and
    ! b=k*eta/(1-k*u). Thus A <= (m+b)/(1-g). All budget arithmetic is outward.
    ! k=n for sums and k=2n conservatively covers products plus dot additions.
    subroutine finish(result, magnitude, k)
        type(reduction_t), intent(inout) :: result
        real(dp), intent(in) :: magnitude
        integer, intent(in) :: k
        real(dp) :: g, floor, absolute_sum

        if (k == 0) return
        if (.not. ieee_is_finite(result%center) .or. &
            .not. ieee_is_finite(magnitude)) then
            result = failed(reduction_range)
            return
        end if
        g = gamma_up(k)
        if (g >= 1.0_dp) then
            result = failed(reduction_range)
            return
        end if
        floor = mul_up(mul_up(real(k, dp), eta), add_up(1.0_dp, g))
        absolute_sum = div_up(add_up(magnitude, floor), sub_down(1.0_dp, g))
        result%radius = add_up(mul_up(g, absolute_sum), floor)
        if (.not. ieee_is_finite(result%radius)) result = failed(reduction_range)
    end subroutine finish

    function interval_sum(x) result(result)
        real(dp), intent(in) :: x(:)
        type(reduction_t) :: result
        type(interval_t) :: enclosure
        integer :: i

        if (.not. all(ieee_is_finite(x))) then
            result = failed(reduction_input)
            return
        end if
        enclosure = interval(0.0_dp)
        do i = 1, size(x)
            enclosure = enclosure + interval(x(i))
        end do
        result = from_interval(enclosure)
    end function interval_sum

    function interval_dot(x, y) result(result)
        real(dp), intent(in) :: x(:), y(:)
        type(reduction_t) :: result
        type(interval_t) :: enclosure
        integer :: i

        if (size(x) /= size(y)) then
            result = failed(reduction_input)
            return
        end if
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
            result = failed(reduction_input)
            return
        end if
        enclosure = interval(0.0_dp)
        do i = 1, size(x)
            enclosure = enclosure + interval(x(i))*interval(y(i))
        end do
        result = from_interval(enclosure)
    end function interval_dot

    function from_interval(enclosure) result(result)
        type(interval_t), intent(in) :: enclosure
        type(reduction_t) :: result

        if (.not. ieee_is_finite(enclosure%lo) .or. &
            .not. ieee_is_finite(enclosure%hi)) then
            result = failed(reduction_range)
            return
        end if
        result%center = 0.5_dp*enclosure%lo + 0.5_dp*enclosure%hi
        result%radius = max(sub_up(result%center, enclosure%lo), &
            sub_up(enclosure%hi, result%center))
        if (.not. ieee_is_finite(result%radius)) result = failed(reduction_range)
    end function from_interval

    function failed(status) result(result)
        integer, intent(in) :: status
        type(reduction_t) :: result

        result%status = status
        result%radius = ieee_value(1.0_dp, ieee_positive_inf)
    end function failed

end module fortnum_certified_reductions
