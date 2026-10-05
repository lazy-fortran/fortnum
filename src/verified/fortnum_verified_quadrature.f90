module fortnum_verified_quadrature
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    use fortnum_interval, only: interval_t, interval, hull, is_empty
    use fortnum_generated_verified_midpoint, only: midpoint_geometry, &
        midpoint_accumulate
    implicit none
    private
    public :: composite_midpoint, midpoint_callback
    public :: composite_midpoint_batch, midpoint_batch_callback

    abstract interface
        subroutine midpoint_callback(cell, midpoint, value, curvature, ok)
            import :: interval_t
            type(interval_t), intent(in) :: cell, midpoint
            type(interval_t), intent(out) :: value, curvature
            logical, intent(out) :: ok
        end subroutine midpoint_callback
        subroutine midpoint_batch_callback(cell, midpoint, value, curvature, ok)
            import :: interval_t
            type(interval_t), intent(in) :: cell, midpoint
            type(interval_t), intent(out) :: value(:), curvature(:)
            logical, intent(out) :: ok
        end subroutine midpoint_batch_callback
    end interface

contains

    ! Each stored edge is an exact real endpoint. On every cell the callback
    ! must enclose f at the mathematical midpoint and a nonnegative bound on
    ! |f''| throughout the cell. f must have an absolutely continuous first
    ! derivative and essentially bounded second derivative on that cell.
    ! This permits piecewise smooth functions with knots among the edges.
    ! Refinement, tails and root-bracket schedules remain caller responsibilities.
    ! No allocations occur here; no callback estimate is promoted to a proof.
    subroutine composite_midpoint(callback, edges, integral, ok)
        procedure(midpoint_callback) :: callback
        real(dp), intent(in) :: edges(:)
        type(interval_t), intent(out) :: integral
        logical, intent(out) :: ok
        type(interval_t) :: midpoint, width, value, curvature, acc, lower, upper
        integer :: j
        logical :: callback_ok

        ok = .false.
        integral%lo = ieee_value(0.0_dp, ieee_quiet_nan)
        integral%hi = integral%lo
        if (size(edges) < 2) return
        if (.not. all(ieee_is_finite(edges))) return
        do j = 1, size(edges) - 1
            if (edges(j + 1) <= edges(j)) return
        end do
        acc = interval(0.0_dp)
        do j = 1, size(edges) - 1
            call midpoint_geometry(interval(edges(j)), &
                interval(edges(j + 1)), midpoint, width)
            if (.not. finite_interval(midpoint)) return
            if (.not. finite_interval(width)) return
            call callback(interval(edges(j), edges(j + 1)), midpoint, &
                          value, curvature, callback_ok)
            if (.not. callback_ok) return
            if (.not. finite_interval(value)) return
            if (.not. finite_interval(curvature)) return
            if (curvature%lo < 0.0_dp) return
            call midpoint_accumulate(acc, width, value, curvature, &
                                             lower, upper)
            if (.not. finite_interval(lower)) return
            if (.not. finite_interval(upper)) return
            acc = hull(lower, upper)
        end do
        integral = acc
        ok = .true.
    end subroutine composite_midpoint

    ! The callback shares integrand work among caller-sized components. Scratch
    ! arrays exist once per invocation; no allocation occurs in the cell loop.
    ! A failure in any component invalidates the complete output packet.
    subroutine composite_midpoint_batch(callback, edges, integral, ok)
        procedure(midpoint_batch_callback) :: callback
        real(dp), intent(in) :: edges(:)
        type(interval_t), intent(out) :: integral(:)
        logical, intent(out) :: ok
        type(interval_t) :: value(size(integral)), curvature(size(integral))
        type(interval_t) :: acc(size(integral)), midpoint, width, lower, upper
        integer :: j, component
        logical :: callback_ok

        ok = .false.
        integral%lo = ieee_value(0.0_dp, ieee_quiet_nan)
        integral%hi = ieee_value(0.0_dp, ieee_quiet_nan)
        if (size(integral) < 1) return
        if (size(edges) < 2) return
        if (.not. all(ieee_is_finite(edges))) return
        do j = 1, size(edges) - 1
            if (edges(j + 1) <= edges(j)) return
        end do
        acc = interval(0.0_dp)
        do j = 1, size(edges) - 1
            call midpoint_geometry(interval(edges(j)), &
                interval(edges(j + 1)), midpoint, width)
            if (.not. finite_interval(midpoint)) return
            if (.not. finite_interval(width)) return
            call callback(interval(edges(j), edges(j + 1)), midpoint, &
                          value, curvature, callback_ok)
            if (.not. callback_ok) return
            do component = 1, size(integral)
                if (.not. finite_interval(value(component))) return
                if (.not. finite_interval(curvature(component))) return
                if (curvature(component)%lo < 0.0_dp) return
                call midpoint_accumulate(acc(component), width, &
                    value(component), curvature(component), lower, upper)
                if (.not. finite_interval(lower)) return
                if (.not. finite_interval(upper)) return
                acc(component) = hull(lower, upper)
            end do
        end do
        integral = acc
        ok = .true.
    end subroutine composite_midpoint_batch

    pure logical function finite_interval(value) result(ok)
        type(interval_t), intent(in) :: value

        ok = .false.
        if (is_empty(value)) return
        if (.not. ieee_is_finite(value%lo)) return
        if (.not. ieee_is_finite(value%hi)) return
        ok = .true.
    end function finite_interval

end module fortnum_verified_quadrature
