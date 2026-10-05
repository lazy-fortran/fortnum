module verified_midpoint_benchmark_callbacks
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, exp, operator(*)
    implicit none
    integer :: active_component = 1
contains
    subroutine scalar_callback(cell, midpoint, value, curvature, ok)
        type(interval_t), intent(in) :: cell, midpoint
        type(interval_t), intent(out) :: value, curvature
        logical, intent(out) :: ok
        type(interval_t) :: derivative

        value = real(active_component, dp)*exp(midpoint)
        derivative = real(active_component, dp)*exp(cell)
        curvature = interval(0.0_dp, derivative%hi)
        ok = .true.
    end subroutine scalar_callback

    subroutine batch_callback(cell, midpoint, value, curvature, ok)
        type(interval_t), intent(in) :: cell, midpoint
        type(interval_t), intent(out) :: value(:), curvature(:)
        logical, intent(out) :: ok
        type(interval_t) :: common, derivative, scaled_derivative
        integer :: component

        common = exp(midpoint)
        derivative = exp(cell)
        do component = 1, size(value)
            value(component) = real(component, dp)*common
            scaled_derivative = real(component, dp)*derivative
            curvature(component) = interval(0.0_dp, scaled_derivative%hi)
        end do
        ok = .true.
    end subroutine batch_callback
end module verified_midpoint_benchmark_callbacks

program bench_verified_midpoint
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, rad
    use fortnum_verified_quadrature, only: composite_midpoint, composite_midpoint_batch
    use verified_midpoint_benchmark_callbacks, only: scalar_callback, batch_callback, &
        active_component
    implicit none
    integer, parameter :: cells = 128, components = 8, samples = 15, repeats = 20
    real(dp) :: edges(cells + 1), scalar_seconds, batch_seconds, checksum
    type(interval_t) :: result(components)
    integer(int64) :: before, after, rate
    integer :: j, sample, repeat
    logical :: ok

    do j = 1, cells + 1
        edges(j) = real(j - 1, dp)/real(cells, dp)
    end do
    call system_clock(count_rate=rate)
    do repeat = 1, 3
        do active_component = 1, components
            call composite_midpoint(scalar_callback, edges, &
                                    result(active_component), ok)
        end do
        call validate(result, ok)
        call composite_midpoint_batch(batch_callback, edges, result, ok)
        call validate(result, ok)
    end do
    checksum = 0.0_dp
    print '(a)', 'sample,cells,components,repeats,scalar_seconds,batch_seconds,checksum'
    do sample = 1, samples
        call system_clock(before)
        do repeat = 1, repeats
            do active_component = 1, components
                call composite_midpoint(scalar_callback, edges, &
                                        result(active_component), ok)
            end do
        end do
        call system_clock(after)
        scalar_seconds = real(after - before, dp)/real(rate, dp)/real(repeats, dp)
        call validate(result, ok)
        checksum = checksum + sum(result%lo) + sum(result%hi)
        call system_clock(before)
        do repeat = 1, repeats
            call composite_midpoint_batch(batch_callback, edges, result, ok)
        end do
        call system_clock(after)
        batch_seconds = real(after - before, dp)/real(rate, dp)/real(repeats, dp)
        call validate(result, ok)
        checksum = checksum + sum(result%lo) + sum(result%hi)
        print '(4(i0,","),3(es24.16,:,","))', sample, cells, components, repeats, &
            scalar_seconds, batch_seconds, checksum
    end do
contains
    subroutine validate(packet, accepted)
        type(interval_t), intent(in) :: packet(:)
        logical, intent(in) :: accepted
        real(qp) :: exact
        integer :: component

        if (.not. accepted) error stop 'benchmark quadrature failure'
        do component = 1, size(packet)
            if (.not. ieee_is_finite(packet(component)%lo)) error stop 'lower nonfinite'
            if (.not. ieee_is_finite(packet(component)%hi)) error stop 'upper nonfinite'
            exact = real(component, qp)*(exp(1.0_qp) - 1.0_qp)
            if (real(packet(component)%lo, qp) > exact .or. &
                real(packet(component)%hi, qp) < exact) error stop 'oracle excluded'
            if (real(rad(packet(component)), qp) > 0.01_qp*exact) &
                error stop 'one percent quadrature budget exceeded'
        end do
    end subroutine validate
end program bench_verified_midpoint
