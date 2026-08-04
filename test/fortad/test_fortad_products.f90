program test_fortad_products
    !! fortad's products for the operators fortsym does not generate a
    !! derivative kernel for.
    !!
    !! There is no reference kernel to compare against here - fortsym emits the
    !! primal and, where it emits a derivative at all, emits it as its own
    !! analytical recurrence rather than as a product of this primal. So the
    !! oracle is the primal itself: central differences pin the tangent, and
    !! the adjoint identity pins the cotangent against the tangent.
    !!
    !! Neither check uses any derivative machinery. Differences of the primal
    !! and a dot product are all that is needed, which is the point.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fortad_sqrt1pm1_stable_jvp, only: fortnum_sqrt1pm1_stable_jvp_fortad
    use fortnum_fortad_sqrt1pm1_stable_vjp, only: fortnum_sqrt1pm1_stable_vjp_fortad
    use fortnum_fortad_sqrt1pm1_raw_jvp, only: fortnum_sqrt1pm1_raw_jvp_fortad
    use fortnum_fortad_sqrt1pm1_raw_vjp, only: fortnum_sqrt1pm1_raw_vjp_fortad
    use fortnum_fortad_dawson_outer_value_jvp, only: fortnum_dawson_outer_value_jvp_fortad
    use fortnum_fortad_dawson_outer_value_vjp, only: fortnum_dawson_outer_value_vjp_fortad
    use fortnum_fortad_toroidal_order_jvp, only: fortnum_toroidal_order_jvp_fortad
    use fortnum_fortad_toroidal_order_vjp, only: fortnum_toroidal_order_vjp_fortad
    use fortnum_fortad_scaled_jacobi_recurrence_jvp, only: fortnum_scaled_jacobi_recurrence_jvp_fortad
    use fortnum_fortad_scaled_jacobi_recurrence_vjp, only: fortnum_scaled_jacobi_recurrence_vjp_fortad
    implicit none

    interface
        subroutine fortnum_sqrt1pm1_stable(x, value)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value
        end subroutine fortnum_sqrt1pm1_stable
        subroutine fortnum_sqrt1pm1_raw(x, value)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value
        end subroutine fortnum_sqrt1pm1_raw
        subroutine fortnum_dawson_outer_value(f, value)
            import :: dp
            real(dp), intent(in) :: f
            real(dp), intent(out) :: value
        end subroutine fortnum_dawson_outer_value
        subroutine fortnum_toroidal_order(degree, order, x, current, next_order, following_order)
            import :: dp
            real(dp), intent(in) :: degree
            real(dp), intent(in) :: order
            real(dp), intent(in) :: x
            real(dp), intent(in) :: current
            real(dp), intent(in) :: next_order
            real(dp), intent(out) :: following_order
        end subroutine fortnum_toroidal_order
        subroutine fortnum_scaled_jacobi_recurrence(degree, alpha, beta, x, scale, previous, current, next)
            import :: dp
            real(dp), intent(in) :: degree
            real(dp), intent(in) :: alpha
            real(dp), intent(in) :: beta
            real(dp), intent(in) :: x
            real(dp), intent(in) :: scale
            real(dp), intent(in) :: previous
            real(dp), intent(in) :: current
            real(dp), intent(out) :: next
        end subroutine fortnum_scaled_jacobi_recurrence
    end interface

    real(dp), parameter :: h = 1.0e-6_dp
    real(dp), parameter :: u = 1.3_dp
    integer :: failures

    failures = 0
    call check(failures)

    if (failures == 0) then
        print *, "test_fortad_products: all cases passed"
    else
        print *, "test_fortad_products: ", failures, " case(s) FAILED"
        error stop 1
    end if

contains

    subroutine check(failures)
        integer, intent(inout) :: failures
        real(dp) :: p(8), d(8), g(8)
        real(dp) :: tangent, plus, minus
        integer :: k

        ! sqrt1pm1_stable
        p(1:1) = [0.35_dp]
        do k = 1, 1
            d(k) = 0.4_dp*k - 0.7_dp
        end do
        call fortnum_sqrt1pm1_stable_jvp_fortad(p(1), d(1), tangent)
        call fortnum_sqrt1pm1_stable(p(1) + h*d(1), plus)
        call fortnum_sqrt1pm1_stable(p(1) - h*d(1), minus)
        call differenced("sqrt1pm1_stable", tangent, (plus - minus)/(2.0_dp*h), failures)
        call fortnum_sqrt1pm1_stable_vjp_fortad(p(1), u, g(1))
        call same("sqrt1pm1_stable adjoint identity", &
                  dot_product(g(1:1), d(1:1)), tangent*u, failures)

        ! sqrt1pm1_raw
        p(1:1) = [0.35_dp]
        do k = 1, 1
            d(k) = 0.4_dp*k - 0.7_dp
        end do
        call fortnum_sqrt1pm1_raw_jvp_fortad(p(1), d(1), tangent)
        call fortnum_sqrt1pm1_raw(p(1) + h*d(1), plus)
        call fortnum_sqrt1pm1_raw(p(1) - h*d(1), minus)
        call differenced("sqrt1pm1_raw", tangent, (plus - minus)/(2.0_dp*h), failures)
        call fortnum_sqrt1pm1_raw_vjp_fortad(p(1), u, g(1))
        call same("sqrt1pm1_raw adjoint identity", &
                  dot_product(g(1:1), d(1:1)), tangent*u, failures)

        ! dawson_outer_value
        p(1:1) = [0.82_dp]
        do k = 1, 1
            d(k) = 0.4_dp*k - 0.7_dp
        end do
        call fortnum_dawson_outer_value_jvp_fortad(p(1), d(1), tangent)
        call fortnum_dawson_outer_value(p(1) + h*d(1), plus)
        call fortnum_dawson_outer_value(p(1) - h*d(1), minus)
        call differenced("dawson_outer_value", tangent, (plus - minus)/(2.0_dp*h), failures)
        call fortnum_dawson_outer_value_vjp_fortad(p(1), u, g(1))
        call same("dawson_outer_value adjoint identity", &
                  dot_product(g(1:1), d(1:1)), tangent*u, failures)

        ! toroidal_order
        p(1:5) = [3.0_dp, 1.0_dp, 1.7_dp, 0.63_dp, -0.41_dp]
        do k = 1, 3
            d(k) = 0.4_dp*k - 0.7_dp
        end do
        call fortnum_toroidal_order_jvp_fortad(p(1), p(2), p(3), d(1), p(4), d(2), p(5), d(3), tangent)
        call fortnum_toroidal_order(p(1), p(2), p(3) + h*d(1), p(4) + h*d(2), p(5) + h*d(3), plus)
        call fortnum_toroidal_order(p(1), p(2), p(3) - h*d(1), p(4) - h*d(2), p(5) - h*d(3), minus)
        call differenced("toroidal_order", tangent, (plus - minus)/(2.0_dp*h), failures)
        call fortnum_toroidal_order_vjp_fortad(p(1), p(2), p(3), p(4), p(5), u, g(1), g(2), g(3))
        call same("toroidal_order adjoint identity", &
                  dot_product(g(1:3), d(1:3)), tangent*u, failures)

        ! scaled_jacobi_recurrence
        p(1:7) = [4.0_dp, 0.5_dp, 0.25_dp, 0.31_dp, 1.2_dp, 0.77_dp, -0.44_dp]
        do k = 1, 4
            d(k) = 0.4_dp*k - 0.7_dp
        end do
        call fortnum_scaled_jacobi_recurrence_jvp_fortad(p(1), p(2), p(3), p(4), d(1), p(5), d(2), p(6), d(3), p(7), d(4), tangent)
        call fortnum_scaled_jacobi_recurrence(p(1), p(2), p(3), p(4) + h*d(1), p(5) + h*d(2), p(6) + h*d(3), p(7) + h*d(4), plus)
        call fortnum_scaled_jacobi_recurrence(p(1), p(2), p(3), p(4) - h*d(1), p(5) - h*d(2), p(6) - h*d(3), p(7) - h*d(4), minus)
        call differenced("scaled_jacobi_recurrence", tangent, (plus - minus)/(2.0_dp*h), failures)
        call fortnum_scaled_jacobi_recurrence_vjp_fortad(p(1), p(2), p(3), p(4), p(5), p(6), p(7), u, g(1), g(2), g(3), g(4))
        call same("scaled_jacobi_recurrence adjoint identity", &
                  dot_product(g(1:4), d(1:4)), tangent*u, failures)
    end subroutine check

    subroutine differenced(label, tangent, difference, failures)
        !! Central differences lose several digits to cancellation, so the
        !! tolerance is set by the step size rather than by the arithmetic.
        character(*), intent(in) :: label
        real(dp), intent(in) :: tangent, difference
        integer, intent(inout) :: failures

        if (abs(tangent - difference) > 1.0e-6_dp*max(1.0_dp, abs(difference))) then
            print *, "FAIL ", label, " against differences: ", difference, &
                " vs ", tangent
            failures = failures + 1
        else
            print *, "pass ", label, " against differences"
        end if
    end subroutine differenced

    subroutine same(label, reference, produced, failures)
        character(*), intent(in) :: label
        real(dp), intent(in) :: reference, produced
        integer, intent(inout) :: failures

        if (abs(produced - reference) > 1.0e-12_dp*max(1.0_dp, abs(reference))) then
            print *, "FAIL ", label, ": ", reference, " vs ", produced
            failures = failures + 1
        else
            print *, "pass ", label
        end if
    end subroutine same

end program test_fortad_products
