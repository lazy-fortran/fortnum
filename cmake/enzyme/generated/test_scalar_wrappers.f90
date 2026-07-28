module scalar_wrapper_primals
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private

    public :: primal_p1, primal_p2, primal_p3, primal_p4
    public :: analytical_p1_jvp

contains

    pure function primal_p1(x1) result(value) &
            bind(c, name="fortnum_scalar_primal_p1")
        real(c_double), value :: x1
        real(c_double) :: value

        value = sin(x1) + 0.5_c_double*x1*x1
    end function primal_p1

    pure function analytical_p1_jvp(x1, tangent1) result(value) &
            bind(c, name="fortnum_scalar_analytical_p1_jvp")
        real(c_double), value :: x1, tangent1
        real(c_double) :: value

        value = (cos(x1) + x1)*tangent1
    end function analytical_p1_jvp

    pure function primal_p2(x1, x2) result(value) &
            bind(c, name="fortnum_scalar_primal_p2")
        real(c_double), value :: x1, x2
        real(c_double) :: value

        value = sin(x1) + x1*x2 + 2.0_c_double*x2*x2
    end function primal_p2

    pure function primal_p3(x1, x2, x3) result(value) &
            bind(c, name="fortnum_scalar_primal_p3")
        real(c_double), value :: x1, x2, x3
        real(c_double) :: value

        value = sin(x1) + x1*x2 + 2.0_c_double*x2*x2 + x2*x3 + &
            3.0_c_double*x3*x3
    end function primal_p3

    pure function primal_p4(x1, x2, x3, x4) result(value) &
            bind(c, name="fortnum_scalar_primal_p4")
        real(c_double), value :: x1, x2, x3, x4
        real(c_double) :: value

        value = sin(x1) + x1*x2 + 2.0_c_double*x2*x2 + x2*x3 + &
            3.0_c_double*x3*x3 + x3*x4 + 4.0_c_double*x4*x4
    end function primal_p4

end module scalar_wrapper_primals

program test_scalar_wrappers
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    use fortnum_generated_enzyme_scalar_p1, only: &
        fortnum_enzyme_scalar_p1_jvp, fortnum_enzyme_scalar_p1_vjp
    use fortnum_generated_enzyme_scalar_p2, only: &
        fortnum_enzyme_scalar_p2_jvp, fortnum_enzyme_scalar_p2_vjp
    use fortnum_generated_enzyme_scalar_p3, only: &
        fortnum_enzyme_scalar_p3_jvp, fortnum_enzyme_scalar_p3_vjp
    use fortnum_generated_enzyme_scalar_p4, only: &
        fortnum_enzyme_scalar_p4_jvp, fortnum_enzyme_scalar_p4_vjp
    implicit none

    integer, parameter :: repetitions = 200000
    real(dp), parameter :: tolerance = 3.0e-13_dp
    real(dp) :: samples(fixture_sample_count)
    real(dp) :: sink
    character(16) :: action, counting, product
    integer :: active_inputs
    logical :: valid

    interface
        subroutine rule_counter_reset() &
                bind(c, name="fortnum_enzyme_rule_counter_reset")
        end subroutine rule_counter_reset

        function rule_counter_calls() result(count) &
                bind(c, name="fortnum_enzyme_rule_counter_calls")
            import c_int64_t
            integer(c_int64_t) :: count
        end function rule_counter_calls

        subroutine rule_counter_disable() &
                bind(c, name="fortnum_enzyme_rule_counter_disable")
        end subroutine rule_counter_disable
    end interface

    call read_fixture_environment("FORTNUM_SCALAR_WRAPPER_ACTION", &
        "validate", action)
    if (trim(action) == "validate") then
        call prove_custom_rule_selection()
        call validate_wrappers()
        stop
    end if
    call read_fixture_environment("FORTNUM_SCALAR_WRAPPER_COUNTING", &
        "disabled", counting)
    if (trim(counting) == "enabled") then
        call rule_counter_reset()
    else if (trim(counting) == "disabled") then
        call rule_counter_disable()
    else
        error stop "counting must be enabled or disabled"
    end if

    call read_fixture_environment("FORTNUM_SCALAR_WRAPPER_PRODUCT", &
        "jvp", product)
    call read_fixture_integer("FORTNUM_SCALAR_WRAPPER_INPUTS", 1, &
        active_inputs, valid)
    if (.not. valid) error stop "invalid active-input count"
    if (active_inputs < 1 .or. active_inputs > 4) then
        error stop "active-input count must be one to four"
    end if
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
        error stop "product must be jvp or vjp"
    end if

    if (trim(action) == "peak-rss") then
        sink = measure_product()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else if (trim(action) == "benchmark") then
        call collect_fixture_samples(measure_product, samples)
        call report_samples()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "scalar-wrapper benchmark produced NaN"

contains

    subroutine prove_custom_rule_selection()
        real(dp) :: derivative

        call rule_counter_reset()
        derivative = fortnum_enzyme_scalar_p1_jvp(0.3_dp, -0.7_dp)
        if (derivative /= derivative) error stop "custom rule produced NaN"
        if (rule_counter_calls() /= 1_c_int64_t) then
            error stop "Enzyme did not select the analytical forward rule"
        end if
        call rule_counter_disable()
    end subroutine prove_custom_rule_selection

    subroutine validate_wrappers()
        real(dp) :: x(4), tangent(4), gradient(4), cotangents(4)
        real(dp) :: jvp, cotangent
        integer :: count, sample

        do count = 1, 4
            do sample = 1, 257
                x = [0.11_dp, -0.23_dp, 0.37_dp, -0.41_dp] + &
                    1.0e-3_dp*real(sample - 129, dp)
                tangent = [0.7_dp, -0.5_dp, 0.3_dp, -0.2_dp] + &
                    2.0e-4_dp*real(sample - 129, dp)
                cotangent = -0.8_dp + 1.0e-3_dp*real(sample, dp)
                call exact_gradient(count, x, gradient)
                call wrapper_jvp(count, x, tangent, jvp)
                if (abs(jvp - dot_product(gradient(:count), &
                    tangent(:count))) > tolerance) then
                    error stop "generated scalar JVP disagrees with formula"
                end if
                call wrapper_vjp(count, x, cotangent, cotangents)
                if (maxval(abs(cotangents(:count) - &
                    cotangent*gradient(:count))) > tolerance) then
                    error stop "generated scalar VJP disagrees with formula"
                end if
                if (abs(jvp*cotangent - dot_product(tangent(:count), &
                    cotangents(:count))) > tolerance) then
                    error stop "generated scalar wrapper fails adjoint identity"
                end if
            end do
        end do
    end subroutine validate_wrappers

    subroutine exact_gradient(count, x, gradient)
        integer, intent(in) :: count
        real(dp), intent(in) :: x(4)
        real(dp), intent(out) :: gradient(4)

        gradient = 0.0_dp
        gradient(1) = cos(x(1)) + x(2)
        if (count == 1) gradient(1) = cos(x(1)) + x(1)
        if (count >= 2) gradient(2) = x(1) + 4.0_dp*x(2)
        if (count >= 3) then
            gradient(2) = gradient(2) + x(3)
            gradient(3) = x(2) + 6.0_dp*x(3)
        end if
        if (count >= 4) then
            gradient(3) = gradient(3) + x(4)
            gradient(4) = x(3) + 8.0_dp*x(4)
        end if
    end subroutine exact_gradient

    subroutine wrapper_jvp(count, x, tangent, value)
        integer, intent(in) :: count
        real(dp), intent(in) :: x(4), tangent(4)
        real(dp), intent(out) :: value

        select case (count)
        case (1)
            value = fortnum_enzyme_scalar_p1_jvp(x(1), tangent(1))
        case (2)
            value = fortnum_enzyme_scalar_p2_jvp( &
                x(1), tangent(1), x(2), tangent(2))
        case (3)
            value = fortnum_enzyme_scalar_p3_jvp( &
                x(1), tangent(1), x(2), tangent(2), x(3), tangent(3))
        case default
            value = fortnum_enzyme_scalar_p4_jvp( &
                x(1), tangent(1), x(2), tangent(2), x(3), tangent(3), &
                x(4), tangent(4))
        end select
    end subroutine wrapper_jvp

    subroutine wrapper_vjp(count, x, cotangent, values)
        integer, intent(in) :: count
        real(dp), intent(in) :: x(4), cotangent
        real(dp), intent(out) :: values(4)

        values = 0.0_dp
        select case (count)
        case (1)
            call fortnum_enzyme_scalar_p1_vjp( &
                x(1), cotangent, values(1))
        case (2)
            call fortnum_enzyme_scalar_p2_vjp( &
                x(1), x(2), cotangent, values(1), values(2))
        case (3)
            call fortnum_enzyme_scalar_p3_vjp( &
                x(1), x(2), x(3), cotangent, values(1), values(2), values(3))
        case default
            call fortnum_enzyme_scalar_p4_vjp( &
                x(1), x(2), x(3), x(4), cotangent, values(1), values(2), &
                values(3), values(4))
        end select
    end subroutine wrapper_vjp

    function measure_product() result(nanoseconds)
        type(fixture_timer_t) :: timer
        real(dp) :: nanoseconds
        real(dp) :: x(4), tangent(4), values(4), value
        integer :: iteration

        x = [0.11_dp, -0.23_dp, 0.37_dp, -0.41_dp]
        tangent = [0.7_dp, -0.5_dp, 0.3_dp, -0.2_dp]
        value = 0.0_dp
        call timer%start()
        if (trim(product) == "jvp") then
            do iteration = 1, repetitions
                x(1) = 0.11_dp + 1.0e-8_dp*real(mod(iteration, 1021), dp)
                call wrapper_jvp(active_inputs, x, tangent, sink)
                value = value + sink
            end do
        else
            do iteration = 1, repetitions
                x(1) = 0.11_dp + 1.0e-8_dp*real(mod(iteration, 1021), dp)
                call wrapper_vjp(active_inputs, x, -0.8_dp, values)
                value = value + sum(values(:active_inputs))
            end do
        end if
        nanoseconds = timer%elapsed_ns()/real(repetitions, dp)
        sink = value
    end function measure_product

    subroutine report_samples()
        real(dp) :: median, mad
        character(32) :: name

        call median_mad(samples, median, mad)
        write (name, "(a,'_p',i0)") trim(product), active_inputs
        call write_fixture_result(trim(name), repetitions, median, mad)
    end subroutine report_samples

end program test_scalar_wrappers
