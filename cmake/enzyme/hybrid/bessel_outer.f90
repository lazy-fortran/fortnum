module bessel_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private

    public :: outer, outer_autodiff

    interface
        function fortnum_bessel_i0_kernel(x) result(value) &
                bind(c, name="fortnum_bessel_i0_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_bessel_i0_kernel

        function fortnum_bessel_i0_kernel_autodiff(x) result(value) &
                bind(c, name="fortnum_bessel_i0_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_bessel_i0_kernel_autodiff
    end interface

contains

    function outer(x) result(value) bind(c, name="fortnum_bessel_outer")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_bessel_i0_kernel(x)
        value = sin(inner) + inner*inner
    end function outer

    function outer_autodiff(x) result(value) &
            bind(c, name="fortnum_bessel_outer_autodiff")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_bessel_i0_kernel_autodiff(x)
        value = sin(inner) + inner*inner
    end function outer_autodiff

end module bessel_outer_kernel

program enzyme_bessel_tournament
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    use fortnum_generated_enzyme_bessel_outer, only: &
        fortnum_enzyme_bessel_outer_jvp
    use fortnum_generated_enzyme_bessel_outer_autodiff, only: &
        fortnum_enzyme_bessel_outer_autodiff_jvp, &
        fortnum_enzyme_bessel_outer_autodiff_vjp_scalar
    use fortnum_generated_bessel_outer_jvp, only: &
        fortnum_bessel_outer_jvp_kernel
    use bessel_outer_kernel, only: outer
    use fortnum_special_bessel, only: bessel_in
    implicit none

    integer, parameter :: region_count = 3
    real(c_double), parameter :: region_centers(region_count) = [ &
        1.0_c_double, 4.0_c_double, 24.0_c_double]
    character(16), parameter :: region_names(region_count) = [ &
        character(16) :: "series", "recurrence", "asymptotic"]
    real(c_double) :: samples(fixture_sample_count), sink, validation_error
    character(32) :: action, candidate, product, region_text
    integer :: region, directions, iterations, candidate_kind, product_kind
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

    call read_fixture_environment("FORTNUM_BESSEL_ACTION", "validate", action)
    if (trim(action) == "validate") then
        call validate_candidates(validation_error)
        write (*, "('PASS max_relative_error=',es12.5)") validation_error
        stop
    end if

    call rule_counter_disable()
    call read_fixture_environment("FORTNUM_BESSEL_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_environment("FORTNUM_BESSEL_PRODUCT", "jvp", product)
    call read_fixture_environment("FORTNUM_BESSEL_REGION", "series", &
        region_text)
    call read_fixture_integer("FORTNUM_BESSEL_DIRECTIONS", 16, directions, &
        valid)
    if (.not. valid) error stop "invalid direction count"
    call read_fixture_integer("FORTNUM_BESSEL_ITERATIONS", 200000, iterations, &
        valid)
    if (.not. valid) error stop "invalid iteration count"
    region = parse_region(trim(region_text))
    call validate_benchmark_configuration()
    select case (trim(candidate))
    case ("analytical")
        candidate_kind = 1
    case ("autodiff")
        candidate_kind = 2
    case default
        candidate_kind = 3
    end select
    if (trim(product) == "jvp") then
        product_kind = 1
    else
        product_kind = 2
    end if

    if (trim(action) == "--benchmark" .or. trim(action) == "benchmark") then
        call collect_fixture_samples(measure_candidate, samples)
        call report_samples()
    else if (trim(action) == "--peak-rss" .or. trim(action) == "peak-rss") then
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "Bessel benchmark produced NaN"

contains

    function analytical_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative, inner, inner_derivative

        inner = bessel_in(0, x)
        inner_derivative = bessel_in(1, x)
        call fortnum_bessel_outer_jvp_kernel( &
            inner, inner_derivative, direction, derivative)
    end function analytical_jvp

    function autodiff_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative

        derivative = fortnum_enzyme_bessel_outer_autodiff_jvp(x, direction)
    end function autodiff_jvp

    function hybrid_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative

        derivative = fortnum_enzyme_bessel_outer_jvp(x, direction)
    end function hybrid_jvp

    function autodiff_vjp(x, cotangent) result(derivative)
        real(c_double), intent(in) :: x, cotangent
        real(c_double) :: derivative

        derivative = fortnum_enzyme_bessel_outer_autodiff_vjp_scalar( &
            x, cotangent)
    end function autodiff_vjp

    subroutine validate_candidates(max_relative_error)
        real(c_double), intent(out) :: max_relative_error
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double), parameter :: direction = -0.4_c_double
        real(c_double) :: x, reference, scale, actual, relative_error
        integer :: i

        max_relative_error = 0.0_c_double
        do i = 1, region_count
            x = region_centers(i)
            reference = (outer(x + h*direction) - outer(x - h*direction)) &
                /(2.0_c_double*h)
            scale = max(1.0_c_double, abs(reference))
            actual = analytical_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-8_c_double) then
                error stop "analytical Bessel JVP mismatch"
            end if
            actual = autodiff_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-8_c_double) then
                error stop "autodiff Bessel JVP mismatch"
            end if
            call rule_counter_reset()
            actual = hybrid_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-8_c_double) then
                error stop "hybrid Bessel JVP mismatch"
            end if
            if (rule_counter_calls() /= 1_c_int64_t) then
                error stop "analytical Bessel rule was not selected"
            end if
            actual = autodiff_vjp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-8_c_double) then
                error stop "autodiff Bessel VJP mismatch"
            end if
        end do
        call rule_counter_disable()
    end subroutine validate_candidates

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: direction, elapsed_ns, local_sink, x
        integer :: direction_index, iteration

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            x = region_centers(region) + 1.0e-6_c_double*real( &
                mod(iteration, 101) - 50, c_double)
            do direction_index = 1, directions
                direction = 0.03_c_double*real( &
                    mod(5*direction_index, 11) - 5, c_double)
                if (product_kind == 1) then
                    select case (candidate_kind)
                    case (1)
                        local_sink = local_sink + analytical_jvp(x, direction)
                    case (2)
                        local_sink = local_sink + autodiff_jvp(x, direction)
                    case default
                        local_sink = local_sink + hybrid_jvp(x, direction)
                    end select
                else
                    select case (candidate_kind)
                    case (1)
                        local_sink = local_sink + analytical_jvp(x, direction)
                    case default
                        local_sink = local_sink + autodiff_vjp(x, direction)
                    end select
                end if
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

    subroutine report_samples()
        real(c_double) :: median, mad
        character(128) :: name

        call median_mad(samples, median, mad)
        write (name, "(a,'_',a,'_',a,'_d',i0)") trim(candidate), &
            trim(product), trim(region_names(region)), directions
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine report_samples

    subroutine validate_benchmark_configuration()
        if (trim(candidate) /= "analytical" .and. &
            trim(candidate) /= "autodiff" .and. &
            trim(candidate) /= "hybrid") then
            error stop "candidate must be analytical, autodiff, or hybrid"
        end if
        if (trim(product) /= "jvp" .and. trim(product) /= "vjp") then
            error stop "product must be jvp or vjp"
        end if
        if (trim(product) == "vjp" .and. trim(candidate) == "hybrid") then
            error stop "hybrid reverse custom rule is not implemented"
        end if
        if (directions < 1 .or. directions > 16) then
            error stop "directions must be 1..16"
        end if
        if (iterations < 1) error stop "iterations must be positive"
    end subroutine validate_benchmark_configuration

    function parse_region(name) result(index)
        character(*), intent(in) :: name
        integer :: index

        select case (name)
        case ("", "series")
            index = 1
        case ("recurrence")
            index = 2
        case ("asymptotic")
            index = 3
        case default
            error stop "region must be series, recurrence, or asymptotic"
        end select
    end function parse_region

end program enzyme_bessel_tournament
