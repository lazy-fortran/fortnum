program enzyme_hyperg_tournament
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    use fortnum_generated_enzyme_hyperg_outer, only: &
        fortnum_enzyme_hyperg_outer_jvp
    use fortnum_generated_enzyme_hyperg_outer_autodiff, only: &
        fortnum_enzyme_hyperg_outer_autodiff_jvp, &
        fortnum_enzyme_hyperg_outer_autodiff_vjp_scalar
    use fortnum_generated_hyperg_asymptotic_outer_jvp, only: &
        fortnum_hyperg_asymptotic_outer_jvp
    use fortnum_generated_hyperg_asymptotic_outer_vjp, only: &
        fortnum_hyperg_asymptotic_outer_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1_a1, &
        hyperg_1f1_a1_jvp, hyperg_1f1_a1_vjp
    use hyperg_outer_kernel, only: outer
    implicit none

    integer, parameter :: region_count = 3
    real(c_double), parameter :: region_x(region_count) = [ &
        -4.0_c_double, 4.0_c_double, 80.0_c_double]
    character(16), parameter :: region_names(region_count) = [ &
        character(16) :: "kummer", "series", "asymptotic"]
    complex(dp), parameter :: fixed_b = cmplx(2.5_dp, 0.0_dp, dp)
    real(c_double) :: samples(fixture_sample_count), sink, validation_error
    character(32) :: action, candidate, product, region_text
    integer :: region, directions, iterations, candidate_kind, product_kind
    logical :: valid

    interface
        function fortnum_hyperg_kernel(x) result(value) &
                bind(c, name="fortnum_hyperg_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_hyperg_kernel
        function fortnum_hyperg_kernel_autodiff(x) result(value) &
                bind(c, name="fortnum_hyperg_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_hyperg_kernel_autodiff
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

    call read_fixture_environment("FORTNUM_HYPERG_ACTION", "validate", action)
    if (trim(action) == "validate") then
        call validate_candidates(validation_error)
        write (*, "('PASS max_relative_error=',es12.5)") validation_error
        stop
    end if
    call rule_counter_disable()
    call read_fixture_environment("FORTNUM_HYPERG_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_environment("FORTNUM_HYPERG_PRODUCT", "jvp", product)
    call read_fixture_environment("FORTNUM_HYPERG_REGION", "series", region_text)
    call read_fixture_integer("FORTNUM_HYPERG_DIRECTIONS", 16, directions, valid)
    if (.not. valid) error stop "invalid direction count"
    call read_fixture_integer("FORTNUM_HYPERG_ITERATIONS", 10000, iterations, valid)
    if (.not. valid) error stop "invalid iteration count"
    region = parse_region(trim(region_text))
    call validate_benchmark_configuration()
    select case (trim(candidate))
    case ("analytical")
        candidate_kind = 1
    case ("autodiff")
        candidate_kind = 2
    case ("analytical_region")
        candidate_kind = 4
    case default
        candidate_kind = 3
    end select
    product_kind = merge(1, 2, trim(product) == "jvp")

    if (trim(action) == "benchmark" .or. trim(action) == "--benchmark") then
        call collect_fixture_samples(measure_candidate, samples)
        call report_samples()
    else if (trim(action) == "peak-rss" .or. trim(action) == "--peak-rss") then
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "hyperg benchmark produced NaN"

contains

    function analytical_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative, z(2), tangent(2), product_value(2), inner
        complex(dp) :: value
        type(fortnum_status_t) :: status

        call hyperg_1f1_a1(fixed_b, cmplx(x, 0.0_dp, dp), value, status)
        if (.not. status_ok(status)) error stop "analytical primal failed"
        inner = real(value, dp)
        z = [x, 0.0_dp]
        tangent = [direction, 0.0_dp]
        call hyperg_1f1_a1_jvp(z, fixed_b, tangent, product_value)
        derivative = 2.0_c_double*inner*product_value(1) &
            /(1.0_c_double + inner*inner)
    end function analytical_jvp

    function analytical_vjp(x, cotangent) result(derivative)
        real(c_double), intent(in) :: x, cotangent
        real(c_double) :: derivative, z(2), inner_bar(2), product_value(2), inner
        complex(dp) :: value
        type(fortnum_status_t) :: status

        call hyperg_1f1_a1(fixed_b, cmplx(x, 0.0_dp, dp), value, status)
        if (.not. status_ok(status)) error stop "analytical primal failed"
        inner = real(value, dp)
        z = [x, 0.0_dp]
        inner_bar = [2.0_c_double*inner*cotangent &
            /(1.0_c_double + inner*inner), 0.0_c_double]
        call hyperg_1f1_a1_vjp(z, fixed_b, inner_bar, product_value)
        derivative = product_value(1)
    end function analytical_vjp

    subroutine validate_candidates(max_relative_error)
        real(c_double), intent(out) :: max_relative_error
        real(c_double), parameter :: h = 2.0e-5_c_double
        real(c_double), parameter :: direction = -0.4_c_double
        real(c_double) :: x, reference, scale, actual, relative_error
        integer :: i

        max_relative_error = 0.0_c_double
        do i = 1, region_count
            x = region_x(i)
            scale = max(1.0_c_double, abs(fortnum_hyperg_kernel(x)))
            relative_error = abs(fortnum_hyperg_kernel_autodiff(x) - &
                fortnum_hyperg_kernel(x))/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 5.0e-12_c_double) &
                error stop "real-specialized hyperg primal mismatch"
            reference = (outer(x + h*direction) - outer(x - h*direction)) &
                /(2.0_c_double*h)
            scale = max(1.0_c_double, abs(reference))
            actual = analytical_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-7_c_double) &
                error stop "analytical hyperg JVP mismatch"
            actual = fortnum_enzyme_hyperg_outer_autodiff_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-7_c_double) &
                error stop "autodiff hyperg JVP mismatch"
            call rule_counter_reset()
            actual = fortnum_enzyme_hyperg_outer_jvp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-7_c_double) &
                error stop "hybrid hyperg JVP mismatch"
            if (rule_counter_calls() /= 1_c_int64_t) &
                error stop "analytical hyperg rule was not selected"
            actual = fortnum_enzyme_hyperg_outer_autodiff_vjp_scalar(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-7_c_double) &
                error stop "autodiff hyperg VJP mismatch"
            actual = analytical_vjp(x, direction)
            relative_error = abs(actual - reference)/scale
            max_relative_error = max(max_relative_error, relative_error)
            if (relative_error > 2.0e-7_c_double) &
                error stop "analytical hyperg VJP mismatch"
            if (i == region_count) then
                call fortnum_hyperg_asymptotic_outer_jvp(x, direction, actual)
                relative_error = abs(actual - reference)/scale
                max_relative_error = max(max_relative_error, relative_error)
                if (relative_error > 2.0e-7_c_double) &
                    error stop "region analytical hyperg JVP mismatch"
                call fortnum_hyperg_asymptotic_outer_vjp(x, direction, actual)
                relative_error = abs(actual - reference)/scale
                max_relative_error = max(max_relative_error, relative_error)
                if (relative_error > 2.0e-7_c_double) &
                    error stop "region analytical hyperg VJP mismatch"
            end if
        end do
        call rule_counter_disable()
    end subroutine validate_candidates

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: x, direction, derivative, elapsed_ns, local_sink
        integer :: direction_index, iteration

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            x = region_x(region) + 1.0e-7_c_double*real(mod(iteration, 101) - 50, c_double)
            do direction_index = 1, directions
                direction = 0.03_c_double*real( &
                    mod(5*direction_index, 11) - 5, c_double)
                if (product_kind == 1) then
                    select case (candidate_kind)
                    case (1)
                        local_sink = local_sink + analytical_jvp(x, direction)
                    case (2)
                        local_sink = local_sink + &
                            fortnum_enzyme_hyperg_outer_autodiff_jvp(x, direction)
                    case (4)
                        call fortnum_hyperg_asymptotic_outer_jvp( &
                            x, direction, derivative)
                        local_sink = local_sink + derivative
                    case default
                        local_sink = local_sink + &
                            fortnum_enzyme_hyperg_outer_jvp(x, direction)
                    end select
                else
                    select case (candidate_kind)
                    case (1)
                        local_sink = local_sink + analytical_vjp(x, direction)
                    case (4)
                        call fortnum_hyperg_asymptotic_outer_vjp( &
                            x, direction, derivative)
                        local_sink = local_sink + derivative
                    case default
                        local_sink = local_sink + &
                            fortnum_enzyme_hyperg_outer_autodiff_vjp_scalar(x, direction)
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
            trim(candidate) /= "autodiff" .and. trim(candidate) /= "hybrid" .and. &
            trim(candidate) /= "analytical_region") &
            error stop "invalid candidate"
        if (trim(product) /= "jvp" .and. trim(product) /= "vjp") &
            error stop "product must be jvp or vjp"
        if (trim(product) == "vjp" .and. trim(candidate) == "hybrid") &
            error stop "hybrid reverse custom rule is not implemented"
        if (trim(candidate) == "analytical_region" .and. region /= region_count) &
            error stop "region analytical candidate requires asymptotic region"
        if (directions < 1 .or. directions > 16) &
            error stop "directions must be 1..16"
        if (iterations < 1) error stop "iterations must be positive"
    end subroutine validate_benchmark_configuration

    function parse_region(name) result(index)
        character(*), intent(in) :: name
        integer :: index

        select case (name)
        case ("kummer")
            index = 1
        case ("", "series")
            index = 2
        case ("asymptotic")
            index = 3
        case default
            error stop "region must be kummer, series, or asymptotic"
        end select
    end function parse_region

end program enzyme_hyperg_tournament
