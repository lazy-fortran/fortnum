module dawson_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private
    public :: outer, outer_autodiff

    interface
        function fortnum_dawson_kernel(x) result(f) &
                bind(c, name="fortnum_dawson_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel

        function fortnum_dawson_kernel_autodiff(x) result(f) &
                bind(c, name="fortnum_dawson_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel_autodiff
    end interface

contains

    function outer(x) result(y) bind(c, name="fortnum_dawson_outer")
        real(c_double), value :: x
        real(c_double) :: y, f

        f = fortnum_dawson_kernel(x)
        y = sin(f) + f*f
    end function outer

    function outer_autodiff(x) result(y) &
            bind(c, name="fortnum_dawson_outer_autodiff")
        real(c_double), value :: x
        real(c_double) :: y, f

        f = fortnum_dawson_kernel_autodiff(x)
        y = sin(f) + f*f
    end function outer_autodiff

end module dawson_outer_kernel

program enzyme_dawson_hybrid
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use dawson_outer_kernel, only: outer
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    use fortnum_generated_enzyme_dawson_outer, only: &
        fortnum_enzyme_dawson_outer_jvp
    use fortnum_generated_enzyme_dawson_outer_autodiff, only: &
        fortnum_enzyme_dawson_outer_autodiff_jvp
    implicit none

    real(c_double), parameter :: direction = -0.4_c_double
    real(c_double), parameter :: h = 1.0e-6_c_double
    real(c_double) :: samples(fixture_sample_count), sink
    character(32) :: action, candidate, argument
    integer :: candidate_kind, repetitions
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

        function fortnum_dawson_kernel(x) result(f) &
                bind(c, name="fortnum_dawson_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel
    end interface

    call get_command_argument(1, argument)
    if (trim(argument) == "--benchmark") then
        action = "benchmark"
    else
        call read_fixture_environment("FORTNUM_DAWSON_ACTION", &
            "validate", action)
    end if
    if (trim(action) == "validate") then
        call validate_candidates()
        write (*, "(a)") "PASS"
        stop
    end if

    call read_fixture_environment("FORTNUM_DAWSON_CANDIDATE", "all", candidate)
    call read_fixture_integer("FORTNUM_DAWSON_ITERATIONS", 200000, &
        repetitions, valid)
    if (.not. valid .or. repetitions < 1) then
        error stop "invalid Dawson iteration count"
    end if
    call rule_counter_disable()

    if (trim(action) == "benchmark" .or. trim(action) == "--benchmark") then
        if (trim(candidate) == "all") then
            do candidate_kind = 1, 3
                call benchmark_candidate()
            end do
        else
            candidate_kind = parse_candidate(trim(candidate))
            call benchmark_candidate()
        end if
    else if (trim(action) == "peak-rss" .or. trim(action) == "--peak-rss") then
        candidate_kind = parse_candidate(trim(candidate))
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "Dawson benchmark produced NaN"

contains

    function analytical_jvp(x, dx) result(dy)
        real(c_double), intent(in) :: x, dx
        real(c_double) :: dy, f, value

        f = fortnum_dawson_kernel(x)
        call fortnum_dawson_outer_kernel(x, f, dx, value, dy)
    end function analytical_jvp

    subroutine validate_candidates()
        real(c_double), parameter :: points(3) = [ &
            0.2_c_double, 0.7_c_double, 6.0_c_double]
        real(c_double) :: got, reference, scale, x
        integer :: i

        do i = 1, size(points)
            x = points(i)
            reference = (outer(x + h*direction) - outer(x - h*direction)) &
                /(2.0_c_double*h)
            scale = max(1.0_c_double, abs(reference))
            got = analytical_jvp(x, direction)
            if (abs(got - reference)/scale > 2.0e-9_c_double) then
                error stop "analytical Dawson JVP mismatch"
            end if
            got = fortnum_enzyme_dawson_outer_autodiff_jvp(x, direction)
            if (abs(got - reference)/scale > 2.0e-9_c_double) then
                error stop "autodiff Dawson JVP mismatch"
            end if
            call rule_counter_reset()
            got = fortnum_enzyme_dawson_outer_jvp(x, direction)
            if (abs(got - reference)/scale > 2.0e-9_c_double) then
                error stop "hybrid Dawson JVP mismatch"
            end if
            if (rule_counter_calls() /= 1_c_int64_t) then
                error stop "analytical Dawson rule was not selected"
            end if
        end do
        call rule_counter_disable()
    end subroutine validate_candidates

    subroutine benchmark_candidate()
        real(c_double) :: median, mad

        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_result(candidate_name(candidate_kind), &
            repetitions, median, mad)
    end subroutine benchmark_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: elapsed_ns, local_sink, xi
        integer :: i

        local_sink = 0.0_c_double
        call timer%start()
        do i = 1, repetitions
            xi = 0.65_c_double + 0.001_c_double*real(mod(i, 101), c_double)
            select case (candidate_kind)
            case (1)
                local_sink = local_sink + analytical_jvp(xi, direction)
            case (2)
                local_sink = local_sink + &
                    fortnum_enzyme_dawson_outer_autodiff_jvp(xi, direction)
            case default
                local_sink = local_sink + &
                    fortnum_enzyme_dawson_outer_jvp(xi, direction)
            end select
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, c_double)
        sink = local_sink
    end function measure_candidate

    function parse_candidate(name) result(kind)
        character(*), intent(in) :: name
        integer :: kind

        select case (name)
        case ("analytical")
            kind = 1
        case ("autodiff")
            kind = 2
        case ("hybrid")
            kind = 3
        case default
            error stop "candidate must be analytical, autodiff, or hybrid"
        end select
    end function parse_candidate

    function candidate_name(kind) result(name)
        integer, intent(in) :: kind
        character(16) :: name

        select case (kind)
        case (1)
            name = "analytical"
        case (2)
            name = "autodiff"
        case default
            name = "hybrid"
        end select
    end function candidate_name

end program enzyme_dawson_hybrid
