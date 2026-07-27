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
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr, &
        c_int, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    use bessel_outer_kernel, only: outer, outer_autodiff
    use fortnum_special_bessel, only: bessel_in
    implicit none

    interface
        function enzyme_fwddiff(f, x, dx) result(derivative) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x, dx
            real(c_double) :: derivative
        end function enzyme_fwddiff

        function enzyme_autodiff(f, x) result(derivative) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x
            real(c_double) :: derivative
        end function enzyme_autodiff

        subroutine rule_reset() bind(c, name="fortnum_bessel_rule_reset")
        end subroutine rule_reset

        function rule_calls() result(count) &
                bind(c, name="fortnum_bessel_rule_calls")
            import :: c_int
            integer(c_int) :: count
        end function rule_calls

        subroutine rule_disable_count() &
                bind(c, name="fortnum_bessel_rule_disable_count")
        end subroutine rule_disable_count

        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    integer, parameter :: region_count = 3
    real(c_double), parameter :: region_centers(region_count) = [ &
        1.0_c_double, 4.0_c_double, 24.0_c_double]
    character(16), parameter :: region_names(region_count) = [ &
        character(16) :: "series", "recurrence", "asymptotic"]
    character(32) :: action, candidate, product, region_text, direction_text
    integer :: region, directions
    integer(int64) :: iterations

    call get_environment_variable("FORTNUM_BESSEL_ACTION", action)
    call get_environment_variable("FORTNUM_BESSEL_CANDIDATE", candidate)
    call get_environment_variable("FORTNUM_BESSEL_PRODUCT", product)
    call get_environment_variable("FORTNUM_BESSEL_REGION", region_text)
    call get_environment_variable("FORTNUM_BESSEL_DIRECTIONS", direction_text)
    call read_int64_env("FORTNUM_BESSEL_ITERATIONS", iterations, 200000_int64)
    region = parse_region(trim(region_text))
    directions = parse_integer(direction_text, 16)

    if (trim(action) == "--benchmark") then
        call rule_disable_count()
        call benchmark_candidate(trim(candidate), trim(product), region, &
            directions, iterations)
    else if (trim(action) == "--peak-rss") then
        call rule_disable_count()
        call benchmark_candidate(trim(candidate), trim(product), region, &
            directions, iterations)
        write (*, "(a,i0)") "peak_rss_bytes=", peak_rss_bytes()
    else
        call validate_candidates()
        write (*, "(a)") "PASS"
    end if

contains

    function analytical_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative, inner, inner_derivative

        inner = bessel_in(0, x)
        inner_derivative = bessel_in(1, x)
        derivative = (cos(inner) + 2.0_c_double*inner)*inner_derivative*direction
    end function analytical_jvp

    function autodiff_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative

        derivative = enzyme_fwddiff(c_funloc(outer_autodiff), x, direction)
    end function autodiff_jvp

    function hybrid_jvp(x, direction) result(derivative)
        real(c_double), intent(in) :: x, direction
        real(c_double) :: derivative

        derivative = enzyme_fwddiff(c_funloc(outer), x, direction)
    end function hybrid_jvp

    function autodiff_vjp(x, cotangent) result(derivative)
        real(c_double), intent(in) :: x, cotangent
        real(c_double) :: derivative

        derivative = cotangent*enzyme_autodiff(c_funloc(outer_autodiff), x)
    end function autodiff_vjp

    subroutine validate_candidates()
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double), parameter :: direction = -0.4_c_double
        real(c_double) :: x, reference, scale
        integer :: i

        do i = 1, region_count
            x = region_centers(i)
            reference = (outer(x + h*direction) - outer(x - h*direction)) &
                /(2.0_c_double*h)
            scale = max(1.0_c_double, abs(reference))
            if (abs(analytical_jvp(x, direction) - reference)/scale &
                > 2.0e-8_c_double) then
                error stop "analytical Bessel JVP mismatch"
            end if
            if (abs(autodiff_jvp(x, direction) - reference)/scale &
                > 2.0e-8_c_double) then
                error stop "autodiff Bessel JVP mismatch"
            end if
            call rule_reset()
            if (abs(hybrid_jvp(x, direction) - reference)/scale &
                > 2.0e-8_c_double) then
                error stop "hybrid Bessel JVP mismatch"
            end if
            if (rule_calls() /= 1_c_int) then
                error stop "analytical Bessel rule was not selected"
            end if
            if (abs(autodiff_vjp(x, direction) - reference)/scale &
                > 2.0e-8_c_double) then
                error stop "autodiff Bessel VJP mismatch"
            end if
        end do
    end subroutine validate_candidates

    subroutine benchmark_candidate(name, derivative_product, region_index, &
            direction_count, count)
        character(*), intent(in) :: name, derivative_product
        integer, intent(in) :: region_index, direction_count
        integer(int64), intent(in) :: count
        integer :: direction_index
        integer(int64) :: iteration, start, finish, rate
        real(c_double) :: direction, elapsed_ns, sink, x

        if (name /= "analytical" .and. name /= "autodiff" .and. &
            name /= "hybrid") then
            error stop "candidate must be analytical, autodiff, or hybrid"
        end if
        if (derivative_product /= "jvp" .and. derivative_product /= "vjp") then
            error stop "product must be jvp or vjp"
        end if
        if (derivative_product == "vjp" .and. name == "hybrid") then
            error stop "hybrid reverse custom rule is not implemented"
        end if
        if (direction_count < 1 .or. direction_count > 16) then
            error stop "directions must be 1..16"
        end if
        sink = 0.0_c_double
        call system_clock(start, rate)
        do iteration = 1, count
            x = region_centers(region_index) + 1.0e-6_c_double*real( &
                mod(iteration, 101_int64) - 50_int64, c_double)
            do direction_index = 1, direction_count
                direction = 0.03_c_double*real( &
                    mod(5*direction_index, 11) - 5, c_double)
                if (derivative_product == "jvp") then
                    select case (name)
                    case ("analytical")
                        sink = sink + analytical_jvp(x, direction)
                    case ("autodiff")
                        sink = sink + autodiff_jvp(x, direction)
                    case ("hybrid")
                        sink = sink + hybrid_jvp(x, direction)
                    end select
                else
                    select case (name)
                    case ("analytical")
                        sink = sink + analytical_jvp(x, direction)
                    case ("autodiff")
                        sink = sink + autodiff_vjp(x, direction)
                    end select
                end if
            end do
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = real(finish - start, c_double)*1.0e9_c_double &
            /(real(rate, c_double)*real(count, c_double))
        write (*, "(a,a,a,a,a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", &
            name, " product=", derivative_product, " region=", &
            trim(region_names(region_index)), &
            " directions=", direction_count, " iterations=", count, &
            " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

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

    function parse_integer(text, default_value) result(value)
        character(*), intent(in) :: text
        integer, intent(in) :: default_value
        integer :: value, ios

        value = default_value
        if (len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid integer environment value"
    end function parse_integer

    subroutine read_int64_env(name, value, default_value)
        character(*), intent(in) :: name
        integer(int64), intent(out) :: value
        integer(int64), intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid int64 environment value"
    end subroutine read_int64_env

end program enzyme_bessel_tournament
