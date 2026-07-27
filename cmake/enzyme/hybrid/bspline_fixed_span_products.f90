module bspline_fixed_span_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    implicit none
    private

    integer, parameter, public :: coefficient_count = 4
    public :: analytical_jvp, analytical_vjp, autodiff_jvp, autodiff_vjp
    public :: spline_value

    interface
        function enzyme_fwddiff(f, x, dx, coefficients, dcoefficients) &
                result(derivative) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x, dx
            real(c_double), intent(in) :: coefficients(*), dcoefficients(*)
            real(c_double) :: derivative
        end function enzyme_fwddiff

        function enzyme_autodiff(f, x, xbar, coefficients, coefficients_bar) &
                result(value) bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x
            real(c_double), intent(inout) :: xbar
            real(c_double), intent(in) :: coefficients(*)
            real(c_double), intent(inout) :: coefficients_bar(*)
            real(c_double) :: value
        end function enzyme_autodiff
    end interface

contains

    pure function spline_value(x, coefficients) result(value) &
            bind(c, name="fortnum_bspline_fixed_span_value")
        real(c_double), intent(in) :: x
        real(c_double), intent(in) :: coefficients(coefficient_count)
        real(c_double) :: value, basis(coefficient_count)

        call cubic_uniform_basis(x, basis)
        value = dot_product(coefficients, basis)
    end function spline_value

    pure subroutine cubic_uniform_basis(x, basis, derivative)
        real(c_double), intent(in) :: x
        real(c_double), intent(out) :: basis(coefficient_count)
        real(c_double), intent(out), optional :: derivative(coefficient_count)
        real(c_double) :: x2, x3

        x2 = x*x
        x3 = x2*x
        basis(1) = (1.0_c_double - 3.0_c_double*x + 3.0_c_double*x2 - x3) &
            /6.0_c_double
        basis(2) = (4.0_c_double - 6.0_c_double*x2 + 3.0_c_double*x3) &
            /6.0_c_double
        basis(3) = (1.0_c_double + 3.0_c_double*x + 3.0_c_double*x2 &
            - 3.0_c_double*x3)/6.0_c_double
        basis(4) = x3/6.0_c_double
        if (present(derivative)) then
            derivative(1) = (-3.0_c_double + 6.0_c_double*x &
                - 3.0_c_double*x2)/6.0_c_double
            derivative(2) = (-12.0_c_double*x + 9.0_c_double*x2) &
                /6.0_c_double
            derivative(3) = (3.0_c_double + 6.0_c_double*x &
                - 9.0_c_double*x2)/6.0_c_double
            derivative(4) = 3.0_c_double*x2/6.0_c_double
        end if
    end subroutine cubic_uniform_basis

    pure function analytical_jvp(x, coefficients, dx, dcoefficients) &
            result(derivative)
        real(c_double), intent(in) :: x, coefficients(coefficient_count), dx
        real(c_double), intent(in) :: dcoefficients(coefficient_count)
        real(c_double) :: derivative, basis(coefficient_count)
        real(c_double) :: basis_derivative(coefficient_count)

        call cubic_uniform_basis(x, basis, basis_derivative)
        derivative = dx*dot_product(coefficients, basis_derivative) &
            + dot_product(dcoefficients, basis)
    end function analytical_jvp

    function autodiff_jvp(x, coefficients, dx, dcoefficients) result(derivative)
        real(c_double), intent(in) :: x, coefficients(coefficient_count), dx
        real(c_double), intent(in) :: dcoefficients(coefficient_count)
        real(c_double) :: derivative

        derivative = enzyme_fwddiff(c_funloc(spline_value), x, dx, &
            coefficients, dcoefficients)
    end function autodiff_jvp

    pure subroutine analytical_vjp(x, coefficients, cotangent, xbar, &
            coefficients_bar)
        real(c_double), intent(in) :: x, coefficients(coefficient_count)
        real(c_double), intent(in) :: cotangent
        real(c_double), intent(out) :: xbar
        real(c_double), intent(out) :: coefficients_bar(coefficient_count)
        real(c_double) :: basis(coefficient_count)
        real(c_double) :: basis_derivative(coefficient_count)

        call cubic_uniform_basis(x, basis, basis_derivative)
        xbar = cotangent*dot_product(coefficients, basis_derivative)
        coefficients_bar = cotangent*basis
    end subroutine analytical_vjp

    subroutine autodiff_vjp(x, coefficients, cotangent, xbar, coefficients_bar)
        real(c_double), intent(in) :: x, coefficients(coefficient_count)
        real(c_double), intent(in) :: cotangent
        real(c_double), intent(out) :: xbar
        real(c_double), intent(out) :: coefficients_bar(coefficient_count)
        real(c_double) :: ignored_value

        xbar = 0.0_c_double
        coefficients_bar = 0.0_c_double
        ignored_value = enzyme_autodiff(c_funloc(spline_value), x, xbar, coefficients, &
            coefficients_bar)
        xbar = cotangent*xbar
        coefficients_bar = cotangent*coefficients_bar
    end subroutine autodiff_vjp

end module bspline_fixed_span_kernel

program enzyme_bspline_fixed_span_products
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: int64
    use bspline_fixed_span_kernel, only: analytical_jvp, analytical_vjp, &
        autodiff_jvp, autodiff_vjp, coefficient_count, spline_value
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(c_double) :: coefficients(coefficient_count), directions(coefficient_count, 16)
    real(c_double) :: x_directions(16), cotangents(16)
    integer :: direction_count
    integer(int64) :: iterations
    character(32) :: action, candidate, product

    call initialize_inputs()
    call get_environment_variable("FORTNUM_BSPLINE_SPAN_ACTION", action)
    call get_environment_variable("FORTNUM_BSPLINE_SPAN_CANDIDATE", candidate)
    call get_environment_variable("FORTNUM_BSPLINE_SPAN_PRODUCT", product)
    call read_integer_env("FORTNUM_BSPLINE_SPAN_DIRECTIONS", direction_count, 16)
    call read_int64_env("FORTNUM_BSPLINE_SPAN_ITERATIONS", iterations, 100000_int64)
    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), trim(product), direction_count, &
            iterations)
    else if (trim(action) == "--peak-rss") then
        call benchmark_candidate(trim(candidate), trim(product), direction_count, &
            iterations)
        write (*, "(a,i0)") "peak_rss_bytes=", peak_rss_bytes()
    else
        call validate_candidates()
        write (*, "(a)") "PASS"
    end if

contains

    subroutine initialize_inputs()
        integer :: i, j

        coefficients = [0.3_c_double, -0.7_c_double, 1.2_c_double, 0.4_c_double]
        do j = 1, 16
            x_directions(j) = 0.03_c_double*real(mod(5*j, 11) - 5, c_double)
            cotangents(j) = 0.1_c_double*real(mod(7*j, 13) - 6, c_double)
            do i = 1, coefficient_count
                directions(i, j) = 0.02_c_double*real( &
                    mod(3*i + 7*j, 17) - 8, c_double)
            end do
        end do
    end subroutine initialize_inputs

    subroutine validate_candidates()
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double) :: analytical, automatic, finite_difference
        real(c_double) :: plus_coefficients(coefficient_count)
        real(c_double) :: minus_coefficients(coefficient_count)
        real(c_double) :: xbar_analytical, xbar_autodiff
        real(c_double) :: coefficient_bar_analytical(coefficient_count)
        real(c_double) :: coefficient_bar_autodiff(coefficient_count)
        real(c_double) :: lhs, rhs, x
        integer :: direction

        x = 0.37_c_double
        do direction = 1, 16
            analytical = analytical_jvp(x, coefficients, x_directions(direction), &
                directions(:, direction))
            automatic = autodiff_jvp(x, coefficients, x_directions(direction), &
                directions(:, direction))
            plus_coefficients = coefficients + h*directions(:, direction)
            minus_coefficients = coefficients - h*directions(:, direction)
            finite_difference = (spline_value(x + h*x_directions(direction), &
                plus_coefficients) - spline_value(x - h*x_directions(direction), &
                minus_coefficients))/(2.0_c_double*h)
            if (abs(analytical - finite_difference) > 2.0e-10_c_double) then
                error stop "analytical B-spline JVP mismatch"
            end if
            if (abs(automatic - finite_difference) > 2.0e-10_c_double) then
                error stop "autodiff B-spline JVP mismatch"
            end if
            call analytical_vjp(x, coefficients, cotangents(direction), &
                xbar_analytical, coefficient_bar_analytical)
            call autodiff_vjp(x, coefficients, cotangents(direction), &
                xbar_autodiff, coefficient_bar_autodiff)
            if (abs(xbar_analytical - xbar_autodiff) > 2.0e-12_c_double .or. &
                maxval(abs(coefficient_bar_analytical &
                - coefficient_bar_autodiff)) > 2.0e-12_c_double) then
                error stop "autodiff B-spline VJP mismatch"
            end if
            lhs = cotangents(direction)*analytical
            rhs = x_directions(direction)*xbar_analytical &
                + dot_product(directions(:, direction), &
                coefficient_bar_analytical)
            if (abs(lhs - rhs) > 2.0e-12_c_double) then
                error stop "B-spline adjoint identity mismatch"
            end if
        end do
    end subroutine validate_candidates

    subroutine benchmark_candidate(name, derivative_product, active_count, count)
        character(*), intent(in) :: name, derivative_product
        integer, intent(in) :: active_count
        integer(int64), intent(in) :: count
        real(c_double) :: derivative, xbar, coefficient_bar(coefficient_count)
        real(c_double) :: elapsed_ns, sink, x
        integer :: direction
        integer(int64) :: iteration, start, finish, rate

        if (active_count < 1 .or. active_count > 16) then
            error stop "directions must be 1..16"
        end if
        if (name /= "analytical" .and. name /= "autodiff") then
            error stop "candidate must be analytical or autodiff"
        end if
        if (derivative_product /= "jvp" .and. derivative_product /= "vjp") then
            error stop "product must be jvp or vjp"
        end if
        sink = 0.0_c_double
        call system_clock(start, rate)
        do iteration = 1, count
            x = 0.37_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024_int64), c_double)
            do direction = 1, active_count
                if (derivative_product == "jvp") then
                    if (name == "analytical") then
                        derivative = analytical_jvp(x, coefficients, &
                            x_directions(direction), directions(:, direction))
                    else
                        derivative = autodiff_jvp(x, coefficients, &
                            x_directions(direction), directions(:, direction))
                    end if
                    sink = sink + derivative
                else
                    if (name == "analytical") then
                        call analytical_vjp(x, coefficients, &
                            cotangents(direction), xbar, coefficient_bar)
                    else
                        call autodiff_vjp(x, coefficients, cotangents(direction), &
                            xbar, coefficient_bar)
                    end if
                    sink = sink + xbar + coefficient_bar(1)
                end if
            end do
        end do
        call system_clock(finish)
        if (sink /= sink) error stop "benchmark produced NaN"
        elapsed_ns = real(finish - start, c_double)*1.0e9_c_double &
            /(real(rate, c_double)*real(count, c_double))
        write (*, "(a,a,a,a,a,i0,a,i0,a,f0.6,a,es12.4)") "candidate=", name, &
            " product=", derivative_product, " directions=", active_count, &
            " iterations=", count, " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

    subroutine read_integer_env(name, value, default_value)
        character(*), intent(in) :: name
        integer, intent(out) :: value
        integer, intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid integer environment value"
    end subroutine read_integer_env

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

end program enzyme_bspline_fixed_span_products
