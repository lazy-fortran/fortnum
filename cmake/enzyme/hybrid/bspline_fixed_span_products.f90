module bspline_fixed_span_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_bspline_fixed_span, only: &
        fortnum_enzyme_bspline_fixed_span_jvp, &
        fortnum_enzyme_bspline_fixed_span_vjp
    implicit none
    private

    integer, parameter, public :: coefficient_count = 4
    public :: analytical_jvp, analytical_vjp, autodiff_jvp, autodiff_vjp
    public :: spline_value

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

        derivative = fortnum_enzyme_bspline_fixed_span_jvp( &
            x, coefficients, dx, dcoefficients)
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

        call fortnum_enzyme_bspline_fixed_span_vjp( &
            x, coefficients, cotangent, xbar, coefficients_bar)
    end subroutine autodiff_vjp

end module bspline_fixed_span_kernel

program enzyme_bspline_fixed_span_products
    use, intrinsic :: iso_c_binding, only: c_double
    use bspline_fixed_span_kernel, only: analytical_jvp, analytical_vjp, &
        autodiff_jvp, autodiff_vjp, coefficient_count, spline_value
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none

    real(c_double) :: coefficients(coefficient_count)
    real(c_double) :: directions(coefficient_count, 16)
    real(c_double) :: x_directions(16), cotangents(16)
    real(c_double) :: samples(fixture_sample_count), sink, validation_error
    character(32) :: action, candidate, product
    integer :: candidate_kind, product_kind, direction_count, iterations
    logical :: valid

    call initialize_inputs()
    call read_fixture_environment("FORTNUM_BSPLINE_SPAN_ACTION", &
        "validate", action)
    if (trim(action) == "validate") then
        call validate_candidates(validation_error)
        write (*, "('PASS max_absolute_error=',es12.5)") validation_error
        stop
    end if

    call read_fixture_environment("FORTNUM_BSPLINE_SPAN_CANDIDATE", &
        "analytical", candidate)
    call read_fixture_environment("FORTNUM_BSPLINE_SPAN_PRODUCT", &
        "jvp", product)
    call read_fixture_integer("FORTNUM_BSPLINE_SPAN_DIRECTIONS", 16, &
        direction_count, valid)
    if (.not. valid) error stop "invalid B-spline direction count"
    call read_fixture_integer("FORTNUM_BSPLINE_SPAN_ITERATIONS", 100000, &
        iterations, valid)
    if (.not. valid) error stop "invalid B-spline iteration count"
    call parse_configuration()

    if (trim(action) == "benchmark" .or. trim(action) == "--benchmark") then
        call benchmark_candidate()
    else if (trim(action) == "peak-rss" .or. trim(action) == "--peak-rss") then
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    else
        error stop "action must be validate, benchmark, or peak-rss"
    end if
    if (sink /= sink) error stop "B-spline benchmark produced NaN"

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

    subroutine validate_candidates(max_absolute_error)
        real(c_double), intent(out) :: max_absolute_error
        real(c_double), parameter :: h = 1.0e-5_c_double
        real(c_double) :: analytical, automatic, finite_difference
        real(c_double) :: plus_coefficients(coefficient_count)
        real(c_double) :: minus_coefficients(coefficient_count)
        real(c_double) :: xbar_analytical, xbar_autodiff
        real(c_double) :: coefficient_bar_analytical(coefficient_count)
        real(c_double) :: coefficient_bar_autodiff(coefficient_count)
        real(c_double) :: lhs, rhs, x
        integer :: direction

        max_absolute_error = 0.0_c_double
        x = 0.37_c_double
        do direction = 1, 16
            analytical = analytical_jvp(x, coefficients, &
                x_directions(direction), directions(:, direction))
            automatic = autodiff_jvp(x, coefficients, &
                x_directions(direction), directions(:, direction))
            plus_coefficients = coefficients + h*directions(:, direction)
            minus_coefficients = coefficients - h*directions(:, direction)
            finite_difference = (spline_value(x + h*x_directions(direction), &
                plus_coefficients) - spline_value(x - h*x_directions(direction), &
                minus_coefficients))/(2.0_c_double*h)
            max_absolute_error = max(max_absolute_error, &
                abs(analytical - finite_difference), &
                abs(automatic - finite_difference))
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
            max_absolute_error = max(max_absolute_error, &
                abs(xbar_analytical - xbar_autodiff), &
                maxval(abs(coefficient_bar_analytical - &
                coefficient_bar_autodiff)))
            if (abs(xbar_analytical - xbar_autodiff) > 2.0e-12_c_double .or. &
                maxval(abs(coefficient_bar_analytical - &
                coefficient_bar_autodiff)) > 2.0e-12_c_double) then
                error stop "autodiff B-spline VJP mismatch"
            end if
            lhs = cotangents(direction)*analytical
            rhs = x_directions(direction)*xbar_analytical &
                + dot_product(directions(:, direction), &
                coefficient_bar_analytical)
            max_absolute_error = max(max_absolute_error, abs(lhs - rhs))
            if (abs(lhs - rhs) > 2.0e-12_c_double) then
                error stop "B-spline adjoint identity mismatch"
            end if
        end do
    end subroutine validate_candidates

    subroutine parse_configuration()
        if (direction_count < 1 .or. direction_count > 16) then
            error stop "directions must be 1..16"
        end if
        if (iterations < 1) error stop "iterations must be positive"
        select case (trim(candidate))
        case ("analytical")
            candidate_kind = 1
        case ("autodiff")
            candidate_kind = 2
        case default
            error stop "candidate must be analytical or autodiff"
        end select
        select case (trim(product))
        case ("jvp")
            product_kind = 1
        case ("vjp")
            product_kind = 2
        case default
            error stop "product must be jvp or vjp"
        end select
    end subroutine parse_configuration

    subroutine benchmark_candidate()
        real(c_double) :: median, mad
        character(96) :: name

        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        write (name, "(a,'_',a,'_d',i0)") trim(candidate), trim(product), &
            direction_count
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine benchmark_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: coefficient_bar(coefficient_count)
        real(c_double) :: derivative, elapsed_ns, local_sink, x, xbar
        integer :: direction, iteration

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            x = 0.37_c_double + 1.0e-12_c_double*real( &
                mod(iteration, 1024), c_double)
            do direction = 1, direction_count
                if (product_kind == 1) then
                    if (candidate_kind == 1) then
                        derivative = analytical_jvp(x, coefficients, &
                            x_directions(direction), directions(:, direction))
                    else
                        derivative = autodiff_jvp(x, coefficients, &
                            x_directions(direction), directions(:, direction))
                    end if
                    local_sink = local_sink + derivative
                else
                    if (candidate_kind == 1) then
                        call analytical_vjp(x, coefficients, &
                            cotangents(direction), xbar, coefficient_bar)
                    else
                        call autodiff_vjp(x, coefficients, &
                            cotangents(direction), xbar, coefficient_bar)
                    end if
                    local_sink = local_sink + xbar + coefficient_bar(1)
                end if
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

end program enzyme_bspline_fixed_span_products
