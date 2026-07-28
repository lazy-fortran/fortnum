module fftw8_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private

    integer, parameter, public :: real_size = 16
    public :: fftw8_outer

    interface
        function fftw8_objective_scalar( &
                x1, x2, x3, x4, x5, x6, x7, x8, &
                x9, x10, x11, x12, x13, x14, x15, x16) result(value) &
                bind(c, name="fortnum_fftw8_objective_scalar")
            import :: c_double
            real(c_double), value :: x1, x2, x3, x4, x5, x6, x7, x8
            real(c_double), value :: x9, x10, x11, x12, x13, x14, x15, x16
            real(c_double) :: value
        end function fftw8_objective_scalar
    end interface

contains

    function fftw8_outer(x) result(value) bind(c, name="fortnum_fftw8_outer")
        real(c_double), intent(in) :: x(real_size)
        real(c_double) :: value, transformed

        transformed = fftw8_objective_scalar( &
            x(1), x(2), x(3), x(4), x(5), x(6), x(7), x(8), &
            x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16))
        value = sin(transformed) + transformed*transformed
    end function fftw8_outer

end module fftw8_outer_kernel

program enzyme_fftw8_custom_rule
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use fftw8_outer_kernel, only: fftw8_outer, real_size
    use fortnum_generated_enzyme_fftw8_outer, only: &
        fortnum_enzyme_fftw8_outer_jvp
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none

    real(c_double), parameter :: h = 1.0e-6_c_double
    real(c_double) :: samples(fixture_sample_count), sink
    character(32) :: action, candidate
    integer :: directions, iterations
    logical :: valid

    interface
        subroutine fftw8_init() bind(c, name="fortnum_fftw8_init")
        end subroutine fftw8_init

        subroutine fftw8_finalize() bind(c, name="fortnum_fftw8_finalize")
        end subroutine fftw8_finalize

        function fftw8_objective(x) result(value) &
                bind(c, name="fortnum_fftw8_objective_array")
            import :: c_double
            real(c_double), intent(in) :: x(*)
            real(c_double) :: value
        end function fftw8_objective

        function fftw8_objective_jvp(x, direction) result(product) &
                bind(c, name="fortnum_fftw8_objective_jvp")
            import :: c_double
            real(c_double), intent(in) :: x(*), direction(*)
            real(c_double) :: product
        end function fftw8_objective_jvp

        subroutine rule_counter_reset() &
                bind(c, name="fortnum_enzyme_rule_counter_reset")
        end subroutine rule_counter_reset

        function rule_counter_calls() result(count) &
                bind(c, name="fortnum_enzyme_rule_counter_calls")
            import :: c_int64_t
            integer(c_int64_t) :: count
        end function rule_counter_calls

        subroutine rule_counter_disable() &
                bind(c, name="fortnum_enzyme_rule_counter_disable")
        end subroutine rule_counter_disable
    end interface

    call fftw8_init()
    call read_fixture_environment("FORTNUM_FFTW_ACTION", "validate", action)
    if (trim(action) == "validate") then
        call validate_candidates()
        call fftw8_finalize()
        write (*, "(a)") "PASS"
        stop
    end if

    call read_fixture_environment( &
        "FORTNUM_FFTW_CANDIDATE", "analytical", candidate)
    call read_fixture_integer( &
        "FORTNUM_FFTW_DIRECTIONS", 16, directions, valid)
    if (.not. valid .or. directions < 1 .or. directions > 16) &
        error stop "directions must be 1..16"
    call read_fixture_integer( &
        "FORTNUM_FFTW_ITERATIONS", 10000, iterations, valid)
    if (.not. valid .or. iterations < 1) &
        error stop "iterations must be positive"
    if (trim(candidate) /= "analytical" .and. trim(candidate) /= "hybrid") &
        error stop "candidate must be analytical or hybrid"
    call rule_counter_disable()

    select case (trim(action))
    case ("benchmark", "--benchmark")
        call collect_fixture_samples(measure_candidate, samples)
        call report_samples()
    case ("peak-rss", "--peak-rss")
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    case default
        error stop "action must be validate, benchmark, or peak-rss"
    end select
    call fftw8_finalize()
    if (sink /= sink) error stop "FFTW benchmark produced NaN"

contains

    function analytical_jvp(x, direction) result(product)
        real(c_double), intent(in) :: x(real_size), direction(real_size)
        real(c_double) :: product, transformed, tangent

        transformed = fftw8_objective(x)
        tangent = fftw8_objective_jvp(x, direction)
        product = (cos(transformed) + 2.0_c_double*transformed)*tangent
    end function analytical_jvp

    subroutine validate_candidates()
        real(c_double), parameter :: tolerance = 4.0e-9_c_double
        real(c_double) :: x(real_size), direction(real_size)
        real(c_double) :: expected, got, objective, reference, scale
        integer :: i

        call fill_inputs(3, 2, x, direction)
        objective = fftw8_objective(x)
        reference = direct_objective(x)
        scale = max(1.0_c_double, abs(reference))
        if (abs(objective - reference) > 5.0e-13_c_double*scale) &
            error stop "FFTW primitive disagrees with direct DFT"

        expected = (fftw8_outer(x + h*direction) - &
            fftw8_outer(x - h*direction))/(2.0_c_double*h)
        got = analytical_jvp(x, direction)
        scale = max(1.0_c_double, abs(expected))
        if (abs(got - expected) > tolerance*scale) &
            error stop "analytical external FFT JVP mismatch"

        call rule_counter_reset()
        got = fortnum_enzyme_fftw8_outer_jvp(x, direction)
        if (abs(got - expected) > tolerance*scale) &
            error stop "hybrid external FFT JVP mismatch"
        if (rule_counter_calls() /= 1_c_int64_t) &
            error stop "external FFT custom rule was not selected"
        call rule_counter_disable()
    end subroutine validate_candidates

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: x(real_size), direction(real_size)
        real(c_double) :: elapsed_ns, local_sink, product
        integer :: direction_index, iteration

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            do direction_index = 1, directions
                call fill_inputs(iteration, direction_index, x, direction)
                if (trim(candidate) == "analytical") then
                    product = analytical_jvp(x, direction)
                else
                    product = fortnum_enzyme_fftw8_outer_jvp(x, direction)
                end if
                local_sink = local_sink + product
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

    subroutine report_samples()
        real(c_double) :: median, mad
        character(96) :: name

        call median_mad(samples, median, mad)
        write (name, "('fftw8_external_',a,'_jvp_d',i0)") &
            trim(candidate), directions
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine report_samples

    pure subroutine fill_inputs(iteration, direction_index, x, direction)
        integer, intent(in) :: iteration, direction_index
        real(c_double), intent(out) :: x(real_size), direction(real_size)
        integer :: i

        do i = 1, real_size
            x(i) = 0.05_c_double*sin(0.013_c_double* &
                real(i + iteration, c_double))
            direction(i) = cos(0.017_c_double* &
                real(i + 3*direction_index, c_double))
        end do
    end subroutine fill_inputs

    pure function direct_objective(x) result(value)
        real(c_double), intent(in) :: x(real_size)
        real(c_double) :: value
        real(c_double), parameter :: pi = acos(-1.0_c_double)
        complex(c_double) :: transformed, z(8)
        real(c_double) :: angle
        integer :: frequency, i, sample

        do i = 1, 8
            z(i) = cmplx(x(2*i - 1), x(2*i), c_double)
        end do
        value = 0.0_c_double
        do frequency = 0, 7
            transformed = cmplx(0.0_c_double, 0.0_c_double, c_double)
            do sample = 0, 7
                angle = -2.0_c_double*pi* &
                    real(frequency*sample, c_double)/8.0_c_double
                transformed = transformed + z(sample + 1)* &
                    cmplx(cos(angle), sin(angle), c_double)
            end do
            value = value + 0.5_c_double*real(transformed, c_double)**2 + &
                0.25_c_double*aimag(transformed)**2
        end do
    end function direct_objective

end program enzyme_fftw8_custom_rule
