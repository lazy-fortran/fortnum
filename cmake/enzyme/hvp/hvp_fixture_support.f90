module fortnum_hvp_fixture_support
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none
    private

    integer, parameter, public :: hvp_size = 4

    abstract interface
        subroutine hvp_candidate(x, direction, target, product)
            import c_double, hvp_size
            real(c_double), intent(in) :: x(hvp_size), direction(hvp_size)
            real(c_double), intent(in) :: target(hvp_size)
            real(c_double), intent(out) :: product(hvp_size)
        end subroutine hvp_candidate
    end interface

    public :: run_hvp_fixture

contains

    subroutine run_hvp_fixture(name, candidate)
        character(*), intent(in) :: name
        procedure(hvp_candidate) :: candidate
        real(c_double) :: samples(fixture_sample_count), sink
        character(32) :: action, selected
        integer :: iterations
        logical :: valid

        call read_fixture_environment("FORTNUM_HVP_ACTION", "validate", action)
        call read_fixture_environment("FORTNUM_HVP_CANDIDATE", "autodiff", &
            selected)
        call read_fixture_integer("FORTNUM_HVP_ITERATIONS", 200000, iterations, &
            valid)
        if (.not. valid .or. iterations < 1) error stop "invalid HVP iterations"
        if (trim(action) == "validate") then
            call validate_candidate()
            write (*, "(a)") "PASS "//trim(name)
        else if (trim(action) == "benchmark") then
            call collect_fixture_samples(measure_candidate, samples)
            call report_samples()
        else if (trim(action) == "peak-rss") then
            sink = measure_candidate()
            write (*, "(i0)") fixture_peak_rss_bytes()
        else
            error stop "HVP action must be validate, benchmark, or peak-rss"
        end if

    contains

        subroutine problem_data(x, direction, target)
            real(c_double), intent(out) :: x(hvp_size), direction(hvp_size)
            real(c_double), intent(out) :: target(hvp_size)

            x = [0.3_c_double, -0.7_c_double, 1.1_c_double, 0.5_c_double]
            direction = [0.2_c_double, 0.4_c_double, -0.3_c_double, 0.8_c_double]
            target = [0.1_c_double, -0.2_c_double, 0.7_c_double, 0.4_c_double]
        end subroutine problem_data

        subroutine validate_candidate()
            real(c_double) :: x(hvp_size), direction(hvp_size), target(hvp_size)
            real(c_double) :: product(hvp_size), expected(hvp_size)

            call problem_data(x, direction, target)
            call evaluate(x, direction, target, product)
            expected = (6.0_c_double*x*x - 2.0_c_double*target)*direction
            if (maxval(abs(product - expected)) > 2.0e-12_c_double) then
                error stop "mixed-mode HVP mismatch"
            end if
        end subroutine validate_candidate

        function measure_candidate() result(elapsed_ns)
            type(fixture_timer_t) :: timer
            real(c_double) :: elapsed_ns, local_sink
            real(c_double) :: x(hvp_size), direction(hvp_size), target(hvp_size)
            real(c_double) :: product(hvp_size)
            integer :: iteration

            call problem_data(x, direction, target)
            local_sink = 0.0_c_double
            call timer%start()
            do iteration = 1, iterations
                x(1) = 0.3_c_double + 1.0e-7_c_double* &
                    real(mod(iteration, 101) - 50, c_double)
                call evaluate(x, direction, target, product)
                local_sink = local_sink + product(1)
            end do
            elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
            sink = local_sink
        end function measure_candidate

        subroutine report_samples()
            real(c_double) :: median, mad

            call median_mad(samples, median, mad)
            call write_fixture_result(trim(name)//"_"//trim(selected), &
                iterations, median, mad)
        end subroutine report_samples

        subroutine evaluate(x, direction, target, product)
            real(c_double), intent(in) :: x(hvp_size), direction(hvp_size)
            real(c_double), intent(in) :: target(hvp_size)
            real(c_double), intent(out) :: product(hvp_size)
            real(c_double), parameter :: step = 1.0e-5_c_double
            real(c_double) :: plus(hvp_size), minus(hvp_size)

            select case (trim(selected))
            case ("autodiff")
                call candidate(x, direction, target, product)
            case ("analytical")
                product = (6.0_c_double*x*x - 2.0_c_double*target)*direction
            case ("diagnostic")
                plus = x + step*direction
                minus = x - step*direction
                product = (2.0_c_double*plus*(plus*plus - target) - &
                    2.0_c_double*minus*(minus*minus - target))/(2.0_c_double*step)
            case default
                error stop "HVP candidate must be autodiff, analytical, or diagnostic"
            end select
        end subroutine evaluate

    end subroutine run_hvp_fixture

end module fortnum_hvp_fixture_support
