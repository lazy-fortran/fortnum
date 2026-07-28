module fixed_quadrature_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_fixed_quadrature_integrand, only: &
        fortnum_enzyme_fixed_quadrature_integrand_jvp
    use fortnum_generated_enzyme_fixed_quadrature_kernel, only: &
        fortnum_enzyme_fixed_quadrature_kernel_jvp
    use fortnum_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_jvp
    implicit none
    private

    integer, parameter :: kernel_rule_size = 32
    real(dp), save :: kernel_nodes(kernel_rule_size)
    real(dp), save :: kernel_weights(kernel_rule_size)

    public :: analytical_integral_jvp, autodiff_integral_jvp, hybrid_integral_jvp
    public :: diagnostic_integral_jvp, exact_integral_jvp
    public :: configure_quadrature_kernel

contains

    pure function integrand(x, p1, p2, p3, p4) result(value) &
            bind(c, name="fortnum_fixed_quadrature_integrand")
        real(c_double), value :: x, p1, p2, p3, p4
        real(c_double) :: value

        value = exp(p1*x) + sin(p2*x) + p3*x*x + p4*x*x*x
    end function integrand

    subroutine configure_quadrature_kernel(nodes, weights)
        real(dp), intent(in) :: nodes(kernel_rule_size), weights(kernel_rule_size)

        kernel_nodes = nodes
        kernel_weights = weights
    end subroutine configure_quadrature_kernel

    function quadrature_kernel(p1, p2, p3, p4) result(value) &
            bind(c, name="fortnum_fixed_quadrature_kernel")
        real(c_double), value :: p1, p2, p3, p4
        real(c_double) :: value
        integer :: i

        value = 0.0_dp
        do i = 1, kernel_rule_size
            value = value + kernel_weights(i)*integrand(kernel_nodes(i), &
                p1, p2, p3, p4)
        end do
    end function quadrature_kernel

    function analytical_integral_jvp(nodes, weights, parameters, direction) &
            result(jvp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), direction(4)
        real(dp) :: jvp, tangent_values(size(nodes)), contracted(1)

        tangent_values = direction(1)*nodes*exp(parameters(1)*nodes) + &
            direction(2)*nodes*cos(parameters(2)*nodes) + &
            direction(3)*nodes*nodes + direction(4)*nodes*nodes*nodes
        call gauss_legendre_jvp(weights, tangent_values, contracted)
        jvp = contracted(1)
    end function analytical_integral_jvp

    function autodiff_integral_jvp(parameters, direction) result(jvp)
        real(dp), intent(in) :: parameters(4), direction(4)
        real(dp) :: jvp

        jvp = fortnum_enzyme_fixed_quadrature_kernel_jvp( &
            parameters(1), direction(1), parameters(2), direction(2), &
            parameters(3), direction(3), parameters(4), direction(4))
    end function autodiff_integral_jvp

    function hybrid_integral_jvp(nodes, weights, parameters, direction) &
            result(jvp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), direction(4)
        real(dp) :: jvp, tangent_values(size(nodes)), contracted(1)
        integer :: i

        do i = 1, size(nodes)
            tangent_values(i) = fortnum_enzyme_fixed_quadrature_integrand_jvp( &
                nodes(i), 0.0_dp, parameters(1), direction(1), &
                parameters(2), direction(2), parameters(3), direction(3), &
                parameters(4), direction(4))
        end do
        call gauss_legendre_jvp(weights, tangent_values, contracted)
        jvp = contracted(1)
    end function hybrid_integral_jvp

    pure function diagnostic_integral_jvp(nodes, weights, parameters, &
            direction) result(jvp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), direction(4)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: jvp, parameters_plus(4), parameters_minus(4)

        parameters_plus = parameters + h*direction
        parameters_minus = parameters - h*direction
        jvp = (quadrature_value(nodes, weights, parameters_plus) - &
            quadrature_value(nodes, weights, parameters_minus))/(2.0_dp*h)
    end function diagnostic_integral_jvp

    pure function quadrature_value(nodes, weights, parameters) result(value)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4)
        real(dp) :: value
        integer :: i

        value = 0.0_dp
        do i = 1, size(nodes)
            value = value + weights(i)*integrand(nodes(i), parameters(1), &
                parameters(2), parameters(3), parameters(4))
        end do
    end function quadrature_value

    pure function exact_integral_jvp(parameters, direction) result(jvp)
        real(dp), intent(in) :: parameters(4), direction(4)
        real(dp) :: jvp, p1, p2

        p1 = parameters(1)
        p2 = parameters(2)
        jvp = direction(1)*(exp(p1)*(p1 - 1.0_dp) + 1.0_dp)/(p1*p1) + &
            direction(2)*(p2*sin(p2) + cos(p2) - 1.0_dp)/(p2*p2) + &
            direction(3)/3.0_dp + direction(4)/4.0_dp
    end function exact_integral_jvp

end module fixed_quadrature_hybrid_kernel

program enzyme_fixed_quadrature_jvp_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_quadrature_hybrid_kernel, only: analytical_integral_jvp, &
        autodiff_integral_jvp, hybrid_integral_jvp, diagnostic_integral_jvp, &
        exact_integral_jvp, configure_quadrature_kernel
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_scaling_result
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: rule_size = 32
    integer, parameter :: max_directions = 4
    real(dp) :: nodes(rule_size), weights(rule_size)
    real(dp) :: parameters(4), directions(4, max_directions), errors(4)
    real(dp) :: samples(fixture_sample_count), sink
    character(32) :: argument, candidate, direction_argument
    integer :: direction_count, candidate_kind, measurement_count
    integer :: repetitions

    call initialize_workload()
    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    call get_command_argument(3, direction_argument)
    if (len_trim(argument) == 0) then
        call get_environment_variable("FORTNUM_BATCH_ACTION", argument)
        call get_environment_variable("FORTNUM_BATCH_CANDIDATE", candidate)
        call get_environment_variable("FORTNUM_BATCH_SIZE", direction_argument)
    end if

    if (trim(argument) == "--benchmark") then
        if (len_trim(candidate) > 0) then
            read (direction_argument, *) direction_count
            call run_single_benchmark(trim(candidate), direction_count)
        else
            call run_benchmark()
        end if
    else if (trim(argument) == "--batch-benchmark") then
        read (direction_argument, *) direction_count
        call run_batch_benchmark(trim(candidate), direction_count)
    else if (trim(argument) == "--batch-peak-rss") then
        read (direction_argument, *) direction_count
        call run_batch_peak_rss(trim(candidate), direction_count)
    else if (trim(argument) == "--peak-rss") then
        read (direction_argument, *) direction_count
        call run_peak_rss(trim(candidate), direction_count)
    else
        call validate_candidates()
        call validate_batch_candidates()
    end if

contains

    subroutine initialize_workload()
        call gauss_legendre_ab(rule_size, 0.0_dp, 1.0_dp, nodes, weights)
        call configure_quadrature_kernel(nodes, weights)
        parameters = [0.7_dp, 1.1_dp, -0.3_dp, 0.2_dp]
        directions(:, 1) = [0.4_dp, -0.6_dp, 0.2_dp, 0.7_dp]
        directions(:, 2) = [-0.3_dp, 0.1_dp, 0.8_dp, -0.2_dp]
        directions(:, 3) = [0.5_dp, 0.4_dp, -0.1_dp, 0.3_dp]
        directions(:, 4) = [-0.2_dp, 0.7_dp, 0.6_dp, -0.5_dp]
    end subroutine initialize_workload

    subroutine validate_candidates()
        real(dp) :: analytical, autodiff, hybrid, diagnostic, reference
        real(dp) :: maximum_error
        integer :: direction

        maximum_error = 0.0_dp
        do direction = 1, max_directions
            reference = exact_integral_jvp(parameters, directions(:, direction))
            analytical = analytical_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            autodiff = autodiff_integral_jvp(parameters, directions(:, direction))
            hybrid = hybrid_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            diagnostic = diagnostic_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            errors(1) = abs(analytical - reference)
            errors(2) = abs(autodiff - reference)
            errors(3) = abs(hybrid - reference)
            errors(4) = abs(diagnostic - reference)
            maximum_error = max(maximum_error, maxval(errors))
            if (maxval(errors(1:3)) > 2.0e-14_dp .or. &
                errors(4) > 2.0e-10_dp) then
                print *, "fixed-quadrature JVP mismatch", direction, errors
                error stop 1
            end if
        end do
        print *, "PASS fixed-quadrature hybrid JVP max_absolute_error", &
            maximum_error
    end subroutine validate_candidates

    subroutine validate_batch_candidates()
        real(dp) :: analytical, autodiff, hybrid, diagnostic, reference
        real(dp) :: batch_parameters(4), batch_errors(4), maximum_error
        integer :: batch, direction

        maximum_error = 0.0_dp
        do batch = 1, 16
            batch_parameters = parameters
            batch_parameters(1) = batch_parameters(1) + &
                1.0e-3_dp*real(batch - 1, dp)
            do direction = 1, 4
                reference = exact_integral_jvp(batch_parameters, &
                    directions(:, direction))
                analytical = analytical_integral_jvp(nodes, weights, &
                    batch_parameters, directions(:, direction))
                autodiff = autodiff_integral_jvp(batch_parameters, &
                    directions(:, direction))
                hybrid = hybrid_integral_jvp(nodes, weights, batch_parameters, &
                    directions(:, direction))
                diagnostic = diagnostic_integral_jvp(nodes, weights, &
                    batch_parameters, directions(:, direction))
                batch_errors(1) = abs(analytical - reference)
                batch_errors(2) = abs(autodiff - reference)
                batch_errors(3) = abs(hybrid - reference)
                batch_errors(4) = abs(diagnostic - reference)
                maximum_error = max(maximum_error, maxval(batch_errors))
                if (maxval(batch_errors(1:3)) > 2.0e-14_dp .or. &
                    batch_errors(4) > 2.0e-10_dp) then
                    print *, "batched quadrature JVP mismatch", batch, &
                        direction, batch_errors
                    error stop 1
                end if
            end do
        end do
        print *, "PASS batched quadrature JVP max_absolute_error", maximum_error
    end subroutine validate_batch_candidates

    subroutine run_benchmark()
        integer, parameter :: counts(3) = [1, 2, 4]
        integer :: count_index

        do count_index = 1, size(counts)
            call run_single_benchmark("analytical", counts(count_index))
            call run_single_benchmark("autodiff", counts(count_index))
            call run_single_benchmark("hybrid", counts(count_index))
            call run_single_benchmark("diagnostic", counts(count_index))
        end do
    end subroutine run_benchmark

    subroutine run_single_benchmark(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        real(dp) :: median, mad

        call validate_request(name, count)
        call select_candidate(name)
        measurement_count = count
        repetitions = 10000
        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_scaling_result( &
            name, count, repetitions, median, mad)
    end subroutine run_single_benchmark

    subroutine run_peak_rss(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count

        call validate_request(name, count)
        call select_candidate(name)
        measurement_count = count
        repetitions = 20000
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine run_batch_benchmark(name, batch_size)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size
        real(dp) :: median, mad

        call validate_batch_request(name, batch_size)
        call select_candidate(name)
        measurement_count = batch_size
        repetitions = 2000
        call collect_fixture_samples(measure_batch_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_scaling_result( &
            name, batch_size, repetitions, median, mad)
    end subroutine run_batch_benchmark

    subroutine run_batch_peak_rss(name, batch_size)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size

        call validate_batch_request(name, batch_size)
        call select_candidate(name)
        measurement_count = batch_size
        repetitions = 5000
        sink = measure_batch_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_batch_peak_rss

    subroutine validate_batch_request(name, batch_size)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size

        call validate_request(name, 4)
        if ((batch_size /= 1) .and. (batch_size /= 4) .and. &
            (batch_size /= 16)) then
            error stop "batch size must be 1, 4, or 16"
        end if
    end subroutine validate_batch_request

    function measure_batch_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        integer :: iteration
        integer :: batch, direction
        real(dp) :: varied_parameters(4), local_sink

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            do batch = 1, measurement_count
                varied_parameters = parameters
                varied_parameters(1) = varied_parameters(1) + &
                    1.0e-6_dp*real(mod(iteration, 101), dp) + &
                    1.0e-3_dp*real(batch - 1, dp)
                do direction = 1, 4
                    select case (candidate_kind)
                    case (1)
                        local_sink = local_sink + analytical_integral_jvp( &
                            nodes, weights, &
                            varied_parameters, directions(:, direction))
                    case (2)
                        local_sink = local_sink + autodiff_integral_jvp( &
                            varied_parameters, directions(:, direction))
                    case (3)
                        local_sink = local_sink + hybrid_integral_jvp( &
                            nodes, weights, &
                            varied_parameters, directions(:, direction))
                    case default
                        local_sink = local_sink + diagnostic_integral_jvp( &
                            nodes, weights, &
                            varied_parameters, directions(:, direction))
                    end select
                end do
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "fixed-quadrature JVP batch benchmark failed"
        end if
        sink = local_sink
    end function measure_batch_candidate

    subroutine validate_request(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count

        if ((name /= "analytical") .and. (name /= "autodiff") .and. &
            (name /= "hybrid") .and. (name /= "diagnostic")) then
            error stop "candidate must be analytical, autodiff, hybrid, or diagnostic"
        end if
        if ((count /= 1) .and. (count /= 2) .and. (count /= 4)) then
            error stop "direction count must be 1, 2, or 4"
        end if
    end subroutine validate_request

    subroutine select_candidate(name)
        character(*), intent(in) :: name

        select case (name)
        case ("analytical")
            candidate_kind = 1
        case ("autodiff")
            candidate_kind = 2
        case ("hybrid")
            candidate_kind = 3
        case default
            candidate_kind = 4
        end select
    end subroutine select_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        integer :: iteration
        integer :: direction
        real(dp) :: varied_parameters(4), local_sink

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            varied_parameters = parameters
            varied_parameters(1) = varied_parameters(1) + &
                1.0e-6_dp*real(mod(iteration, 101), dp)
            do direction = 1, measurement_count
                select case (candidate_kind)
                case (1)
                    local_sink = local_sink + analytical_integral_jvp( &
                        nodes, weights, &
                        varied_parameters, directions(:, direction))
                case (2)
                    local_sink = local_sink + autodiff_integral_jvp( &
                        varied_parameters, directions(:, direction))
                case (3)
                    local_sink = local_sink + hybrid_integral_jvp(nodes, weights, &
                        varied_parameters, directions(:, direction))
                case default
                    local_sink = local_sink + diagnostic_integral_jvp( &
                        nodes, weights, &
                        varied_parameters, directions(:, direction))
                end select
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "fixed-quadrature JVP benchmark failed"
        end if
        sink = local_sink
    end function measure_candidate

end program enzyme_fixed_quadrature_jvp_hybrid
