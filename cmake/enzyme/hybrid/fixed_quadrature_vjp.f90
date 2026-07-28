module fixed_quadrature_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_fixed_quadrature_integrand, only: &
        fortnum_enzyme_fixed_quadrature_integrand_vjp
    use fortnum_generated_enzyme_fixed_quadrature_kernel, only: &
        fortnum_enzyme_fixed_quadrature_kernel_vjp
    use fortnum_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_vjp
    implicit none
    private

    integer, parameter :: kernel_rule_size = 32
    real(dp), save :: kernel_nodes(kernel_rule_size)
    real(dp), save :: kernel_weights(kernel_rule_size)

    public :: analytical_integral_vjp, autodiff_integral_vjp, hybrid_integral_vjp
    public :: diagnostic_integral_vjp, exact_integral_vjp
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

    pure function analytical_integral_vjp(nodes, weights, parameters, &
            cotangent) result(vjp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), cotangent
        real(dp) :: vjp(4), node_cotangents(size(nodes)), output_cotangent(1)
        integer :: i

        output_cotangent(1) = cotangent
        call gauss_legendre_vjp(weights, output_cotangent, node_cotangents)
        vjp = 0.0_dp
        do i = 1, size(nodes)
            vjp(1) = vjp(1) + node_cotangents(i)*nodes(i)* &
                exp(parameters(1)*nodes(i))
            vjp(2) = vjp(2) + node_cotangents(i)*nodes(i)* &
                cos(parameters(2)*nodes(i))
            vjp(3) = vjp(3) + node_cotangents(i)*nodes(i)*nodes(i)
            vjp(4) = vjp(4) + node_cotangents(i)*nodes(i)*nodes(i)*nodes(i)
        end do
    end function analytical_integral_vjp

    function autodiff_integral_vjp(parameters, cotangent) result(vjp)
        real(dp), intent(in) :: parameters(4), cotangent
        real(dp) :: vjp(4)

        call fortnum_enzyme_fixed_quadrature_kernel_vjp( &
            parameters(1), parameters(2), parameters(3), parameters(4), &
            cotangent, vjp(1), vjp(2), vjp(3), vjp(4))
    end function autodiff_integral_vjp

    function hybrid_integral_vjp(nodes, weights, parameters, cotangent) &
            result(vjp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), cotangent
        real(dp) :: vjp(4), node_cotangents(size(nodes)), output_cotangent(1)
        real(dp) :: ignored_xbar, parameter_bars(4)
        integer :: i

        output_cotangent(1) = cotangent
        call gauss_legendre_vjp(weights, output_cotangent, node_cotangents)
        vjp = 0.0_dp
        do i = 1, size(nodes)
            call fortnum_enzyme_fixed_quadrature_integrand_vjp( &
                nodes(i), parameters(1), parameters(2), parameters(3), &
                parameters(4), node_cotangents(i), ignored_xbar, &
                parameter_bars(1), parameter_bars(2), parameter_bars(3), &
                parameter_bars(4))
            vjp = vjp + parameter_bars
        end do
    end function hybrid_integral_vjp

    pure function diagnostic_integral_vjp(nodes, weights, parameters, &
            cotangent) result(vjp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), cotangent
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: vjp(4), parameters_plus(4), parameters_minus(4)
        integer :: parameter

        do parameter = 1, 4
            parameters_plus = parameters
            parameters_minus = parameters
            parameters_plus(parameter) = parameters_plus(parameter) + h
            parameters_minus(parameter) = parameters_minus(parameter) - h
            vjp(parameter) = cotangent*( &
                quadrature_value(nodes, weights, parameters_plus) - &
                quadrature_value(nodes, weights, parameters_minus))/(2.0_dp*h)
        end do
    end function diagnostic_integral_vjp

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

    pure function exact_integral_vjp(parameters, cotangent) result(vjp)
        real(dp), intent(in) :: parameters(4), cotangent
        real(dp) :: vjp(4), p1, p2

        p1 = parameters(1)
        p2 = parameters(2)
        vjp(1) = cotangent*(exp(p1)*(p1 - 1.0_dp) + 1.0_dp)/(p1*p1)
        vjp(2) = cotangent*(p2*sin(p2) + cos(p2) - 1.0_dp)/(p2*p2)
        vjp(3) = cotangent/3.0_dp
        vjp(4) = cotangent/4.0_dp
    end function exact_integral_vjp

end module fixed_quadrature_vjp_hybrid_kernel

program enzyme_fixed_quadrature_vjp_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fixed_quadrature_vjp_hybrid_kernel, only: analytical_integral_vjp, &
        autodiff_integral_vjp, hybrid_integral_vjp, diagnostic_integral_vjp, &
        exact_integral_vjp, configure_quadrature_kernel
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_scaling_result
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: rule_size = 32
    integer, parameter :: max_cotangents = 4
    real(dp) :: nodes(rule_size), weights(rule_size), parameters(4)
    real(dp) :: cotangents(max_cotangents)
    real(dp) :: samples(fixture_sample_count), sink
    character(32) :: argument, candidate, count_argument
    integer :: cotangent_count, candidate_kind, measurement_count
    integer :: repetitions

    call initialize_workload()
    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    call get_command_argument(3, count_argument)
    if (len_trim(argument) == 0) then
        call get_environment_variable("FORTNUM_BATCH_ACTION", argument)
        call get_environment_variable("FORTNUM_BATCH_CANDIDATE", candidate)
        call get_environment_variable("FORTNUM_BATCH_SIZE", count_argument)
    end if
    if (trim(argument) == "--benchmark") then
        if (len_trim(candidate) > 0) then
            read (count_argument, *) cotangent_count
            call run_single_benchmark(trim(candidate), cotangent_count)
        else
            call run_benchmark()
        end if
    else if (trim(argument) == "--batch-benchmark") then
        read (count_argument, *) cotangent_count
        call run_batch_benchmark(trim(candidate), cotangent_count)
    else if (trim(argument) == "--batch-peak-rss") then
        read (count_argument, *) cotangent_count
        call run_batch_peak_rss(trim(candidate), cotangent_count)
    else if (trim(argument) == "--peak-rss") then
        read (count_argument, *) cotangent_count
        call run_peak_rss(trim(candidate), cotangent_count)
    else
        call validate_candidates()
        call validate_batch_candidates()
    end if

contains

    subroutine initialize_workload()
        call gauss_legendre_ab(rule_size, 0.0_dp, 1.0_dp, nodes, weights)
        call configure_quadrature_kernel(nodes, weights)
        parameters = [0.7_dp, 1.1_dp, -0.3_dp, 0.2_dp]
        cotangents = [1.3_dp, -0.4_dp, 0.7_dp, -1.1_dp]
    end subroutine initialize_workload

    subroutine validate_candidates()
        real(dp) :: analytical(4), autodiff(4), hybrid(4), diagnostic(4)
        real(dp) :: reference(4), errors(4), maximum_error
        integer :: cotangent

        maximum_error = 0.0_dp
        do cotangent = 1, max_cotangents
            reference = exact_integral_vjp(parameters, cotangents(cotangent))
            analytical = analytical_integral_vjp(nodes, weights, parameters, &
                cotangents(cotangent))
            autodiff = autodiff_integral_vjp(parameters, cotangents(cotangent))
            hybrid = hybrid_integral_vjp(nodes, weights, parameters, &
                cotangents(cotangent))
            diagnostic = diagnostic_integral_vjp(nodes, weights, parameters, &
                cotangents(cotangent))
            errors(1) = maxval(abs(analytical - reference))
            errors(2) = maxval(abs(autodiff - reference))
            errors(3) = maxval(abs(hybrid - reference))
            errors(4) = maxval(abs(diagnostic - reference))
            maximum_error = max(maximum_error, maxval(errors))
            if (maxval(errors(1:3)) > 2.0e-14_dp .or. &
                errors(4) > 2.0e-10_dp) then
                print *, "fixed-quadrature VJP mismatch", cotangent, errors
                error stop 1
            end if
        end do
        print *, "PASS fixed-quadrature hybrid VJP max_absolute_error", &
            maximum_error
    end subroutine validate_candidates

    subroutine validate_batch_candidates()
        real(dp) :: analytical(4), autodiff(4), hybrid(4), diagnostic(4)
        real(dp) :: reference(4), batch_parameters(4), batch_errors(4)
        real(dp) :: maximum_error
        integer :: batch

        maximum_error = 0.0_dp
        do batch = 1, 16
            batch_parameters = parameters
            batch_parameters(1) = batch_parameters(1) + &
                1.0e-3_dp*real(batch - 1, dp)
            reference = exact_integral_vjp(batch_parameters, cotangents(1))
            analytical = analytical_integral_vjp(nodes, weights, &
                batch_parameters, cotangents(1))
            autodiff = autodiff_integral_vjp(batch_parameters, cotangents(1))
            hybrid = hybrid_integral_vjp(nodes, weights, batch_parameters, &
                cotangents(1))
            diagnostic = diagnostic_integral_vjp(nodes, weights, &
                batch_parameters, cotangents(1))
            batch_errors(1) = maxval(abs(analytical - reference))
            batch_errors(2) = maxval(abs(autodiff - reference))
            batch_errors(3) = maxval(abs(hybrid - reference))
            batch_errors(4) = maxval(abs(diagnostic - reference))
            maximum_error = max(maximum_error, maxval(batch_errors))
            if (maxval(batch_errors(1:3)) > 2.0e-14_dp .or. &
                batch_errors(4) > 2.0e-10_dp) then
                print *, "batched quadrature VJP mismatch", batch, batch_errors
                error stop 1
            end if
        end do
        print *, "PASS batched quadrature VJP max_absolute_error", maximum_error
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

        call validate_request(name, 1)
        if ((batch_size /= 1) .and. (batch_size /= 4) .and. &
            (batch_size /= 16)) then
            error stop "batch size must be 1, 4, or 16"
        end if
    end subroutine validate_batch_request

    function measure_batch_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        integer :: iteration
        integer :: batch
        real(dp) :: varied_parameters(4), result(4), local_sink

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            do batch = 1, measurement_count
                varied_parameters = parameters
                varied_parameters(1) = varied_parameters(1) + &
                    1.0e-6_dp*real(mod(iteration, 101), dp) + &
                    1.0e-3_dp*real(batch - 1, dp)
                select case (candidate_kind)
                case (1)
                    result = analytical_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                case (2)
                    result = autodiff_integral_vjp(varied_parameters, &
                        cotangents(1))
                case (3)
                    result = hybrid_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                case default
                    result = diagnostic_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                end select
                local_sink = local_sink + sum(result)
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "fixed-quadrature VJP batch benchmark failed"
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
            error stop "cotangent count must be 1, 2, or 4"
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
        integer :: cotangent
        real(dp) :: varied_parameters(4), result(4), local_sink

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            varied_parameters = parameters
            varied_parameters(1) = varied_parameters(1) + &
                1.0e-6_dp*real(mod(iteration, 101), dp)
            do cotangent = 1, measurement_count
                select case (candidate_kind)
                case (1)
                    result = analytical_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                case (2)
                    result = autodiff_integral_vjp(varied_parameters, &
                        cotangents(cotangent))
                case (3)
                    result = hybrid_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                case default
                    result = diagnostic_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                end select
                local_sink = local_sink + sum(result)
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "fixed-quadrature VJP benchmark failed"
        end if
        sink = local_sink
    end function measure_candidate

end program enzyme_fixed_quadrature_vjp_hybrid
