module fixed_quadrature_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fixed_quadrature_full_vjp_autodiff, only: full_quadrature_vjp
    use fortnum_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_vjp
    implicit none
    private

    integer, parameter :: kernel_rule_size = 32
    real(dp), save :: kernel_nodes(kernel_rule_size)
    real(dp), save :: kernel_weights(kernel_rule_size)

    type, bind(c) :: integrand_gradient_t
        real(c_double) :: values(5)
    end type integrand_gradient_t

    public :: analytical_integral_vjp, autodiff_integral_vjp, hybrid_integral_vjp
    public :: diagnostic_integral_vjp, exact_integral_vjp
    public :: configure_quadrature_kernel

    interface
        function enzyme_autodiff(f, x, p1, p2, p3, p4) result(gradient) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr, integrand_gradient_t
            type(c_funptr), value :: f
            real(c_double), value :: x, p1, p2, p3, p4
            type(integrand_gradient_t) :: gradient
        end function enzyme_autodiff
    end interface

contains

    pure function integrand(x, p1, p2, p3, p4) result(value) &
            bind(c, name="fortnum_fixed_quadrature_vjp_integrand")
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
            bind(c, name="fortnum_fixed_quadrature_vjp_kernel")
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

        vjp = full_quadrature_vjp(parameters, cotangent)
    end function autodiff_integral_vjp

    function hybrid_integral_vjp(nodes, weights, parameters, cotangent) &
            result(vjp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), cotangent
        real(dp) :: vjp(4), node_cotangents(size(nodes)), output_cotangent(1)
        type(integrand_gradient_t) :: gradient
        integer :: i, parameter

        output_cotangent(1) = cotangent
        call gauss_legendre_vjp(weights, output_cotangent, node_cotangents)
        vjp = 0.0_dp
        do i = 1, size(nodes)
            gradient = enzyme_autodiff(c_funloc(integrand), nodes(i), &
                parameters(1), parameters(2), parameters(3), parameters(4))
            do parameter = 1, 4
                vjp(parameter) = vjp(parameter) + &
                    node_cotangents(i)*gradient%values(parameter + 1)
            end do
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
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fixed_quadrature_vjp_hybrid_kernel, only: analytical_integral_vjp, &
        autodiff_integral_vjp, hybrid_integral_vjp, diagnostic_integral_vjp, &
        exact_integral_vjp, configure_quadrature_kernel
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: rule_size = 32
    integer, parameter :: max_cotangents = 4
    real(dp) :: nodes(rule_size), weights(rule_size), parameters(4)
    real(dp) :: cotangents(max_cotangents)
    character(32) :: argument, candidate, count_argument
    integer :: cotangent_count

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

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
        real(dp) :: reference(4), errors(4)
        integer :: cotangent

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
            if (maxval(errors(1:3)) > 2.0e-14_dp .or. &
                errors(4) > 2.0e-10_dp) then
                print *, "fixed-quadrature VJP mismatch", cotangent, errors
                error stop 1
            end if
        end do
        print *, "PASS fixed-quadrature hybrid VJP"
    end subroutine validate_candidates

    subroutine validate_batch_candidates()
        real(dp) :: analytical(4), autodiff(4), hybrid(4), diagnostic(4)
        real(dp) :: reference(4), batch_parameters(4), batch_errors(4)
        integer :: batch

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
            if (maxval(batch_errors(1:3)) > 2.0e-14_dp .or. &
                batch_errors(4) > 2.0e-10_dp) then
                print *, "batched quadrature VJP mismatch", batch, batch_errors
                error stop 1
            end if
        end do
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
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 10000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_request(name, count)
        do sample = 1, 3
            call time_candidate(name, count, reps/100_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, count, reps, elapsed(sample))
        end do
        call report(name, count, elapsed, reps)
    end subroutine run_single_benchmark

    subroutine run_peak_rss(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), parameter :: reps = 20000_int64
        real(dp) :: elapsed

        call validate_request(name, count)
        call time_candidate(name, count, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine run_batch_benchmark(name, batch_size)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size
        integer, parameter :: samples = 31
        integer(int64), parameter :: reps = 2000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_batch_request(name, batch_size)
        do sample = 1, 3
            call time_batch_candidate(name, batch_size, reps/20_int64, sink)
        end do
        do sample = 1, samples
            call time_batch_candidate(name, batch_size, reps, elapsed(sample))
        end do
        call report(name, batch_size, elapsed, reps)
    end subroutine run_batch_benchmark

    subroutine run_batch_peak_rss(name, batch_size)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size
        integer(int64), parameter :: reps = 5000_int64
        real(dp) :: elapsed

        call validate_batch_request(name, batch_size)
        call time_batch_candidate(name, batch_size, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
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

    subroutine time_batch_candidate(name, batch_size, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: batch_size
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        integer :: batch
        real(dp) :: varied_parameters(4), result(4), sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            do batch = 1, batch_size
                varied_parameters = parameters
                varied_parameters(1) = varied_parameters(1) + &
                    1.0e-6_dp*real(mod(iteration, 101_int64), dp) + &
                    1.0e-3_dp*real(batch - 1, dp)
                select case (name)
                case ("analytical")
                    result = analytical_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                case ("autodiff")
                    result = autodiff_integral_vjp(varied_parameters, &
                        cotangents(1))
                case ("hybrid")
                    result = hybrid_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                case ("diagnostic")
                    result = diagnostic_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(1))
                end select
                sink = sink + sum(result)
            end do
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
        if (sink == huge(sink)) print *, sink
    end subroutine time_batch_candidate

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

    subroutine time_candidate(name, count, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        integer :: cotangent
        real(dp) :: varied_parameters(4), result(4), sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            varied_parameters = parameters
            varied_parameters(1) = varied_parameters(1) + &
                1.0e-6_dp*real(mod(iteration, 101_int64), dp)
            do cotangent = 1, count
                select case (name)
                case ("analytical")
                    result = analytical_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                case ("autodiff")
                    result = autodiff_integral_vjp(varied_parameters, &
                        cotangents(cotangent))
                case ("hybrid")
                    result = hybrid_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                case ("diagnostic")
                    result = diagnostic_integral_vjp(nodes, weights, &
                        varied_parameters, cotangents(cotangent))
                end select
                sink = sink + sum(result)
            end do
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
        if (sink == huge(sink)) print *, sink
    end subroutine time_candidate

    subroutine report(name, count, values, reps)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        real(dp), intent(in) :: values(:)
        integer(int64), intent(in) :: reps
        real(dp) :: ordered(size(values)), deviations(size(values))
        real(dp) :: median, mad

        ordered = values
        call sort_values(ordered)
        median = ordered((size(ordered) + 1)/2)
        deviations = abs(values - median)
        call sort_values(deviations)
        mad = deviations((size(deviations) + 1)/2)
        write (*, "(a,',',i0,',',i0,',',f12.4,',',f12.4)") &
            name, count, reps, median, mad
    end subroutine report

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: temporary
        integer :: i, j

        do i = 1, size(values) - 1
            do j = i + 1, size(values)
                if (values(j) < values(i)) then
                    temporary = values(i)
                    values(i) = values(j)
                    values(j) = temporary
                end if
            end do
        end do
    end subroutine sort_values

end program enzyme_fixed_quadrature_vjp_hybrid
