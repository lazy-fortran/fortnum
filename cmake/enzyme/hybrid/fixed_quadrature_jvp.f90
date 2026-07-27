module fixed_quadrature_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortnum_kinds, only: dp
    use fortnum_quadrature, only: gauss_legendre_jvp
    implicit none
    private

    public :: analytical_integral_jvp, hybrid_integral_jvp
    public :: diagnostic_integral_jvp, exact_integral_jvp

    interface
        function enzyme_fwddiff(f, x, dx, p1, dp1, p2, dp2, p3, dp3, &
                p4, dp4) result(df) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x, dx, p1, dp1, p2, dp2
            real(c_double), value :: p3, dp3, p4, dp4
            real(c_double) :: df
        end function enzyme_fwddiff
    end interface

contains

    pure function integrand(x, p1, p2, p3, p4) result(value) &
            bind(c, name="fortnum_fixed_quadrature_integrand")
        real(c_double), value :: x, p1, p2, p3, p4
        real(c_double) :: value

        value = exp(p1*x) + sin(p2*x) + p3*x*x + p4*x*x*x
    end function integrand

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

    function hybrid_integral_jvp(nodes, weights, parameters, direction) &
            result(jvp)
        real(dp), intent(in) :: nodes(:), weights(:), parameters(4), direction(4)
        real(dp) :: jvp, tangent_values(size(nodes)), contracted(1)
        integer :: i

        do i = 1, size(nodes)
            tangent_values(i) = enzyme_fwddiff(c_funloc(integrand), &
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
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fixed_quadrature_hybrid_kernel, only: analytical_integral_jvp, &
        hybrid_integral_jvp, diagnostic_integral_jvp, exact_integral_jvp
    use fortnum_quadrature, only: gauss_legendre_ab
    implicit none

    integer, parameter :: rule_size = 32
    integer, parameter :: max_directions = 4
    real(dp) :: nodes(rule_size), weights(rule_size)
    real(dp) :: parameters(4), directions(4, max_directions), errors(3)
    character(32) :: argument, candidate, direction_argument
    integer :: direction_count

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
    call get_command_argument(3, direction_argument)

    if (trim(argument) == "--benchmark") then
        if (len_trim(candidate) > 0) then
            read (direction_argument, *) direction_count
            call run_single_benchmark(trim(candidate), direction_count)
        else
            call run_benchmark()
        end if
    else if (trim(argument) == "--peak-rss") then
        read (direction_argument, *) direction_count
        call run_peak_rss(trim(candidate), direction_count)
    else
        call validate_candidates()
    end if

contains

    subroutine initialize_workload()
        call gauss_legendre_ab(rule_size, 0.0_dp, 1.0_dp, nodes, weights)
        parameters = [0.7_dp, 1.1_dp, -0.3_dp, 0.2_dp]
        directions(:, 1) = [0.4_dp, -0.6_dp, 0.2_dp, 0.7_dp]
        directions(:, 2) = [-0.3_dp, 0.1_dp, 0.8_dp, -0.2_dp]
        directions(:, 3) = [0.5_dp, 0.4_dp, -0.1_dp, 0.3_dp]
        directions(:, 4) = [-0.2_dp, 0.7_dp, 0.6_dp, -0.5_dp]
    end subroutine initialize_workload

    subroutine validate_candidates()
        real(dp) :: analytical, hybrid, diagnostic, reference
        integer :: direction

        do direction = 1, max_directions
            reference = exact_integral_jvp(parameters, directions(:, direction))
            analytical = analytical_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            hybrid = hybrid_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            diagnostic = diagnostic_integral_jvp(nodes, weights, parameters, &
                directions(:, direction))
            errors(1) = abs(analytical - reference)
            errors(2) = abs(hybrid - reference)
            errors(3) = abs(diagnostic - reference)
            if (maxval(errors(1:2)) > 2.0e-14_dp .or. &
                errors(3) > 2.0e-10_dp) then
                print *, "fixed-quadrature JVP mismatch", direction, errors
                error stop 1
            end if
        end do
        print *, "PASS fixed-quadrature hybrid JVP"
    end subroutine validate_candidates

    subroutine run_benchmark()
        integer, parameter :: counts(3) = [1, 2, 4]
        integer :: count_index

        do count_index = 1, size(counts)
            call run_single_benchmark("analytical", counts(count_index))
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

    subroutine validate_request(name, count)
        character(*), intent(in) :: name
        integer, intent(in) :: count

        if ((name /= "analytical") .and. (name /= "hybrid") .and. &
            (name /= "diagnostic")) then
            error stop "candidate must be analytical, hybrid, or diagnostic"
        end if
        if ((count /= 1) .and. (count /= 2) .and. (count /= 4)) then
            error stop "direction count must be 1, 2, or 4"
        end if
    end subroutine validate_request

    subroutine time_candidate(name, count, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer, intent(in) :: count
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        integer :: direction
        real(dp) :: varied_parameters(4), sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            varied_parameters = parameters
            varied_parameters(1) = varied_parameters(1) + &
                1.0e-6_dp*real(mod(iteration, 101_int64), dp)
            do direction = 1, count
                select case (name)
                case ("analytical")
                    sink = sink + analytical_integral_jvp(nodes, weights, &
                        varied_parameters, directions(:, direction))
                case ("hybrid")
                    sink = sink + hybrid_integral_jvp(nodes, weights, &
                        varied_parameters, directions(:, direction))
                case ("diagnostic")
                    sink = sink + diagnostic_integral_jvp(nodes, weights, &
                        varied_parameters, directions(:, direction))
                end select
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

end program enzyme_fixed_quadrature_jvp_hybrid
