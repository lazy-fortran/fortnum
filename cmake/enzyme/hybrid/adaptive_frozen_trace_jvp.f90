module adaptive_frozen_trace_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use adaptive_integrand_autodiff, only: integrand_jvp
    use fortnum_integrate, only: integrate_workspace_t, integrate_result_t, &
        integrate_qag, integrate_qag_jvp
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    real(dp), parameter :: xgk(8) = [ &
        0.9914553711208126392068546975263_dp, &
        0.9491079123427585245261896840479_dp, &
        0.8648644233597690727897127886409_dp, &
        0.7415311855993944398638647732808_dp, &
        0.5860872354676911302941448382587_dp, &
        0.4058451513773971669066064120770_dp, &
        0.2077849550078984676006894037732_dp, 0.0_dp]
    real(dp), parameter :: wgk(8) = [ &
        0.02293532201052922496373200805897_dp, &
        0.06309209262997855329070066318920_dp, &
        0.1047900103222501838398763225415_dp, &
        0.1406532597155259187451895905102_dp, &
        0.1690047266392679028265834265986_dp, &
        0.1903505780647854099132564024210_dp, &
        0.2044329400752988924141619992346_dp, &
        0.2094821410847278280129991748917_dp]

    type :: parameter_t
        real(dp) :: value
    end type parameter_t

    public :: analytical_value_jvp, compact_analytical_value_jvp
    public :: autodiff_value_jvp, hybrid_value_jvp, compact_hybrid_value_jvp
    public :: diagnostic_value_jvp, exact_jvp

    interface
        function enzyme_fwddiff(f, p, dp_seed) result(derivative) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: p, dp_seed
            real(c_double) :: derivative
        end function enzyme_fwddiff

    end interface

contains

    pure function frozen_trace_value(p) result(value) &
            bind(c, name="fortnum_adaptive_frozen_trace_value")
        real(c_double), value :: p
        real(c_double) :: value

        value = panel_value(p, 0.0_dp, 0.5_dp) + &
            panel_value(p, 0.5_dp, 0.75_dp) + &
            panel_value(p, 0.75_dp, 1.0_dp)
    end function frozen_trace_value

    pure function panel_value(p, a, b) result(value)
        real(dp), intent(in) :: p, a, b
        real(dp) :: value, center, half_length, offset
        integer :: node

        center = 0.5_dp*(a + b)
        half_length = 0.5_dp*(b - a)
        value = wgk(8)*exp(p*center)
        do node = 1, 7
            offset = half_length*xgk(node)
            value = value + wgk(node)*( &
                exp(p*(center - offset)) + exp(p*(center + offset)))
        end do
        value = half_length*value
    end function panel_value

    pure function panel_tangent_value(p, a, b) result(value)
        real(dp), intent(in) :: p, a, b
        real(dp) :: value, center, half_length, offset, x_minus, x_plus
        integer :: node

        center = 0.5_dp*(a + b)
        half_length = 0.5_dp*(b - a)
        value = wgk(8)*center*exp(p*center)
        do node = 1, 7
            offset = half_length*xgk(node)
            x_minus = center - offset
            x_plus = center + offset
            value = value + wgk(node)*( &
                x_minus*exp(p*x_minus) + x_plus*exp(p*x_plus))
        end do
        value = half_length*value
    end function panel_tangent_value

    function panel_hybrid_tangent_value(p, a, b) result(value)
        real(dp), intent(in) :: p, a, b
        real(dp) :: value, center, half_length, offset
        integer :: node

        center = 0.5_dp*(a + b)
        half_length = 0.5_dp*(b - a)
        value = wgk(8)*integrand_jvp(center, p)
        do node = 1, 7
            offset = half_length*xgk(node)
            value = value + wgk(node)*( &
                integrand_jvp(center - offset, p) + &
                integrand_jvp(center + offset, p))
        end do
        value = half_length*value
    end function panel_hybrid_tangent_value

    function analytical_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        call integrate_qag_jvp(tangent_integrand, result, derivative, status, &
            ctx=parameter)
        if (.not. status_ok(status)) error stop "analytical trace JVP failed"
    end function analytical_value_jvp

    function compact_analytical_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        derivative = panel_tangent_value(p, 0.0_dp, 0.5_dp) + &
            panel_tangent_value(p, 0.5_dp, 0.75_dp) + &
            panel_tangent_value(p, 0.75_dp, 1.0_dp)
    end function compact_analytical_value_jvp

    function autodiff_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        derivative = enzyme_fwddiff(c_funloc(frozen_trace_value), p, 1.0_dp)
    end function autodiff_value_jvp

    function hybrid_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        call integrate_qag_jvp(autodiff_tangent_integrand, result, derivative, &
            status, ctx=parameter)
        if (.not. status_ok(status)) error stop "hybrid trace JVP failed"
    end function hybrid_value_jvp

    function compact_hybrid_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        derivative = panel_hybrid_tangent_value(p, 0.0_dp, 0.5_dp) + &
            panel_hybrid_tangent_value(p, 0.5_dp, 0.75_dp) + &
            panel_hybrid_tangent_value(p, 0.75_dp, 1.0_dp)
    end function compact_hybrid_value_jvp

    function diagnostic_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_trace(parameter, workspace, result, status)
        call require_trace(result, status)
        derivative = (frozen_trace_value(p + h) - frozen_trace_value(p - h)) &
            /(2.0_dp*h)
    end function diagnostic_value_jvp

    subroutine build_trace(parameter, workspace, result, status)
        type(parameter_t), intent(in) :: parameter
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_result_t), intent(inout) :: result
        type(fortnum_status_t), intent(out) :: status

        call integrate_qag(primal_integrand, 0.0_dp, 1.0_dp, 0.0_dp, &
            1.0e-12_dp, workspace, result, status, key=15, ctx=parameter)
    end subroutine build_trace

    subroutine require_trace(result, status)
        type(integrate_result_t), intent(in) :: result
        type(fortnum_status_t), intent(in) :: status
        logical :: found(3)
        integer :: panel

        if ((.not. status_ok(status)) .or. result%nsub /= 3) then
            error stop "unexpected adaptive trace"
        end if
        found = .false.
        do panel = 1, result%nsub
            if (result%sub_a(panel) == 0.0_dp .and. &
                result%sub_b(panel) == 0.5_dp) found(1) = .true.
            if (result%sub_a(panel) == 0.5_dp .and. &
                result%sub_b(panel) == 0.75_dp) found(2) = .true.
            if (result%sub_a(panel) == 0.75_dp .and. &
                result%sub_b(panel) == 1.0_dp) found(3) = .true.
        end do
        if (.not. all(found)) error stop "adaptive trace panels changed"
    end subroutine require_trace

    pure function exact_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative

        derivative = (exp(p)*(p - 1.0_dp) + 1.0_dp)/(p*p)
    end function exact_jvp

    function primal_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = exp(context%value*x)
        class default
            error stop "missing adaptive parameter"
        end select
    end function primal_integrand

    function tangent_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = x*exp(context%value*x)
        class default
            error stop "missing adaptive parameter"
        end select
    end function tangent_integrand

    function autodiff_tangent_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = integrand_jvp(x, context%value)
        class default
            error stop "missing adaptive parameter"
        end select
    end function autodiff_tangent_integrand

end module adaptive_frozen_trace_kernel

program enzyme_adaptive_frozen_trace_jvp
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use adaptive_frozen_trace_kernel, only: analytical_value_jvp, &
        compact_analytical_value_jvp, autodiff_value_jvp, hybrid_value_jvp, &
        compact_hybrid_value_jvp, diagnostic_value_jvp, exact_jvp
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    character(32) :: argument, candidate
    real(dp) :: errors(6), reference

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--tournament") then
        call run_tournament()
    else if (trim(argument) == "--benchmark") then
        call run_benchmark(trim(candidate))
    else if (trim(argument) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        reference = exact_jvp(12.0_dp)
        errors(1) = abs(analytical_value_jvp(12.0_dp) - reference)
        errors(2) = abs(compact_analytical_value_jvp(12.0_dp) - reference)
        errors(3) = abs(autodiff_value_jvp(12.0_dp) - reference)
        errors(4) = abs(hybrid_value_jvp(12.0_dp) - reference)
        errors(5) = abs(compact_hybrid_value_jvp(12.0_dp) - reference)
        errors(6) = abs(diagnostic_value_jvp(12.0_dp) - reference)
        if (maxval(errors(1:5)) > 1.0e-7_dp .or. &
            errors(6) > 1.0e-4_dp) then
            print *, "adaptive frozen-trace JVP mismatch", errors
            error stop 1
        end if
        print *, "PASS adaptive frozen-trace autodiff JVP", errors
    end if

contains

    subroutine run_tournament()
        integer, parameter :: candidate_count = 6, samples = 31
        integer(int64), parameter :: reps = 2000_int64
        character(20), parameter :: names(candidate_count) = [ &
            character(20) :: "analytical", "analytical_compact", "autodiff", &
            "hybrid", "hybrid_compact", "diagnostic"]
        real(dp) :: elapsed(candidate_count, samples), sink
        integer :: candidate_index, order_index, sample

        do candidate_index = 1, candidate_count
            call time_candidate(names(candidate_index), reps/10_int64, sink)
        end do
        do sample = 1, samples
            do order_index = 1, candidate_count
                candidate_index = 1 + mod(sample + order_index - 2, &
                    candidate_count)
                call time_candidate(names(candidate_index), reps, &
                    elapsed(candidate_index, sample))
            end do
        end do
        do candidate_index = 1, candidate_count
            call report(names(candidate_index), elapsed(candidate_index, :), &
                reps)
        end do
    end subroutine run_tournament

    subroutine run_benchmark(name)
        character(*), intent(in) :: name
        integer, parameter :: samples = 15
        integer(int64), parameter :: reps = 10000_int64
        real(dp) :: elapsed(samples), sink
        integer :: sample

        call validate_candidate(name)
        do sample = 1, 3
            call time_candidate(name, reps/10_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate(name, reps, elapsed(sample))
        end do
        call report(name, elapsed, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name
        integer(int64), parameter :: reps = 10000_int64
        real(dp) :: elapsed

        call validate_candidate(name)
        call time_candidate(name, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_candidate(name)
        character(*), intent(in) :: name

        if ((name /= "analytical") .and. (name /= "analytical_compact") .and. &
            (name /= "autodiff") .and. (name /= "hybrid") .and. &
            (name /= "hybrid_compact") .and. (name /= "diagnostic")) then
            error stop "unknown smooth adaptive candidate"
        end if
    end subroutine validate_candidate

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: iteration, start, finish, rate
        real(dp) :: p, sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do iteration = 1, reps
            p = 12.0_dp + 0.001_dp*real(mod(iteration, 101_int64), dp)
            select case (name)
            case ("analytical")
                sink = sink + analytical_value_jvp(p)
            case ("analytical_compact")
                sink = sink + compact_analytical_value_jvp(p)
            case ("autodiff")
                sink = sink + autodiff_value_jvp(p)
            case ("hybrid")
                sink = sink + hybrid_value_jvp(p)
            case ("hybrid_compact")
                sink = sink + compact_hybrid_value_jvp(p)
            case ("diagnostic")
                sink = sink + diagnostic_value_jvp(p)
            end select
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp)/ &
            (real(rate, dp)*real(reps, dp))
        if (sink == huge(sink)) print *, sink
    end subroutine time_candidate

    subroutine report(name, values, reps)
        character(*), intent(in) :: name
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
        write (*, "(a,',',i0,',',f12.4,',',f12.4)") &
            name, reps, median, mad
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

end program enzyme_adaptive_frozen_trace_jvp
