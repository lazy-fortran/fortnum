module adaptive_frozen_trace_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use adaptive_integrand_autodiff, only: integrand_jvp, singular_integrand_jvp
    use fortnum_generated_enzyme_adaptive_frozen_trace, only: &
        fortnum_enzyme_adaptive_frozen_trace_jvp
    use fortnum_generated_enzyme_singular_frozen_trace, only: &
        fortnum_enzyme_singular_frozen_trace_jvp
    use fortnum_integrate, only: integrate_workspace_t, integrate_epstab_t, &
        integrate_result_t, integrate_qag, integrate_qag_jvp, integrate_qags, &
        integrate_qags_jvp
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
    real(dp), parameter :: xgk21(11) = [ &
        0.9956571630258080807355272806890_dp, 0.9739065285171717200779640120845_dp, &
        0.9301574913557082260012071800595_dp, 0.8650633666889845107320966884235_dp, &
        0.7808177265864168970637175783450_dp, 0.6794095682990244062343273651149_dp, &
        0.5627571346686046833390000992727_dp, 0.4333953941292471907992659431658_dp, &
        0.2943928627014601981311266031039_dp, 0.1488743389816312108848260011297_dp, &
        0.0_dp]
    real(dp), parameter :: wgk21(11) = [ &
        0.01169463886737187427806439606219_dp, 0.03255816230796472747881897245939_dp, &
        0.05475589657435199603138130024458_dp, 0.07503967481091995276704314091619_dp, &
        0.09312545458369760553506546508337_dp, 0.1093871588022976418992105903258_dp, &
        0.1234919762620658510779581098311_dp, 0.1347092173114733259280540017717_dp, &
        0.1427759385770600807970942731387_dp, 0.1477391049013384913748415159721_dp, &
        0.1494455540029169056649364683898_dp]
    integer, parameter :: singular_trace_n = 6
    real(dp), parameter :: singular_trace_a(singular_trace_n) = [ &
        0.0_dp, 0.03125_dp, 0.0625_dp, 0.125_dp, 0.25_dp, 0.5_dp]
    real(dp), parameter :: singular_trace_b(singular_trace_n) = [ &
        0.03125_dp, 0.0625_dp, 0.125_dp, 0.25_dp, 0.5_dp, 1.0_dp]

    type :: parameter_t
        real(dp) :: value
    end type parameter_t

    public :: analytical_value_jvp, compact_analytical_value_jvp
    public :: autodiff_value_jvp, hybrid_value_jvp, compact_hybrid_value_jvp
    public :: diagnostic_value_jvp, exact_jvp
    public :: singular_analytical_value_jvp, singular_compact_analytical_value_jvp
    public :: singular_autodiff_value_jvp, singular_hybrid_value_jvp
    public :: singular_compact_hybrid_value_jvp, singular_diagnostic_value_jvp
    public :: singular_exact_jvp

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

    function singular_frozen_trace_value(p) result(value) &
            bind(c, name="fortnum_singular_frozen_trace_value")
        real(c_double), value :: p
        real(c_double) :: value
        integer :: panel

        value = 0.0_dp
        do panel = 1, singular_trace_n
            value = value + singular_panel_value(p, singular_trace_a(panel), &
                singular_trace_b(panel))
        end do
    end function singular_frozen_trace_value

    pure function singular_panel_value(p, a, b) result(value)
        real(dp), intent(in) :: p, a, b
        real(dp) :: value, center, half_length, offset
        integer :: node

        center = 0.5_dp*(a + b)
        half_length = 0.5_dp*(b - a)
        value = wgk21(11)*exp(p)/sqrt(center)
        do node = 1, 10
            offset = half_length*xgk21(node)
            value = value + wgk21(node)*exp(p)*( &
                1.0_dp/sqrt(center - offset) + 1.0_dp/sqrt(center + offset))
        end do
        value = half_length*value
    end function singular_panel_value

    function singular_panel_hybrid_tangent(p, a, b) result(value)
        real(dp), intent(in) :: p, a, b
        real(dp) :: value, center, half_length, offset
        integer :: node

        center = 0.5_dp*(a + b)
        half_length = 0.5_dp*(b - a)
        value = wgk21(11)*singular_integrand_jvp(center, p)
        do node = 1, 10
            offset = half_length*xgk21(node)
            value = value + wgk21(node)*( &
                singular_integrand_jvp(center - offset, p) + &
                singular_integrand_jvp(center + offset, p))
        end do
        value = half_length*value
    end function singular_panel_hybrid_tangent

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
        derivative = fortnum_enzyme_adaptive_frozen_trace_jvp(p, 1.0_dp)
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

    function singular_analytical_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        call integrate_qags_jvp(singular_tangent_integrand, result, derivative, &
            status, ctx=parameter)
        if (.not. status_ok(status)) error stop "singular analytical JVP failed"
    end function singular_analytical_value_jvp

    function singular_compact_analytical_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        derivative = singular_frozen_trace_value(p)
    end function singular_compact_analytical_value_jvp

    function singular_autodiff_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        derivative = fortnum_enzyme_singular_frozen_trace_jvp(p, 1.0_dp)
    end function singular_autodiff_value_jvp

    function singular_hybrid_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        call integrate_qags_jvp(singular_autodiff_tangent_integrand, result, &
            derivative, status, ctx=parameter)
        if (.not. status_ok(status)) error stop "singular hybrid JVP failed"
    end function singular_hybrid_value_jvp

    function singular_compact_hybrid_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status
        integer :: panel

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        derivative = 0.0_dp
        do panel = 1, singular_trace_n
            derivative = derivative + singular_panel_hybrid_tangent(p, &
                singular_trace_a(panel), singular_trace_b(panel))
        end do
    end function singular_compact_hybrid_value_jvp

    function singular_diagnostic_value_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp) :: derivative
        type(parameter_t) :: parameter
        type(integrate_workspace_t) :: workspace
        type(integrate_result_t) :: result
        type(fortnum_status_t) :: status

        parameter%value = p
        call build_singular_trace(parameter, workspace, result, status)
        call configure_singular_trace(result, status)
        derivative = (singular_frozen_trace_value(p + h) - &
            singular_frozen_trace_value(p - h))/(2.0_dp*h)
    end function singular_diagnostic_value_jvp

    subroutine build_singular_trace(parameter, workspace, result, status)
        type(parameter_t), intent(in) :: parameter
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_result_t), intent(inout) :: result
        type(fortnum_status_t), intent(out) :: status
        type(integrate_epstab_t) :: epstab

        call integrate_qags(singular_primal_integrand, 0.0_dp, 1.0_dp, &
            0.0_dp, 1.0e-10_dp, workspace, epstab, result, status, &
            ctx=parameter)
    end subroutine build_singular_trace

    subroutine configure_singular_trace(result, status)
        type(integrate_result_t), intent(in) :: result
        type(fortnum_status_t), intent(in) :: status

        if (.not. status_ok(status)) error stop "singular QAGS trace failed"
        if (result%nsub /= singular_trace_n) &
            error stop "unexpected singular QAGS trace size"
        if (any(result%sub_a(1:singular_trace_n) /= singular_trace_a) .or. &
            any(result%sub_b(1:singular_trace_n) /= singular_trace_b)) &
            error stop "singular QAGS trace panels changed"
    end subroutine configure_singular_trace

    pure function singular_exact_jvp(p) result(derivative)
        real(dp), intent(in) :: p
        real(dp) :: derivative

        derivative = 2.0_dp*exp(p)
    end function singular_exact_jvp

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

    function singular_primal_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = exp(context%value)/sqrt(x)
        class default
            error stop "missing singular parameter"
        end select
    end function singular_primal_integrand

    function singular_tangent_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = exp(context%value)/sqrt(x)
        class default
            error stop "missing singular parameter"
        end select
    end function singular_tangent_integrand

    function singular_autodiff_tangent_integrand(x, context) result(value)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: context
        real(dp) :: value

        select type (context)
            type is (parameter_t)
            value = singular_integrand_jvp(x, context%value)
        class default
            error stop "missing singular parameter"
        end select
    end function singular_autodiff_tangent_integrand

end module adaptive_frozen_trace_kernel

program enzyme_adaptive_frozen_trace_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use adaptive_frozen_trace_kernel, only: analytical_value_jvp, &
        compact_analytical_value_jvp, autodiff_value_jvp, hybrid_value_jvp, &
        compact_hybrid_value_jvp, diagnostic_value_jvp, exact_jvp
    use adaptive_frozen_trace_kernel, only: singular_analytical_value_jvp, &
        singular_compact_analytical_value_jvp, singular_autodiff_value_jvp, &
        singular_hybrid_value_jvp, singular_compact_hybrid_value_jvp, &
        singular_diagnostic_value_jvp, singular_exact_jvp
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_result
    implicit none

    character(32) :: argument, candidate
    real(dp) :: errors(6), reference, samples(fixture_sample_count), sink
    integer :: candidate_kind, repetitions
    logical :: singular_mode

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--tournament") then
        call run_tournament()
    else if (trim(argument) == "--singular-tournament") then
        call run_singular_tournament()
    else if (trim(argument) == "--singular-benchmark") then
        call run_singular_benchmark(trim(candidate))
    else if (trim(argument) == "--singular-peak-rss") then
        call run_singular_peak_rss(trim(candidate))
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
        call validate_singular_candidates()
        print *, "PASS adaptive frozen-trace autodiff JVP", errors
    end if

contains

    subroutine run_tournament()
        singular_mode = .false.
        call run_interleaved_tournament()
    end subroutine run_tournament

    subroutine run_singular_tournament()
        call validate_singular_candidates()
        singular_mode = .true.
        call run_interleaved_tournament()
    end subroutine run_singular_tournament

    subroutine run_interleaved_tournament()
        integer, parameter :: candidate_count = 6, tournament_samples = 31
        character(20), parameter :: names(candidate_count) = [ &
            character(20) :: "analytical", "analytical_compact", "autodiff", &
            "hybrid", "hybrid_compact", "diagnostic"]
        real(dp) :: elapsed(candidate_count, tournament_samples), median, mad
        integer :: candidate_index, order_index, sample

        repetitions = 200
        do candidate_index = 1, candidate_count
            call select_candidate(names(candidate_index))
            sink = measure_candidate()
        end do
        repetitions = 2000
        do sample = 1, tournament_samples
            do order_index = 1, candidate_count
                candidate_index = 1 + mod(sample + order_index - 2, &
                    candidate_count)
                call select_candidate(names(candidate_index))
                elapsed(candidate_index, sample) = measure_candidate()
            end do
        end do
        do candidate_index = 1, candidate_count
            call median_mad(elapsed(candidate_index, :), median, mad)
            call write_fixture_result( &
                names(candidate_index), repetitions, median, mad)
        end do
    end subroutine run_interleaved_tournament

    subroutine run_singular_benchmark(name)
        character(*), intent(in) :: name

        call validate_singular_name(name)
        singular_mode = .true.
        repetitions = 5000
        call benchmark_one(name)
    end subroutine run_singular_benchmark

    subroutine run_singular_peak_rss(name)
        character(*), intent(in) :: name

        call validate_singular_name(name)
        singular_mode = .true.
        repetitions = 5000
        call select_candidate(name)
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_singular_peak_rss

    subroutine validate_singular_name(name)
        character(*), intent(in) :: name

        if ((name /= "analytical") .and. (name /= "analytical_compact") .and. &
            (name /= "autodiff") .and. (name /= "hybrid") .and. &
            (name /= "hybrid_compact") .and. (name /= "diagnostic")) then
            error stop "unknown singular adaptive candidate"
        end if
    end subroutine validate_singular_name

    subroutine validate_singular_candidates()
        real(dp) :: values(6), reference

        reference = singular_diagnostic_value_jvp(0.7_dp)
        values(1) = singular_analytical_value_jvp(0.7_dp)
        values(2) = singular_compact_analytical_value_jvp(0.7_dp)
        values(3) = singular_autodiff_value_jvp(0.7_dp)
        values(4) = singular_hybrid_value_jvp(0.7_dp)
        values(5) = singular_compact_hybrid_value_jvp(0.7_dp)
        values(6) = reference
        if (maxval(abs(values(1:5) - reference)) > 2.0e-8_dp .or. &
            abs(reference - singular_exact_jvp(0.7_dp)) > 2.0e-2_dp) then
            print *, "singular frozen-trace JVP mismatch", values - reference
            error stop 1
        end if
    end subroutine validate_singular_candidates

    subroutine run_benchmark(name)
        character(*), intent(in) :: name

        call validate_candidate(name)
        singular_mode = .false.
        repetitions = 10000
        call benchmark_one(name)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name

        call validate_candidate(name)
        singular_mode = .false.
        repetitions = 10000
        call select_candidate(name)
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine validate_candidate(name)
        character(*), intent(in) :: name

        if ((name /= "analytical") .and. (name /= "analytical_compact") .and. &
            (name /= "autodiff") .and. (name /= "hybrid") .and. &
            (name /= "hybrid_compact") .and. (name /= "diagnostic")) then
            error stop "unknown smooth adaptive candidate"
        end if
    end subroutine validate_candidate

    subroutine benchmark_one(name)
        character(*), intent(in) :: name
        real(dp) :: median, mad

        call select_candidate(name)
        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_result(name, repetitions, median, mad)
    end subroutine benchmark_one

    subroutine select_candidate(name)
        character(*), intent(in) :: name

        select case (name)
        case ("analytical")
            candidate_kind = 1
        case ("analytical_compact")
            candidate_kind = 2
        case ("autodiff")
            candidate_kind = 3
        case ("hybrid")
            candidate_kind = 4
        case ("hybrid_compact")
            candidate_kind = 5
        case default
            candidate_kind = 6
        end select
    end subroutine select_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        real(dp) :: p, local_sink
        integer :: iteration

        local_sink = 0.0_dp
        call timer%start()
        do iteration = 1, repetitions
            if (singular_mode) then
                p = 0.7_dp
                select case (candidate_kind)
                case (1)
                    local_sink = local_sink + singular_analytical_value_jvp(p)
                case (2)
                    local_sink = local_sink + &
                        singular_compact_analytical_value_jvp(p)
                case (3)
                    local_sink = local_sink + singular_autodiff_value_jvp(p)
                case (4)
                    local_sink = local_sink + singular_hybrid_value_jvp(p)
                case (5)
                    local_sink = local_sink + singular_compact_hybrid_value_jvp(p)
                case default
                    local_sink = local_sink + singular_diagnostic_value_jvp(p)
                end select
            else
                p = 12.0_dp + 0.001_dp*real(mod(iteration, 101), dp)
                select case (candidate_kind)
                case (1)
                    local_sink = local_sink + analytical_value_jvp(p)
                case (2)
                    local_sink = local_sink + compact_analytical_value_jvp(p)
                case (3)
                    local_sink = local_sink + autodiff_value_jvp(p)
                case (4)
                    local_sink = local_sink + hybrid_value_jvp(p)
                case (5)
                    local_sink = local_sink + compact_hybrid_value_jvp(p)
                case default
                    local_sink = local_sink + diagnostic_value_jvp(p)
                end select
            end if
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) then
            error stop "adaptive frozen-trace benchmark failed"
        end if
        sink = local_sink
    end function measure_candidate

end program enzyme_adaptive_frozen_trace_jvp
