module minimize_fixture
    use fortnum_kinds, only: dp
    use fortnum_ad_interfaces
    use fortnum_status, only: FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_NOT_IMPLEMENTED
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    type :: refinement_fixture_t
        integer :: mode = 1, calls = 0
    end type
contains
    subroutine refinement_objective(n, x, f, g, context, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: f, g(n)
        class(*), intent(inout) :: context
        type(fortnum_ad_status_t), intent(out) :: status
        real(dp) :: delta(n)
        call ad_status_set(status, FORTNUM_OK, '', &
            FORTNUM_AD_BACKEND_ANALYTICAL, FORTNUM_AD_QUALITY_EXACT)
        select type (context)
            type is (refinement_fixture_t)
            context%calls = context%calls + 1
            delta = x - [0.3_dp, -0.2_dp]
            g = [4*delta(1) + delta(2), delta(1) + 3*delta(2)]
            f = dot_product(delta, g)/2
            select case (context%mode)
            case (1)
                f = 1.0e20_dp + f
            case (2)
                f = f/137
                g = g/137
            case (3)
                f = (delta(1)**2 - delta(2)**2)/2
                g = [delta(1), -delta(2)]
            case (4, 6)
                if (context%calls > 1) then
                    status%status%code = FORTNUM_DOMAIN_ERROR
                    if (context%mode == 6) status%status%code = FORTNUM_NOT_IMPLEMENTED
                end if
            case (5)
                delta = x - [1.0e16_dp, 0.0_dp]
                g = delta + [1.0_dp, 0.0_dp]
                f = sum(delta**2)/2 + delta(1)
            case (7)
                f = ieee_value(1.0_dp, ieee_quiet_nan)
            case (8, 9)
                if (context%calls >= 19) then
                    status%status%code = FORTNUM_DOMAIN_ERROR
                    if (context%mode == 8) status%status%code = FORTNUM_NOT_IMPLEMENTED
                end if
            case (10)
                f = sum(delta**2)/2 + 1e-8_dp*sum(cos(1e4_dp*delta))
                g = delta - 1e-4_dp*sin(1e4_dp*delta)
            case (11)
                f = sum(delta**2)/2 + 10*sum(delta**3)/3 + 25*sum(delta**4)
                g = delta + 10*delta**2 + 100*delta**3
            end select
        end select
    end subroutine

    subroutine objective(n, x, f, g, context, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: f, g(n)
        class(*), intent(inout) :: context
        type(fortnum_ad_status_t), intent(out) :: status
        real(dp) :: r(3), jac(3, 2), a(3), b(3)
        call ad_status_set(status, FORTNUM_OK, '', &
            FORTNUM_AD_BACKEND_ANALYTICAL, FORTNUM_AD_QUALITY_EXACT)
        select type (context)
            type is (integer)
            select case (context)
            case (1, 8:13)
                g = [4*x(1) + x(2) - 3, x(1) + 3*x(2) + 2]
                f = 2*x(1)**2 + x(1)*x(2) + 1.5_dp*x(2)**2 - 3*x(1) + 2*x(2)
                if (context == 12) then
                    if (abs(x(1) - 4.0_dp) > 0) status%status%code = FORTNUM_DOMAIN_ERROR
                else if (x(1) <= 0) then
                    select case (context)
                    case (8)
                        status%status%code = FORTNUM_DOMAIN_ERROR
                    case (9)
                        f = ieee_value(1.0_dp, ieee_quiet_nan)
                    case (10)
                        g = ieee_value(1.0_dp, ieee_quiet_nan)
                    case (11)
                        status%status%code = FORTNUM_NOT_IMPLEMENTED
                    case (13)
                        status%quality = FORTNUM_AD_QUALITY_NONSMOOTH
                    end select
                end if
            case (2)
                f = (1 - x(1))**2 + 100*(x(2) - x(1)**2)**2
                g = [-2*(1 - x(1)) - 400*x(1)*(x(2) - x(1)**2), &
                    200*(x(2) - x(1)**2)]
            case (3)
                ! Synthetic exponential observation residual with a Gaussian prior.
                a = [1.0_dp, 0.5_dp, -0.5_dp]
                b = [0.2_dp, -0.4_dp, 0.3_dp]
                r = exp(a*x(1) + b*x(2)) - exp(a*0.3_dp - b*0.2_dp)
                jac(:, 1) = a*exp(a*x(1) + b*x(2))
                jac(:, 2) = b*exp(a*x(1) + b*x(2))
                f = sum(r*r)/2 + 0.05_dp*sum((x - [0.3_dp, -0.2_dp])**2)
                g = matmul(transpose(jac), r) + 0.1_dp*(x - [0.3_dp, -0.2_dp])
            case (4)
                f = ieee_value(1.0_dp, ieee_quiet_nan)
                g = 1
            case (6)
                f = sum(x*x)
                g = 2*x
                status%status%code = FORTNUM_DOMAIN_ERROR
            case (7)
                f = sum(x*x)
                g = ieee_value(1.0_dp, ieee_quiet_nan)
            case (5)
                f = 1
                g = 1
            case default
                f = sum(x*x)/2
                g = x
            end select
        end select
    end subroutine
end module
program test_minimize
    use minimize_fixture
    use fortnum_minimize
    use fortnum_status, only: FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    real(dp) :: x(2), xp(2), xm(2), g(2), dummy(2), f, fp, fm, fd(2), step
    type(refinement_fixture_t) :: fixture
    type(stationary_result_t) :: refined
    type(stationary_options_t) :: refine_options
    real(dp), allocatable :: empty(:)
    type(bfgs_result_t) :: result
    type(bfgs_options_t) :: options
    type(fortnum_ad_status_t) :: status
    integer :: context, i
    context = 1
    x = [4.0_dp, -3.0_dp]
    call bfgs_minimize(objective, x, context, result)
    call require(result%converged, 'quadratic convergence')
    call require(maxval(abs(x - [1.0_dp, -1.0_dp])) < 1e-7_dp, &
        'exact quadratic minimizer')
    call require(result%function_calls > result%iterations, 'evaluation accounting')
    do i = 8, 10
        context = i
        x = [4.0_dp, -3.0_dp]
        call bfgs_minimize(objective, x, context, result)
        call require(result%converged, 'recoverable trial convergence')
        call require(maxval(abs(x - [1.0_dp, -1.0_dp])) < 1e-7_dp, &
            'recoverable trial exact minimizer')
        call require(result%function_calls > result%iterations + 2, &
            'rejected trial evaluation accounting')
    end do
    context = 8
    x = [-1.0_dp, -3.0_dp]
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, &
        'invalid initial point still fails')
    call require(result%function_calls == 1, 'invalid initial point evaluated once')
    call require(maxval(abs(x - [-1.0_dp, -3.0_dp])) <= 0, 'invalid initial point unchanged')
    context = 11
    x = [4.0_dp, -3.0_dp]
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, &
        'non-domain callback failure still terminates')
    call require(result%function_calls == 2, 'non-domain trial is not retried')
    call require(maxval(abs(x - [4.0_dp, -3.0_dp])) <= 0, 'non-domain failure preserves point')
    context = 13
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, &
        'nonsmooth trial still terminates')
    call require(result%function_calls == 2, 'nonsmooth trial is not retried')
    context = 12
    x = [4.0_dp, -3.0_dp]
    options%max_line_search = 4
    call bfgs_minimize(objective, x, context, result, options)
    call require(result%status%code == FORTNUM_CONVERGENCE_ERROR, &
        'domain trials exhaust finite line-search budget')
    call require(.not. result%converged, 'exhausted domain trials are not success')
    call require(result%function_calls == 5, 'domain trial budget accounting')
    call require(result%iterations == 0, 'domain trials never accepted')
    call require(maxval(abs(x - [4.0_dp, -3.0_dp])) <= 0, 'domain exhaustion preserves point')
    call require(abs(result%value - 15.5_dp) <= 0, 'domain exhaustion preserves value')
    call require(abs(result%gradient_norm - sqrt(109.0_dp)) < 1e-13_dp, &
        'domain exhaustion preserves gradient norm')
    options = bfgs_options_t()
    context = 2
    x = [-1.2_dp, 1.0_dp]
    call bfgs_minimize(objective, x, context, result)
    call require(result%converged, 'Rosenbrock convergence')
    call require(maxval(abs(x - 1)) < 1e-6_dp, 'Rosenbrock minimizer')
    context = 3
    x = [0.1_dp, 0.2_dp]
    step = 1e-5_dp
    call objective(2, x, f, g, context, status)
    do i = 1, 2
        xp = x
        xm = x
        xp(i) = xp(i) + step
        xm(i) = xm(i) - step
        call objective(2, xp, fp, dummy, context, status)
        call objective(2, xm, fm, dummy, context, status)
        fd(i) = (fp - fm)/(2*step)
    end do
    call require(maxval(abs(g - fd)) < 1e-9_dp, 'residual gradient finite difference')
    call bfgs_minimize(objective, x, context, result)
    call require(result%converged, 'residual convergence')
    call require(maxval(abs(x - [0.3_dp, -0.2_dp])) < 1e-6_dp, 'residual minimizer')
    context = 4
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, 'nonfinite callback')
    context = 6
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, 'callback status failure')
    context = 7
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, 'nonfinite gradient')
    context = 5
    call bfgs_minimize(objective, x, context, result)
    call require(result%status%code == FORTNUM_CONVERGENCE_ERROR, 'line search failure')
    call require(.not. result%converged, 'line search not success')
    context = 1
    options%max_iterations = 0
    call bfgs_minimize(objective, x, context, result, options)
    call require(result%status%code == FORTNUM_CONVERGENCE_ERROR, 'iteration limit')
    options = bfgs_options_t()
    options%contraction = 1
    call bfgs_minimize(objective, x, context, result, options)
    call require(result%status%code == FORTNUM_DOMAIN_ERROR, 'invalid options')
    options = bfgs_options_t()
    options%step_tolerance = 1e5_dp
    call bfgs_minimize(objective, x, context, result, options)
    call require(.not. result%converged, 'stagnation not success')

    fixture = refinement_fixture_t()
    x = [0.3_dp + 2e-5_dp, -0.2_dp - 3e-5_dp]
    xp = x
    call refinement_objective(2, x, fp, g, fixture, status)
    xm = [0.3_dp, -0.2_dp]
    call refinement_objective(2, xm, fm, dummy, fixture, status)
    call require(fp == fm, 'scalar objective cannot resolve the local quadratic improvement')
    fixture%calls = 0
    refine_options%gradient_tolerance = 1e-12_dp
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(refined%converged, 'large-offset quadratic gradient refinement')
    call require(norm2(x - [0.3_dp, -0.2_dp]) < 1e-12_dp, 'independent known stationary point')
    call require(refined%accepted_steps == 1, 'one local root step')
    call require(refined%curvature_valid, 'two-resolution positive curvature')
    call require(refined%minimum_eigenvalue > 2, 'known positive quadratic curvature')
    call require(refined%gradient_norm <= refine_options%gradient_tolerance, 'explicit Euclidean gradient tolerance')
    call require(refined%initial_gradient_norm > 1e-6_dp, 'nonstationary initial point')
    call require(refined%function_calls == fixture%calls, 'mutable callback accounting')
    call require(refined%function_calls <= 70, 'bounded one-step callback cost')

    fixture = refinement_fixture_t(mode=2)
    x = xp
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(refined%converged, 'normalized Gaussian objective stationarity')
    call require(norm2(x - [0.3_dp, -0.2_dp]) < 1e-10_dp, 'normalized Gaussian known mode')
    fixture = refinement_fixture_t(mode=11)
    x = xp
    refine_options%max_steps = 1
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(refined%accepted_steps == 1, 'useful one-iteration root trial accepted')
    call require(refined%last_root_status == FORTNUM_CONVERGENCE_ERROR, 'root budget status is diagnostic')
    call require(.not. refined%converged, 'useful step alone is not convergence')
    call require(refined%gradient_norm < refined%initial_gradient_norm, 'nonlinear gradient merit decreases')
    refine_options%max_steps = 8
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(refined%converged, 'nonlinear gradient root reaches explicit tolerance')
    call require(norm2(x - [0.3_dp, -0.2_dp]) < 1e-12_dp, 'independent nonlinear stationary point')
    fixture = refinement_fixture_t(mode=3)
    x = xp
    call stationary_refine(refinement_objective, x, fixture, refined)
    call require(.not. refined%converged, 'indefinite stationary curvature rejected')
    call require(refined%accepted_steps == 0 .and. all(x == xp), 'indefinite point unchanged')
    fixture = refinement_fixture_t(mode=5)
    x = [1.0e16_dp, 0.0_dp]
    call stationary_refine(refinement_objective, x, fixture, refined)
    call require(.not. refined%converged, 'unresolvable root step is not convergence')
    call require(refined%accepted_steps == 0, 'nondecreasing gradient root trial rejected')
    call require(all(x == [1.0e16_dp, 0.0_dp]), 'failed root leaves point unchanged')
    call require(refined%gradient_norm == 1, 'failed root preserves gradient receipt')
    do i = 4, 6, 2
        fixture = refinement_fixture_t(mode=i)
        x = xp
        call stationary_refine(refinement_objective, x, fixture, refined)
        context = FORTNUM_DOMAIN_ERROR
        if (i == 6) context = FORTNUM_NOT_IMPLEMENTED
        call require(refined%status%code == context, 'refinement preserves true callback error')
        call require(refined%function_calls == 2 .and. fixture%calls == 2, 'failed Hessian callback counted once')
        call require(all(x == xp), 'failed Hessian leaves point unchanged')
    end do
    do i = 8, 9
        fixture = refinement_fixture_t(mode=i)
        x = xp
        call stationary_refine(refinement_objective, x, fixture, refined)
        context = FORTNUM_DOMAIN_ERROR
        if (i == 8) context = FORTNUM_NOT_IMPLEMENTED
        call require(refined%status%code == context, 'root callback preserves failure status')
        call require(refined%function_calls == 19 .and. fixture%calls == 19, 'root callback failure bounded')
        call require(all(x == xp), 'root callback failure preserves accepted point')
    end do
    fixture = refinement_fixture_t(mode=10)
    x = [0.3_dp, -0.2_dp]
    call stationary_refine(refinement_objective, x, fixture, refined)
    call require(.not. refined%converged, 'unstable Hessian rejects zero gradient')
    call require(refined%hessian_difference > 0.1_dp, 'two-resolution discrepancy detected')
    fixture = refinement_fixture_t(mode=7)
    call stationary_refine(refinement_objective, x, fixture, refined)
    call require(refined%status%code == FORTNUM_DOMAIN_ERROR, 'nonfinite refinement callback rejected')
    call require(refined%function_calls == 1, 'nonfinite callback bounded')
    fixture = refinement_fixture_t()
    x = [1.0_dp, 1.0_dp]
    call stationary_refine(refinement_objective, x, fixture, refined)
    call require(.not. refined%converged .and. refined%accepted_steps == 0, 'nonlocal Newton step rejected')
    call require(all(x == 1), 'nonlocal step preserves point')
    allocate (empty(0))
    call stationary_refine(refinement_objective, empty, fixture, refined)
    call require(refined%status%code == FORTNUM_DOMAIN_ERROR, 'empty refinement input rejected')
    call require(refined%function_calls == 0, 'invalid dimension has no callback')
    refine_options%gradient_tolerance = 0
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(refined%status%code == FORTNUM_DOMAIN_ERROR, 'zero gradient tolerance rejected')
    call require(refined%function_calls == 0, 'invalid option has no callback')
    refine_options = stationary_options_t(max_steps=0)
    x = xp
    call stationary_refine(refinement_objective, x, fixture, refined, refine_options)
    call require(.not. refined%converged .and. all(x == xp), 'zero-step budget preserves nonstationary point')
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, message
            stop 1
        end if
    end subroutine
end program
