module minimize_fixture
    use fortnum_kinds, only: dp
    use fortnum_ad_interfaces
    use fortnum_status, only: FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
contains
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
            case (1)
                g = [4*x(1) + x(2) - 3, x(1) + 3*x(2) + 2]
                f = 2*x(1)**2 + x(1)*x(2) + 1.5_dp*x(2)**2 - 3*x(1) + 2*x(2)
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
