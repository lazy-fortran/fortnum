module fortnum_minimize
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ad_interfaces, only: grad_fn, fortnum_ad_status_t, ad_status_ok
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_symmetric_eigen, only: symmetric_eigen
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    type, public :: bfgs_options_t
        integer :: max_iterations = 500
        integer :: max_line_search = 50
        real(dp) :: gradient_tolerance = 1.0e-8_dp
        real(dp) :: step_tolerance = 1.0e-14_dp
        real(dp) :: armijo = 1.0e-4_dp
        real(dp) :: contraction = 0.5_dp
    end type
    type, public :: bfgs_result_t
        type(fortnum_status_t) :: status
        logical :: converged = .false.
        integer :: iterations = 0
        integer :: function_calls = 0
        real(dp) :: value = huge(1.0_dp)
        real(dp) :: gradient_norm = huge(1.0_dp)
    end type
    type, public :: stationary_options_t
        integer :: max_steps = 8
        real(dp) :: gradient_tolerance = 1.0e-6_dp
        real(dp) :: hessian_step = 1.0e-4_dp
        real(dp) :: curvature_tolerance = 5.0e-4_dp
        real(dp) :: relative_radius = 1.0e-3_dp
    end type
    type, public :: stationary_result_t
        type(fortnum_status_t) :: status
        logical :: converged = .false.
        logical :: curvature_valid = .false.
        integer :: accepted_steps = 0, function_calls = 0
        integer :: last_root_status = FORTNUM_OK
        real(dp) :: value = huge(1.0_dp)
        real(dp) :: initial_gradient_norm = huge(1.0_dp)
        real(dp) :: gradient_norm = huge(1.0_dp)
        real(dp) :: hessian_difference = huge(1.0_dp)
        real(dp) :: minimum_eigenvalue = -huge(1.0_dp)
        real(dp) :: last_step_norm = 0.0_dp
    end type
    type :: stationary_state_t
        procedure(grad_fn), pointer, nopass :: objective => null()
        class(*), pointer :: context => null()
        type(stationary_options_t) :: options
        type(stationary_result_t) :: result
        logical :: failed = .false.
    end type
    type :: stationary_adapter_t
        type(stationary_state_t), pointer :: state => null()
    end type
    public :: bfgs_minimize, stationary_refine
contains
    subroutine bfgs_minimize(objective, x, context, result, options)
        procedure(grad_fn) :: objective
        real(dp), intent(inout) :: x(:)
        class(*), intent(inout) :: context
        type(bfgs_result_t), intent(out) :: result
        type(bfgs_options_t), intent(in), optional :: options
        type(bfgs_options_t) :: opt
        type(fortnum_ad_status_t) :: callback_status
        real(dp), allocatable :: h(:, :), g(:), gt(:), p(:), xt(:), s(:), y(:), hy(:)
        real(dp) :: f, ft, alpha, slope, ys, yhy, threshold
        integer :: n, iteration, trial, i, j
        logical :: accepted

        result = bfgs_result_t()
        opt = bfgs_options_t()
        if (present(options)) opt = options
        call status_set(result%status, FORTNUM_DOMAIN_ERROR, 'Invalid BFGS input')
        n = size(x)
        if (n < 1) return
        if (.not. all(ieee_is_finite(x))) return
        if (opt%max_iterations < 0 .or. opt%max_line_search < 1) return
        if (.not. ieee_is_finite(opt%gradient_tolerance)) return
        if (.not. ieee_is_finite(opt%step_tolerance)) return
        if (.not. ieee_is_finite(opt%armijo)) return
        if (.not. ieee_is_finite(opt%contraction)) return
        if (opt%gradient_tolerance < 0 .or. opt%step_tolerance < 0) return
        if (opt%armijo <= 0 .or. opt%armijo >= 1) return
        if (opt%contraction <= 0 .or. opt%contraction >= 1) return
        allocate (h(n, n), g(n), gt(n), p(n), xt(n), s(n), y(n), hy(n))
        call objective(n, x, f, g, context, callback_status)
        result%function_calls = 1
        if (.not. ad_status_ok(callback_status)) then
            call status_set(result%status, FORTNUM_DOMAIN_ERROR, &
                'Objective callback failed')
            return
        end if
        if (.not. ieee_is_finite(f) .or. .not. all(ieee_is_finite(g))) then
            call status_set(result%status, FORTNUM_DOMAIN_ERROR, &
                'Nonfinite objective or gradient')
            return
        end if
        result%value = f
        result%gradient_norm = norm2(g)
        h = 0
        do i = 1, n
            h(i, i) = 1
        end do
        do iteration = 0, opt%max_iterations
            if (result%gradient_norm <= opt%gradient_tolerance) then
                result%converged = .true.
                call status_set(result%status, FORTNUM_OK, 'Gradient tolerance reached')
                return
            end if
            if (iteration == opt%max_iterations) exit
            p = -matmul(h, g)
            slope = dot_product(g, p)
            if (.not. ieee_is_finite(slope) .or. slope >= 0) then
                h = 0
                do i = 1, n
                    h(i, i) = 1
                end do
                p = -g
                slope = -dot_product(g, g)
            end if
            if (.not. ieee_is_finite(slope)) then
                call status_set(result%status, FORTNUM_CONVERGENCE_ERROR, &
                    'Nonfinite search direction')
                return
            end if
            alpha = 1
            accepted = .false.
            do trial = 1, opt%max_line_search
                xt = x + alpha*p
                if (.not. all(ieee_is_finite(xt))) then
                    alpha = alpha*opt%contraction
                    cycle
                end if
                call objective(n, xt, ft, gt, context, callback_status)
                result%function_calls = result%function_calls + 1
                if (.not. ad_status_ok(callback_status)) then
                    if (callback_status%status%code == FORTNUM_DOMAIN_ERROR) then
                        alpha = alpha*opt%contraction
                        cycle
                    end if
                    call status_set(result%status, FORTNUM_DOMAIN_ERROR, &
                        'Objective callback failed during line search')
                    return
                end if
                if (.not. ieee_is_finite(ft) .or. .not. all(ieee_is_finite(gt))) then
                    alpha = alpha*opt%contraction
                    cycle
                end if
                if (ft <= f + opt%armijo*alpha*slope) then
                    accepted = .true.
                    exit
                end if
                alpha = alpha*opt%contraction
            end do
            if (.not. accepted) then
                call status_set(result%status, FORTNUM_CONVERGENCE_ERROR, &
                    'Line search exhausted')
                return
            end if
            s = xt - x
            y = gt - g
            x = xt
            g = gt
            f = ft
            result%iterations = iteration + 1
            result%value = f
            result%gradient_norm = norm2(g)
            if (result%gradient_norm > opt%gradient_tolerance) then
                if (norm2(s) <= opt%step_tolerance*max(1.0_dp, norm2(x))) then
                    call status_set(result%status, FORTNUM_CONVERGENCE_ERROR, &
                        'Step stagnated before gradient convergence')
                    return
                end if
            end if
            ys = dot_product(y, s)
            threshold = sqrt(epsilon(1.0_dp))*norm2(y)*norm2(s)
            ! Skip unreliable curvature rather than corrupt the positive metric.
            if (ieee_is_finite(ys) .and. ys > threshold) then
                hy = matmul(h, y)
                yhy = dot_product(y, hy)
                do j = 1, n
                    do i = 1, n
                        h(i, j) = h(i, j) + (1 + yhy/ys)*s(i)*s(j)/ys &
                            - (hy(i)*s(j) + s(i)*hy(j))/ys
                    end do
                end do
            end if
        end do
        call status_set(result%status, FORTNUM_CONVERGENCE_ERROR, &
            'Iteration limit reached')
    end subroutine
    subroutine stationary_refine(objective, x, context, result, options)
        !! Opt-in local gradient-root refinement; this does not establish a global minimum.
        !! Scalar values remain diagnostic: acceptance requires a smaller Euclidean gradient
        !! and stable positive curvature, within a local relative step radius.
        procedure(grad_fn) :: objective
        real(dp), intent(inout) :: x(:)
        class(*), target, intent(inout) :: context
        type(stationary_result_t), intent(out) :: result
        type(stationary_options_t), optional, intent(in) :: options
        type(stationary_state_t), target :: state

        state%objective => objective
        state%context => context
        if (present(options)) state%options = options
        call stationary_run(x, state)
        result = state%result
    end subroutine

    subroutine stationary_run(x, state)
        real(dp), intent(inout) :: x(:)
        type(stationary_state_t), target, intent(inout) :: state
        type(stationary_adapter_t) :: adapter
        type(fortnum_status_t) :: root_status
        real(dp), allocatable :: g(:), gt(:), h(:, :), values(:), vectors(:, :), step(:), trial(:)
        real(dp) :: f, ft, radius, trial_norm
        integer :: n, iteration
        logical :: valid

        call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Invalid stationary refinement input')
        n = size(x)
        if (n < 1 .or. .not. all(ieee_is_finite(x))) return
        associate (opt => state%options)
            if (opt%max_steps < 0 .or. opt%max_steps > 8) return
            if (.not. ieee_is_finite(opt%gradient_tolerance) .or. &
                .not. ieee_is_finite(opt%hessian_step) .or. &
                .not. ieee_is_finite(opt%curvature_tolerance) .or. &
                .not. ieee_is_finite(opt%relative_radius)) return
            if (opt%gradient_tolerance <= 0 .or. opt%hessian_step <= 0 .or. &
                opt%curvature_tolerance <= 0 .or. opt%relative_radius <= 0) return
        end associate
        allocate (g(n), gt(n), h(n, n), values(n), vectors(n, n), step(n), trial(n))
        call stationary_evaluate(x, f, g, state)
        if (state%failed) return
        state%result%value = f
        state%result%initial_gradient_norm = norm2(g)
        state%result%gradient_norm = norm2(g)
        if (.not. ieee_is_finite(state%result%gradient_norm)) then
            call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Nonfinite gradient norm')
            return
        end if
        call status_set(state%result%status, FORTNUM_CONVERGENCE_ERROR, &
            'Stationary refinement did not reach gradient tolerance')
        adapter%state => state
        do iteration = 0, state%options%max_steps
            call stationary_curvature(x, h, values, vectors, valid, state)
            if (.not. valid) return
            if (state%result%gradient_norm <= state%options%gradient_tolerance) then
                state%result%converged = .true.
                call status_set(state%result%status, FORTNUM_OK, 'Gradient tolerance and curvature checks reached')
                return
            end if
            if (iteration == state%options%max_steps) return
            radius = state%options%relative_radius*(1 + norm2(x))
            if (.not. ieee_is_finite(radius)) return
            step = -matmul(vectors, matmul(transpose(vectors), g)/values)
            if (.not. all(ieee_is_finite(step))) return
            if (norm2(step) > radius) return
            call multiroot_hybrid(stationary_root, n, x, trial, root_status, &
                xtol=epsilon(1.0_dp), ftol=state%options%gradient_tolerance/sqrt(real(n, dp)), &
                max_iter=1, ctx=adapter)
            state%result%last_root_status = root_status%code
            if (state%failed) return
            if (.not. all(ieee_is_finite(trial))) return
            state%result%last_step_norm = norm2(trial - x)
            if (.not. ieee_is_finite(state%result%last_step_norm)) return
            if (state%result%last_step_norm > radius) return
            call stationary_evaluate(trial, ft, gt, state)
            if (state%failed) return
            trial_norm = norm2(gt)
            if (.not. (trial_norm < state%result%gradient_norm)) return
            call stationary_curvature(trial, h, values, vectors, valid, state)
            if (.not. valid) return
            x = trial
            g = gt
            state%result%value = ft
            state%result%gradient_norm = trial_norm
            state%result%accepted_steps = state%result%accepted_steps + 1
        end do
    end subroutine

    subroutine stationary_evaluate(x, f, g, state)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        type(stationary_state_t), intent(inout) :: state
        type(fortnum_ad_status_t) :: callback_status

        f = ieee_value(1.0_dp, ieee_quiet_nan)
        g = f
        if (state%failed) return
        if (.not. all(ieee_is_finite(x))) then
            state%failed = .true.
            call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Nonfinite callback input')
            return
        end if
        call state%objective(size(x), x, f, g, state%context, callback_status)
        state%result%function_calls = state%result%function_calls + 1
        if (.not. ad_status_ok(callback_status)) then
            state%failed = .true.
            state%result%status = callback_status%status
            if (state%result%status%code == FORTNUM_OK) &
                call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Invalid gradient quality')
        else if (.not. ieee_is_finite(f) .or. .not. all(ieee_is_finite(g))) then
            state%failed = .true.
            call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Nonfinite objective or gradient')
        end if
    end subroutine

    subroutine stationary_curvature(x, h, values, vectors, valid, state)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: h(:, :), values(:), vectors(:, :)
        logical, intent(out) :: valid
        type(stationary_state_t), intent(inout) :: state
        real(dp) :: coarse(size(x), size(x)), ev(size(x)), vec(size(x), size(x))
        real(dp) :: scale, symmetry
        integer :: info

        valid = .false.
        state%result%curvature_valid = .false.
        state%result%hessian_difference = huge(1.0_dp)
        state%result%minimum_eigenvalue = -huge(1.0_dp)
        call stationary_hessian(x, state%options%hessian_step, coarse, state)
        if (state%failed) return
        call stationary_hessian(x, state%options%hessian_step/2, h, state)
        if (state%failed) return
        if (.not. all(ieee_is_finite(coarse)) .or. .not. all(ieee_is_finite(h))) return
        scale = max(tiny(1.0_dp), maxval(abs(coarse)), maxval(abs(h)))
        state%result%hessian_difference = maxval(abs(h - coarse))/scale
        symmetry = max(maxval(abs(h - transpose(h))), maxval(abs(coarse - transpose(coarse))))/scale
        if (max(symmetry, state%result%hessian_difference) > state%options%curvature_tolerance) return
        coarse = (coarse + transpose(coarse))/2
        h = (h + transpose(h))/2
        call symmetric_eigen(coarse, ev, vec, info)
        if (info /= 0) return
        if (.not. all(ieee_is_finite(ev)) .or. .not. all(ieee_is_finite(vec))) return
        state%result%minimum_eigenvalue = minval(ev)
        if (minval(ev) <= 0) return
        call symmetric_eigen(h, values, vectors, info)
        if (info /= 0) return
        if (.not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(vectors))) return
        state%result%minimum_eigenvalue = min(minval(ev), minval(values))
        if (state%result%minimum_eigenvalue <= 0) return
        valid = .true.
        state%result%curvature_valid = .true.
    end subroutine

    subroutine stationary_hessian(x, relative_step, h, state)
        real(dp), intent(in) :: x(:), relative_step
        real(dp), intent(out) :: h(:, :)
        type(stationary_state_t), intent(inout) :: state
        real(dp) :: plus(size(x)), minus(size(x)), gp(size(x)), gm(size(x)), f, step
        integer :: j

        do j = 1, size(x)
            step = relative_step*max(1.0_dp, abs(x(j)))
            plus = x
            minus = x
            plus(j) = x(j) + step
            minus(j) = x(j) - step
            if (.not. ieee_is_finite(step) .or. plus(j) == minus(j)) then
                state%failed = .true.
                call status_set(state%result%status, FORTNUM_DOMAIN_ERROR, 'Unresolvable Hessian difference step')
                return
            end if
            call stationary_evaluate(plus, f, gp, state)
            if (state%failed) return
            call stationary_evaluate(minus, f, gm, state)
            if (state%failed) return
            h(:, j) = (gp - gm)/(2*step)
        end do
    end subroutine

    subroutine stationary_root(x, residual, jacobian, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: residual(:), jacobian(:, :)
        class(*), optional, intent(in) :: ctx
        real(dp) :: f, values(size(x)), vectors(size(x), size(x))
        logical :: valid

        residual = ieee_value(1.0_dp, ieee_quiet_nan)
        jacobian = residual(1)
        if (.not. present(ctx)) return
        select type (ctx)
            type is (stationary_adapter_t)
            call stationary_evaluate(x, f, residual, ctx%state)
            if (ctx%state%failed) return
            call stationary_curvature(x, jacobian, values, vectors, valid, ctx%state)
            if (.not. valid) then
                residual = ieee_value(1.0_dp, ieee_quiet_nan)
                jacobian = residual(1)
            end if
        end select
    end subroutine
end module
