module fortnum_minimize
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ad_interfaces, only: grad_fn, fortnum_ad_status_t, ad_status_ok
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
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
    public :: bfgs_minimize
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
                    call status_set(result%status, FORTNUM_DOMAIN_ERROR, &
                        'Objective callback failed during line search')
                    return
                end if
                if (.not. ieee_is_finite(ft) .or. .not. all(ieee_is_finite(gt))) then
                    call status_set(result%status, FORTNUM_DOMAIN_ERROR, &
                        'Nonfinite objective or gradient during line search')
                    return
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
end module
