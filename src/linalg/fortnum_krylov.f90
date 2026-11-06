module fortnum_krylov
    !! Restarted matrix-free GMRES and flexible GMRES for real and complex host-side linear operators,
    !! conjugate gradients and preconditioned MINRES for real symmetric operators.
    !!
    !! Modified Gram--Schmidt is followed by a reorthogonalization pass.
    !! Complex Givens rotations update the Hessenberg least-squares problem
    !! without forming normal equations.
    use fortnum_kinds, only: dp
    use fortnum_parallel_vector, only: parallel_vector_axpy, &
        parallel_vector_dot, parallel_vector_gemv_n, parallel_vector_gemv_t, &
        parallel_vector_norm2, parallel_vector_scale, parallel_vector_sub, &
        parallel_vector_xpay
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    private

    integer, parameter, public :: KRYLOV_OK = 0
    integer, parameter, public :: KRYLOV_MAX_ITERATIONS = 1
    integer, parameter, public :: KRYLOV_BREAKDOWN = 2
    integer, parameter, public :: KRYLOV_STAGNATION = 3
    integer, parameter, public :: KRYLOV_NONPOSITIVE_CURVATURE = 4
    integer, parameter, public :: KRYLOV_INVALID_ARGUMENT = -1

    abstract interface
        subroutine real_matvec_t(input, output)
            import :: dp
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
        end subroutine real_matvec_t

        subroutine real_matmat_t(input, output)
            import :: dp
            real(dp), intent(in) :: input(:, :)
            real(dp), intent(out) :: output(:, :)
        end subroutine real_matmat_t

        subroutine complex_matvec_t(input, output)
            import :: dp
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
        end subroutine complex_matvec_t

        subroutine real_preconditioner_t(input, output)
            import :: dp
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
        end subroutine real_preconditioner_t

        subroutine complex_global_sum_t(values)
            !! Sum values over all processes that share a distributed vector.
            import :: dp
            complex(dp), intent(inout) :: values(:)
        end subroutine complex_global_sum_t

        subroutine complex_preconditioner_t(input, output)
            import :: dp
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
        end subroutine complex_preconditioner_t
    end interface

    type, abstract, public :: real_linear_operator_t
        !! Real linear map with state.  Extend it and implement `apply` to pass
        !! an operator (or a preconditioner, M^{-1}) to `real_conjugate_gradient`.
    contains
        procedure(real_linear_operator_apply_t), deferred :: apply
    end type real_linear_operator_t

    type, extends(real_linear_operator_t) :: procedure_operator_t
        !! Private adapter that wraps a stateless callback.
        procedure(real_matvec_t), pointer, nopass :: apply_procedure => null()
    contains
        procedure :: apply => apply_procedure_operator
    end type procedure_operator_t

    abstract interface
        subroutine real_linear_operator_apply_t(self, input, output)
            import :: dp, real_linear_operator_t
            class(real_linear_operator_t), intent(inout) :: self
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
        end subroutine real_linear_operator_apply_t
    end interface

    public :: complex_fgmres_operator
    public :: complex_gmres_operator
    public :: complex_preconditioner_t
    public :: complex_global_sum_t
    public :: real_fgmres_operator
    public :: complex_matvec_t
    public :: real_conjugate_gradient
    public :: real_conjugate_gradient_operator
    public :: real_conjugate_gradient_matmat_operator
    public :: real_gmres_operator
    public :: real_matmat_t
    public :: real_matvec_t
    public :: real_minres
    public :: real_minres_operator
    public :: real_preconditioner_t

    type, abstract, public :: complex_linear_operator_t
    contains
        procedure(complex_linear_apply_t), deferred :: apply
    end type

    type, extends(complex_linear_operator_t) :: complex_procedure_operator_t
        procedure(complex_matvec_t), pointer, nopass :: apply_procedure => null()
    contains
        procedure :: apply => apply_complex_procedure_operator
    end type

    abstract interface
        subroutine complex_linear_apply_t(self, input, output)
            import :: dp, complex_linear_operator_t
            class(complex_linear_operator_t), intent(inout) :: self
            complex(dp), intent(in) :: input(:)
            complex(dp), intent(out) :: output(:)
        end subroutine
    end interface

    public :: real_gmres, complex_gmres, real_fgmres, complex_fgmres

contains

    subroutine real_conjugate_gradient( &
            operator, right_hand_side, solution, tolerance, max_iterations, &
            info, iterations, residual_norm, preconditioner)
        !! Solve an SPD system given as an operator object with (preconditioned) CG.
        !!
        !! The operator must represent a symmetric positive-definite map.  It is
        !! an object so that its state (matrices, factors, work arrays) travels
        !! with the apply, and `apply` may update that state (intent(inout)).
        !! An optional preconditioner object applies M^{-1} to a residual.  The
        !! solution is used as the initial guess and is updated in place.
        !!
        !! The true residual b - A x is checked every 100 iterations.  If the
        !! iteration stalls or reaches max_iterations, the best true-residual
        !! iterate is returned with KRYLOV_MAX_ITERATIONS.
        class(real_linear_operator_t), intent(inout) :: operator
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        class(real_linear_operator_t), intent(inout), optional :: preconditioner

        real(dp), allocatable :: residual(:), direction(:), preconditioned(:)
        real(dp), allocatable :: operator_direction(:), work(:), best_solution(:)
        real(dp) :: right_hand_side_norm, target, rho, next_rho
        real(dp) :: denominator, step, beta, best_residual_norm
        integer :: last_gain, stall_iterations
        integer, parameter :: check_every = 100
        real(dp), parameter :: gain_factor = 0.9_dp

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        if (size(right_hand_side) < 1) return
        if (size(solution) /= size(right_hand_side)) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return

        allocate( &
            residual(size(solution)), direction(size(solution)), &
            preconditioned(size(solution)), operator_direction(size(solution)), &
            work(size(solution)), best_solution(size(solution)))
        right_hand_side_norm = parallel_vector_norm2(right_hand_side)
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = 0.0_dp
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        call operator%apply(solution, work)
        call parallel_vector_sub(right_hand_side, work, residual)
        residual_norm = parallel_vector_norm2(residual)
        if (residual_norm <= target) then
            info = KRYLOV_OK
            return
        end if
        best_solution = solution
        best_residual_norm = residual_norm
        last_gain = 0
        stall_iterations = max(1000, max_iterations/10)

        if (present(preconditioner)) then
            call preconditioner%apply(residual, preconditioned)
        else
            call parallel_vector_scale(1.0_dp, residual, preconditioned)
        end if
        rho = parallel_vector_dot(residual, preconditioned)
        if (rho <= 0.0_dp) then
            info = KRYLOV_BREAKDOWN
            return
        end if
        if (.not. ieee_is_finite(rho)) then
            info = KRYLOV_BREAKDOWN
            return
        end if
        call parallel_vector_scale(1.0_dp, preconditioned, direction)

        do while (iterations < max_iterations)
            call operator%apply(direction, operator_direction)
            denominator = parallel_vector_dot(direction, operator_direction)
            if (denominator <= 0.0_dp) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            if (.not. ieee_is_finite(denominator)) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            step = rho/denominator
            call parallel_vector_axpy(step, direction, solution)
            call parallel_vector_axpy(-step, operator_direction, residual)
            iterations = iterations + 1
            residual_norm = parallel_vector_norm2(residual)
            if (residual_norm <= target) then
                ! Recompute the residual from the operator before accepting
                ! recurrence-based convergence; this catches loss of
                ! orthogonality and inaccurate device callbacks.
                call operator%apply(solution, work)
                call parallel_vector_sub(right_hand_side, work, residual)
                residual_norm = parallel_vector_norm2(residual)
                if (residual_norm <= target) then
                    info = KRYLOV_OK
                    return
                end if
                if (residual_norm < gain_factor*best_residual_norm) &
                    last_gain = iterations
                if (residual_norm < best_residual_norm) then
                    best_solution = solution
                    best_residual_norm = residual_norm
                end if
                call parallel_vector_sub(right_hand_side, work, residual)
            end if

            if (mod(iterations, check_every) == 0 .or. &
                iterations == max_iterations) then
                ! Periodic true-residual check: replace the recurrence
                ! residual, keep the best iterate, and stop on convergence,
                ! on a stall, or at the iteration cap.
                call operator%apply(solution, work)
                call parallel_vector_sub(right_hand_side, work, residual)
                residual_norm = parallel_vector_norm2(residual)
                if (residual_norm <= target) then
                    info = KRYLOV_OK
                    return
                end if
                if (residual_norm < gain_factor*best_residual_norm) &
                    last_gain = iterations
                if (residual_norm < best_residual_norm) then
                    best_solution = solution
                    best_residual_norm = residual_norm
                end if
                if (iterations >= max_iterations .or. &
                    iterations - last_gain >= stall_iterations) exit
            end if

            if (present(preconditioner)) then
                call preconditioner%apply(residual, preconditioned)
            else
                call parallel_vector_scale(1.0_dp, residual, preconditioned)
            end if
            next_rho = parallel_vector_dot(residual, preconditioned)
            if (next_rho <= 0.0_dp) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            if (.not. ieee_is_finite(next_rho)) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            beta = next_rho/rho
            call parallel_vector_xpay(preconditioned, beta, direction)
            rho = next_rho
        end do

        ! Stalled or capped: return the best true-residual iterate.
        solution = best_solution
        residual_norm = best_residual_norm
        info = KRYLOV_MAX_ITERATIONS
    end subroutine real_conjugate_gradient

    subroutine real_conjugate_gradient_operator( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            info, iterations, residual_norm, preconditioner)
        !! Solve an SPD matrix-free system with (preconditioned) CG.
        !!
        !! Procedure-argument form of `real_conjugate_gradient`: the callbacks
        !! are wrapped in operator objects and the same CG core runs, so both
        !! entry points give identical results.  The callbacks carry no
        !! context; use `real_conjugate_gradient` for an operator with state.
        procedure(real_matvec_t) :: matvec
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        procedure(real_preconditioner_t), optional :: preconditioner

        type(procedure_operator_t) :: operator_adapter, preconditioner_adapter

        operator_adapter%apply_procedure => matvec
        if (present(preconditioner)) then
            preconditioner_adapter%apply_procedure => preconditioner
            call real_conjugate_gradient( &
                operator_adapter, right_hand_side, solution, tolerance, &
                max_iterations, info, iterations, residual_norm, &
                preconditioner_adapter)
        else
            call real_conjugate_gradient( &
                operator_adapter, right_hand_side, solution, tolerance, &
                max_iterations, info, iterations, residual_norm)
        end if
    end subroutine real_conjugate_gradient_operator

    subroutine real_minres( &
            operator, right_hand_side, solution, tolerance, max_iterations, &
            info, iterations, residual_norm, preconditioner, curvature_direction)
        !! Preconditioned MINRES (Paige and Saunders 1975) for a symmetric,
        !! possibly indefinite or singular operator A.
        !!
        !! The optional preconditioner applies an SPD M ~ A^{-1}; the iterates
        !! minimize the M-norm of the residual over the Krylov space, so
        !! |b - A x|_M does not increase.  The solution is the initial guess
        !! and is updated in place.  The residual r_k = b - A x_k is carried by
        !! the recurrence of Liu and Roosta ("MINRES: from negative curvature
        !! detection to monotonicity properties", SIAM J. Optim. 32, 2022),
        !! the Lanczos step subtracts the previous vector before forming the
        !! diagonal coefficient (modified Gram--Schmidt).
        !!
        !! Converged (KRYLOV_OK) when |b - A x| <= tolerance |b| (Euclidean,
        !! relative to the right-hand side without the max(|b|, 1) floor of
        !! the CG and GMRES solvers, as inner solves of inexact Newton
        !! methods need).  residual_norm is the explicit |b - A x| of the
        !! returned solution; a recurrence that met the tolerance while the
        !! explicit residual does not is KRYLOV_STAGNATION.  At the iteration
        !! cap the last iterate is returned with KRYLOV_MAX_ITERATIONS.  A
        !! singular A with consistent b is admitted: the Krylov space stays in
        !! the range of A.
        !!
        !! With curvature_direction present, the first residual r_{k-1} of
        !! nonpositive curvature, (M r)^T A (M r) <= 0, ends the solve with
        !! KRYLOV_NONPOSITIVE_CURVATURE (the test c_{k-1} gamma_k >= 0 of the
        !! QR recurrence, Liu and Roosta): solution is x_{k-1},
        !! curvature_direction is M r_{k-1} and residual_norm is |r_{k-1}| from
        !! the recurrence.  For A x = -g at x = 0, M r_{k-1} is a descent
        !! direction of nonpositive curvature (Newton-MR, fortnum_newton_mr).
        !!
        !! A preconditioner with r^T M r < 0, or a non-finite operator or
        !! preconditioner output (the way a callback reports failure), is
        !! KRYLOV_BREAKDOWN with the last iterate in solution.
        class(real_linear_operator_t), intent(inout) :: operator
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        class(real_linear_operator_t), intent(inout), optional :: preconditioner
        real(dp), intent(out), optional :: curvature_direction(:)

        real(dp), allocatable :: v(:), y(:), r1(:), r2(:), w(:), w1(:), w2(:)
        real(dp), allocatable :: residual(:), correction(:)
        real(dp) :: beta, previous_beta, alpha, delta, gamma_bar, epsilon_k
        real(dp) :: previous_epsilon, delta_bar, gamma, cosine, sine
        real(dp) :: phi, phi_bar, right_hand_side_norm, target, recurrence_norm
        integer :: n, k
        logical :: detect

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        n = size(right_hand_side)
        detect = present(curvature_direction)
        if (n < 1) return
        if (size(solution) /= n) return
        if (detect) then
            if (size(curvature_direction) /= n) return
        end if
        if (.not. tolerance > 0.0_dp .or. max_iterations < 1) return

        right_hand_side_norm = parallel_vector_norm2(right_hand_side)
        if (right_hand_side_norm == 0.0_dp) then
            solution = 0.0_dp
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if
        target = tolerance*right_hand_side_norm
        allocate (v(n), y(n), r1(n), r2(n), w(n), w1(n), w2(n), residual(n), &
            correction(n))
        if (any(solution /= 0.0_dp)) then
            call operator%apply(solution, y)
            call parallel_vector_sub(right_hand_side, y, r2)
        else
            r2 = right_hand_side
        end if
        residual = r2
        info = KRYLOV_BREAKDOWN
        if (.not. all(ieee_is_finite(r2))) return
        residual_norm = parallel_vector_norm2(residual)
        if (residual_norm <= target) then
            info = KRYLOV_OK
            return
        end if
        call apply_preconditioner(r2, y)
        beta = parallel_vector_dot(r2, y)
        if (.not. ieee_is_finite(beta) .or. .not. beta > 0.0_dp) return
        beta = sqrt(beta)
        r1 = 0.0_dp
        w = 0.0_dp
        w1 = 0.0_dp
        correction = 0.0_dp
        previous_beta = 0.0_dp
        delta_bar = 0.0_dp
        epsilon_k = 0.0_dp
        phi_bar = beta
        cosine = -1.0_dp
        sine = 0.0_dp

        info = KRYLOV_MAX_ITERATIONS
        do k = 1, max_iterations
            call parallel_vector_scale(1.0_dp/beta, y, v)
            call operator%apply(v, y)
            if (k > 1) call parallel_vector_axpy(-beta/previous_beta, r1, y)
            alpha = parallel_vector_dot(v, y)
            call parallel_vector_axpy(-alpha/beta, r2, y)
            r1 = r2
            r2 = y
            call apply_preconditioner(r2, y)
            previous_beta = beta
            beta = parallel_vector_dot(r2, y)
            if (.not. ieee_is_finite(beta) .or. .not. beta >= 0.0_dp .or. &
                .not. ieee_is_finite(alpha)) then
                info = KRYLOV_BREAKDOWN
                exit
            end if
            beta = sqrt(beta)
            previous_epsilon = epsilon_k
            delta = cosine*delta_bar + sine*alpha
            gamma_bar = sine*delta_bar - cosine*alpha
            epsilon_k = sine*beta
            delta_bar = -cosine*beta
            if (detect .and. cosine*gamma_bar >= 0.0_dp) then
                ! The residual of x_(k-1) has nonpositive curvature.
                call apply_preconditioner(residual, curvature_direction)
                call parallel_vector_axpy(1.0_dp, correction, solution)
                residual_norm = parallel_vector_norm2(residual)
                info = KRYLOV_NONPOSITIVE_CURVATURE
                if (.not. all(ieee_is_finite(curvature_direction))) &
                    info = KRYLOV_BREAKDOWN
                return
            end if
            gamma = sqrt(gamma_bar**2 + beta**2)
            if (.not. gamma > 0.0_dp) then
                info = KRYLOV_BREAKDOWN
                exit
            end if
            cosine = gamma_bar/gamma
            sine = beta/gamma
            phi = cosine*phi_bar
            phi_bar = sine*phi_bar
            w2 = w1
            w1 = w
            w = (v - previous_epsilon*w2 - delta*w1)/gamma
            call parallel_vector_axpy(phi, w, correction)
            iterations = k
            if (beta > 0.0_dp) then
                residual = sine**2*residual - (cosine*phi_bar/beta)*r2
            else
                ! Invariant Krylov space: the iterate solves the system.
                residual = 0.0_dp
            end if
            recurrence_norm = parallel_vector_norm2(residual)
            if (.not. ieee_is_finite(recurrence_norm)) then
                info = KRYLOV_BREAKDOWN
                exit
            end if
            if (recurrence_norm <= target) then
                info = KRYLOV_OK
                exit
            end if
        end do

        call parallel_vector_axpy(1.0_dp, correction, solution)
        call operator%apply(solution, y)
        call parallel_vector_sub(right_hand_side, y, residual)
        residual_norm = parallel_vector_norm2(residual)
        if (info == KRYLOV_OK .and. .not. residual_norm <= target) &
            info = KRYLOV_STAGNATION
    contains
        subroutine apply_preconditioner(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)

            if (present(preconditioner)) then
                call preconditioner%apply(input, output)
            else
                output = input
            end if
        end subroutine apply_preconditioner
    end subroutine real_minres

    subroutine real_minres_operator( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            info, iterations, residual_norm, preconditioner, curvature_direction)
        !! Procedure-argument form of `real_minres` for stateless callbacks.
        procedure(real_matvec_t) :: matvec
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        procedure(real_preconditioner_t), optional :: preconditioner
        real(dp), intent(out), optional :: curvature_direction(:)

        type(procedure_operator_t) :: operator_adapter, preconditioner_adapter

        operator_adapter%apply_procedure => matvec
        if (present(preconditioner)) then
            preconditioner_adapter%apply_procedure => preconditioner
            call real_minres(operator_adapter, right_hand_side, solution, &
                tolerance, max_iterations, info, iterations, residual_norm, &
                preconditioner_adapter, curvature_direction)
        else
            call real_minres(operator_adapter, right_hand_side, solution, &
                tolerance, max_iterations, info, iterations, residual_norm, &
                curvature_direction=curvature_direction)
        end if
    end subroutine real_minres_operator

    subroutine apply_procedure_operator(self, input, output)
        class(procedure_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        call self%apply_procedure(input, output)
    end subroutine apply_procedure_operator

    subroutine real_conjugate_gradient_matmat_operator( &
            matmat, right_hand_side, solution, tolerance, max_iterations, &
            info, iterations, residual_norm, preconditioner)
        !! Solve independent SPD systems while batching operator products.
        !!
        !! Each right-hand side follows its own CG recurrence. Active search
        !! directions are presented to one matrix-matrix callback per
        !! iteration, allowing a kernel backend to reuse pairwise work.
        procedure(real_matmat_t) :: matmat
        real(dp), intent(in) :: right_hand_side(:, :)
        real(dp), intent(inout) :: solution(:, :)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations
        integer, intent(out) :: info(:), iterations(:)
        real(dp), intent(out) :: residual_norm(:)
        procedure(real_matmat_t), optional :: preconditioner

        logical, allocatable :: active(:), candidate(:)
        real(dp), allocatable :: residual(:, :), direction(:, :)
        real(dp), allocatable :: preconditioned(:, :), operator_direction(:, :)
        real(dp), allocatable :: work(:, :)
        real(dp), allocatable :: right_hand_side_norm(:), target_norm(:)
        real(dp), allocatable :: rho(:), next_rho(:), denominator(:)
        real(dp) :: step, beta
        integer :: n_samples, n_rhs, column

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        n_samples = size(right_hand_side, 1)
        n_rhs = size(right_hand_side, 2)
        if (n_samples < 1 .or. n_rhs < 1) return
        if (any(shape(solution) /= [n_samples, n_rhs])) return
        if (size(info) /= n_rhs .or. size(iterations) /= n_rhs .or. &
            size(residual_norm) /= n_rhs) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return

        allocate( &
            residual(n_samples, n_rhs), direction(n_samples, n_rhs), &
            preconditioned(n_samples, n_rhs), &
            operator_direction(n_samples, n_rhs), work(n_samples, n_rhs), &
            right_hand_side_norm(n_rhs), target_norm(n_rhs), rho(n_rhs), &
            next_rho(n_rhs), &
            denominator(n_rhs), active(n_rhs), candidate(n_rhs))

        info = KRYLOV_MAX_ITERATIONS
        iterations = 0
        residual_norm = huge(1.0_dp)
        active = .false.
        candidate = .false.
        do column = 1, n_rhs
            right_hand_side_norm(column) = &
                parallel_vector_norm2(right_hand_side(:, column))
            target_norm(column) = tolerance*max( &
                right_hand_side_norm(column), 1.0_dp)
            if (right_hand_side_norm(column) == 0.0_dp) then
                solution(:, column) = 0.0_dp
            end if
        end do

        call matmat(solution, work)
        residual = right_hand_side - work
        if (present(preconditioner)) then
            call preconditioner(residual, preconditioned)
        else
            preconditioned = residual
        end if
        do column = 1, n_rhs
            residual_norm(column) = parallel_vector_norm2(residual(:, column))
            if (residual_norm(column) <= target_norm(column)) then
                info(column) = KRYLOV_OK
            else
                rho(column) = parallel_vector_dot( &
                    residual(:, column), preconditioned(:, column))
                if (rho(column) <= 0.0_dp .or. &
                    .not. ieee_is_finite(rho(column))) then
                    info(column) = KRYLOV_BREAKDOWN
                else
                    direction(:, column) = preconditioned(:, column)
                    active(column) = .true.
                end if
            end if
        end do

        do while (any(active))
            call matmat(direction, operator_direction)
            candidate = .false.
            do column = 1, n_rhs
                if (active(column)) then
                    denominator(column) = parallel_vector_dot( &
                        direction(:, column), operator_direction(:, column))
                    if (denominator(column) <= 0.0_dp .or. &
                        .not. ieee_is_finite(denominator(column))) then
                        info(column) = KRYLOV_BREAKDOWN
                        active(column) = .false.
                    else
                        step = rho(column)/denominator(column)
                        call parallel_vector_axpy( &
                            step, direction(:, column), solution(:, column))
                        call parallel_vector_axpy( &
                            -step, operator_direction(:, column), &
                            residual(:, column))
                        iterations(column) = iterations(column) + 1
                        residual_norm(column) = &
                            parallel_vector_norm2(residual(:, column))
                        if (residual_norm(column) <= target_norm(column) .or. &
                            iterations(column) >= max_iterations) then
                            candidate(column) = .true.
                        end if
                    end if
                end if
            end do

            if (any(candidate)) then
                call matmat(solution, work)
                do column = 1, n_rhs
                    if (candidate(column)) then
                        call parallel_vector_sub( &
                            right_hand_side(:, column), work(:, column), &
                            residual(:, column))
                        residual_norm(column) = &
                            parallel_vector_norm2(residual(:, column))
                        if (residual_norm(column) <= target_norm(column)) then
                            info(column) = KRYLOV_OK
                            active(column) = .false.
                            direction(:, column) = 0.0_dp
                        else if (iterations(column) >= max_iterations) then
                            info(column) = KRYLOV_MAX_ITERATIONS
                            active(column) = .false.
                            direction(:, column) = 0.0_dp
                        else
                            candidate(column) = .false.
                        end if
                    end if
                end do
            end if

            if (present(preconditioner)) then
                call preconditioner(residual, preconditioned)
            else
                preconditioned = residual
            end if
            do column = 1, n_rhs
                if (active(column)) then
                    next_rho(column) = parallel_vector_dot( &
                        residual(:, column), preconditioned(:, column))
                    if (next_rho(column) <= 0.0_dp .or. &
                        .not. ieee_is_finite(next_rho(column))) then
                        info(column) = KRYLOV_BREAKDOWN
                        active(column) = .false.
                    else
                        beta = next_rho(column)/rho(column)
                        call parallel_vector_xpay( &
                            preconditioned(:, column), beta, direction(:, column))
                        rho(column) = next_rho(column)
                    end if
                end if
            end do
        end do

        call matmat(solution, work)
        do column = 1, n_rhs
            call parallel_vector_sub( &
                right_hand_side(:, column), work(:, column), residual(:, column))
            residual_norm(column) = parallel_vector_norm2(residual(:, column))
            if (residual_norm(column) <= target_norm(column)) then
                info(column) = KRYLOV_OK
            end if
        end do

    end subroutine real_conjugate_gradient_matmat_operator

    subroutine real_gmres( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            restart, info, iterations, residual_norm)
        class(real_linear_operator_t), intent(inout) :: matvec
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm

        real(dp), allocatable :: basis(:, :), givens_cosine(:)
        real(dp), allocatable :: givens_sine(:), hessenberg(:, :)
        real(dp), allocatable :: residual(:), rotated_rhs(:), work(:), y(:)
        real(dp), allocatable :: correction(:)
        real(dp) :: temporary
        real(dp) :: basis_norm, right_hand_side_norm, target
        integer :: cycle_size, inner, previous, used
        logical :: arnoldi_breakdown

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        if (size(right_hand_side) < 1) return
        if (size(solution) /= size(right_hand_side)) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return
        if (restart < 1 .or. restart > max_iterations) return

        allocate( &
            basis(size(solution), restart + 1), &
            hessenberg(restart + 1, restart), residual(size(solution)), &
            work(size(solution)), y(restart), rotated_rhs(restart + 1), &
            correction(restart + 1), &
            givens_cosine(restart), givens_sine(restart))
        right_hand_side_norm = parallel_vector_norm2(right_hand_side)
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = 0.0_dp
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        do while (iterations < max_iterations)
            call matvec%apply(solution, work)
            call parallel_vector_sub(right_hand_side, work, residual)
            residual_norm = parallel_vector_norm2(residual)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if

            basis = 0.0_dp
            hessenberg = 0.0_dp
            rotated_rhs = 0.0_dp
            rotated_rhs(1) = residual_norm
            call parallel_vector_scale(1.0_dp/residual_norm, residual, basis(:, 1))
            cycle_size = min(restart, max_iterations - iterations)
            used = 0
            arnoldi_breakdown = .false.

            do inner = 1, cycle_size
                call matvec%apply(basis(:, inner), work)
                ! Classical Gram-Schmidt with one reorthogonalisation (CGS2).
                call parallel_vector_gemv_t( &
                    basis, inner, work, hessenberg(:, inner))
                call parallel_vector_gemv_n( &
                    -1.0_dp, basis, inner, hessenberg(:, inner), work)
                call parallel_vector_gemv_t(basis, inner, work, correction)
                call parallel_vector_gemv_n( &
                    -1.0_dp, basis, inner, correction, work)
                hessenberg(1:inner, inner) = &
                    hessenberg(1:inner, inner) + correction(1:inner)
                basis_norm = parallel_vector_norm2(work)
                hessenberg(inner + 1, inner) = basis_norm
                if (basis_norm > 64.0_dp*epsilon(1.0_dp)) then
                    call parallel_vector_scale( &
                        1.0_dp/basis_norm, work, basis(:, inner + 1))
                else
                    arnoldi_breakdown = .true.
                end if

                do previous = 1, inner - 1
                    temporary = &
                        givens_cosine(previous)*hessenberg(previous, inner) + &
                        givens_sine(previous)*hessenberg(previous + 1, inner)
                    hessenberg(previous + 1, inner) = &
                        -givens_sine(previous)*hessenberg(previous, inner) + &
                        givens_cosine(previous)* &
                        hessenberg(previous + 1, inner)
                    hessenberg(previous, inner) = temporary
                end do
                call real_givens( &
                    hessenberg(inner, inner), &
                    hessenberg(inner + 1, inner), givens_cosine(inner), &
                    givens_sine(inner), hessenberg(inner, inner))
                hessenberg(inner + 1, inner) = 0.0_dp
                temporary = &
                    givens_cosine(inner)*rotated_rhs(inner) + &
                    givens_sine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner + 1) = &
                    -givens_sine(inner)*rotated_rhs(inner) + &
                    givens_cosine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner) = temporary
                iterations = iterations + 1
                used = inner
                if (abs(rotated_rhs(inner + 1)) <= target) exit
                if (arnoldi_breakdown) exit
            end do

            call real_upper_triangular_solve( &
                hessenberg, rotated_rhs, used, y, info)
            if (info /= KRYLOV_OK) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            call parallel_vector_gemv_n(1.0_dp, basis, used, y, solution)
            call matvec%apply(solution, work)
            call parallel_vector_sub(right_hand_side, work, residual)
            residual_norm = parallel_vector_norm2(residual)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            if (arnoldi_breakdown) then
                info = KRYLOV_BREAKDOWN
                return
            end if
        end do
        info = KRYLOV_MAX_ITERATIONS
    end subroutine real_gmres

    subroutine real_gmres_operator( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            restart, info, iterations, residual_norm)
        procedure(real_matvec_t) :: matvec
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        type(procedure_operator_t) :: operator_adapter

        operator_adapter%apply_procedure => matvec
        call real_gmres(operator_adapter, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm)
    end subroutine real_gmres_operator

    pure subroutine real_givens(first, second, cosine, sine, rotated)
        real(dp), intent(in) :: first, second
        real(dp), intent(out) :: cosine, sine, rotated
        real(dp) :: magnitude

        magnitude = hypot(first, second)
        if (magnitude == 0.0_dp) then
            cosine = 1.0_dp
            sine = 0.0_dp
            rotated = 0.0_dp
        else
            cosine = first/magnitude
            sine = second/magnitude
            rotated = magnitude
        end if
    end subroutine real_givens

    pure subroutine real_upper_triangular_solve( &
            matrix, right_hand_side, size_used, solution, info)
        real(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        integer, intent(in) :: size_used
        real(dp), intent(out) :: solution(:)
        integer, intent(out) :: info
        integer :: column

        solution = 0.0_dp
        info = KRYLOV_BREAKDOWN
        do column = size_used, 1, -1
            if (abs(matrix(column, column)) <= &
                64.0_dp*epsilon(1.0_dp)) return
            solution(column) = ( &
                right_hand_side(column) - &
                dot_product( &
                matrix(column, column + 1:size_used), &
                solution(column + 1:size_used)))/matrix(column, column)
        end do
        info = KRYLOV_OK
    end subroutine real_upper_triangular_solve

    subroutine complex_gmres( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            restart, info, iterations, residual_norm)
        class(complex_linear_operator_t), intent(inout) :: matvec
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm

        complex(dp), allocatable :: basis(:, :), givens_sine(:)
        complex(dp), allocatable :: hessenberg(:, :), residual(:)
        complex(dp), allocatable :: rotated_rhs(:), work(:), y(:)
        real(dp), allocatable :: givens_cosine(:)
        complex(dp) :: coefficient, temporary
        real(dp) :: basis_norm, right_hand_side_norm, target
        integer :: cycle_size, inner, previous, used
        logical :: arnoldi_breakdown

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        if (size(right_hand_side) < 1) return
        if (size(solution) /= size(right_hand_side)) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return
        if (restart < 1 .or. restart > max_iterations) return

        allocate( &
            basis(size(solution), restart + 1), &
            hessenberg(restart + 1, restart), residual(size(solution)), &
            work(size(solution)), y(restart), rotated_rhs(restart + 1), &
            givens_cosine(restart), givens_sine(restart))
        right_hand_side_norm = norm2(abs(right_hand_side))
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = cmplx(0.0_dp, 0.0_dp, dp)
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        do while (iterations < max_iterations)
            call matvec%apply(solution, work)
            residual = right_hand_side - work
            residual_norm = norm2(abs(residual))
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if

            basis = cmplx(0.0_dp, 0.0_dp, dp)
            hessenberg = cmplx(0.0_dp, 0.0_dp, dp)
            rotated_rhs = cmplx(0.0_dp, 0.0_dp, dp)
            rotated_rhs(1) = cmplx(residual_norm, 0.0_dp, dp)
            basis(:, 1) = residual/residual_norm
            cycle_size = min(restart, max_iterations - iterations)
            used = 0
            arnoldi_breakdown = .false.

            do inner = 1, cycle_size
                call matvec%apply(basis(:, inner), work)
                do previous = 1, inner
                    coefficient = dot_product(basis(:, previous), work)
                    hessenberg(previous, inner) = coefficient
                    work = work - coefficient*basis(:, previous)
                end do
                do previous = 1, inner
                    coefficient = dot_product(basis(:, previous), work)
                    hessenberg(previous, inner) = &
                        hessenberg(previous, inner) + coefficient
                    work = work - coefficient*basis(:, previous)
                end do
                basis_norm = norm2(abs(work))
                hessenberg(inner + 1, inner) = &
                    cmplx(basis_norm, 0.0_dp, dp)
                if (basis_norm > 64.0_dp*epsilon(1.0_dp)) then
                    basis(:, inner + 1) = work/basis_norm
                else
                    arnoldi_breakdown = .true.
                end if

                do previous = 1, inner - 1
                    temporary = &
                        givens_cosine(previous)*hessenberg(previous, inner) + &
                        givens_sine(previous)*hessenberg(previous + 1, inner)
                    hessenberg(previous + 1, inner) = &
                        -conjg(givens_sine(previous))* &
                        hessenberg(previous, inner) + &
                        givens_cosine(previous)* &
                        hessenberg(previous + 1, inner)
                    hessenberg(previous, inner) = temporary
                end do
                call complex_givens( &
                    hessenberg(inner, inner), &
                    hessenberg(inner + 1, inner), givens_cosine(inner), &
                    givens_sine(inner), hessenberg(inner, inner))
                hessenberg(inner + 1, inner) = cmplx(0.0_dp, 0.0_dp, dp)
                temporary = &
                    givens_cosine(inner)*rotated_rhs(inner) + &
                    givens_sine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner + 1) = &
                    -conjg(givens_sine(inner))*rotated_rhs(inner) + &
                    givens_cosine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner) = temporary
                iterations = iterations + 1
                used = inner
                if (abs(rotated_rhs(inner + 1)) <= target) exit
                if (arnoldi_breakdown) exit
            end do

            call upper_triangular_solve( &
                hessenberg, rotated_rhs, used, y, info)
            if (info /= KRYLOV_OK) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            solution = solution + matmul(basis(:, :used), y(:used))
            call matvec%apply(solution, work)
            residual_norm = norm2(abs(right_hand_side - work))
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            if (arnoldi_breakdown) then
                info = KRYLOV_BREAKDOWN
                return
            end if
        end do
        info = KRYLOV_MAX_ITERATIONS
    end subroutine complex_gmres

    subroutine complex_gmres_operator( &
            matvec, right_hand_side, solution, tolerance, max_iterations, &
            restart, info, iterations, residual_norm)
        procedure(complex_matvec_t) :: matvec
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        type(complex_procedure_operator_t) :: operator_adapter

        operator_adapter%apply_procedure => matvec
        call complex_gmres(operator_adapter, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm)
    end subroutine complex_gmres_operator

    subroutine real_fgmres( &
            matvec, preconditioner, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio)
        !! Restarted flexible GMRES with right preconditioning.
        !!
        !! Solves A x = b through A M_j^{-1} u_j = v_j, x += M_j^{-1} u_j,
        !! where the callback may apply a different M_j^{-1} each iteration
        !! (for example an inner Krylov solve).  The preconditioned basis is
        !! stored, so the cost is one extra vector per iteration.
        !!
        !! solution is the initial guess and is updated in place.  After each
        !! restart cycle the true residual is compared with the one at the
        !! start of the cycle; if it fell by less than the factor
        !! stagnation_ratio (default 0.99) the routine returns
        !! KRYLOV_STAGNATION instead of looping, keeping the improved iterate.
        class(real_linear_operator_t), intent(inout) :: matvec
        class(real_linear_operator_t), intent(inout) :: preconditioner
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        real(dp), intent(in), optional :: stagnation_ratio

        real(dp), allocatable :: basis(:, :), preconditioned(:, :)
        real(dp), allocatable :: givens_cosine(:), givens_sine(:)
        real(dp), allocatable :: hessenberg(:, :), residual(:)
        real(dp), allocatable :: rotated_rhs(:), work(:), y(:), correction(:)
        real(dp) :: temporary, basis_norm, right_hand_side_norm, target
        real(dp) :: previous_norm, ratio
        integer :: cycle_size, inner, previous, used
        logical :: arnoldi_breakdown

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        if (size(right_hand_side) < 1) return
        if (size(solution) /= size(right_hand_side)) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return
        if (restart < 1 .or. restart > max_iterations) return
        ratio = 0.99_dp
        if (present(stagnation_ratio)) ratio = stagnation_ratio
        if (ratio <= 0.0_dp) return

        allocate( &
            basis(size(solution), restart + 1), &
            preconditioned(size(solution), restart), &
            hessenberg(restart + 1, restart), residual(size(solution)), &
            work(size(solution)), y(restart), rotated_rhs(restart + 1), &
            correction(restart + 1), &
            givens_cosine(restart), givens_sine(restart))
        right_hand_side_norm = parallel_vector_norm2(right_hand_side)
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = 0.0_dp
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        call matvec%apply(solution, work)
        call parallel_vector_sub(right_hand_side, work, residual)
        residual_norm = parallel_vector_norm2(residual)
        do while (iterations < max_iterations)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            previous_norm = residual_norm

            basis = 0.0_dp
            hessenberg = 0.0_dp
            rotated_rhs = 0.0_dp
            rotated_rhs(1) = residual_norm
            call parallel_vector_scale(1.0_dp/residual_norm, residual, basis(:, 1))
            cycle_size = min(restart, max_iterations - iterations)
            used = 0
            arnoldi_breakdown = .false.

            do inner = 1, cycle_size
                call preconditioner%apply(basis(:, inner), preconditioned(:, inner))
                call matvec%apply(preconditioned(:, inner), work)
                ! Classical Gram-Schmidt with one reorthogonalisation (CGS2).
                call parallel_vector_gemv_t( &
                    basis, inner, work, hessenberg(:, inner))
                call parallel_vector_gemv_n( &
                    -1.0_dp, basis, inner, hessenberg(:, inner), work)
                call parallel_vector_gemv_t(basis, inner, work, correction)
                call parallel_vector_gemv_n( &
                    -1.0_dp, basis, inner, correction, work)
                hessenberg(1:inner, inner) = &
                    hessenberg(1:inner, inner) + correction(1:inner)
                basis_norm = parallel_vector_norm2(work)
                hessenberg(inner + 1, inner) = basis_norm
                if (basis_norm > 64.0_dp*epsilon(1.0_dp)) then
                    call parallel_vector_scale( &
                        1.0_dp/basis_norm, work, basis(:, inner + 1))
                else
                    arnoldi_breakdown = .true.
                end if

                do previous = 1, inner - 1
                    temporary = &
                        givens_cosine(previous)*hessenberg(previous, inner) + &
                        givens_sine(previous)*hessenberg(previous + 1, inner)
                    hessenberg(previous + 1, inner) = &
                        -givens_sine(previous)*hessenberg(previous, inner) + &
                        givens_cosine(previous)* &
                        hessenberg(previous + 1, inner)
                    hessenberg(previous, inner) = temporary
                end do
                call real_givens( &
                    hessenberg(inner, inner), &
                    hessenberg(inner + 1, inner), givens_cosine(inner), &
                    givens_sine(inner), hessenberg(inner, inner))
                hessenberg(inner + 1, inner) = 0.0_dp
                temporary = &
                    givens_cosine(inner)*rotated_rhs(inner) + &
                    givens_sine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner + 1) = &
                    -givens_sine(inner)*rotated_rhs(inner) + &
                    givens_cosine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner) = temporary
                iterations = iterations + 1
                used = inner
                if (abs(rotated_rhs(inner + 1)) <= target) exit
                if (arnoldi_breakdown) exit
            end do

            call real_upper_triangular_solve( &
                hessenberg, rotated_rhs, used, y, info)
            if (info /= KRYLOV_OK) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            call parallel_vector_gemv_n(1.0_dp, preconditioned, used, y, solution)
            call matvec%apply(solution, work)
            call parallel_vector_sub(right_hand_side, work, residual)
            residual_norm = parallel_vector_norm2(residual)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            if (arnoldi_breakdown) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            if (residual_norm > ratio*previous_norm) then
                info = KRYLOV_STAGNATION
                return
            end if
        end do
        info = KRYLOV_MAX_ITERATIONS
    end subroutine real_fgmres

    subroutine real_fgmres_operator( &
            matvec, preconditioner, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio)
        !! Restarted flexible GMRES with right preconditioning.
        !!
        !! Solves A x = b through A M_j^{-1} u_j = v_j, x += M_j^{-1} u_j,
        !! where the callback may apply a different M_j^{-1} each iteration
        !! (for example an inner Krylov solve).  The preconditioned basis is
        !! stored, so the cost is one extra vector per iteration.
        !!
        !! solution is the initial guess and is updated in place.  After each
        !! restart cycle the true residual is compared with the one at the
        !! start of the cycle; if it fell by less than the factor
        !! stagnation_ratio (default 0.99) the routine returns
        !! KRYLOV_STAGNATION instead of looping, keeping the improved iterate.
        procedure(real_matvec_t) :: matvec
        procedure(real_preconditioner_t) :: preconditioner
        real(dp), intent(in) :: right_hand_side(:)
        real(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        real(dp), intent(in), optional :: stagnation_ratio
        type(procedure_operator_t) :: operator_adapter, preconditioner_adapter

        operator_adapter%apply_procedure => matvec
        preconditioner_adapter%apply_procedure => preconditioner
        call real_fgmres(operator_adapter, preconditioner_adapter, &
            right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio)
    end subroutine real_fgmres_operator

    subroutine complex_fgmres( &
            matvec, preconditioner, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio, global_sum)
        !! Complex counterpart of real_fgmres_operator (same contract).
        !! global_sum, when present, sums partial inner products of vectors
        !! distributed over processes (e.g. an MPI allreduce); every norm and
        !! inner product of the method passes through it.
        class(complex_linear_operator_t), intent(inout) :: matvec
        class(complex_linear_operator_t), intent(inout) :: preconditioner
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        real(dp), intent(in), optional :: stagnation_ratio
        procedure(complex_global_sum_t), optional :: global_sum

        complex(dp), allocatable :: basis(:, :), preconditioned(:, :)
        complex(dp), allocatable :: givens_sine(:), hessenberg(:, :)
        complex(dp), allocatable :: residual(:), rotated_rhs(:), work(:), y(:)
        complex(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: givens_cosine(:)
        complex(dp) :: temporary
        real(dp) :: basis_norm, right_hand_side_norm, target
        real(dp) :: previous_norm, ratio
        integer :: cycle_size, inner, previous, used, row, pass
        logical :: arnoldi_breakdown

        info = KRYLOV_INVALID_ARGUMENT
        iterations = 0
        residual_norm = huge(1.0_dp)
        if (size(right_hand_side) < 1) return
        if (size(solution) /= size(right_hand_side)) return
        if (tolerance <= 0.0_dp .or. max_iterations < 1) return
        if (restart < 1 .or. restart > max_iterations) return
        ratio = 0.99_dp
        if (present(stagnation_ratio)) ratio = stagnation_ratio
        if (ratio <= 0.0_dp) return

        allocate( &
            basis(size(solution), restart + 1), &
            preconditioned(size(solution), restart), &
            hessenberg(restart + 1, restart), residual(size(solution)), &
            work(size(solution)), y(restart), rotated_rhs(restart + 1), &
            coefficients(restart + 1), &
            givens_cosine(restart), givens_sine(restart))
        right_hand_side_norm = global_norm(right_hand_side)
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = cmplx(0.0_dp, 0.0_dp, dp)
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        call matvec%apply(solution, work)
        residual = right_hand_side - work
        residual_norm = global_norm(residual)
        do while (iterations < max_iterations)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            previous_norm = residual_norm

            basis = cmplx(0.0_dp, 0.0_dp, dp)
            hessenberg = cmplx(0.0_dp, 0.0_dp, dp)
            rotated_rhs = cmplx(0.0_dp, 0.0_dp, dp)
            rotated_rhs(1) = cmplx(residual_norm, 0.0_dp, dp)
            basis(:, 1) = residual/residual_norm
            cycle_size = min(restart, max_iterations - iterations)
            used = 0
            arnoldi_breakdown = .false.

            do inner = 1, cycle_size
                call preconditioner%apply(basis(:, inner), preconditioned(:, inner))
                call matvec%apply(preconditioned(:, inner), work)
                ! Classical Gram-Schmidt with one reorthogonalisation (CGS2),
                ! parallel over basis columns (dots) and rows (updates).
                do pass = 1, 2
                    !$omp parallel do default(shared) private(previous) &
                    !$omp& if(size(work)*inner > 20000)
                    do previous = 1, inner
                        coefficients(previous) = &
                            dot_product(basis(:, previous), work)
                    end do
                    !$omp end parallel do
                    if (present(global_sum)) call global_sum(coefficients(1:inner))
                    !$omp parallel do default(shared) private(row, previous) &
                    !$omp& if(size(work)*inner > 20000)
                    do row = 1, size(work)
                        do previous = 1, inner
                            work(row) = work(row) - &
                                coefficients(previous)*basis(row, previous)
                        end do
                    end do
                    !$omp end parallel do
                    if (pass == 1) then
                        hessenberg(1:inner, inner) = coefficients(1:inner)
                    else
                        hessenberg(1:inner, inner) = &
                            hessenberg(1:inner, inner) + coefficients(1:inner)
                    end if
                end do
                basis_norm = global_norm(work)
                hessenberg(inner + 1, inner) = &
                    cmplx(basis_norm, 0.0_dp, dp)
                if (basis_norm > 64.0_dp*epsilon(1.0_dp)) then
                    basis(:, inner + 1) = work/basis_norm
                else
                    arnoldi_breakdown = .true.
                end if

                do previous = 1, inner - 1
                    temporary = &
                        givens_cosine(previous)*hessenberg(previous, inner) + &
                        givens_sine(previous)*hessenberg(previous + 1, inner)
                    hessenberg(previous + 1, inner) = &
                        -conjg(givens_sine(previous))* &
                        hessenberg(previous, inner) + &
                        givens_cosine(previous)* &
                        hessenberg(previous + 1, inner)
                    hessenberg(previous, inner) = temporary
                end do
                call complex_givens( &
                    hessenberg(inner, inner), &
                    hessenberg(inner + 1, inner), givens_cosine(inner), &
                    givens_sine(inner), hessenberg(inner, inner))
                hessenberg(inner + 1, inner) = cmplx(0.0_dp, 0.0_dp, dp)
                temporary = &
                    givens_cosine(inner)*rotated_rhs(inner) + &
                    givens_sine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner + 1) = &
                    -conjg(givens_sine(inner))*rotated_rhs(inner) + &
                    givens_cosine(inner)*rotated_rhs(inner + 1)
                rotated_rhs(inner) = temporary
                iterations = iterations + 1
                used = inner
                if (abs(rotated_rhs(inner + 1)) <= target) exit
                if (arnoldi_breakdown) exit
            end do

            call upper_triangular_solve( &
                hessenberg, rotated_rhs, used, y, info)
            if (info /= KRYLOV_OK) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            solution = solution + matmul(preconditioned(:, :used), y(:used))
            call matvec%apply(solution, work)
            residual = right_hand_side - work
            residual_norm = global_norm(residual)
            if (residual_norm <= target) then
                info = KRYLOV_OK
                return
            end if
            if (arnoldi_breakdown) then
                info = KRYLOV_BREAKDOWN
                return
            end if
            if (residual_norm > ratio*previous_norm) then
                info = KRYLOV_STAGNATION
                return
            end if
        end do
        info = KRYLOV_MAX_ITERATIONS
    contains
        function global_norm(vector) result(norm)
            complex(dp), intent(in) :: vector(:)
            real(dp) :: norm
            complex(dp) :: partial(1)
            if (.not. present(global_sum)) then
                norm = norm2(abs(vector))
                return
            end if
            partial(1) = cmplx(sum(abs(vector)**2), 0.0_dp, dp)
            call global_sum(partial)
            norm = sqrt(real(partial(1), dp))
        end function global_norm
    end subroutine complex_fgmres

    subroutine complex_fgmres_operator( &
            matvec, preconditioner, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio, global_sum)
        !! Complex counterpart of real_fgmres_operator (same contract).
        procedure(complex_matvec_t) :: matvec
        procedure(complex_preconditioner_t) :: preconditioner
        complex(dp), intent(in) :: right_hand_side(:)
        complex(dp), intent(inout) :: solution(:)
        real(dp), intent(in) :: tolerance
        integer, intent(in) :: max_iterations, restart
        integer, intent(out) :: info, iterations
        real(dp), intent(out) :: residual_norm
        real(dp), intent(in), optional :: stagnation_ratio
        procedure(complex_global_sum_t), optional :: global_sum
        type(complex_procedure_operator_t) :: operator_adapter, preconditioner_adapter

        operator_adapter%apply_procedure => matvec
        preconditioner_adapter%apply_procedure => preconditioner
        call complex_fgmres(operator_adapter, preconditioner_adapter, &
            right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio, global_sum)
    end subroutine complex_fgmres_operator

    pure subroutine complex_givens(first, second, cosine, sine, rotated)
        complex(dp), intent(in) :: first, second
        real(dp), intent(out) :: cosine
        complex(dp), intent(out) :: sine, rotated

        complex(dp) :: phase
        real(dp) :: magnitude, scale

        if (abs(second) == 0.0_dp) then
            cosine = 1.0_dp
            sine = cmplx(0.0_dp, 0.0_dp, dp)
            rotated = first
        else if (abs(first) == 0.0_dp) then
            cosine = 0.0_dp
            sine = conjg(second)/abs(second)
            rotated = cmplx(abs(second), 0.0_dp, dp)
        else
            scale = abs(first) + abs(second)
            magnitude = scale*sqrt( &
                (abs(first)/scale)**2 + (abs(second)/scale)**2)
            phase = first/abs(first)
            cosine = abs(first)/magnitude
            sine = phase*conjg(second)/magnitude
            rotated = phase*magnitude
        end if
    end subroutine complex_givens

    pure subroutine upper_triangular_solve( &
            matrix, right_hand_side, size_used, solution, info)
        complex(dp), intent(in) :: matrix(:, :), right_hand_side(:)
        integer, intent(in) :: size_used
        complex(dp), intent(out) :: solution(:)
        integer, intent(out) :: info

        integer :: column

        solution = cmplx(0.0_dp, 0.0_dp, dp)
        info = KRYLOV_BREAKDOWN
        do column = size_used, 1, -1
            if (abs(matrix(column, column)) <= &
                64.0_dp*epsilon(1.0_dp)) return
            solution(column) = ( &
                right_hand_side(column) - &
                dot_product( &
                conjg(matrix(column, column + 1:size_used)), &
                solution(column + 1:size_used)))/matrix(column, column)
        end do
        info = KRYLOV_OK
    end subroutine upper_triangular_solve

    subroutine apply_complex_procedure_operator(self, input, output)
        class(complex_procedure_operator_t), intent(inout) :: self
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)
        call self%apply_procedure(input, output)
    end subroutine
end module fortnum_krylov
