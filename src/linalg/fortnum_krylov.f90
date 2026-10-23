module fortnum_krylov
    !! Restarted matrix-free GMRES and flexible GMRES for real and complex host-side linear operators.
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
    public :: real_fgmres_operator
    public :: complex_matvec_t
    public :: real_conjugate_gradient
    public :: real_conjugate_gradient_operator
    public :: real_conjugate_gradient_matmat_operator
    public :: real_gmres_operator
    public :: real_matmat_t
    public :: real_matvec_t
    public :: real_preconditioner_t

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
            call matvec(solution, work)
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
                call matvec(basis(:, inner), work)
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
            call matvec(solution, work)
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
            call matvec(solution, work)
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
                call matvec(basis(:, inner), work)
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
            call matvec(solution, work)
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
    end subroutine complex_gmres_operator

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

        call matvec(solution, work)
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
                call preconditioner(basis(:, inner), preconditioned(:, inner))
                call matvec(preconditioned(:, inner), work)
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
            call matvec(solution, work)
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
    end subroutine real_fgmres_operator

    subroutine complex_fgmres_operator( &
            matvec, preconditioner, right_hand_side, solution, tolerance, &
            max_iterations, restart, info, iterations, residual_norm, &
            stagnation_ratio)
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
        right_hand_side_norm = norm2(abs(right_hand_side))
        target = tolerance*max(right_hand_side_norm, 1.0_dp)
        if (right_hand_side_norm == 0.0_dp) then
            solution = cmplx(0.0_dp, 0.0_dp, dp)
            residual_norm = 0.0_dp
            info = KRYLOV_OK
            return
        end if

        call matvec(solution, work)
        residual = right_hand_side - work
        residual_norm = norm2(abs(residual))
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
                call preconditioner(basis(:, inner), preconditioned(:, inner))
                call matvec(preconditioned(:, inner), work)
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
            solution = solution + matmul(preconditioned(:, :used), y(:used))
            call matvec(solution, work)
            residual = right_hand_side - work
            residual_norm = norm2(abs(residual))
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

end module fortnum_krylov
