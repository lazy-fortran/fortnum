module fortnum_phi
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t
    use fortnum_matrix_exponential, only: matrix_exponential
    use fortnum_parallel_vector, only: parallel_vector_norm2, &
        parallel_vector_gemv_t, parallel_vector_gemv_n, &
        PARALLEL_VECTOR_CHUNK
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: PHI_OK = 0, PHI_INVALID_ARGUMENT = -1
    integer, parameter, public :: PHI_NO_CONVERGENCE = 1, PHI_NONFINITE = 2

    type, public :: phi_options_t
        real(dp) :: rtol = 1.0e-12_dp, atol = 0.0_dp
        integer :: min_dimension = 8, max_dimension = 64
        integer :: max_steps = 10000, max_rejections = 100
    end type

    type, public :: phi_report_t
        integer :: matvecs = 0, max_dimension = 0
        integer :: substeps = 0, rejections = 0
        real(dp) :: estimated_error = 0.0_dp, reached_fraction = 0.0_dp
    end type

    type, extends(real_linear_operator_t) :: augmented_operator_t
        class(real_linear_operator_t), pointer :: operator => null()
        real(dp), pointer :: vectors(:, :) => null()
        integer :: n = 0, p = 0
        real(dp) :: tau = 0.0_dp, balance = 1.0_dp
    contains
        procedure :: apply => apply_augmented
    end type

    public :: phi_action

contains
    subroutine phi_action(operator, tau, vectors, result, info, options, report)
        !! Columns 0:p represent sum_k phi_k(tau*A)*vectors(:,k).
        class(real_linear_operator_t), intent(inout), target :: operator
        real(dp), intent(in) :: tau
        real(dp), intent(in), target :: vectors(:, 0:)
        real(dp), intent(out) :: result(:)
        integer, intent(out) :: info
        type(phi_options_t), intent(in), optional :: options
        type(phi_report_t), intent(out), optional :: report
        type(phi_options_t) :: controls
        type(phi_report_t) :: diagnostics
        type(augmented_operator_t) :: augmented
        real(dp), allocatable :: state(:), candidate(:), basis(:, :), hess(:, :)
        real(dp), allocatable :: work(:), correction(:), projection(:, :), e(:, :)
        real(dp) :: beta, original_norm, residual, step, time, error, target
        real(dp) :: initial_norm, factor, local_scale, input_scale
        integer :: n, p, dimension, j, k, pass, dense_info, attempts
        logical :: breakdown, accepted

        info = PHI_INVALID_ARGUMENT
        result = 0.0_dp
        if (present(report)) report = diagnostics
        if (present(options)) controls = options
        n = size(vectors, 1)
        p = size(vectors, 2) - 1
        if (n < 1 .or. p < 0 .or. size(result) /= n) return
        if (.not. ieee_is_finite(tau)) return
        if (.not. all(ieee_is_finite(vectors))) return
        if (.not. ieee_is_finite(controls%rtol)) return
        if (.not. ieee_is_finite(controls%atol)) return
        if (min(controls%rtol, controls%atol) < 0.0_dp) return
        if (max(controls%rtol, controls%atol) <= 0.0_dp) return
        if (controls%min_dimension < 1) return
        if (controls%max_dimension < controls%min_dimension) return
        if (controls%max_steps < 1 .or. controls%max_rejections < 0) return
        dimension = min(n + p, controls%max_dimension)
        allocate (state(n + p), candidate(n + p), basis(n + p, dimension + 1))
        allocate (hess(dimension + 1, dimension), work(n + p))
        allocate (correction(dimension), projection(dimension + 2, dimension + 2))
        allocate (e(dimension + 2, dimension + 2))
        augmented%operator => operator
        augmented%vectors => vectors
        augmented%n = n
        augmented%p = p
        augmented%tau = tau
        input_scale = maxval(abs(vectors))
        if (input_scale <= 0.0_dp) then
            info = PHI_OK
            diagnostics%reached_fraction = 1.0_dp
            if (present(report)) report = diagnostics
            return
        end if
        ! Normalize the augmented state to preserve very small or large inputs.
        augmented%balance = input_scale
        state = 0.0_dp
        state(1:n) = vectors(:, 0)/input_scale
        if (p > 0) state(n + p) = 1.0_dp
        initial_norm = parallel_vector_norm2(state(1:n))
        do k = 1, p
            work(1:n) = vectors(:, k)/input_scale
            initial_norm = max(initial_norm, parallel_vector_norm2(work(1:n)))
        end do
        time = 0.0_dp
        step = 1.0_dp
        attempts = 0
        info = PHI_NO_CONVERGENCE
        do while (time < 1.0_dp)
            if (diagnostics%substeps >= controls%max_steps) exit
            step = min(step, 1.0_dp - time)
            beta = parallel_vector_norm2(state)
            if (.not. ieee_is_finite(beta)) then
                info = PHI_NONFINITE
                exit
            end if
            if (beta <= 0.0_dp) then
                time = 1.0_dp
                exit
            end if
            basis(:, 1) = state/beta
            hess = 0.0_dp
            accepted = .false.
            do j = 1, dimension
                call augmented%apply(basis(:, j), work)
                diagnostics%matvecs = diagnostics%matvecs + 1
                if (.not. all(ieee_is_finite(work))) then
                    info = PHI_NONFINITE
                    exit
                end if
                original_norm = parallel_vector_norm2(work)
                if (.not. ieee_is_finite(original_norm)) then
                    info = PHI_NONFINITE
                    exit
                end if
                ! Reuse the Krylov CGS2 kernels to avoid loss of orthogonality.
                do pass = 1, 2
                    call parallel_vector_gemv_t(basis, j, work, correction)
                    hess(1:j, j) = hess(1:j, j) + correction(1:j)
                    call parallel_vector_gemv_n(-1.0_dp, basis, j, correction, work)
                end do
                residual = parallel_vector_norm2(work)
                hess(j + 1, j) = residual
                breakdown = residual <= 16.0_dp*epsilon(1.0_dp)*original_norm
                if (j == n + p) breakdown = .true.
                if (.not. breakdown) basis(:, j + 1) = work/residual
                if (.not. breakdown .and. j < dimension) then
                    if (j < controls%min_dimension) cycle
                    if (mod(j - controls%min_dimension, 4) /= 0) cycle
                end if
                projection = 0.0_dp
                projection(1:j, 1:j) = step*hess(1:j, 1:j)
                ! Two integrated defect terms guard cancellation of the leading term.
                projection(j + 1, j) = step*residual
                projection(j + 2, j + 1) = step
                call matrix_exponential(projection(1:j + 2, 1:j + 2), &
                    e(1:j + 2, 1:j + 2), dense_info)
                if (dense_info /= 0) then
                    info = PHI_NONFINITE
                    exit
                end if
                candidate = 0.0_dp
                call parallel_vector_gemv_n(beta, basis, j, e(:, 1), candidate)
                error = beta*max(abs(e(j + 1, 1)), abs(e(j + 2, 1)))
                if (breakdown) error = 0.0_dp
                local_scale = max(initial_norm, parallel_vector_norm2(candidate(1:n)))
                target = step*(controls%atol/input_scale + controls%rtol*local_scale)
                diagnostics%max_dimension = max(diagnostics%max_dimension, j)
                if (.not. all(ieee_is_finite(candidate))) then
                    info = PHI_NONFINITE
                    exit
                end if
                if (error <= 0.1_dp*target) then
                    accepted = .true.
                    exit
                end if
            end do
            if (info == PHI_NONFINITE) exit
            if (accepted) then
                state = candidate
                time = min(1.0_dp, time + step)
                diagnostics%substeps = diagnostics%substeps + 1
                diagnostics%estimated_error = &
                    diagnostics%estimated_error + input_scale*error
                diagnostics%reached_fraction = time
                factor = 2.0_dp
                if (error > 0.0_dp) then
                    factor = 0.8_dp*(0.1_dp*target/error)**(1.0_dp/max(1, j))
                    factor = min(2.0_dp, max(1.0_dp, factor))
                end if
                step = step*factor
            else
                diagnostics%rejections = diagnostics%rejections + 1
                attempts = attempts + 1
                if (attempts > controls%max_rejections) exit
                step = step*0.5_dp
                if (time + step <= time) exit
            end if
        end do
        diagnostics%reached_fraction = time
        if (time >= 1.0_dp) info = PHI_OK
        ! A failed evaluation must never masquerade as a value at tau.
        if (info == PHI_OK) then
            if (all(ieee_is_finite(input_scale*state(1:n)))) then
                result = input_scale*state(1:n)
            else
                info = PHI_NONFINITE
            end if
        end if
        if (present(report)) report = diagnostics
    end subroutine

    subroutine apply_augmented(self, input, output)
        class(augmented_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: i, k, n, p
        n = self%n
        p = self%p
        call self%operator%apply(input(1:n), output(1:n))
        !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
        do i = 1, n
            output(i) = self%tau*output(i)
        end do
        !$omp end parallel do simd
        do k = 1, p
            !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
            do i = 1, n
                output(i) = output(i) + &
                    (self%vectors(i, p - k + 1)/self%balance)*input(n + k)
            end do
            !$omp end parallel do simd
        end do
        if (p > 0) then
            output(n + p) = 0.0_dp
            output(n + 1:n + p - 1) = input(n + 2:n + p)
        end if
    end subroutine
end module
