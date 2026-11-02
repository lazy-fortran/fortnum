module phi_test_operators
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t
    implicit none
    type, extends(real_linear_operator_t) :: matrix_t
        real(dp), allocatable :: a(:, :)
        integer :: calls = 0
    contains
        procedure :: apply => multiply
    end type
contains
    subroutine multiply(self, input, output)
        class(matrix_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        output = matmul(self%a, input)
        self%calls = self%calls + 1
    end subroutine
end module

program test_fortnum_phi
    use fortnum_kinds, only: dp
    use fortnum_phi, only: phi_action, phi_options_t, phi_report_t, &
        PHI_OK, PHI_INVALID_ARGUMENT, PHI_NO_CONVERGENCE, &
        PHI_NONFINITE
    use phi_test_operators, only: matrix_t
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: qp = selected_real_kind(30)
    type(matrix_t) :: operator
    type(phi_options_t) :: options
    type(phi_report_t) :: report
    real(dp), allocatable :: v(:, :), actual(:), expected(:)
    real(dp) :: tau, error, started, finished
    integer :: n, k, i, regime, info

    call cpu_time(started)
    options%rtol = 2.0e-14_dp
    ! The full Dirichlet spectrum is independent of the Arnoldi projection.
    n = 96
    allocate (operator%a(n, n), v(n, 0:3), actual(n), expected(n))
    operator%a = 0.0_dp
    do i = 1, n
        operator%a(i, i) = -2.0_dp
        if (i < n) then
            operator%a(i + 1, i) = 1.0_dp
            operator%a(i, i + 1) = 1.0_dp
        end if
    end do
    do regime = 1, 3
        tau = 0.000025_dp
        if (regime == 2) tau = 0.25_dp
        if (regime == 3) tau = 2500.0_dp
        do k = 0, 3
            v = 0.0_dp
            do i = 1, n
                v(i, k) = 1.0_dp + cos(0.37_dp*i) + sin(0.83_dp*i)
            end do
            call dirichlet_reference(tau, v, expected)
            call phi_action(operator, tau, v, actual, info, options, report)
            call require(info == PHI_OK, 'Dirichlet phi converges')
            error = norm2(actual - expected)/norm2(expected)
            write (*, '(a,i1,a,es10.2,a,es10.2,a,i0,a,i0)') &
                'Dirichlet k=', k, ' norm1(tau*A)=', 4*tau, &
                ' relative error=', error, ' m=', report%max_dimension, &
                ' substeps=', report%substeps
            call require(error <= 1.0e-10_dp, 'Dirichlet relative error')
        end do
    end do
    do k = 0, 3
        do i = 1, n
            v(i, k) = cos(real((k + 1)*i, dp))
        end do
    end do
    call dirichlet_reference(100.0_dp, v, expected)
    call phi_action(operator, 100.0_dp, v, actual, info, options, report)
    call require(info == PHI_OK, 'mixed phi combination converges')
    call require(norm2(actual - expected)/norm2(expected) <= 1.0e-10_dp, &
        'mixed phi combination matches eigenpairs')
    call require(report%max_dimension > options%min_dimension, &
        'Krylov dimension increases for a difficult action')
    call require(report%matvecs > 0, 'operator work is reported')

    ! Force substepping with a small basis, retaining the same spectral oracle.
    options%max_dimension = 12
    options%max_rejections = 1000
    call phi_action(operator, 100.0_dp, v, actual, info, options, report)
    call require(info == PHI_OK, 'small basis converges with substeps')
    call require(report%substeps > 1 .and. report%rejections > 0, &
        'small basis exercises rejection and substepping')
    call require(report%max_dimension <= 12, 'dimension budget is respected')
    call require(norm2(actual - expected)/norm2(expected) <= 1.0e-10_dp, &
        'substepped combination matches eigenpairs')
    write (*, '(a,i0,a,i0,a,es10.2)') 'forced substeps=', report%substeps, &
        ' rejections=', report%rejections, ' estimated error=', &
        report%estimated_error
    options = phi_options_t(rtol=2.0e-14_dp)
    deallocate (operator%a, v, actual, expected)

    ! Upwind convection-diffusion, driven by a slowly decaying boundary state.
    n = 18
    allocate (operator%a(n, n), v(n, 0:3), actual(n), expected(n))
    operator%a = 0.0_dp
    operator%a(1, 1) = -1.0e-5_dp
    do i = 2, n
        operator%a(i, i) = -2.1_dp
        operator%a(i, i - 1) = 1.8_dp
        if (i < n) operator%a(i, i + 1) = 0.2_dp
    end do
    do regime = 1, 3
        tau = 10.0_dp**(2*regime - 3)
        if (regime == 3) tau = 1.0e4_dp/maxval(sum(abs(operator%a), dim=1))
        do k = 0, 3
            v = 0.0_dp
            do i = 1, n
                v(i, k) = sin(0.4_dp*i) + 0.5_dp
            end do
            call dense_reference(operator%a, tau, v, expected)
            call phi_action(operator, tau, v, actual, info, options, report)
            call require(info == PHI_OK, 'nonnormal phi converges')
            error = norm2(actual - expected)/norm2(expected)
            write (*, '(a,i1,a,es10.2,a,es10.2,a,i0)') 'nonnormal k=', k, &
                ' norm1(tau*A)=', tau*maxval(sum(abs(operator%a), dim=1)), &
                ' relative error=', error, ' m=', report%max_dimension
            call require(error <= 1.0e-10_dp, 'nonnormal relative error')
        end do
    end do
    do k = 0, 3
        do i = 1, n
            v(i, k) = sin(real(i + k, dp))
        end do
    end do
    call dense_reference(operator%a, -0.2_dp, v, expected)
    call phi_action(operator, -0.2_dp, v, actual, info, options, report)
    call require(info == PHI_OK, 'negative tau converges')
    call require(norm2(actual - expected)/norm2(expected) < 1.0e-10_dp, &
        'negative tau matches the independent dense reference')

    ! A rotating block catches cancellation in defect-based error estimates.
    operator%a = 0.0_dp
    do i = 1, n, 2
        operator%a(i, i + 1) = 3.0_dp*i
        operator%a(i + 1, i) = -3.0_dp*i
    end do
    call dense_reference(operator%a, 1.0_dp, v, expected)
    call phi_action(operator, 1.0_dp, v, actual, info, options, report)
    call require(info == PHI_OK, 'oscillatory phi converges')
    call require(norm2(actual - expected)/norm2(expected) < 1.0e-10_dp, &
        'oscillatory phi matches the dense Taylor oracle')

    ! Invalid arguments and bounded failure return no purported endpoint.
    options%min_dimension = 2
    options%max_dimension = 2
    options%max_rejections = 0
    call phi_action(operator, 10.0_dp, v, actual, info, options, report)
    call require(info == PHI_NO_CONVERGENCE, 'budget exhaustion is visible')
    call require(maxval(abs(actual)) <= 0.0_dp, 'failed action does not return a value')
    call require(report%reached_fraction < 1.0_dp, 'failure reports partial time')
    options%max_rejections = 100
    options%max_steps = 1
    call phi_action(operator, 0.1_dp, v, actual, info, options, report)
    call require(info == PHI_NO_CONVERGENCE, 'substep budget is enforced')
    options%rtol = -1.0_dp
    call phi_action(operator, 1.0_dp, v, actual, info, options)
    call require(info == PHI_INVALID_ARGUMENT, 'negative tolerance is rejected')
    options = phi_options_t()
    options%rtol = ieee_value(0.0_dp, ieee_quiet_nan)
    call phi_action(operator, 1.0_dp, v, actual, info, options)
    call require(info == PHI_INVALID_ARGUMENT, 'NaN tolerance is rejected')
    options = phi_options_t()
    call phi_action(operator, 1.0_dp, v, actual(:n - 1), info, options)
    call require(info == PHI_INVALID_ARGUMENT, 'shape mismatch is rejected')
    operator%a(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call phi_action(operator, 1.0_dp, v, actual, info, options)
    call require(info == PHI_NONFINITE, 'nonfinite operator action is rejected')
    v(1, 0) = ieee_value(0.0_dp, ieee_quiet_nan)
    call phi_action(operator, 1.0_dp, v, actual, info, options)
    call require(info == PHI_INVALID_ARGUMENT, 'nonfinite input is rejected')
    operator%a = 0.0_dp
    v = 0.0_dp
    call phi_action(operator, 1.0_dp, v, actual, info, report=report)
    call require(info == PHI_OK .and. maxval(abs(actual)) <= 0.0_dp, 'zero combination')
    call require(report%matvecs == 0, 'zero combination needs no operator call')
    do k = 0, 3
        v(:, k) = real(k + 1, dp)
    end do
    expected = v(:, 0) + v(:, 1) + v(:, 2)/2.0_dp + v(:, 3)/6.0_dp
    call phi_action(operator, 1.0_dp, v, actual, info, options)
    call require(info == PHI_OK, 'zero operator converges')
    call require(norm2(actual - expected) < 1.0e-12_dp, 'phi_k(0)=1/k!')
    call phi_action(operator, 1.0_dp, v(:, 0:0), actual, info, options)
    call require(info == PHI_OK, 'phi0-only input converges')
    call require(norm2(actual - v(:, 0)) < 1.0e-12_dp, 'exp(0) is identity')
    ! Normalizing the augmented state must retain the relative accuracy contract.
    operator%a = 0.0_dp
    do i = 1, n
        operator%a(i, i) = -1.0_dp
    end do
    do regime = 1, 2
        tau = 1.0e-200_dp
        if (regime == 2) tau = 1.0e200_dp
        v = tau
        expected = tau*(exp(-1.0_dp) + (1.0_dp - exp(-1.0_dp)) + &
            exp(-1.0_dp) + (0.5_dp - exp(-1.0_dp)))
        call phi_action(operator, 1.0_dp, v, actual, info, options)
        call require(info == PHI_OK, 'scaled input converges')
        call require(maxval(abs((actual - expected)/tau)) < 1.0e-12_dp, &
            'scaled input retains relative precision')
    end do
    v = 1.0_dp
    call phi_action(operator, 0.0_dp, v, actual, info, options)
    call require(info == PHI_OK, 'zero time with a nonzero operator')
    call require(maxval(abs(actual - 8.0_dp/3.0_dp)) < 1.0e-12_dp, &
        'zero time retains every phi term')
    call phi_action(operator, 0.5_dp, v(:, 0:0), actual, info, options)
    call require(info == PHI_OK, 'phi0 without augmentation converges')
    call require(maxval(abs(actual - exp(-0.5_dp))) < 1.0e-12_dp, &
        'phi0 without augmentation matches an exact exponential')
    ! Strided columns exercise the public assumed-shape input contract.
    call phi_action(operator, 0.5_dp, v(:, 0:3:2), actual, info, options)
    call require(info == PHI_OK, 'strided columns converge')
    call require(maxval(abs(actual - (2.0_dp - exp(-0.5_dp)))) < &
        1.0e-12_dp, 'strided columns match phi0 plus phi1')
    call cpu_time(finished)
    write (*, '(a,f9.3,a)') 'phi tests CPU time: ', finished - started, ' s'
contains
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            write (*, '(a)') 'FAIL: '//message
            error stop 1
        end if
    end subroutine

    subroutine dirichlet_reference(tau, vectors, value)
        real(dp), intent(in) :: tau, vectors(:, 0:)
        real(dp), intent(out) :: value(:)
        real(qp) :: eigenvector(size(value)), answer(size(value))
        real(qp) :: eigenvalue, coefficient
        integer :: mode, k, i, n
        n = size(value)
        answer = 0.0_qp
        do mode = 1, n
            do i = 1, n
                eigenvector(i) = sqrt(2.0_qp/(n + 1))* &
                    sin(acos(-1.0_qp)*i*mode/(n + 1))
            end do
            eigenvalue = -4.0_qp*sin(acos(-1.0_qp)*mode/(2*(n + 1)))**2
            do k = 0, ubound(vectors, 2)
                coefficient = sum(eigenvector*real(vectors(:, k), qp))
                answer = answer + eigenvector*coefficient* &
                    scalar_phi(real(tau, qp)*eigenvalue, k)
            end do
        end do
        value = real(answer, dp)
    end subroutine

    function scalar_phi(z, k) result(value)
        real(qp), intent(in) :: z
        integer, intent(in) :: k
        real(qp) :: value, term, factorial
        integer :: j
        factorial = 1.0_qp
        do j = 2, k
            factorial = factorial*j
        end do
        if (abs(z) <= 1.0_qp) then
            term = 1.0_qp/factorial
            value = term
            do j = 1, 100
                term = term*z/(j + k)
                value = value + term
                if (abs(term) < 1.0e-32_qp*abs(value)) exit
            end do
        else
            value = exp(z)
            factorial = 1.0_qp
            do j = 1, k
                value = (value - 1.0_qp/factorial)/z
                factorial = factorial*j
            end do
        end if
    end function

    subroutine dense_reference(a, tau, vectors, value)
        real(dp), intent(in) :: a(:, :), tau, vectors(:, 0:)
        real(dp), intent(out) :: value(:)
        real(qp), allocatable :: augmented(:, :), exponential(:, :), term(:, :)
        real(qp), allocatable :: initial(:), answer(:), scratch(:, :)
        real(qp) :: anorm
        integer :: n, p, order, i, j, squarings
        n = size(a, 1)
        p = ubound(vectors, 2)
        order = n + p
        allocate (augmented(order, order), exponential(order, order))
        allocate (term(order, order), initial(order), answer(order))
        allocate (scratch(order, order))
        augmented = 0.0_qp
        augmented(1:n, 1:n) = real(tau, qp)*real(a, qp)
        do j = 1, p
            augmented(1:n, n + j) = real(vectors(:, p - j + 1), qp)
            if (j < p) augmented(n + j, n + j + 1) = 1.0_qp
        end do
        initial = 0.0_qp
        initial(1:n) = real(vectors(:, 0), qp)
        if (p > 0) initial(order) = 1.0_qp
        anorm = maxval(sum(abs(augmented), dim=1))
        squarings = max(0, ceiling(log(max(anorm, 0.125_qp)/0.125_qp)/log(2.0_qp)))
        augmented = scale(augmented, -squarings)
        term = 0.0_qp
        do i = 1, order
            term(i, i) = 1.0_qp
        end do
        exponential = term
        do j = 1, 100
            scratch = matmul(augmented, term)
            term = scratch/j
            exponential = exponential + term
            if (maxval(abs(term)) < 1.0e-33_qp) exit
        end do
        call require(j <= 100, 'quadruple Taylor reference converges')
        do j = 1, squarings
            scratch = matmul(exponential, exponential)
            exponential = scratch
        end do
        answer = matmul(exponential, initial)
        value = real(answer(1:n), dp)
    end subroutine
end program
