module minres_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_krylov, only: real_linear_operator_t
    implicit none
    private

    public :: dense_operator_t, diagonal_preconditioner_t

    type, extends(real_linear_operator_t) :: dense_operator_t
        real(dp), allocatable :: matrix(:, :)
        integer :: applications = 0
    contains
        procedure :: apply => apply_dense
    end type dense_operator_t

    type, extends(real_linear_operator_t) :: diagonal_preconditioner_t
        real(dp), allocatable :: inverse(:)
        logical :: fail = .false.
    contains
        procedure :: apply => apply_diagonal
    end type diagonal_preconditioner_t

contains

    subroutine apply_dense(self, input, output)
        class(dense_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        self%applications = self%applications + 1
        output = matmul(self%matrix, input)
    end subroutine apply_dense

    subroutine apply_diagonal(self, input, output)
        class(diagonal_preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = self%inverse*input
        if (self%fail) output = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine apply_diagonal
end module minres_fixtures

program test_fortnum_minres
    !! Preconditioned MINRES against independent oracles: a dense LU solve
    !! (SPD and indefinite), the minimum-norm solution of a consistent
    !! singular system built from known eigenpairs, the curvature of the
    !! returned nonpositive-curvature direction, and the explicit residual.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_krylov, only: real_minres, real_minres_operator, KRYLOV_OK, &
        KRYLOV_BREAKDOWN, KRYLOV_INVALID_ARGUMENT, KRYLOV_MAX_ITERATIONS, &
        KRYLOV_NONPOSITIVE_CURVATURE
    use fortnum_linalg, only: dense_solve, LINALG_OK
    use minres_fixtures, only: dense_operator_t, diagonal_preconditioner_t
    implicit none

    integer, parameter :: n = 60
    type(dense_operator_t) :: operator
    type(diagonal_preconditioner_t) :: preconditioner
    real(dp) :: q(n, n), lambda(n), b(n), x(n), xd(n), p(n), u(n)
    real(dp) :: residual_norm, residual_norm2, scale
    integer :: info, iterations, info2, iterations2, i, j

    ! Orthogonal Q = I - 2 u u^T / |u|^2 (Householder): the eigenvectors.
    do i = 1, n
        u(i) = sin(1.3_dp*i) + 0.2_dp
        b(i) = cos(0.7_dp*i) + 0.5_dp
    end do
    q = 0.0_dp
    do i = 1, n
        q(i, i) = 1.0_dp
    end do
    do j = 1, n
        q(:, j) = q(:, j) - 2.0_dp*u*u(j)/dot_product(u, u)
    end do

    ! 1. SPD with a wide spectrum and a diagonal (Jacobi) preconditioner.
    do i = 1, n
        lambda(i) = 10.0_dp**(4.0_dp*(i - 1)/(n - 1))
    end do
    call build(lambda)
    do i = 1, n
        operator%matrix(i, i) = operator%matrix(i, i) + 5.0_dp*i
    end do
    call dense_solve(operator%matrix, b, xd, info)
    call require(info == LINALG_OK, "dense oracle solves the SPD system")
    call jacobi()
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-11_dp, 500, info, iterations, &
        residual_norm, preconditioner)
    call require(info == KRYLOV_OK, "preconditioned MINRES converges (SPD)")
    call require(residual_norm <= 1.0e-11_dp*norm2(b), &
        "explicit residual meets the relative tolerance (SPD)")
    call require(norm2(b - matmul(operator%matrix, x)) == residual_norm, &
        "residual_norm is the explicit residual of the solution")
    call require(maxval(abs(x - xd)) < 1.0e-8_dp*maxval(abs(xd)), &
        "MINRES matches the dense oracle (SPD)")

    ! Procedure form runs the same core: identical results.
    xd = 0.0_dp
    call real_minres_operator(apply_matrix, b, xd, 1.0e-11_dp, 500, info2, &
        iterations2, residual_norm2, apply_jacobi)
    call require(info2 == info .and. iterations2 == iterations .and. &
        all(xd == x) .and. residual_norm2 == residual_norm, &
        "procedure MINRES matches operator MINRES bitwise")

    ! Initial guess: the converged solution is returned at once.
    call real_minres(operator, b, x, 1.0e-11_dp, 500, info, iterations, &
        residual_norm, preconditioner)
    call require(info == KRYLOV_OK .and. iterations == 0, &
        "a converged initial guess needs no iteration")

    ! 2. Indefinite nonsingular: eigenvalues of both signs.
    do i = 1, n
        lambda(i) = merge(-1.0_dp, 1.0_dp, mod(i, 3) == 0)*(1.0_dp + 0.5_dp*i)
    end do
    call build(lambda)
    call dense_solve(operator%matrix, b, xd, info)
    call require(info == LINALG_OK, "dense oracle solves the indefinite system")
    call jacobi()
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm, preconditioner)
    call require(info == KRYLOV_OK, "preconditioned MINRES converges (indefinite)")
    call require(maxval(abs(x - xd)) < 1.0e-9_dp*maxval(abs(xd)), &
        "MINRES matches the dense oracle (indefinite)")
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm)
    call require(info == KRYLOV_OK .and. &
        maxval(abs(x - xd)) < 1.0e-9_dp*maxval(abs(xd)), &
        "unpreconditioned MINRES matches the dense oracle (indefinite)")

    ! Monotone residual: a capped solve has |r| below |b| and the cap status.
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 5, info, iterations, &
        residual_norm)
    call require(info == KRYLOV_MAX_ITERATIONS .and. iterations == 5 .and. &
        residual_norm < norm2(b), "capped MINRES returns a reduced residual")

    ! 3. Nonpositive curvature: the direction M r has p^T A p <= 0 and
    !    descends for A x = b at x = 0 (b . p > 0).
    p = 0.0_dp
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm, preconditioner, curvature_direction=p)
    call require(info == KRYLOV_NONPOSITIVE_CURVATURE, &
        "an indefinite operator triggers the curvature exit")
    call require(dot_product(p, matmul(operator%matrix, p)) <= &
        1.0e-12_dp*norm2(p)**2*maxval(abs(lambda)), &
        "curvature direction has nonpositive curvature")
    call require(dot_product(b, p) > 0.0_dp, "curvature direction descends")
    call require(abs(residual_norm - norm2(b - matmul(operator%matrix, x))) &
        <= 1.0e-10_dp*norm2(b), "recurred residual of x_(k-1) is the true one")
    ! An SPD operator never triggers it.
    do i = 1, n
        lambda(i) = 1.0_dp + i
    end do
    call build(lambda)
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm, curvature_direction=p)
    call require(info == KRYLOV_OK, "SPD operator has no curvature exit")

    ! 4. Singular consistent: b in the range, x the minimum-norm solution.
    do i = 1, n
        lambda(i) = merge(0.0_dp, (1.0_dp + i)*merge(-1.0_dp, 1.0_dp, mod(i, 4) == 0), &
            mod(i, 5) == 0)
    end do
    call build(lambda)
    xd = 0.0_dp
    b = 0.0_dp
    do i = 1, n
        if (lambda(i) /= 0.0_dp) then
            scale = cos(0.3_dp*i)
            b = b + scale*q(:, i)
            xd = xd + (scale/lambda(i))*q(:, i)
        end if
    end do
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm)
    call require(info == KRYLOV_OK, "MINRES converges on a consistent singular system")
    call require(maxval(abs(x - xd)) < 1.0e-9_dp*maxval(abs(xd)), &
        "MINRES returns the minimum-norm solution")

    ! 5. Failure and argument checks.
    call jacobi()
    preconditioner%fail = .true.
    x = 0.0_dp
    call real_minres(operator, b, x, 1.0e-12_dp, 500, info, iterations, &
        residual_norm, preconditioner)
    call require(info == KRYLOV_BREAKDOWN, "a failed preconditioner is a breakdown")
    preconditioner%fail = .false.
    call real_minres(operator, b, x, 0.0_dp, 500, info, iterations, residual_norm)
    call require(info == KRYLOV_INVALID_ARGUMENT, "zero tolerance is rejected")
    call real_minres(operator, b(1:n - 1), x, 1.0e-8_dp, 500, info, iterations, &
        residual_norm)
    call require(info == KRYLOV_INVALID_ARGUMENT, "size mismatch is rejected")
    b = 0.0_dp
    x = 1.0_dp
    call real_minres(operator, b, x, 1.0e-8_dp, 500, info, iterations, residual_norm)
    call require(info == KRYLOV_OK .and. all(x == 0.0_dp), &
        "zero right-hand side gives the zero solution")

    print '(a)', "PASS"

contains

    subroutine build(values)
        real(dp), intent(in) :: values(:)
        integer :: k

        if (.not. allocated(operator%matrix)) allocate (operator%matrix(n, n))
        operator%matrix = 0.0_dp
        do k = 1, n
            operator%matrix = operator%matrix + values(k)* &
                spread(q(:, k), 2, n)*spread(q(:, k), 1, n)
        end do
        operator%matrix = 0.5_dp*(operator%matrix + transpose(operator%matrix))
    end subroutine build

    subroutine jacobi()
        integer :: k

        if (allocated(preconditioner%inverse)) deallocate (preconditioner%inverse)
        allocate (preconditioner%inverse(n))
        do k = 1, n
            preconditioner%inverse(k) = 1.0_dp/abs(operator%matrix(k, k))
        end do
    end subroutine jacobi

    subroutine apply_matrix(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(operator%matrix, input)
    end subroutine apply_matrix

    subroutine apply_jacobi(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = preconditioner%inverse*input
    end subroutine apply_jacobi

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//message
            error stop 1
        end if
    end subroutine require
end program test_fortnum_minres
