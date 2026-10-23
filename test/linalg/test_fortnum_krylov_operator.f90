module krylov_operator_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_krylov, only: real_linear_operator_t
    implicit none
    private

    public :: stateful_spd_t, diagonal_inverse_t

    type, extends(real_linear_operator_t) :: stateful_spd_t
        real(dp), allocatable :: matrix(:, :)
        integer :: applications = 0
    contains
        procedure :: apply => apply_stateful_spd
    end type stateful_spd_t

    type, extends(real_linear_operator_t) :: diagonal_inverse_t
        real(dp), allocatable :: inverse(:)
    contains
        procedure :: apply => apply_diagonal_inverse
    end type diagonal_inverse_t

contains

    subroutine apply_stateful_spd(self, input, output)
        class(stateful_spd_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        self%applications = self%applications + 1
        output = matmul(self%matrix, input)
    end subroutine apply_stateful_spd

    subroutine apply_diagonal_inverse(self, input, output)
        class(diagonal_inverse_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = self%inverse*input
    end subroutine apply_diagonal_inverse
end module krylov_operator_fixtures

program test_fortnum_krylov_operator
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_krylov, only: real_conjugate_gradient, &
        real_conjugate_gradient_operator, KRYLOV_BREAKDOWN, KRYLOV_OK, &
        KRYLOV_INVALID_ARGUMENT
    use fortnum_linalg, only: dense_solve, LINALG_OK
    use krylov_operator_fixtures, only: stateful_spd_t, diagonal_inverse_t
    implicit none

    integer, parameter :: n = 40
    type(stateful_spd_t) :: operator
    type(diagonal_inverse_t) :: preconditioner
    real(dp) :: right_hand_side(n), solution(n), legacy_solution(n)
    real(dp) :: dense_solution(n), residual_norm, legacy_residual_norm
    integer :: info, iterations, legacy_info, legacy_iterations, row, column

    ! Diagonally dominant SPD matrix with a widely varying diagonal, so the
    ! Jacobi preconditioner changes the iteration count.
    allocate (operator%matrix(n, n))
    operator%matrix = 0.0_dp
    do row = 1, n
        operator%matrix(row, row) = 2.0_dp + 10.0_dp**(3.0_dp*row/n)
        do column = max(1, row - 2), row - 1
            operator%matrix(row, column) = -0.3_dp/(row - column)
            operator%matrix(column, row) = operator%matrix(row, column)
        end do
        right_hand_side(row) = sin(0.7_dp*row) + 0.5_dp
    end do
    allocate (preconditioner%inverse(n))
    do row = 1, n
        preconditioner%inverse(row) = 1.0_dp/operator%matrix(row, row)
    end do
    call dense_solve(operator%matrix, right_hand_side, dense_solution, info)
    call require(info == LINALG_OK, "dense oracle solves the SPD system")

    solution = 0.0_dp
    call real_conjugate_gradient( &
        operator, right_hand_side, solution, 1.0e-12_dp, 200, info, &
        iterations, residual_norm)
    call require(info == KRYLOV_OK, "operator CG converges")
    call require(maxval(abs(solution - dense_solution)) < 1.0e-10_dp, &
        "operator CG matches the dense oracle")
    call require(operator%applications >= iterations, &
        "operator state is updated by every apply")

    ! The procedure-argument entry point runs the same core, so the results
    ! agree bitwise.
    legacy_solution = 0.0_dp
    call real_conjugate_gradient_operator( &
        apply_matrix, right_hand_side, legacy_solution, 1.0e-12_dp, 200, &
        legacy_info, legacy_iterations, legacy_residual_norm)
    call require(legacy_info == info .and. legacy_iterations == iterations, &
        "procedure CG reports the same status and iteration count")
    call require(all(legacy_solution == solution) .and. &
        legacy_residual_norm == residual_norm, &
        "procedure CG matches operator CG bitwise")

    solution = 0.0_dp
    operator%applications = 0
    call real_conjugate_gradient( &
        operator, right_hand_side, solution, 1.0e-12_dp, 200, info, &
        iterations, residual_norm, preconditioner)
    call require(info == KRYLOV_OK, "preconditioned operator CG converges")
    call require(maxval(abs(solution - dense_solution)) < 1.0e-10_dp, &
        "preconditioned operator CG matches the dense oracle")
    call require(residual_norm < 1.0e-10_dp, &
        "preconditioned operator CG reports the true residual")

    legacy_solution = 0.0_dp
    call real_conjugate_gradient_operator( &
        apply_matrix, right_hand_side, legacy_solution, 1.0e-12_dp, 200, &
        legacy_info, legacy_iterations, legacy_residual_norm, &
        apply_inverse_diagonal)
    call require(legacy_info == info .and. legacy_iterations == iterations, &
        "preconditioned procedure CG reports the same iteration count")
    call require(all(legacy_solution == solution) .and. &
        legacy_residual_norm == residual_norm, &
        "preconditioned procedure CG matches operator CG bitwise")

    operator%matrix = -operator%matrix
    solution = 0.0_dp
    call real_conjugate_gradient( &
        operator, right_hand_side, solution, 1.0e-12_dp, 200, info, &
        iterations, residual_norm)
    call require(info == KRYLOV_BREAKDOWN, &
        "operator CG rejects a non-positive-definite operator")
    operator%matrix = -operator%matrix

    call real_conjugate_gradient( &
        operator, right_hand_side, solution, -1.0_dp, 200, info, &
        iterations, residual_norm)
    call require(info == KRYLOV_INVALID_ARGUMENT, &
        "operator CG rejects a negative tolerance")

contains

    subroutine apply_matrix(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(operator%matrix, input)
    end subroutine apply_matrix

    subroutine apply_inverse_diagonal(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = preconditioner%inverse*input
    end subroutine apply_inverse_diagonal

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//message
            error stop 1
        end if
    end subroutine require
end program test_fortnum_krylov_operator
