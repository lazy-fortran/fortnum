program test_fortnum_krylov
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_krylov, only: complex_gmres_operator, KRYLOV_OK, &
        KRYLOV_INVALID_ARGUMENT
    use fortnum_linalg, only: dense_solve, LINALG_OK
    implicit none

    complex(dp) :: matrix(5, 5), right_hand_side(5), solution(5)
    complex(dp) :: dense_solution(5), zero_solution(5)
    real(dp) :: residual_norm
    integer :: column, info, iterations, row

    matrix = cmplx(0.0_dp, 0.0_dp, dp)
    do column = 1, 5
        matrix(column, column) = cmplx(2.5_dp + 0.2_dp*column, 0.1_dp, dp)
        if (column < 5) then
            matrix(column, column + 1) = cmplx(-0.7_dp, 0.3_dp, dp)
            matrix(column + 1, column) = cmplx(0.15_dp, -0.05_dp, dp)
        end if
    end do
    matrix(1, 5) = cmplx(0.4_dp, -0.2_dp, dp)
    right_hand_side = [ &
        (cmplx(real(row, dp), (-1.0_dp)**row*0.25_dp, dp), row=1, 5)]

    call dense_solve(matrix, right_hand_side, dense_solution, info)
    call require(info == LINALG_OK, "dense oracle solves the test system")
    solution = cmplx(0.0_dp, 0.0_dp, dp)
    call complex_gmres_operator( &
        apply_matrix, right_hand_side, solution, 1.0e-12_dp, 30, 3, &
        info, iterations, residual_norm)
    call require(info == KRYLOV_OK, "restarted complex GMRES converges")
    call require(iterations > 3 .and. iterations <= 30, &
        "restart path performs a bounded number of iterations")
    call require(maxval(abs(solution - dense_solution)) < 2.0e-11_dp, &
        "complex GMRES matches the independent dense LU oracle")
    call require(residual_norm < 1.0e-11_dp, &
        "complex GMRES reports the true small residual")

    zero_solution = cmplx(1.0_dp, -1.0_dp, dp)
    call complex_gmres_operator( &
        apply_matrix, cmplx(0.0_dp, 0.0_dp, dp)*right_hand_side, &
        zero_solution, 1.0e-12_dp, 30, 3, info, iterations, residual_norm)
    call require(info == KRYLOV_OK .and. iterations == 0, &
        "zero right-hand side converges without an iteration")
    call require(maxval(abs(zero_solution)) < 1.0e-14_dp, &
        "zero right-hand side returns the zero solution")

    call complex_gmres_operator( &
        apply_matrix, right_hand_side, solution, -1.0_dp, 30, 3, &
        info, iterations, residual_norm)
    call require(info == KRYLOV_INVALID_ARGUMENT, &
        "negative tolerance is rejected")

contains

    subroutine apply_matrix(input, output)
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)

        output = matmul(matrix, input)
    end subroutine apply_matrix

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            error stop 1
        end if
    end subroutine require

end program test_fortnum_krylov
