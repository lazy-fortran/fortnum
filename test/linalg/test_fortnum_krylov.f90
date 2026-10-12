program test_fortnum_krylov
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_krylov, only: complex_gmres_operator, real_gmres_operator, &
        real_conjugate_gradient_operator, &
        real_conjugate_gradient_matmat_operator, KRYLOV_BREAKDOWN, &
        KRYLOV_OK, KRYLOV_INVALID_ARGUMENT, KRYLOV_MAX_ITERATIONS
    use fortnum_linalg, only: dense_solve, LINALG_OK
    implicit none

    complex(dp) :: matrix(5, 5), right_hand_side(5), solution(5)
    complex(dp) :: dense_solution(5), zero_solution(5)
    real(dp) :: real_matrix(5, 5), real_rhs(5), real_solution(5)
    real(dp) :: real_dense_solution(5)
    real(dp) :: diagonal(5), indefinite_matrix(5, 5)
    real(dp) :: multi_rhs(5, 2), multi_solution(5, 2)
    real(dp) :: multi_dense_solution(5, 2), multi_residual_norm(2)
    real(dp) :: residual_norm
    integer, parameter :: stall_size = 2000
    real(dp), allocatable :: stall_coupling(:), stall_diagonal(:)
    real(dp), allocatable :: stall_rhs(:), stall_solution(:)
    real(dp), allocatable :: stall_matrix(:, :), stall_residual(:)
    real(dp), allocatable :: stall_applied(:)
    integer :: column, info, iterations, row, site
    integer :: multi_info(2), multi_iterations(2)

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

    real_matrix = real(matrix, dp)
    real_rhs = real(right_hand_side, dp)
    call dense_solve(real_matrix, real_rhs, real_dense_solution, info)
    call require(info == LINALG_OK, "dense real oracle solves the test system")
    real_solution = 0.0_dp
    call real_gmres_operator( &
        apply_real_matrix, real_rhs, real_solution, 1.0e-12_dp, 30, 3, &
        info, iterations, residual_norm)
    call require(info == KRYLOV_OK, "restarted real GMRES converges")
    call require(maxval(abs(real_solution - real_dense_solution)) < &
        2.0e-11_dp, "real GMRES matches the independent dense LU oracle")
    call require(residual_norm < 1.0e-11_dp, &
        "real GMRES reports the true small residual")

    real_matrix = 0.0_dp
    do row = 1, 5
        real_matrix(row, row) = 2.0_dp + 0.2_dp*row
        if (row < 5) then
            real_matrix(row, row + 1) = -0.4_dp
            real_matrix(row + 1, row) = -0.4_dp
        end if
    end do
    real_rhs = [1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -1.0_dp]
    diagonal = [(real_matrix(row, row), row=1, 5)]
    call dense_solve(real_matrix, real_rhs, real_dense_solution, info)
    call require(info == LINALG_OK, "dense SPD oracle solves the test system")
    real_solution = 0.0_dp
    call real_conjugate_gradient_operator( &
        apply_spd_matrix, real_rhs, real_solution, 1.0e-12_dp, 20, info, &
        iterations, residual_norm, diagonal_preconditioner)
    call require(info == KRYLOV_OK, "preconditioned CG converges")
    call require(iterations > 0 .and. iterations <= 5, &
        "CG uses a bounded number of SPD iterations")
    call require(maxval(abs(real_solution - real_dense_solution)) < 2.0e-11_dp, &
        "CG matches the independent dense LU oracle")
    call require(residual_norm < 1.0e-11_dp, "CG reports the true residual")

    multi_rhs(:, 1) = real_rhs
    multi_rhs(:, 2) = [-0.7_dp, 0.2_dp, 1.5_dp, -1.1_dp, 0.4_dp]
    call dense_solve(real_matrix, multi_rhs, multi_dense_solution, info)
    call require(info == LINALG_OK, &
        "dense multi-RHS oracle solves the test system")
    multi_solution = 0.0_dp
    call real_conjugate_gradient_matmat_operator( &
        apply_spd_matrix_matmat, multi_rhs, multi_solution, 1.0e-12_dp, 20, &
        multi_info, multi_iterations, multi_residual_norm, &
        diagonal_matmat_preconditioner)
    do column = 1, 2
        call require(multi_info(column) == KRYLOV_OK, &
            "batched CG converges for every right-hand side")
        call require(multi_iterations(column) > 0 .and. &
            multi_iterations(column) <= 5, &
            "batched CG reports bounded iterations")
    end do
    call require(maxval(abs(multi_solution - multi_dense_solution)) < &
        2.0e-11_dp, "batched CG matches the independent dense oracle")
    call require(maxval(multi_residual_norm) < 1.0e-11_dp, &
        "batched CG reports true residuals")

    indefinite_matrix = -real_matrix
    real_solution = 0.0_dp
    call real_conjugate_gradient_operator( &
        apply_indefinite_matrix, real_rhs, real_solution, 1.0e-12_dp, 20, &
        info, iterations, residual_norm)
    call require(info == KRYLOV_BREAKDOWN, &
        "CG rejects a non-positive-definite operator")

    call real_conjugate_gradient_operator( &
        apply_spd_matrix, real_rhs, real_solution, -1.0_dp, 20, info, &
        iterations, residual_norm)
    call require(info == KRYLOV_INVALID_ARGUMENT, "CG rejects negative tolerance")

    allocate(stall_coupling(0:stall_size), stall_diagonal(stall_size))
    allocate(stall_rhs(stall_size), stall_solution(stall_size))
    allocate(stall_matrix(stall_size, stall_size), stall_residual(stall_size))
    allocate(stall_applied(stall_size))
    do site = 0, stall_size
        stall_coupling(site) = exp(2.0_dp*sin(0.37_dp*site)**2)
    end do
    do site = 1, stall_size
        stall_diagonal(site) = stall_coupling(site - 1) + &
            stall_coupling(site) + 1.0e-6_dp
        stall_rhs(site) = 1.0_dp + sin(real(site, dp))
    end do
    stall_matrix = 0.0_dp
    do site = 1, stall_size
        stall_matrix(site, site) = stall_diagonal(site)
        if (site < stall_size) then
            stall_matrix(site, site + 1) = -stall_coupling(site)
            stall_matrix(site + 1, site) = -stall_coupling(site)
        end if
    end do
    stall_solution = 0.0_dp
    call real_conjugate_gradient_operator( &
        apply_stall_operator, stall_rhs, stall_solution, 1.0e-14_dp, &
        50000, info, iterations, residual_norm)
    ! The dense matrix is an independent oracle for the true residual.  At the
    ! stall floor the residual is at the roundoff level of the matvec, so the
    ! norm agrees with the operator's own evaluation to rounding only.  The
    ! reported norm is therefore compared with the same operator's recomputed
    ! residual, and the dense oracle checks the relative residual.
    call apply_stall_operator(stall_solution, stall_applied)
    stall_residual = stall_rhs - matmul(stall_matrix, stall_solution)
    call require(info == KRYLOV_MAX_ITERATIONS, &
        "CG below the roundoff floor reports the iteration cap or stall")
    call require(iterations < 50000/2, &
        "CG stops well before the cap once the true residual stalls")
    call require(abs(residual_norm - norm2(stall_rhs - stall_applied)) <= &
        1.0e-10_dp*norm2(stall_rhs - stall_applied), &
        "CG reports the true residual of the returned iterate")
    call require(norm2(stall_residual)/norm2(stall_rhs) < 1.0e-9_dp, &
        "CG returns the best iterate rather than the drifted final iterate")

contains

    subroutine apply_matrix(input, output)
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)

        output = matmul(matrix, input)
    end subroutine apply_matrix

    subroutine apply_real_matrix(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(real_matrix, input)
    end subroutine apply_real_matrix

    subroutine apply_spd_matrix(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(real_matrix, input)
    end subroutine apply_spd_matrix

    subroutine apply_spd_matrix_matmat(input, output)
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)

        output = matmul(real_matrix, input)
    end subroutine apply_spd_matrix_matmat

    subroutine apply_indefinite_matrix(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(indefinite_matrix, input)
    end subroutine apply_indefinite_matrix

    subroutine apply_stall_operator(input, output)
        !! Matrix-free action of the SPD tridiagonal stall_matrix.
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: last

        last = size(input)
        output = stall_diagonal*input
        output(1:last - 1) = output(1:last - 1) - &
            stall_coupling(1:last - 1)*input(2:last)
        output(2:last) = output(2:last) - &
            stall_coupling(1:last - 1)*input(1:last - 1)
    end subroutine apply_stall_operator

    subroutine diagonal_preconditioner(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = input/diagonal
    end subroutine diagonal_preconditioner

    subroutine diagonal_matmat_preconditioner(input, output)
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)

        output = input/spread(diagonal, dim=2, ncopies=size(input, 2))
    end subroutine diagonal_matmat_preconditioner

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            error stop 1
        end if
    end subroutine require

end program test_fortnum_krylov
