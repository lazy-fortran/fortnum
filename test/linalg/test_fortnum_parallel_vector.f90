program test_fortnum_parallel_vector
    !! Thread-count invariance, determinism and accuracy of the Krylov vector
    !! kernels, and a CG/GMRES solve that must be bitwise thread-independent.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    !$ use omp_lib, only: omp_set_num_threads
    use fortnum_krylov, only: real_conjugate_gradient_operator, &
        real_gmres_operator, KRYLOV_OK
    use fortnum_parallel_vector, only: parallel_vector_axpy, &
        parallel_vector_dot, parallel_vector_norm2, parallel_vector_gemv_t, &
        parallel_vector_gemv_n, parallel_vector_set_fast
    implicit none

    integer, parameter :: n = 100003, columns = 7
    integer, parameter :: thread_counts(4) = [1, 2, 3, 8]
    integer, parameter :: qp = selected_real_kind(30)
    real(dp), allocatable :: x(:), y(:), basis(:, :), reference_y(:)
    real(dp), allocatable :: solution(:, :), rhs(:)
    real(dp) :: reference_dot, reference_coeff(columns), coeff(columns)
    real(dp) :: value, residual_norm, absolute_sum
    real(qp) :: exact
    integer :: k, i, j, info, iterations, reference_iterations(2)

    allocate (x(n), y(n), reference_y(n), basis(n, columns), rhs(n))
    allocate (solution(n, 2))
    call random_number(x)
    call random_number(y)
    call random_number(basis)
    x = x - 0.5_dp
    y = y - 0.5_dp
    basis = basis - 0.5_dp

    exact = 0.0_qp
    absolute_sum = 0.0_dp
    do i = 1, n
        exact = exact + real(x(i), qp)*real(y(i), qp)
        absolute_sum = absolute_sum + abs(x(i)*y(i))
    end do

    do k = 1, size(thread_counts)
        !$ call omp_set_num_threads(thread_counts(k))
        value = parallel_vector_dot(x, y)
        call parallel_vector_gemv_t(basis, columns, y, coeff)
        reference_y = y
        call parallel_vector_axpy(0.37_dp, x, reference_y)
        call parallel_vector_gemv_n(-1.0_dp, basis, columns, coeff, reference_y)
        if (k == 1) then
            reference_dot = value
            reference_coeff = coeff
            call solve(rhs, solution(:, 1), reference_iterations(1))
            if (abs(real(value, qp) - exact) > 1.0e-13_qp*absolute_sum) &
                call fail("deterministic dot matches the quad-precision sum")
            if (abs(parallel_vector_norm2(x)**2 - sum(real(x, qp)**2)) > &
                1.0e-13_qp*sum(real(x, qp)**2)) &
                call fail("norm2 matches the quad-precision sum")
        else
            if (value /= reference_dot) &
                call fail("dot is bitwise thread-count independent")
            if (any(coeff /= reference_coeff)) &
                call fail("gemv_t is bitwise thread-count independent")
            call solve(rhs, solution(:, 2), reference_iterations(2))
            if (reference_iterations(2) /= reference_iterations(1)) &
                call fail("CG iteration count is thread-count independent")
            if (any(solution(:, 2) /= solution(:, 1))) &
                call fail("CG solution is bitwise thread-count independent")
        end if
    end do

    ! Fast mode agrees to roundoff only.
    call parallel_vector_set_fast(.true.)
    !$ call omp_set_num_threads(8)
    value = parallel_vector_dot(x, y)
    if (abs(value - reference_dot) > 1.0e-12_dp*absolute_sum) &
        call fail("fast dot agrees to roundoff")
    call parallel_vector_set_fast(.false.)

    ! GMRES on the same SPD operator converges and is thread-independent.
    !$ call omp_set_num_threads(1)
    call gmres_solve(solution(:, 1), iterations)
    !$ call omp_set_num_threads(8)
    call gmres_solve(solution(:, 2), reference_iterations(2))
    if (reference_iterations(2) /= iterations) call fail("GMRES iterations thread-independent")
    if (any(solution(:, 2) /= solution(:, 1))) &
        call fail("GMRES solution is bitwise thread-count independent")

    print '(a)', "fortnum_parallel_vector: all checks passed"

contains

    subroutine apply_laplacian(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: m

        m = size(input)
        output = 2.1_dp*input
        output(1:m - 1) = output(1:m - 1) - input(2:m)
        output(2:m) = output(2:m) - input(1:m - 1)
    end subroutine apply_laplacian

    subroutine jacobi(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = input/2.1_dp
    end subroutine jacobi

    subroutine solve(b, u, count)
        real(dp), intent(inout) :: b(:)
        real(dp), intent(inout) :: u(:)
        integer, intent(out) :: count

        b = 1.0_dp
        u = 0.0_dp
        call real_conjugate_gradient_operator( &
            apply_laplacian, b, u, 1.0e-10_dp, 5000, info, count, &
            residual_norm, jacobi)
        if (info /= KRYLOV_OK) call fail("CG converges")
    end subroutine solve

    subroutine gmres_solve(u, count)
        real(dp), intent(inout) :: u(:)
        integer, intent(out) :: count

        rhs = 1.0_dp
        u = 0.0_dp
        call real_gmres_operator( &
            apply_laplacian, rhs, u, 1.0e-10_dp, 600, 40, info, count, &
            residual_norm)
        if (info /= KRYLOV_OK) call fail("GMRES converges")
    end subroutine gmres_solve

    subroutine fail(message)
        character(*), intent(in) :: message

        write (error_unit, '(2a)') "FAIL: ", message
        error stop 1
    end subroutine fail

end program test_fortnum_parallel_vector
