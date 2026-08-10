program test_fortnum_bspline_lsq
    ! Behavioral oracle for matrix-free B-spline CGLS fitting. A cubic
    ! polynomial is represented exactly by an order-four spline, so the
    ! expected values do not depend on the fitting implementation.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_basis
    use fortnum_bspline_lsq, only: bspline_1d_lsq_cgls
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none

    integer, parameter :: n_data = 41
    real(dp), parameter :: atol = 2.0e-10_dp
    type(bspline_workspace_t) :: ws
    type(fortnum_status_t) :: status
    real(dp) :: breakpts(5), x_data(n_data), f_data(n_data)
    real(dp), allocatable :: coeff(:), basis(:), residual(:)
    real(dp) :: x, fitted, worst, normal_residual
    integer :: i, iterations, nfail
    real(dp), parameter :: probe(6) = [0.013_dp, 0.193_dp, 0.347_dp, &
        0.517_dp, 0.887_dp, 1.0_dp]

    nfail = 0
    breakpts = [0.0_dp, 0.22_dp, 0.51_dp, 0.78_dp, 1.0_dp]
    call bspline_init(ws, 4, size(breakpts), status)
    call check_status("init", status, nfail)
    call bspline_set_knots(ws, breakpts, status)
    call check_status("set_knots", status, nfail)

    do i = 1, n_data
        x_data(i) = (real(i, dp) - 0.37_dp)/(real(n_data, dp) - 0.37_dp)
        f_data(i) = cubic(x_data(i))
    end do

    allocate (coeff(ws%ncoef), basis(ws%ncoef), residual(n_data))
    call bspline_1d_lsq_cgls(ws, x_data, f_data, coeff, 500, 1.0e-12_dp, &
        status, iterations, normal_residual)
    if (status%code /= FORTNUM_OK) then
        write (error_unit, "(a,a)") "FAIL [fit status] ", trim(status%msg)
        nfail = nfail + 1
    end if
    if (iterations < 1 .or. iterations > 500) then
        write (error_unit, "(a,i0)") "FAIL [iteration count] ", iterations
        nfail = nfail + 1
    end if
    if (normal_residual > atol) then
        write (error_unit, "(a,es12.4)") &
            "FAIL [residual norm] ", normal_residual
        nfail = nfail + 1
    end if

    worst = 0.0_dp
    do i = 1, size(probe)
        x = probe(i)
        call bspline_eval_basis(ws, x, basis, status)
        call check_status("probe basis", status, nfail)
        fitted = dot_product(coeff, basis)
        worst = max(worst, abs(fitted - cubic(x)))
    end do
    if (worst > atol) then
        write (error_unit, "(a,es12.4)") "FAIL [polynomial oracle] ", worst
        nfail = nfail + 1
    end if

    ! Independently form A^T(f-A*c) to verify the least-squares normal
    ! residual, rather than only checking the routine's reported norm.
    do i = 1, n_data
        call bspline_eval_basis(ws, x_data(i), basis, status)
        residual(i) = f_data(i) - dot_product(coeff, basis)
    end do
    normal_residual = 0.0_dp
    do i = 1, ws%ncoef
        normal_residual = max(normal_residual, abs(sum_normal_column(ws, x_data, residual, i)))
    end do
    if (normal_residual > 5.0e-9_dp) then
        write (error_unit, "(a,es12.4)") &
            "FAIL [normal residual oracle] ", normal_residual
        nfail = nfail + 1
    end if

    call bspline_1d_lsq_cgls(ws, x_data, f_data, coeff, 500, 0.0_dp, status)
    if (status%code /= FORTNUM_DOMAIN_ERROR) then
        write (error_unit, "(a)") "FAIL [invalid tolerance status]"
        nfail = nfail + 1
    end if

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    pure function cubic(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = 1.0_dp + 2.0_dp*x - 0.5_dp*x*x + 0.25_dp*x*x*x
    end function cubic

    function sum_normal_column(workspace, x, values, column) result(value)
        type(bspline_workspace_t), intent(in) :: workspace
        real(dp), intent(in) :: x(:), values(:)
        integer, intent(in) :: column
        real(dp) :: value
        real(dp) :: local_basis(workspace%ncoef)
        type(fortnum_status_t) :: local_status
        integer :: k

        value = 0.0_dp
        do k = 1, size(x)
            call bspline_eval_basis(workspace, x(k), local_basis, local_status)
            value = value + local_basis(column)*values(k)
        end do
    end function sum_normal_column

    subroutine check_status(label, current, failures)
        character(*), intent(in) :: label
        type(fortnum_status_t), intent(in) :: current
        integer, intent(inout) :: failures

        if (current%code /= FORTNUM_OK) then
            write (error_unit, "(a,a,a)") "FAIL [", label, "]"
            failures = failures + 1
        end if
    end subroutine check_status

end program test_fortnum_bspline_lsq
