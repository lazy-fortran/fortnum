program test_finite_difference
    ! Proves the central-FD JVP helper on analytic reference kernels whose
    ! Jacobian is known in closed form: a linear map A x (J = A, exact) and a
    ! componentwise quadratic (J diagonal, linear in x). The reference kernels
    ! live here, not in fortnum, so the test validates the harness itself.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, fd_jvp, fd_jvp_step, &
        rel_err
    implicit none

    integer :: nfail
    nfail = 0

    call test_linear_map(nfail)
    call test_quadratic(nfail)
    call test_step_safeguard(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! 3x3 linear map y = A x. Its analytic JVP is J v = A v, identical for any
    ! x; central FD must reproduce A v.
    subroutine test_linear_map(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(3), v(3)
        x = [1.0_dp, -2.0_dp, 0.5_dp]
        v = [0.3_dp, 1.1_dp, -0.7_dp]
        if (.not. check_jvp_vs_fd("linear_map", f_linear, jvp_linear, &
                x, v, tol=1.0e-7_dp)) nfail = nfail + 1
    end subroutine test_linear_map

    ! Componentwise quadratic y_i = x_i^2. J = diag(2 x_i), so JVP = 2 x .* v.
    subroutine test_quadratic(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(4), v(4)
        x = [0.5_dp, 1.5_dp, -3.0_dp, 2.0_dp]
        v = [1.0_dp, -1.0_dp, 0.25_dp, 0.5_dp]
        if (.not. check_jvp_vs_fd("quadratic", f_quad, jvp_quad, &
                x, v, tol=1.0e-6_dp)) nfail = nfail + 1
    end subroutine test_quadratic

    ! The safeguarded step must stay strictly positive and finite even for a
    ! zero base point or a tiny direction.
    subroutine test_step_safeguard(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: h
        h = fd_jvp_step([0.0_dp, 0.0_dp], [1.0e-30_dp, 0.0_dp])
        if (.not. (h > 0.0_dp .and. h == h .and. h < huge(1.0_dp))) then
            write (error_unit, '(a,es12.4)') &
                "FAIL [step_safeguard] degenerate step h=", h
            nfail = nfail + 1
        end if
    end subroutine test_step_safeguard

    ! ----------------------------------------------------- reference kernels

    subroutine f_linear(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = 2.0_dp*x(1) - x(2) + 3.0_dp*x(3)
        y(2) = x(1) + 4.0_dp*x(2)
        y(3) = -x(1) + 0.5_dp*x(2) + 2.0_dp*x(3)
    end subroutine f_linear

    subroutine jvp_linear(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = 2.0_dp*v(1) - v(2) + 3.0_dp*v(3)
        jv(2) = v(1) + 4.0_dp*v(2)
        jv(3) = -v(1) + 0.5_dp*v(2) + 2.0_dp*v(3)
    end subroutine jvp_linear

    subroutine f_quad(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = x*x
    end subroutine f_quad

    subroutine jvp_quad(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv = 2.0_dp*x*v
    end subroutine jvp_quad

end program test_finite_difference
