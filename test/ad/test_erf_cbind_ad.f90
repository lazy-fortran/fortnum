program test_erf_cbind_ad
    ! Derivative tests for fortnum_special_erf_cbind (transparent policy).
    !
    !   d/dx erf(x)  =  2/sqrt(pi) exp(-x^2)
    !   d/dx erfc(x) = -2/sqrt(pi) exp(-x^2)
    !
    ! Each analytic JVP is checked against central finite differences, and the
    ! grad (VJP) is checked against the JVP via the dot-product identity, using
    ! the fortnum_ad_test_utils harness.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    use fortnum_special_erf_cbind, only: fortnum_erf, fortnum_erfc, &
        fortnum_erf_jvp, fortnum_erfc_jvp, fortnum_erf_grad, fortnum_erfc_grad
    implicit none

    real(dp), parameter :: tol_fd  = 1.0e-7_dp   ! central-FD tolerance (h ~ eps^1/3)
    real(dp), parameter :: tol_adj = 1.0e-13_dp  ! adjoint identity tolerance

    integer :: nfail
    nfail = 0

    call test_erf_jvp(nfail)
    call test_erfc_jvp(nfail)
    call test_grad_adjoint(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    subroutine test_erf_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: xpts(6), v(1), x(1)
        integer  :: i
        logical  :: ok

        xpts = [-3.0_dp, -0.7_dp, 0.0_dp, 0.5_dp, 1.5_dp, 4.0_dp]
        v    = [1.0_dp]
        do i = 1, size(xpts)
            x(1) = xpts(i)
            ok = check_jvp_vs_fd("fortnum_erf_jvp", f_erf, fortnum_erf_jvp, &
                x, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_erf_jvp

    subroutine test_erfc_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: xpts(6), v(1), x(1)
        integer  :: i
        logical  :: ok

        xpts = [-3.0_dp, -0.7_dp, 0.0_dp, 0.5_dp, 1.5_dp, 4.0_dp]
        v    = [1.0_dp]
        do i = 1, size(xpts)
            x(1) = xpts(i)
            ok = check_jvp_vs_fd("fortnum_erfc_jvp", f_erfc, fortnum_erfc_jvp, &
                x, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_erfc_jvp

    subroutine test_grad_adjoint(nfail)
        ! Scalar output: dot-product identity u*(Jv) = v*(J^T u).
        integer, intent(inout) :: nfail
        real(dp) :: x(1), u(1), v(1)
        x = [1.2_dp]
        u = [0.7_dp]
        v = [1.0_dp]
        if (.not. dot_product_identity("fortnum_erf_grad_adjoint", &
                fortnum_erf_jvp, fortnum_erf_grad, x, u, v, tol_adj)) &
            nfail = nfail + 1
        x = [-0.8_dp]
        u = [-1.3_dp]
        if (.not. dot_product_identity("fortnum_erfc_grad_adjoint", &
                fortnum_erfc_jvp, fortnum_erfc_grad, x, u, v, tol_adj)) &
            nfail = nfail + 1
    end subroutine test_grad_adjoint

    subroutine f_erf(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = fortnum_erf(x(1))
    end subroutine f_erf

    subroutine f_erfc(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = fortnum_erfc(x(1))
    end subroutine f_erfc

end program test_erf_cbind_ad
