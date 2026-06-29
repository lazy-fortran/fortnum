program test_hyperg_ad
    ! Derivative tests for fortnum_special_hypergeometric_1f1 (analytic_rule).
    !
    ! The a=1 primal is the analytic map z |-> M(1,b,z), viewed as R^2 -> R^2
    ! over (Re z, Im z).  Its real Jacobian is the Cauchy-Riemann matrix built
    ! from dM/dz = (1/b) M(2,b+1,z) (DLMF 13.3.15).  The harness checks the
    ! forward product against central finite differences and the dot-product
    ! adjoint identity u.(Jv) = v.(J^T u) between the JVP and VJP.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1_a1, &
        hyperg_1f1_a1_jvp, hyperg_1f1_a1_vjp
    implicit none

    real(dp), parameter :: tol_fd  = 1.0e-7_dp ! central-FD (h ~ eps^1/3)
    real(dp), parameter :: tol_adj = 1.0e-13_dp ! adjoint identity

    ! b shared with the harness wrappers (inactive parameter of the map).
    complex(dp), save :: b_shared = (1.0_dp, 0.0_dp)

    integer :: nfail
    nfail = 0

    call test_jvp_vs_fd(nfail)
    call test_adjoint(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    subroutine test_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: z(2), v(2)
        complex(dp) :: bvals(4), zvals(4)
        integer  :: i, j
        logical  :: ok

        bvals = [cmplx(2.0_dp, 0.0_dp, dp), cmplx(1.5_dp, 0.5_dp, dp), &
            cmplx(3.0_dp, -1.0_dp, dp), cmplx(2.0_dp, 1.0_dp, dp)]
        zvals = [cmplx(0.4_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.5_dp, dp), &
            cmplx(2.0_dp, -0.3_dp, dp), cmplx(0.8_dp, 0.8_dp, dp)]

        ! Perturb the real and the imaginary part separately.
        do j = 1, size(bvals)
            b_shared = bvals(j)
            do i = 1, size(zvals)
                z = [real(zvals(i), dp), aimag(zvals(i))]
                v = [1.0_dp, 0.0_dp]
                ok = check_jvp_vs_fd("hyperg_jvp_re", f_primal, jvp_wrap, &
                    z, v, tol_fd)
                if (.not. ok) nfail = nfail + 1
                v = [0.0_dp, 1.0_dp]
                ok = check_jvp_vs_fd("hyperg_jvp_im", f_primal, jvp_wrap, &
                    z, v, tol_fd)
                if (.not. ok) nfail = nfail + 1
            end do
        end do
    end subroutine test_jvp_vs_fd

    subroutine test_adjoint(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: z(2), u(2), v(2)
        b_shared = cmplx(2.0_dp, 0.5_dp, dp)
        z = [1.0_dp, 0.3_dp]
        u = [0.7_dp, -1.1_dp]
        v = [1.3_dp, 0.4_dp]
        if (.not. dot_product_identity("hyperg_adjoint", jvp_wrap, vjp_wrap, &
            z, u, v, tol_adj)) nfail = nfail + 1
        b_shared = cmplx(1.5_dp, -0.4_dp, dp)
        z = [2.0_dp, -0.6_dp]
        u = [-0.5_dp, 0.9_dp]
        v = [0.2_dp, 1.7_dp]
        if (.not. dot_product_identity("hyperg_adjoint_2", jvp_wrap, vjp_wrap, &
            z, u, v, tol_adj)) nfail = nfail + 1
    end subroutine test_adjoint

    ! Primal map (Re z, Im z) |-> (Re M, Im M) for a = 1, b = b_shared.
    subroutine f_primal(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        complex(dp) :: m
        type(fortnum_status_t) :: st
        call hyperg_1f1_a1(b_shared, cmplx(x(1), x(2), dp), m, st)
        y(1) = real(m, dp)
        y(2) = aimag(m)
    end subroutine f_primal

    subroutine jvp_wrap(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        call hyperg_1f1_a1_jvp(x, b_shared, v, jv)
    end subroutine jvp_wrap

    subroutine vjp_wrap(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        call hyperg_1f1_a1_vjp(x, b_shared, u, jtu)
    end subroutine vjp_wrap

end program test_hyperg_ad
