program test_complex_bessel_ad
    ! Derivative tests for fortnum_special_complex_bessel (analytic_rule).
    !
    ! The active argument is complex z.  Each function is parametrized along a
    ! complex direction dir by a real step s:  g(s) = B_n(z0 + s*dir), mapped to
    ! R^2 as (Re g, Im g).  The analytic JVP is B_n'(z)*dir from the recurrence
    !   J_n' =  (J_{n-1} - J_{n+1})/2   (DLMF 10.6.1)
    !   I_n' =  (I_{n-1} + I_{n+1})/2   (DLMF 10.29.2)
    !   K_n' = -(K_{n-1} + K_{n+1})/2   (DLMF 10.29.4)
    ! checked against central finite differences, and the JVP/VJP adjoint
    ! identity (the VJP is the transpose of the 2x1 Jacobian).  Scaling flag is
    ! .false.: the unscaled functions are analytic in z (the KODE=1 case the
    ! KiLCA derivatives dIm/dKm use).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    use fortnum_status, only: fortnum_status_t
    use fortnum_special_complex_bessel, only: bessel_j_complex, &
        bessel_i_complex, bessel_k_complex, &
        bessel_j_complex_jvp, bessel_i_complex_jvp, bessel_k_complex_jvp

    implicit none

    real(dp), parameter :: tol_fd  = 1.0e-7_dp
    real(dp), parameter :: tol_adj = 1.0e-12_dp

    ! Base point, direction and order shared with the harness wrappers.
    complex(dp), save :: z0   = (0.0_dp, 0.0_dp)
    complex(dp), save :: dir  = (1.0_dp, 0.0_dp)
    integer,     save :: nord = 0
    integer,     save :: kind_sel = 0   ! 0 J, 1 I, 2 K

    integer :: nfail
    nfail = 0

    call sweep("J", 0, nfail)
    call sweep("I", 1, nfail)
    call sweep("K", 2, nfail)
    call adjoint_checks(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine sweep(label, sel, nfail)
        character(*), intent(in)    :: label
        integer,      intent(in)    :: sel
        integer,      intent(inout) :: nfail
        complex(dp) :: pts(4)
        integer     :: orders(3), i, j
        real(dp)    :: x(1), v(1)
        logical     :: ok

        ! Moderate complex z with Re z > 0 (K needs it); covers KiLCA's domain.
        pts = [ (1.0_dp, 1.0_dp), (2.5_dp, -0.8_dp), &
                (4.0_dp, 1.5_dp), (6.0_dp, -2.0_dp) ]
        orders = [0, 1, 3]
        x = [0.0_dp]
        v = [1.0_dp]
        kind_sel = sel
        do j = 1, size(orders)
            nord = orders(j)
            do i = 1, size(pts)
                z0  = pts(i)
                dir = (0.7_dp, -0.4_dp)
                ok = check_jvp_vs_fd(label//"_jvp", f_primal, f_jvp, x, v, tol_fd)
                if (.not. ok) nfail = nfail + 1
            end do
        end do
    end subroutine sweep

    subroutine adjoint_checks(nfail)
        ! u.(Jv) = v.(J^T u): JVP maps R^1 -> R^2, VJP maps R^2 -> R^1.
        integer, intent(inout) :: nfail
        real(dp) :: x(1), u(2), v(1)
        x = [0.0_dp]
        u = [0.7_dp, -1.3_dp]
        v = [1.0_dp]
        z0 = (3.0_dp, 1.0_dp)
        dir = (0.5_dp, 0.6_dp)
        kind_sel = 1
        nord = 2
        if (.not. dot_product_identity("I_adjoint", f_jvp, f_vjp, x, u, v, &
                tol_adj)) nfail = nfail + 1
        kind_sel = 2
        if (.not. dot_product_identity("K_adjoint", f_jvp, f_vjp, x, u, v, &
                tol_adj)) nfail = nfail + 1
    end subroutine adjoint_checks

    ! Primal g(s) = B_nord(z0 + s*dir) mapped to R^2 = (Re, Im).
    subroutine f_primal(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        complex(dp) :: z, val
        type(fortnum_status_t) :: s
        z = z0 + x(1)*dir
        call eval(z, val, s)
        y(1) = real(val, dp)
        y(2) = aimag(val)
    end subroutine f_primal

    ! Forward product: jv = d/ds (Re,Im) g = (Re,Im)(B'(z)*dir) * v.
    subroutine f_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: z, prod
        type(fortnum_status_t) :: s
        z = z0 + x(1)*dir
        call eval_jvp(z, dir, prod, s)
        jv(1) = real(prod, dp)*v(1)
        jv(2) = aimag(prod)*v(1)
    end subroutine f_jvp

    ! Reverse product: jtu = J^T u, J = [Re(d), Im(d)]^T (2x1).
    subroutine f_vjp(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        complex(dp) :: z, prod
        type(fortnum_status_t) :: s
        z = z0 + x(1)*dir
        call eval_jvp(z, dir, prod, s)
        jtu(1) = u(1)*real(prod, dp) + u(2)*aimag(prod)
    end subroutine f_vjp

    subroutine eval(z, val, s)
        complex(dp),            intent(in)  :: z
        complex(dp),            intent(out) :: val
        type(fortnum_status_t), intent(out) :: s
        select case (kind_sel)
        case (0)
            call bessel_j_complex(nord, z, val, s)
        case (1)
            call bessel_i_complex(nord, z, .false., val, s)
        case default
            call bessel_k_complex(nord, z, .false., val, s)
        end select
    end subroutine eval

    subroutine eval_jvp(z, d, prod, s)
        complex(dp),            intent(in)  :: z, d
        complex(dp),            intent(out) :: prod
        type(fortnum_status_t), intent(out) :: s
        select case (kind_sel)
        case (0)
            call bessel_j_complex_jvp(nord, z, d, prod, s)
        case (1)
            call bessel_i_complex_jvp(nord, z, .false., d, prod, s)
        case default
            call bessel_k_complex_jvp(nord, z, .false., d, prod, s)
        end select
    end subroutine eval_jvp

end program test_complex_bessel_ad
