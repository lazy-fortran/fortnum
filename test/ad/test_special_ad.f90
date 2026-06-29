program test_special_ad
    ! Derivative tests for fortnum_special (issue #40, analytic_rule).
    !
    ! Checks each analytic derivative against central finite differences using
    ! the fortnum_ad_test_utils harness. For scalar-output functions the
    ! dot-product identity reduces to u*(Jv) = v*(J^T u) with scalars, which
    ! is trivially satisfied when JVP and grad agree; the harness still runs it
    ! via dawson_grad acting as VJP.
    !
    ! gamma_lower_jvp packs x = [x_val, a_val] so the harness can perturb x(1)
    ! while a remains fixed; the primal wrapper unpacks accordingly.
    !
    ! bessel_in_jvp / bessel_kn_jvp take an explicit integer order n; harness-
    ! compatible wrappers close over n via internal procedures.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity, &
                                     fd_jvp, rel_err
    use fortnum_special_dawson, only: dawson, dawson_jvp, dawson_grad
    use fortnum_special_gamma,  only: gamma_lower, gamma_lower_jvp, &
                                      gamma_lower_jvp_da, gamma_reg_p, &
                                      gamma_reg_p_jvp, gamma_reg_p_grad
    use fortnum_special_bessel, only: bessel_in, bessel_kn, bessel_in_array, &
                                      bessel_in_jvp, bessel_kn_jvp, &
                                      bessel_in_array_jvp, bessel_in_array_vjp
    implicit none

    real(dp), parameter :: tol_fd = 1.0e-7_dp   ! central-FD tolerance (h ~ eps^1/3)
    real(dp), parameter :: tol_adj = 1.0e-13_dp  ! adjoint identity tolerance

    ! Shared order variables closed over by the Bessel harness wrappers.
    integer, save :: n_in = 0
    integer, save :: n_kn = 0
    integer, save :: nmax_arr = 0

    integer :: nfail
    nfail = 0

    call test_dawson_jvp(nfail)
    call test_dawson_grad_adjoint(nfail)
    call test_gamma_lower_jvp(nfail)
    call test_gamma_lower_jvp_da(nfail)
    call test_gamma_reg_p_jvp(nfail)
    call test_gamma_reg_p_grad_adjoint(nfail)
    call test_bessel_in_jvp(nfail)
    call test_bessel_kn_jvp(nfail)
    call test_bessel_in_array_jvp(nfail)
    call test_bessel_in_array_adjoint(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! ---------------------------------------------------------------- Dawson

    subroutine test_dawson_jvp(nfail)
        integer, intent(inout) :: nfail
        ! Test F'(x) = 1 - 2*x*F(x) via JVP vs central FD at several points.
        real(dp) :: xpts(5), v(1), x(1)
        integer  :: i
        logical  :: ok

        xpts = [0.1_dp, 0.5_dp, 1.5_dp, 3.0_dp, 11.0_dp]
        v    = [1.0_dp]
        do i = 1, size(xpts)
            x(1) = xpts(i)
            ok = check_jvp_vs_fd("dawson_jvp", f_dawson, dawson_jvp, x, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
        ! negative x
        x(1) = -2.0_dp
        ok = check_jvp_vs_fd("dawson_jvp_neg", f_dawson, dawson_jvp, x, v, tol_fd)
        if (.not. ok) nfail = nfail + 1
    end subroutine test_dawson_jvp

    subroutine test_dawson_grad_adjoint(nfail)
        ! For scalar output dawson, dot-product identity: u*(Jv) = v*(J^T u).
        ! dawson_jvp plays the role of JVP; dawson_grad plays the role of VJP.
        integer, intent(inout) :: nfail
        real(dp) :: x(1), u(1), v(1)
        x = [1.2_dp]
        u = [0.7_dp]
        v = [1.0_dp]
        if (.not. dot_product_identity("dawson_grad_adjoint", &
                dawson_jvp, dawson_grad, x, u, v, tol_adj)) nfail = nfail + 1
        x = [-0.8_dp]
        u = [-1.3_dp]
        if (.not. dot_product_identity("dawson_grad_adjoint_neg", &
                dawson_jvp, dawson_grad, x, u, v, tol_adj)) nfail = nfail + 1
    end subroutine test_dawson_grad_adjoint

    ! ---------------------------------------------------------------- gamma_lower

    subroutine test_gamma_lower_jvp(nfail)
        ! gamma_lower_jvp tests at several (a, x) pairs.
        ! Packing: x_arr(1) = x limit (active), x_arr(2) = a shape (inactive).
        integer, intent(inout) :: nfail
        real(dp) :: x_arr(2), v(2), xpts(4,2)
        integer  :: i
        logical  :: ok

        ! [x_val, a_val] pairs
        xpts(1,:) = [0.5_dp, 1.5_dp]
        xpts(2,:) = [1.0_dp, 2.0_dp]
        xpts(3,:) = [3.0_dp, 2.0_dp]
        xpts(4,:) = [5.0_dp, 3.5_dp]

        ! Perturb only x(1) = x; x(2)=a is inactive so v(2)=0.
        v = [1.0_dp, 0.0_dp]
        do i = 1, size(xpts, 1)
            x_arr = xpts(i,:)
            ok = check_jvp_vs_fd("gamma_lower_jvp", f_gamma_lower, &
                gamma_lower_jvp_wrap, x_arr, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_gamma_lower_jvp

    ! ------------------------------------------------------ gamma_lower d/da

    subroutine test_gamma_lower_jvp_da(nfail)
        ! Full forward product of gamma_lower in [x, a]. Checks d/dx (v=[1,0]),
        ! d/da (v=[0,1]), and a mixed direction (v=[1,1]) against central FD.
        integer, intent(inout) :: nfail
        real(dp) :: x_arr(2), pts(4,2)
        integer  :: i
        logical  :: ok

        pts(1,:) = [0.5_dp, 1.5_dp]
        pts(2,:) = [1.0_dp, 2.0_dp]
        pts(3,:) = [3.0_dp, 2.5_dp]
        pts(4,:) = [5.0_dp, 3.5_dp]

        do i = 1, size(pts, 1)
            x_arr = pts(i,:)
            ok = check_jvp_vs_fd("gamma_lower_jvp_da_dx", f_gamma_lower_pad, &
                gamma_lower_jvp_da, x_arr, [1.0_dp, 0.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
            ok = check_jvp_vs_fd("gamma_lower_jvp_da_da", f_gamma_lower_pad, &
                gamma_lower_jvp_da, x_arr, [0.0_dp, 1.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
            ok = check_jvp_vs_fd("gamma_lower_jvp_da_mix", f_gamma_lower_pad, &
                gamma_lower_jvp_da, x_arr, [1.0_dp, 1.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_gamma_lower_jvp_da

    ! ---------------------------------------------------------- gamma_reg_p

    subroutine test_gamma_reg_p_jvp(nfail)
        ! Forward product of P(a,x) in [x, a]: d/dx, d/da, mixed vs central FD.
        integer, intent(inout) :: nfail
        real(dp) :: x_arr(2), pts(4,2)
        integer  :: i
        logical  :: ok

        pts(1,:) = [0.5_dp, 1.5_dp]
        pts(2,:) = [1.0_dp, 2.0_dp]
        pts(3,:) = [3.0_dp, 2.5_dp]
        pts(4,:) = [5.0_dp, 3.5_dp]

        do i = 1, size(pts, 1)
            x_arr = pts(i,:)
            ok = check_jvp_vs_fd("gamma_reg_p_jvp_dx", f_gamma_reg_p, &
                gamma_reg_p_jvp, x_arr, [1.0_dp, 0.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
            ok = check_jvp_vs_fd("gamma_reg_p_jvp_da", f_gamma_reg_p, &
                gamma_reg_p_jvp, x_arr, [0.0_dp, 1.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
            ok = check_jvp_vs_fd("gamma_reg_p_jvp_mix", f_gamma_reg_p, &
                gamma_reg_p_jvp, x_arr, [1.0_dp, 1.0_dp], tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_gamma_reg_p_jvp

    subroutine test_gamma_reg_p_grad_adjoint(nfail)
        ! Adjoint identity u.(Jv) = v.(J^T u) for the scalar map P(a,x), with
        ! gamma_reg_p_jvp as JVP and gamma_reg_p_grad as VJP. n=2 inputs, m=1.
        integer, intent(inout) :: nfail
        real(dp) :: x_arr(2), u(1), v(2)

        x_arr = [3.0_dp, 2.5_dp]
        u = [0.7_dp]
        v = [0.9_dp, -1.3_dp]
        if (.not. dot_product_identity("gamma_reg_p_grad_adjoint", &
                gamma_reg_p_jvp, gamma_reg_p_grad, x_arr, u, v, tol_adj)) &
            nfail = nfail + 1

        x_arr = [0.5_dp, 1.5_dp]
        u = [-1.1_dp]
        v = [1.4_dp, 0.6_dp]
        if (.not. dot_product_identity("gamma_reg_p_grad_adjoint2", &
                gamma_reg_p_jvp, gamma_reg_p_grad, x_arr, u, v, tol_adj)) &
            nfail = nfail + 1
    end subroutine test_gamma_reg_p_grad_adjoint

    ! ---------------------------------------------------------------- bessel_in

    subroutine test_bessel_in_jvp(nfail)
        integer, intent(inout) :: nfail
        ! Test at several (n, x) combinations.
        real(dp) :: xpts(4), v(1), x(1)
        integer  :: orders(3), i, j
        logical  :: ok

        xpts  = [0.5_dp, 1.5_dp, 4.0_dp, 8.0_dp]
        orders = [0, 1, 3]
        v = [1.0_dp]

        do j = 1, size(orders)
            n_in = orders(j)
            do i = 1, size(xpts)
                x(1) = xpts(i)
                ok = check_jvp_vs_fd("bessel_in_jvp_n" // char(48 + n_in), &
                    f_bessel_in, jvp_bessel_in_wrap, x, v, tol_fd)
                if (.not. ok) nfail = nfail + 1
            end do
        end do
    end subroutine test_bessel_in_jvp

    ! ---------------------------------------------------------------- bessel_kn

    subroutine test_bessel_kn_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: xpts(4), v(1), x(1)
        integer  :: orders(3), i, j
        logical  :: ok

        xpts  = [0.5_dp, 1.5_dp, 4.0_dp, 8.0_dp]
        orders = [0, 1, 3]
        v = [1.0_dp]

        do j = 1, size(orders)
            n_kn = orders(j)
            do i = 1, size(xpts)
                x(1) = xpts(i)
                ok = check_jvp_vs_fd("bessel_kn_jvp_n" // char(48 + n_kn), &
                    f_bessel_kn, jvp_bessel_kn_wrap, x, v, tol_fd)
                if (.not. ok) nfail = nfail + 1
            end do
        end do
    end subroutine test_bessel_kn_jvp

    ! ------------------------------------------------------------ bessel_in_array

    subroutine test_bessel_in_array_jvp(nfail)
        integer, intent(inout) :: nfail
        ! Whole-array column Jacobian vs central FD across all orders 0..nmax.
        real(dp) :: xpts(4), x(1), v(1)
        real(dp), allocatable :: jv_ad(:), jv_fd(:)
        real(dp) :: worst
        integer  :: i, k

        xpts = [0.5_dp, 1.5_dp, 4.0_dp, 8.0_dp]
        nmax_arr = 5
        v = [1.0_dp]
        allocate (jv_ad(nmax_arr + 1), jv_fd(nmax_arr + 1))
        do i = 1, size(xpts)
            x(1) = xpts(i)
            call jvp_bessel_in_array_wrap(x, v, jv_ad)
            call fd_jvp(f_bessel_in_array, x, v, jv_fd)
            worst = 0.0_dp
            do k = 1, nmax_arr + 1
                worst = max(worst, rel_err(jv_ad(k), jv_fd(k)))
            end do
            if (worst > tol_fd) then
                write (error_unit, '(a,es12.4,a,es12.4)') &
                    "FAIL [bessel_in_array_jvp] worst rel_err=", worst, &
                    " tol=", tol_fd
                nfail = nfail + 1
            end if
        end do
    end subroutine test_bessel_in_array_jvp

    subroutine test_bessel_in_array_adjoint(nfail)
        ! Adjoint identity u.(J v) = v.(J^T u) for the array column Jacobian.
        integer, intent(inout) :: nfail
        real(dp) :: x(1), v(1)
        real(dp), allocatable :: u(:)
        integer :: i

        nmax_arr = 5
        allocate (u(nmax_arr + 1))
        v = [1.0_dp]
        do i = 1, nmax_arr + 1
            u(i) = 0.3_dp*real(i, dp) - 0.5_dp
        end do
        x = [1.7_dp]
        if (.not. dot_product_identity("bessel_in_array_adjoint", &
                jvp_bessel_in_array_wrap, vjp_bessel_in_array_wrap, &
                x, u, v, tol_adj)) nfail = nfail + 1
        x = [4.5_dp]
        if (.not. dot_product_identity("bessel_in_array_adjoint2", &
                jvp_bessel_in_array_wrap, vjp_bessel_in_array_wrap, &
                x, u, v, tol_adj)) nfail = nfail + 1
    end subroutine test_bessel_in_array_adjoint

    ! ---------------------------------------------------------------- primal wrappers

    subroutine f_dawson(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = dawson(x(1))
    end subroutine f_dawson

    ! gamma_lower primal: x(1)=x limit, x(2)=a shape.
    subroutine f_gamma_lower(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = gamma_lower(x(2), x(1))
    end subroutine f_gamma_lower

    ! gamma_lower_jvp wrapper: x(1)=x limit (active), x(2)=a (inactive).
    ! Passes directly to the module routine; v(2) must be 0 in the caller.
    subroutine gamma_lower_jvp_wrap(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        call gamma_lower_jvp(x, v, jv)
    end subroutine gamma_lower_jvp_wrap

    ! gamma_lower primal padded into the square harness: y(2)=0 makes the
    ! padded output component a deterministic match for the padded JVP.
    subroutine f_gamma_lower_pad(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = 0.0_dp
        y(1) = gamma_lower(x(2), x(1))
    end subroutine f_gamma_lower_pad

    ! gamma_reg_p primal: x(1)=x limit, x(2)=a shape; padded as above.
    subroutine f_gamma_reg_p(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y = 0.0_dp
        y(1) = gamma_reg_p(x(2), x(1))
    end subroutine f_gamma_reg_p

    subroutine f_bessel_in(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = bessel_in(n_in, x(1))
    end subroutine f_bessel_in

    subroutine jvp_bessel_in_wrap(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        call bessel_in_jvp(n_in, x(1), v(1), jv(1))
    end subroutine jvp_bessel_in_wrap

    subroutine f_bessel_kn(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = bessel_kn(n_kn, x(1))
    end subroutine f_bessel_kn

    subroutine jvp_bessel_kn_wrap(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        call bessel_kn_jvp(n_kn, x(1), v(1), jv(1))
    end subroutine jvp_bessel_kn_wrap

    subroutine f_bessel_in_array(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        real(dp) :: vals(0:nmax_arr)
        call bessel_in_array(nmax_arr, x(1), vals)
        y(1:nmax_arr + 1) = vals(0:nmax_arr)
    end subroutine f_bessel_in_array

    subroutine jvp_bessel_in_array_wrap(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: jvloc(0:nmax_arr)
        call bessel_in_array_jvp(nmax_arr, x(1), v(1), jvloc)
        jv(1:nmax_arr + 1) = jvloc(0:nmax_arr)
    end subroutine jvp_bessel_in_array_wrap

    subroutine vjp_bessel_in_array_wrap(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        real(dp) :: jtuloc
        call bessel_in_array_vjp(nmax_arr, x(1), u(1:nmax_arr + 1), jtuloc)
        jtu(1) = jtuloc
    end subroutine vjp_bessel_in_array_wrap

end program test_special_ad
