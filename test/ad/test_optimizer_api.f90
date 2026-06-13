program test_optimizer_api
    ! Optimizer-facing API tests (issue #41). Exercises the backend-opaque
    ! derivative interfaces of fortnum_ad_interfaces and the flat active-vector
    ! layout of fortnum_active_vector, then verifies the products with the
    ! shared harness. Backend-agnostic: builds and passes WITHOUT Enzyme.
    !
    ! (a) Scalar objective gradient through a grad_fn kernel, checked vs central
    !     finite difference. The objective is a Rosenbrock-style smooth function
    !     of a flat active vector whose two named blocks ("a", "b") are unpacked
    !     inside the kernel context.
    ! (b) Vector residual JVP/VJP through jvp_fn/vjp_fn kernels, checked by the
    !     mandatory dot-product adjoint identity u.(Jv) = v.(J^T u).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_kinds,  only: dp
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_OK
    use fortnum_active_vector, only: fortnum_active_layout_t, layout_init, &
        layout_add, unpack_block
    use fortnum_ad_interfaces, only: fortnum_ad_status_t, ad_status_set, &
        ad_status_ok, FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    implicit none

    integer :: nfail
    nfail = 0

    call test_layout_roundtrip(nfail)
    call test_scalar_grad(nfail)
    call test_residual_adjoint(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! Scaled relative error helper.
    pure real(dp) function rel_err(got, want)
        real(dp), intent(in) :: got, want
        rel_err = abs(got - want) / max(abs(want), 1.0_dp)
    end function rel_err

    ! Builds the two-block layout used by both kernels: scalar "a", scalar "b".
    subroutine build_layout(layout)
        type(fortnum_active_layout_t), intent(out) :: layout
        type(fortnum_status_t) :: st
        call layout_init(layout, 2)
        call layout_add(layout, "a", 1, st)
        call layout_add(layout, "b", 1, st)
    end subroutine build_layout

    ! Test: pack/unpack round-trips through named blocks and reports n.
    subroutine test_layout_roundtrip(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_active_layout_t) :: layout
        type(fortnum_status_t) :: st
        real(dp) :: x(2), va(1), vb(1)
        logical  :: ok

        call build_layout(layout)
        x = [3.0_dp, -7.0_dp]
        call unpack_block(layout, x, "a", va, st)
        ok = status_ok(st) .and. (va(1) == 3.0_dp)
        call unpack_block(layout, x, "b", vb, st)
        ok = ok .and. status_ok(st) .and. (vb(1) == -7.0_dp) .and. layout%n == 2
        if (.not. ok) then
            write (error_unit, '(a)') "FAIL [layout] roundtrip/n"
            nfail = nfail + 1
        end if
    end subroutine test_layout_roundtrip

    ! (a) Scalar objective gradient via a grad_fn kernel, verified vs central
    !     FD. Objective f(a,b) = (1-a)^2 + 100 (b - a^2)^2 (Rosenbrock), with a
    !     and b unpacked from the flat vector inside the kernel context.
    subroutine test_scalar_grad(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_active_layout_t) :: layout
        type(fortnum_ad_status_t) :: st
        real(dp) :: x(2), f, g(2), gfd(2), h, worst
        integer  :: i
        logical  :: ok

        call build_layout(layout)
        x = [0.7_dp, 0.4_dp]
        call rosen_grad(2, x, f, g, layout, st)
        ok = ad_status_ok(st) .and. st%backend == FORTNUM_AD_BACKEND_ANALYTIC &
            .and. st%quality == FORTNUM_AD_QUALITY_EXACT

        ! Central FD of each gradient component, kernel as a black box.
        h = 1.0e-6_dp
        do i = 1, 2
            gfd(i) = fd_component(layout, x, i, h)
        end do
        worst = 0.0_dp
        do i = 1, 2
            worst = max(worst, rel_err(g(i), gfd(i)))
        end do
        if (worst > 1.0e-6_dp) then
            write (error_unit, '(a,es12.4)') &
                "FAIL [grad] vs FD worst rel_err=", worst
            ok = .false.
        end if
        if (.not. ok) nfail = nfail + 1
    end subroutine test_scalar_grad

    ! Central FD of the objective along the i-th coordinate using the value
    ! produced by the grad_fn kernel (f output), treating it as opaque.
    real(dp) function fd_component(layout, x, i, h) result(d)
        type(fortnum_active_layout_t), intent(in) :: layout
        real(dp), intent(in) :: x(:)
        integer,  intent(in) :: i
        real(dp), intent(in) :: h
        type(fortnum_ad_status_t) :: st
        real(dp) :: xp(2), xm(2), fp, fm, g(2)
        xp = x; xm = x
        xp(i) = xp(i) + h
        xm(i) = xm(i) - h
        call rosen_grad(2, xp, fp, g, layout, st)
        call rosen_grad(2, xm, fm, g, layout, st)
        d = (fp - fm) / (2.0_dp*h)
    end function fd_component

    ! grad_fn kernel: Rosenbrock objective over the flat active vector. Unpacks
    ! a and b through the layout context, computes f and the analytic gradient.
    subroutine rosen_grad(n, x, f, g, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(n)
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        type(fortnum_status_t) :: us
        real(dp) :: va(1), vb(1), a, b
        call unpack_block(layout, x, "a", va, us)
        call unpack_block(layout, x, "b", vb, us)
        a = va(1); b = vb(1)
        f = (1.0_dp - a)**2 + 100.0_dp*(b - a*a)**2
        g(1) = -2.0_dp*(1.0_dp - a) - 400.0_dp*a*(b - a*a)
        g(2) = 200.0_dp*(b - a*a)
        call ad_status_set(st, FORTNUM_OK, "", &
            FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT)
    end subroutine rosen_grad

    ! (b) Vector residual JVP/VJP through jvp_fn/vjp_fn kernels, verified by the
    !     dot-product adjoint identity. Residual r: R^3 -> R^3,
    !       r1 = x1^2 + x2,  r2 = sin(x2) x3,  r3 = x1 x3 - x2.
    !     The harness check_jvp_vs_fd also confirms the JVP vs central FD.
    subroutine test_residual_adjoint(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(3), u(3), v(3)
        logical  :: ok

        x = [0.5_dp, 1.3_dp, -0.8_dp]
        v = [0.2_dp, -0.4_dp, 1.1_dp]
        u = [1.0_dp, -0.5_dp, 0.3_dp]

        ok = check_jvp_vs_fd("optimizer-jvp", res_primal, res_jvp, x, v, &
            1.0e-6_dp)
        ok = dot_product_identity("optimizer-adjoint", res_jvp, res_vjp, &
            x, u, v, 1.0e-12_dp) .and. ok
        if (.not. ok) nfail = nfail + 1
    end subroutine test_residual_adjoint

    ! Primal residual, harness vector_fn_i shape.
    subroutine res_primal(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = x(1)**2 + x(2)
        y(2) = sin(x(2))*x(3)
        y(3) = x(1)*x(3) - x(2)
    end subroutine res_primal

    ! Forward product J(x) v, harness jvp_fn_i shape. This is the body a
    ! jvp_fn optimizer kernel would call after unpacking its context.
    subroutine res_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = 2.0_dp*x(1)*v(1) + v(2)
        jv(2) = cos(x(2))*x(3)*v(2) + sin(x(2))*v(3)
        jv(3) = x(3)*v(1) - v(2) + x(1)*v(3)
    end subroutine res_jvp

    ! Reverse product J(x)^T u, harness vjp_fn_i shape; transpose of res_jvp.
    subroutine res_vjp(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        jtu(1) = 2.0_dp*x(1)*u(1) + x(3)*u(3)
        jtu(2) = u(1) + cos(x(2))*x(3)*u(2) - u(3)
        jtu(3) = sin(x(2))*u(2) + x(1)*u(3)
    end subroutine res_vjp

end program test_optimizer_api
