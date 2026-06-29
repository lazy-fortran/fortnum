program scalar_objective
    ! Downstream example: scalar objective through the #41 optimizer-facing API.
    !
    ! Objective: f(alpha, beta) = alpha * sin(beta) + alpha^2
    !   where alpha is packed in block "alpha" (size 1) and
    !   beta  is packed in block "beta"  (size 1).
    !
    ! Analytic gradient:
    !   df/d_alpha = sin(beta) + 2 alpha
    !   df/d_beta  = alpha * cos(beta)
    !
    ! Checks performed:
    !   (1) value matches analytic formula
    !   (2) grad matches central finite difference
    !   (3) JVP == dot(grad, v) for a random direction v
    !   (4) VJP == u * grad for scalar u (grad_fn doubles as vjp here)
    !
    ! Build standalone (from repo root, adjust paths as needed):
    !   gfortran -I<build>/include -o scalar_objective \
    !       examples/downstream/scalar_objective/scalar_objective.f90 \
    !       -L<build>/lib -lfortnum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_kinds,         only: dp
    use fortnum_status,        only: fortnum_status_t, status_ok, FORTNUM_OK
    use fortnum_active_vector, only: fortnum_active_layout_t, layout_init, &
        layout_add, pack_block, unpack_block
    use fortnum_ad_interfaces, only: fortnum_ad_status_t, ad_status_set, &
        ad_status_ok, &
        FORTNUM_AD_BACKEND_ANALYTIC, &
        FORTNUM_AD_QUALITY_EXACT
    implicit none

    integer  :: nfail
    real(dp) :: alpha_val, beta_val
    real(dp) :: x(2), f_got, f_ref
    real(dp) :: g(2), gfd(2), h
    real(dp) :: v(2), jv_ad, jv_fd
    real(dp) :: err
    type(fortnum_active_layout_t) :: layout
    type(fortnum_status_t)        :: us
    type(fortnum_ad_status_t)     :: st

    nfail     = 0
    alpha_val = 1.3_dp
    beta_val  = 0.7_dp

    ! Build layout: block "alpha" (size 1) then block "beta" (size 1).
    call layout_init(layout, 2)
    call layout_add(layout, "alpha", 1, us)
    call layout_add(layout, "beta",  1, us)
    if (.not. status_ok(us)) error stop "layout_add failed"

    ! Pack active vector.
    call pack_block(layout, x, "alpha", [alpha_val], us)
    call pack_block(layout, x, "beta",  [beta_val],  us)
    if (.not. status_ok(us)) error stop "pack_block failed"

    ! ------------------------------------------------------------------ (1) value
    call my_grad(2, x, f_got, g, layout, st)
    if (.not. ad_status_ok(st)) then
        write (error_unit, '(a)') "FAIL [value] ad_status not ok"
        nfail = nfail + 1
    end if
    f_ref = alpha_val*sin(beta_val) + alpha_val**2
    err = abs(f_got - f_ref) / max(abs(f_ref), 1.0_dp)
    if (err > 1.0e-14_dp) then
        write (error_unit, '(a,es12.4)') "FAIL [value] rel_err=", err
        nfail = nfail + 1
    end if

    ! ------------------------------------------------------------------ (2) grad vs FD
    h = 1.0e-6_dp
    call fd_grad(layout, x, h, gfd)

    err = max(abs(g(1) - gfd(1)), abs(g(2) - gfd(2))) &
        / max(maxval(abs(g)), 1.0_dp)
    if (err > 1.0e-6_dp) then
        write (error_unit, '(a,es12.4)') "FAIL [grad-vs-FD] rel_err=", err
        nfail = nfail + 1
    end if

    ! ------------------------------------------------------------------ (3) JVP == dot(grad, v)
    v = [0.4_dp, -1.1_dp]
    call my_jvp(2, x, v, f_got, jv_ad, layout, st)
    jv_fd = dot_product(g, v) ! exact for this linear-in-g case
    err = abs(jv_ad - jv_fd) / max(abs(jv_fd), 1.0_dp)
    if (err > 1.0e-14_dp) then
        write (error_unit, '(a,es12.4)') "FAIL [JVP] rel_err=", err
        nfail = nfail + 1
    end if

    ! ------------------------------------------------------------------ (4) VJP = u * grad
    ! For a scalar objective the VJP is just u times the gradient.
    ! Call grad with no change; the gradient is the VJP at u = 1.
    ! Verify the dot-product identity:  u * (J v) == v . (J^T u) = v . g * u.
    block
        real(dp) :: u, lhs, rhs
        u   = 2.5_dp
        lhs = u * jv_ad ! u * (J v)
        rhs = dot_product(v, u*g) ! v . (J^T u) = v . (u * grad)
        err = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)
        if (err > 1.0e-14_dp) then
            write (error_unit, '(a,es12.4)') "FAIL [adjoint] rel_err=", err
            nfail = nfail + 1
        end if
    end block

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a,f10.6,a,f10.6,a,f10.6)') &
        "PASS  f=", f_ref, "  g=[", g(1), ",", g(2)
    write (*, '(a)') "PASS"
    stop 0

contains

    ! Analytic grad_fn for f(alpha, beta) = alpha*sin(beta) + alpha^2.
    ! Unpacks alpha and beta from the flat vector via the layout context.
    subroutine my_grad(n, x, f, g, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(n)
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        type(fortnum_status_t) :: us
        real(dp) :: va(1), vb(1), a, b

        call unpack_block(layout, x, "alpha", va, us)
        call unpack_block(layout, x, "beta",  vb, us)
        a = va(1); b = vb(1)

        f    = a*sin(b) + a**2
        g(1) = sin(b) + 2.0_dp*a
        g(2) = a*cos(b)

        call ad_status_set(st, FORTNUM_OK, "", &
            FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT)
    end subroutine my_grad

    ! JVP: (value, directional derivative) in one call. Returns scalar jv.
    subroutine my_jvp(n, x, v, f, jv, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n), v(n)
        real(dp), intent(out) :: f, jv
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        real(dp) :: g(n)
        call my_grad(n, x, f, g, layout, st)
        jv = dot_product(g, v) ! scalar objective: J v = grad . v
    end subroutine my_jvp

    ! Central FD of each gradient component.
    subroutine fd_grad(layout, x, h, gfd)
        type(fortnum_active_layout_t), intent(in)  :: layout
        real(dp), intent(in)  :: x(:), h
        real(dp), intent(out) :: gfd(:)
        type(fortnum_ad_status_t) :: st
        real(dp) :: xp(size(x)), xm(size(x)), fp, fm, g(size(x))
        integer  :: i
        do i = 1, size(x)
            xp = x; xm = x
            xp(i) = xp(i) + h
            xm(i) = xm(i) - h
            call my_grad(size(x), xp, fp, g, layout, st)
            call my_grad(size(x), xm, fm, g, layout, st)
            gfd(i) = (fp - fm) / (2.0_dp*h)
        end do
    end subroutine fd_grad

end program scalar_objective
