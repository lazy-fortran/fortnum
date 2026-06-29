program test_bspline_ad
    ! Verify the derivative routines in fortnum_bspline:
    !   bspline_eval_jvp      -- d/dx [sum_i c_i B_i(x)] . vx   (x active)
    !   bspline_eval_vjp      -- (d/dx [sum_i c_i B_i(x)])^T . u (x active)
    !   bspline_eval_coef_jvp -- d/dc [sum_i c_i B_i(x)] . vc   (coef active)
    !   bspline_eval_coef_vjp -- (d/dc [sum_i c_i B_i(x)])^T . u (coef active)
    !
    ! Checks:
    !   1. JVP vs central finite difference of the spline value (inside a span).
    !   2. JVP vs the analytic derivative computed directly from
    !      bspline_eval_deriv row 1 (the analytic spline-derivative routine).
    !   3. Dot-product (adjoint) identity u.(Jv) = v.(J^T u).
    !   4. Span-boundary non-smoothness: crossing an interior knot changes the
    !      primal_only span index and is flagged via check_smoothness.
    !   5. Coefficient products: coef JVP vs FD, coef VJP (u=1) equals the basis
    !      vector B(x), and the coef adjoint identity.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_deriv, bspline_eval_jvp, &
        bspline_eval_vjp, bspline_span_index, bspline_eval_basis, &
        bspline_eval_coef_jvp, bspline_eval_coef_vjp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ad_test_utils, only: rel_err, fd_jvp_step, &
        check_smoothness, ad_status_t, AD_SMOOTH, AD_NONSMOOTH
    implicit none

    integer :: nfail
    nfail = 0

    call test_jvp_vs_fd(nfail)
    call test_jvp_vs_analytic(nfail)
    call test_dotprod_identity(nfail)
    call test_span_boundary_nonsmooth(nfail)
    call test_coef_products(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! Build a representative order-4 (cubic) clamped B-spline on 6 breakpoints.
    subroutine setup(ws, coef)
        type(bspline_workspace_t), intent(out)              :: ws
        real(dp), allocatable,     intent(out)              :: coef(:)
        type(fortnum_status_t) :: s
        real(dp), parameter :: brk(6) = &
            [0.0_dp, 0.2_dp, 0.45_dp, 0.7_dp, 0.9_dp, 1.0_dp]
        integer :: i

        call bspline_init(ws, 4, 6, s)
        if (s%code /= FORTNUM_OK) error stop "setup init failed"
        call bspline_set_knots(ws, brk, s)
        if (s%code /= FORTNUM_OK) error stop "setup set_knots failed"
        allocate (coef(ws%ncoef))
        do i = 1, ws%ncoef
            coef(i) = sin(1.3_dp*real(i, dp)) + 0.5_dp*real(i, dp)
        end do
    end subroutine setup

    subroutine test_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t)    :: s
        real(dp), allocatable     :: coef(:)
        real(dp) :: x, vx, jv, h, sp, sm
        real(dp), allocatable :: v0(:)

        call setup(ws, coef)
        allocate (v0(ws%ncoef))
        x  = 0.53_dp ! interior, inside one span
        vx = 1.0_dp

        call bspline_eval_jvp(ws, x, coef, vx, jv, s)
        if (s%code /= FORTNUM_OK) then
            write (error_unit, '(a)') "FAIL [jvp_fd] jvp status"
            nfail = nfail + 1
            return
        end if

        h = fd_jvp_step([x], [vx])
        call eval_spline(ws, x + h*vx, coef, sp)
        call eval_spline(ws, x - h*vx, coef, sm)
        if (rel_err(jv, (sp - sm)/(2.0_dp*h)) > 1.0e-7_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [jvp_vs_fd] ad=", jv, " fd=", (sp - sm)/(2.0_dp*h), &
                " rel_err=", rel_err(jv, (sp - sm)/(2.0_dp*h))
            nfail = nfail + 1
        end if
    end subroutine test_jvp_vs_fd

    subroutine test_jvp_vs_analytic(nfail)
        integer, intent(inout) :: nfail
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t)    :: s
        real(dp), allocatable     :: coef(:), dvals(:, :)
        real(dp) :: x, vx, jv, dsdx

        call setup(ws, coef)
        x  = 0.37_dp
        vx = 2.5_dp

        call bspline_eval_jvp(ws, x, coef, vx, jv, s)
        allocate (dvals(0:1, ws%ncoef))
        call bspline_eval_deriv(ws, x, 1, dvals, s)
        dsdx = dot_product(coef, dvals(1, :))

        if (rel_err(jv, dsdx*vx) > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [jvp_analytic] jv=", jv, " ds/dx*vx=", dsdx*vx, &
                " rel_err=", rel_err(jv, dsdx*vx)
            nfail = nfail + 1
        end if
    end subroutine test_jvp_vs_analytic

    subroutine test_dotprod_identity(nfail)
        integer, intent(inout) :: nfail
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t)    :: s
        real(dp), allocatable     :: coef(:)
        real(dp) :: x, u, vx, jv, jtu

        call setup(ws, coef)
        x  = 0.61_dp
        u  = 1.7_dp
        vx = -0.9_dp

        call bspline_eval_jvp(ws, x, coef, vx, jv, s)
        call bspline_eval_vjp(ws, x, coef, u, jtu, s)

        ! scalar->scalar map: u.(J vx) = vx.(J^T u)
        if (rel_err(u*jv, vx*jtu) > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [dotprod] u.(Jv)=", u*jv, " v.(J^Tu)=", vx*jtu, &
                " rel_err=", rel_err(u*jv, vx*jtu)
            nfail = nfail + 1
        end if
    end subroutine test_dotprod_identity

    subroutine test_span_boundary_nonsmooth(nfail)
        integer, intent(inout) :: nfail
        type(bspline_workspace_t) :: ws
        real(dp), allocatable     :: coef(:)
        real(dp)          :: eps, xk
        integer           :: tag_base, tag_pert
        type(ad_status_t) :: st

        call setup(ws, coef)
        eps = 1.0e-9_dp
        xk  = 0.45_dp ! an interior breakpoint / knot

        tag_base = bspline_span_index(ws, xk - eps)
        tag_pert = bspline_span_index(ws, xk + eps)
        st = check_smoothness(tag_base, tag_pert, expect=AD_NONSMOOTH)
        if (.not. st%ok) then
            write (error_unit, '(a)') &
                "FAIL [span_boundary] crossing a knot not flagged non-smooth"
            nfail = nfail + 1
        end if

        tag_pert = bspline_span_index(ws, xk - 2.0_dp*eps)
        st = check_smoothness(tag_base, tag_pert, expect=AD_SMOOTH)
        if (.not. st%ok) then
            write (error_unit, '(a)') &
                "FAIL [span_boundary] same span incorrectly flagged non-smooth"
            nfail = nfail + 1
        end if
    end subroutine test_span_boundary_nonsmooth

    ! (5) Coefficient-active products: s = sum_i c_i B_i(x) is linear in c with
    ! Jacobian B(x). Checks coef JVP vs central FD, coef VJP (u=1) == basis
    ! vector, and the adjoint identity u.(J vc) = vc.(J^T u).
    subroutine test_coef_products(nfail)
        integer, intent(inout) :: nfail
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t)    :: s
        real(dp), allocatable     :: coef(:), vc(:), jtu(:), basis(:), cp(:), cm(:)
        real(dp) :: x, u, jv, fd, sp, sm, h, e
        integer  :: i
        logical  :: ok

        call setup(ws, coef)
        allocate (vc(ws%ncoef), jtu(ws%ncoef), basis(ws%ncoef))
        allocate (cp(ws%ncoef), cm(ws%ncoef))
        x = 0.53_dp
        do i = 1, ws%ncoef
            vc(i) = cos(0.7_dp*real(i, dp))
        end do

        call bspline_eval_coef_jvp(ws, x, vc, jv, s)
        if (s%code /= FORTNUM_OK) then
            write (error_unit, '(a)') "FAIL [coef_jvp] status"
            nfail = nfail + 1
            return
        end if

        h  = 1.0e-6_dp
        cp = coef + h*vc
        cm = coef - h*vc
        call eval_spline(ws, x, cp, sp)
        call eval_spline(ws, x, cm, sm)
        fd = (sp - sm)/(2.0_dp*h)
        if (rel_err(jv, fd) > 1.0e-7_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [coef_jvp_vs_fd] jv=", jv, " fd=", fd, &
                " rel_err=", rel_err(jv, fd)
            nfail = nfail + 1
        end if

        call bspline_eval_coef_vjp(ws, x, 1.0_dp, jtu, s)
        call bspline_eval_basis(ws, x, basis, s)
        ok = .true.
        do i = 1, ws%ncoef
            if (jtu(i) /= basis(i)) ok = .false.
        end do
        if (.not. ok) then
            write (error_unit, '(a)') "FAIL [coef_vjp grad /= basis vector]"
            nfail = nfail + 1
        end if

        u = 1.7_dp
        call bspline_eval_coef_vjp(ws, x, u, jtu, s)
        e = abs(u*jv - dot_product(vc, jtu))/max(abs(u*jv), 1.0_dp)
        if (e > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [coef_dotprod] u.(Jv)=", u*jv, &
                " v.(J^Tu)=", dot_product(vc, jtu), " rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_coef_products

    ! Spline value s(x) = sum_i c_i B_i(x).
    subroutine eval_spline(ws, x, coef, val)
        type(bspline_workspace_t), intent(in)  :: ws
        real(dp),                  intent(in)  :: x, coef(:)
        real(dp),                  intent(out) :: val
        type(fortnum_status_t) :: s
        real(dp), allocatable  :: b(:)
        allocate (b(ws%ncoef))
        call bspline_eval_basis(ws, x, b, s)
        val = dot_product(coef, b)
    end subroutine eval_spline

end program test_bspline_ad
