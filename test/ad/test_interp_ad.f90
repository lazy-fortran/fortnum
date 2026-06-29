program test_interp_ad
    ! Verify derivative routines in fortnum_polynomial:
    !   lagrange_weights_jvp  -- d p(x)/dx . vx    (x active)
    !   lagrange_weights_vjp  -- (d p(x)/dx)^T . u (x active)
    !   lagrange_fval_jvp     -- d p/df . vf        (f values active)
    !   lagrange_fval_vjp     -- (d p/df)^T . u     (f values active)
    !
    ! Checks per derivative product:
    !   1. JVP vs central finite difference.
    !   2. Analytic derivative of a known polynomial (d/dx of cubic = quadratic).
    !   3. Dot-product identity u.(Jv) = v.(J^T u).
    !
    ! Grid-boundary non-smoothness is confirmed by check_smoothness.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_polynomial, only: lagrange_weights_jvp, lagrange_weights_vjp, &
        lagrange_fval_jvp, lagrange_fval_vjp, lagrange_weights
    use fortnum_ad_test_utils, only: fd_jvp_step, rel_err, &
        check_smoothness, ad_status_t, AD_SMOOTH, AD_NONSMOOTH
    implicit none

    integer :: nfail
    nfail = 0

    call test_xjvp_vs_analytic(nfail)
    call test_xjvp_vs_fd(nfail)
    call test_xvjp_identity(nfail)
    call test_fval_jvp_vs_fd(nfail)
    call test_fval_dotprod_identity(nfail)
    call test_cell_boundary_nonsmooth(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! ------------------------------------------------------------------
    ! Cubic polynomial p(x) = x^3 - 2 x^2 + x - 0.5, analytic derivative
    ! dp/dx = 3 x^2 - 4 x + 1.  Interpolated on 4 equidistant nodes in [0,1];
    ! the interpolant is exact for a degree-3 polynomial, so the derivative
    ! weights recover the exact analytic derivative.
    subroutine test_xjvp_vs_analytic(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 4
        real(dp), parameter :: xp(n) = [ &
            0.0_dp, 0.25_dp, 0.75_dp, 1.0_dp]
        real(dp) :: f(n), x, vx, jv, dpda
        integer  :: i

        do i = 1, n
            f(i) = xp(i)**3 - 2.0_dp*xp(i)**2 + xp(i) - 0.5_dp
        end do

        x   = 0.4_dp
        vx  = 1.0_dp

        call lagrange_weights_jvp(n, x, xp, f, vx, jv)
        dpda = 3.0_dp*x**2 - 4.0_dp*x + 1.0_dp

        if (rel_err(jv, dpda) > 1.0e-10_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [xjvp_analytic] got=", jv, " want=", dpda, &
                " rel_err=", rel_err(jv, dpda)
            nfail = nfail + 1
        end if
    end subroutine test_xjvp_vs_analytic


    ! ------------------------------------------------------------------
    ! JVP w.r.t. x vs central FD.  Use the interpolant of a smooth function
    ! on 6 nodes; the FD step is safeguarded so relative error is ~1e-8 and
    ! we test at 1e-5 to give room for interpolation error.
    subroutine test_xjvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 6
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        real(dp) :: xp(n), f(n), x, vx, jv_ad, jv_fd, h, pp, pm
        real(dp) :: coef_p(n), coef_m(n)
        integer  :: i

        do i = 1, n
            xp(i) = 0.1_dp + 0.8_dp * 0.5_dp * &
                (1.0_dp - cos(pi*real(i-1, dp)/real(n-1, dp)))
        end do
        do i = 1, n
            f(i) = sin(pi*xp(i))
        end do

        x  = 0.5_dp
        vx = 1.0_dp

        call lagrange_weights_jvp(n, x, xp, f, vx, jv_ad)

        h = fd_jvp_step([x], [vx])
        call lagrange_weights(n, x + h*vx, xp, coef_p)
        call lagrange_weights(n, x - h*vx, xp, coef_m)
        pp = dot_product(f, coef_p)
        pm = dot_product(f, coef_m)
        jv_fd = (pp - pm) / (2.0_dp*h)

        if (rel_err(jv_ad, jv_fd) > 1.0e-5_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [xjvp_vs_fd] ad=", jv_ad, " fd=", jv_fd, &
                " rel_err=", rel_err(jv_ad, jv_fd)
            nfail = nfail + 1
        end if
    end subroutine test_xjvp_vs_fd


    ! ------------------------------------------------------------------
    ! Dot-product identity for the x-direction map (scalar->scalar):
    !   u . (J vx) = vx . (J^T u)
    ! reduces to  u*(dp/dx)*vx == vx*(dp/dx)*u  (trivially true), but both
    ! code paths are exercised.
    subroutine test_xvjp_identity(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 4
        real(dp), parameter :: xp(n) = [0.0_dp, 0.3_dp, 0.7_dp, 1.0_dp]
        real(dp) :: f(n), x, u, vx, jv, jtu
        integer  :: i

        do i = 1, n
            f(i) = exp(-xp(i))
        end do
        x  = 0.45_dp
        u  = 2.3_dp
        vx = -1.1_dp

        call lagrange_weights_jvp(n, x, xp, f, vx, jv)
        call lagrange_weights_vjp(n, x, xp, f, u, jtu)

        if (rel_err(u*jv, vx*jtu) > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [x_dotprod] u.(Jv)=", u*jv, " v.(J^Tu)=", vx*jtu, &
                " rel_err=", rel_err(u*jv, vx*jtu)
            nfail = nfail + 1
        end if
    end subroutine test_xvjp_identity


    ! ------------------------------------------------------------------
    ! JVP w.r.t. f values: map f -> p(x) = sum_i coef_i f_i is linear, so
    ! JVP = dot(coef, vf).  Verified against central FD on that linear map.
    subroutine test_fval_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 5
        real(dp), parameter :: xp(n) = [ &
            0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp) :: f(n), vf(n), jv_ad, jv_fd, h
        real(dp) :: coef(n)
        real(dp) :: x
        integer  :: i

        x = 0.35_dp

        do i = 1, n
            f(i) = cos(real(i, dp))
        end do
        vf = [0.3_dp, -0.7_dp, 1.1_dp, -0.2_dp, 0.5_dp]

        call lagrange_fval_jvp(n, x, xp, vf, jv_ad)

        ! Central FD on the linear map: only the weights matter, not f itself.
        ! p(f + h vf) - p(f - h vf) = 2h * dot(coef, vf).
        h = fd_jvp_step(f, vf)
        call lagrange_weights(n, x, xp, coef)
        jv_fd = dot_product(coef, (f + h*vf)) / (2.0_dp*h) - &
            dot_product(coef, (f - h*vf)) / (2.0_dp*h)

        if (rel_err(jv_ad, jv_fd) > 1.0e-10_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [fval_jvp_fd] ad=", jv_ad, " fd=", jv_fd, &
                " rel_err=", rel_err(jv_ad, jv_fd)
            nfail = nfail + 1
        end if
    end subroutine test_fval_jvp_vs_fd


    ! ------------------------------------------------------------------
    ! Dot-product identity for the f-values map f -> [p(x)].
    ! Map: R^n -> R^1.  Identity: u . (J vf) = vf . (J^T u).
    subroutine test_fval_dotprod_identity(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 5
        real(dp), parameter :: xp(n) = [ &
            0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp) :: vf(n), jtu(n), jv, u, x, lhs, rhs
        integer  :: i

        x  = 0.6_dp
        u  = 3.7_dp
        do i = 1, n
            vf(i) = sin(real(i, dp)*0.4_dp)
        end do

        call lagrange_fval_jvp(n, x, xp, vf, jv)
        call lagrange_fval_vjp(n, x, xp, u, jtu)

        lhs = u * jv ! u . (J vf)   [scalar dot scalar]
        rhs = dot_product(vf, jtu) ! vf . (J^T u)

        if (rel_err(lhs, rhs) > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [fval_dotprod] u.(Jv)=", lhs, " v.(J^Tu)=", rhs, &
                " rel_err=", rel_err(lhs, rhs)
            nfail = nfail + 1
        end if
    end subroutine test_fval_dotprod_identity


    ! ------------------------------------------------------------------
    ! Confirm that a point crossing a grid node is flagged as a non-smooth
    ! event via check_smoothness.  The cell tag is the index of the largest
    ! node <= x; it changes discontinuously when x crosses a node.
    subroutine test_cell_boundary_nonsmooth(nfail)
        integer, intent(inout) :: nfail
        integer,  parameter :: n = 5
        real(dp), parameter :: xp(n) = [ &
            0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
        real(dp)          :: eps
        integer           :: tag_base, tag_pert
        type(ad_status_t) :: s

        eps = 1.0e-8_dp

        ! Just below 0.5 vs just above 0.5: different cells -> non-smooth.
        tag_base = cell_tag(xp, 0.5_dp - eps)
        tag_pert = cell_tag(xp, 0.5_dp + eps)
        s = check_smoothness(tag_base, tag_pert, expect=AD_NONSMOOTH)
        if (.not. s%ok) then
            write (error_unit, '(a)') &
                "FAIL [cell_boundary] crossing a node not flagged non-smooth"
            nfail = nfail + 1
        end if

        ! Both strictly below 0.5: same cell -> smooth.
        tag_pert = cell_tag(xp, 0.5_dp - 2.0_dp*eps)
        s = check_smoothness(tag_base, tag_pert, expect=AD_SMOOTH)
        if (.not. s%ok) then
            write (error_unit, '(a)') &
                "FAIL [cell_boundary] same cell incorrectly flagged non-smooth"
            nfail = nfail + 1
        end if
    end subroutine test_cell_boundary_nonsmooth


    ! ------------------------------------------------------------------
    ! Returns the index of the largest node <= x.
    pure function cell_tag(xp, x) result(tag)
        real(dp), intent(in) :: xp(:)
        real(dp), intent(in) :: x
        integer :: tag, i
        tag = 1
        do i = 1, size(xp)
            if (xp(i) <= x) tag = i
        end do
    end function cell_tag

end program test_interp_ad
