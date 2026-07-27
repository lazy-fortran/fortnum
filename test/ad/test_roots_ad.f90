program test_roots_ad
    ! Derivative tests for fortnum_roots implicit_rule (issue #40).
    !
    ! Test problem: f(x, p) = x^2 - p = 0, positive root x*(p) = sqrt(p).
    !   f_x = 2 x* = 2 sqrt(p)
    !   f_p = -1
    !   dx*/dp = -f_p/f_x = 1/(2 sqrt(p))  (analytic)
    !
    ! Tests:
    !   1. root_grad vs analytic 1/(2 sqrt(p)).
    !   2. root_grad vs central FD: re-solve at p+h and p-h with root_brent.
    !   3. root_jvp (scalar p as 1-vector) vs analytic.
    !   4. dot-product identity: root_jvp and root_vjp satisfy u.(Jv)=v.(J^T u).
    !   5. Near-multiple-root guard: |f_x| ~ 0 -> FORTNUM_DOMAIN_ERROR.
    !   6. Vector-p case: f(x,p1,p2) = x^2 - p1 - p2 = 0; grad = [1,1]/(2x*).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_kinds,  only: dp
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_DOMAIN_ERROR
    use fortnum_roots,  only: root_brent, root_implicit_jvp, &
        root_implicit_vjp, root_jvp, root_vjp, root_grad
    use fortnum_multiroot, only: multiroot_hybrid, multiroot_implicit_jvp, &
        multiroot_implicit_vjp, multiroot_jvp, multiroot_vjp, multiroot_grad
    implicit none

    type :: diagonal_preconditioner_t
        real(dp) :: inverse_diagonal(2)
        integer :: calls = 0
    end type diagonal_preconditioner_t

    ! p_solve is shared between solve_x2mp and f_x2mp via host association.
    real(dp) :: p_solve
    ! p_mr is shared between solve_mr and fdf_mr via host association.
    real(dp) :: p_mr(2)
    ! p_tangent_solve is used only by the independent complete-solve oracle.
    real(dp) :: p_tangent_solve(2)

    integer :: nfail
    nfail = 0

    call test_grad_analytic(nfail)
    call test_grad_vs_fd(nfail)
    call test_jvp_scalar_p(nfail)
    call test_dot_product_id(nfail)
    call test_near_multiple_root(nfail)
    call test_vector_p(nfail)
    call test_implicit_tangent_boundary(nfail)
    call test_implicit_adjoint_boundary(nfail)
    call test_multiroot_grad_vs_fd(nfail)
    call test_multiroot_jvp_vs_fd(nfail)
    call test_multiroot_implicit_tangent_boundary(nfail)
    call test_multiroot_implicit_adjoint_boundary(nfail)
    call test_multiroot_dot_product_id(nfail)
    call test_multiroot_singular(nfail)
    call test_multiroot_preconditioner_hook(nfail)
    call test_multiroot_condition_diagnostic(nfail)

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

    ! f(x) = x^2 - p used as root_fn_t callback (module-level p via entry).
    ! Use a wrapper that closes over a module variable to avoid global state;
    ! since fortnum_roots requires a procedure(root_fn_t), we use a fixed p
    ! passed through a local module variable approach -- but Fortran does not
    ! permit closures.  Use an internal procedure that references a host variable.

    ! Test 1: root_grad vs analytic 1/(2 sqrt(p)).
    subroutine test_grad_analytic(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p, xstar, f_x, f_p, dxdp, dxdp_exact
        type(fortnum_status_t) :: st

        p      = 4.0_dp
        xstar  = sqrt(p) ! 2.0
        f_x    = 2.0_dp * xstar ! 4.0
        f_p    = -1.0_dp
        dxdp_exact = 1.0_dp / (2.0_dp * xstar) ! 0.25

        call root_grad(f_x, f_p, dxdp, st)
        if (.not. status_ok(st)) then
            write (error_unit, '(a)') "FAIL [grad_analytic] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (rel_err(dxdp, dxdp_exact) > 1.0e-14_dp) then
            write (error_unit, '(a,es24.16,a,es24.16)') &
                "FAIL [grad_analytic] got=", dxdp, " want=", dxdp_exact
            nfail = nfail + 1
        end if
    end subroutine test_grad_analytic

    ! Test 2: root_grad vs central FD (re-solve at p+h and p-h).
    subroutine test_grad_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p, xstar, f_x, f_p, dxdp, dxdp_fd
        real(dp) :: xp, xm, h
        type(fortnum_status_t) :: st

        p  = 3.0_dp
        h  = 1.0e-5_dp

        ! Solve at p, p+h, p-h via root_brent.
        call solve_x2mp(p,   xstar)
        call solve_x2mp(p+h, xp)
        call solve_x2mp(p-h, xm)

        f_x = 2.0_dp * xstar
        f_p = -1.0_dp
        call root_grad(f_x, f_p, dxdp, st)

        dxdp_fd = (xp - xm) / (2.0_dp * h)

        if (rel_err(dxdp, dxdp_fd) > 1.0e-8_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [grad_vs_fd] analytic=", dxdp, " fd=", dxdp_fd, &
                " rel_err=", rel_err(dxdp, dxdp_fd)
            nfail = nfail + 1
        end if
    end subroutine test_grad_vs_fd

    ! Test 3: root_jvp (1-vector p) vs analytic for p=2, dp=0.7.
    subroutine test_jvp_scalar_p(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p, xstar, f_x, f_p(1), dp_vec(1), dx, dx_exact
        type(fortnum_status_t) :: st

        p      = 2.0_dp
        xstar  = sqrt(p)
        f_x    = 2.0_dp * xstar
        f_p    = [-1.0_dp]
        dp_vec = [0.7_dp]
        dx_exact = (-f_p(1) * dp_vec(1)) / f_x ! = 0.7 / (2 sqrt(2))

        call root_jvp(f_x, f_p, dp_vec, dx, st)
        if (.not. status_ok(st)) then
            write (error_unit, '(a)') "FAIL [jvp_scalar_p] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (rel_err(dx, dx_exact) > 1.0e-14_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [jvp_scalar_p] got=", dx, " want=", dx_exact, &
                " rel_err=", rel_err(dx, dx_exact)
            nfail = nfail + 1
        end if
    end subroutine test_jvp_scalar_p

    ! Test 4: dot-product identity u.(J v) = v.(J^T u) for vector-p case.
    ! Map: p (2-vector) -> x*(p) (scalar). J is 1x2 row vector [-fp1/fx, -fp2/fx].
    ! JVP: dx = -(f_p . dp)/f_x  (scalar).
    ! VJP: jtu_i = -(f_p_i/f_x)*u  (2-vector).
    ! Identity: u * (J v) == v . (J^T u)
    !   lhs = u * dx
    !   rhs = dot(v, jtu)
    subroutine test_dot_product_id(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p1, p2, xstar, f_x, f_p(2)
        real(dp) :: dp_vec(2), u, jtu(2), dx, lhs, rhs, e
        type(fortnum_status_t) :: st

        p1    = 3.0_dp
        p2    = 1.0_dp
        ! f(x,p1,p2) = x^2 - p1 - p2 -> root = sqrt(p1+p2) = 2
        xstar = sqrt(p1 + p2)
        f_x   = 2.0_dp * xstar
        f_p   = [-1.0_dp, -1.0_dp]

        dp_vec = [0.3_dp, -0.5_dp]
        u      = 1.7_dp

        call root_jvp(f_x, f_p, dp_vec, dx, st)
        call root_vjp(f_x, f_p, u, jtu, st)

        lhs = u * dx
        rhs = dot_product(dp_vec, jtu)
        e   = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)

        if (e > 1.0e-14_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [dot_product_id] u.(Jv)=", lhs, &
                " v.(J^T u)=", rhs, " rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_dot_product_id

    ! Test 5: near-multiple root guard returns FORTNUM_DOMAIN_ERROR.
    subroutine test_near_multiple_root(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: dxdp, dx, jtu(1)
        type(fortnum_status_t) :: st

        ! f_x nearly zero: double root at x=0 for f(x)=x^2.
        call root_grad(1.0e-20_dp, -1.0_dp, dxdp, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [near_multiple_grad] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if

        call root_jvp(1.0e-20_dp, [-1.0_dp], [1.0_dp], dx, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [near_multiple_jvp] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if

        call root_vjp(1.0e-20_dp, [-1.0_dp], 1.0_dp, jtu, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [near_multiple_vjp] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if
    end subroutine test_near_multiple_root

    ! Test 6: vector-p, root_grad not applicable; use root_jvp/root_vjp.
    ! f(x,p1,p2) = x^2 - p1 - 2*p2 = 0, p1=3, p2=0.5 -> x*=2.
    ! f_p = [-1, -2], f_x = 4.
    ! Sensitivity in direction dp=[1,0]: dx* = 1/4.
    ! Sensitivity in direction dp=[0,1]: dx* = 2/4 = 0.5.
    subroutine test_vector_p(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: f_x, f_p(2), dp_vec(2), dx
        type(fortnum_status_t) :: st

        f_x = 4.0_dp
        f_p = [-1.0_dp, -2.0_dp]

        dp_vec = [1.0_dp, 0.0_dp]
        call root_jvp(f_x, f_p, dp_vec, dx, st)
        if (rel_err(dx, 0.25_dp) > 1.0e-14_dp) then
            write (error_unit, '(a,es24.16,a,es24.16)') &
                "FAIL [vector_p dir1] got=", dx, " want=", 0.25_dp
            nfail = nfail + 1
        end if

        dp_vec = [0.0_dp, 1.0_dp]
        call root_jvp(f_x, f_p, dp_vec, dx, st)
        if (rel_err(dx, 0.5_dp) > 1.0e-14_dp) then
            write (error_unit, '(a,es24.16,a,es24.16)') &
                "FAIL [vector_p dir2] got=", dx, " want=", 0.5_dp
            nfail = nfail + 1
        end if
    end subroutine test_vector_p

    ! Generic boundary vs central difference of complete root solves.
    ! f(x,p) = x^3 + p1*x - p2 = 0 has one real root for p1 > 0.
    subroutine test_implicit_tangent_boundary(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p(2), pp(2), pm(2), tp(2), xstar, xp, xm, dx, dx_fd, h
        type(fortnum_status_t) :: st

        p = [0.7_dp, 3.0_dp]
        tp = [0.4_dp, -0.6_dp]
        h = 1.0e-5_dp
        call solve_tangent_root(p, xstar)
        pp = p + h*tp
        pm = p - h*tp
        call solve_tangent_root(pp, xp)
        call solve_tangent_root(pm, xm)

        call root_implicit_jvp(tangent_residual_jvp, xstar, p, tp, dx, st)
        dx_fd = (xp - xm)/(2.0_dp*h)

        if (.not. status_ok(st)) then
            write (error_unit, '(a)') &
                "FAIL [implicit_tangent_boundary] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (rel_err(dx, dx_fd) > 1.0e-8_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [implicit_tangent_boundary] implicit=", dx, &
                " fd=", dx_fd, " rel_err=", rel_err(dx, dx_fd)
            nfail = nfail + 1
        end if
    end subroutine test_implicit_tangent_boundary

    subroutine tangent_residual_jvp(x, p, tp, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, p(:), tp(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        f_x = 3.0_dp*x*x + p(1)
        f_p_tp = x*tp(1) - tp(2)
    end subroutine tangent_residual_jvp

    ! Scalar-objective gradient from the adjoint boundary vs complete re-solves.
    subroutine test_implicit_adjoint_boundary(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p(2), pp(2), pm(2), xstar, xp, xm
        real(dp) :: u, jtu(2), jtu_fd(2), h
        type(fortnum_status_t) :: st
        integer :: i

        p = [0.7_dp, 3.0_dp]
        u = 1.3_dp
        h = 1.0e-5_dp
        call solve_tangent_root(p, xstar)
        call root_implicit_vjp(adjoint_residual_vjp, xstar, p, u, jtu, st)

        do i = 1, 2
            pp = p
            pm = p
            pp(i) = pp(i) + h
            pm(i) = pm(i) - h
            call solve_tangent_root(pp, xp)
            call solve_tangent_root(pm, xm)
            jtu_fd(i) = u*(xp - xm)/(2.0_dp*h)
        end do

        if (.not. status_ok(st)) then
            write (error_unit, '(a)') &
                "FAIL [implicit_adjoint_boundary] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (maxval(abs(jtu - jtu_fd)) > 1.0e-8_dp) then
            write (error_unit, '(a,2es24.16)') &
                "FAIL [implicit_adjoint_boundary] implicit=", jtu
            write (error_unit, '(a,2es24.16)') &
                "FAIL [implicit_adjoint_boundary] fd=", jtu_fd
            nfail = nfail + 1
        end if
    end subroutine test_implicit_adjoint_boundary

    subroutine adjoint_residual_vjp(x, p, u, f_x, f_p_t_u, context)
        real(dp), intent(in) :: x, p(:), u
        real(dp), intent(out) :: f_x, f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        f_x = 3.0_dp*x*x + p(1)
        f_p_t_u(1) = x*u
        f_p_t_u(2) = -u
    end subroutine adjoint_residual_vjp

    subroutine solve_tangent_root(p, x)
        real(dp), intent(in) :: p(2)
        real(dp), intent(out) :: x
        type(fortnum_status_t) :: st

        p_tangent_solve = p
        call root_brent(tangent_residual, 0.0_dp, 3.0_dp, x, st, &
            xtol=1.0e-13_dp, ftol=1.0e-14_dp)
        if (.not. status_ok(st)) error stop "tangent root solve failed"
    end subroutine solve_tangent_root

    pure function tangent_residual(x) result(f)
        real(dp), intent(in) :: x
        real(dp) :: f

        f = x*x*x + p_tangent_solve(1)*x - p_tangent_solve(2)
    end function tangent_residual

    ! n-dim test system F(x, p) = 0 with p in R^2:
    !   F1 = x1^2 + x2 - p1
    !   F2 = x1 + x2^2 - p2
    ! At p = (2, 2) the root is x* = (1, 1).
    !   J_x = [[2 x1, 1], [1, 2 x2]] = [[2, 1], [1, 2]] at the root, det = 3.
    !   J_p = dF/dp = [[-1, 0], [0, -1]] (constant).
    !   dx*/dp = -J_x^{-1} J_p = J_x^{-1} = (1/3) [[2, -1], [-1, 2]].

    ! Test 7: multiroot_grad (scalar-p, only p1 active) vs central FD.
    ! f_p = dF/dp1 = [-1, 0]; sensitivity is column 1 of dx*/dp = (1/3)[2, -1].
    subroutine test_multiroot_grad_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: jac_x(2, 2), f_p(2), dxdp(2)
        real(dp) :: pbase(2), xp(2), xm(2), dxdp_fd(2), h
        type(fortnum_status_t) :: st
        integer :: i

        pbase = [2.0_dp, 2.0_dp]
        jac_x = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        f_p = [-1.0_dp, 0.0_dp]

        call multiroot_grad(jac_x, f_p, dxdp, st)
        if (.not. status_ok(st)) then
            write (error_unit, '(a)') "FAIL [mr_grad_fd] unexpected status error"
            nfail = nfail + 1
            return
        end if

        h = 1.0e-5_dp
        call solve_mr([pbase(1) + h, pbase(2)], xp)
        call solve_mr([pbase(1) - h, pbase(2)], xm)
        dxdp_fd = (xp - xm) / (2.0_dp * h)

        do i = 1, 2
            if (rel_err(dxdp(i), dxdp_fd(i)) > 1.0e-7_dp) then
                write (error_unit, '(a,i0,a,es24.16,a,es24.16)') &
                    "FAIL [mr_grad_fd] comp ", i, " analytic=", dxdp(i), &
                    " fd=", dxdp_fd(i)
                nfail = nfail + 1
            end if
        end do
    end subroutine test_multiroot_grad_vs_fd

    ! Test 8: multiroot_jvp vs central FD along a parameter tangent.
    subroutine test_multiroot_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: jac_x(2, 2), f_p(2, 2), tp(2), dx(2)
        real(dp) :: pbase(2), xp(2), xm(2), dx_fd(2), h
        type(fortnum_status_t) :: st
        integer :: i

        pbase = [2.0_dp, 2.0_dp]
        jac_x = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        f_p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
        tp = [0.3_dp, -0.5_dp]

        call multiroot_jvp(jac_x, f_p, tp, dx, st)
        if (.not. status_ok(st)) then
            write (error_unit, '(a)') "FAIL [mr_jvp_fd] unexpected status error"
            nfail = nfail + 1
            return
        end if

        h = 1.0e-5_dp
        call solve_mr(pbase + h * tp, xp)
        call solve_mr(pbase - h * tp, xm)
        dx_fd = (xp - xm) / (2.0_dp * h)

        do i = 1, 2
            if (rel_err(dx(i), dx_fd(i)) > 1.0e-7_dp) then
                write (error_unit, '(a,i0,a,es24.16,a,es24.16)') &
                    "FAIL [mr_jvp_fd] comp ", i, " jvp=", dx(i), &
                    " fd=", dx_fd(i)
                nfail = nfail + 1
            end if
        end do
    end subroutine test_multiroot_jvp_vs_fd

    subroutine test_multiroot_implicit_tangent_boundary(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p(2), pp(2), pm(2), tp(2)
        real(dp) :: xstar(2), xp(2), xm(2), dx(2), dx_fd(2), h
        type(fortnum_status_t) :: st

        p = [2.0_dp, 2.0_dp]
        tp = [0.3_dp, -0.5_dp]
        h = 1.0e-5_dp
        call solve_mr(p, xstar)
        pp = p + h*tp
        pm = p - h*tp
        call solve_mr(pp, xp)
        call solve_mr(pm, xm)

        call multiroot_implicit_jvp(vector_residual_jvp, xstar, p, tp, &
            dx, st)
        dx_fd = (xp - xm)/(2.0_dp*h)

        if (.not. status_ok(st)) then
            write (error_unit, '(a)') &
                "FAIL [multiroot_implicit_tangent] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (maxval(abs(dx - dx_fd)) > 1.0e-7_dp) then
            write (error_unit, '(a,2es24.16)') &
                "FAIL [multiroot_implicit_tangent] implicit=", dx
            write (error_unit, '(a,2es24.16)') &
                "FAIL [multiroot_implicit_tangent] fd=", dx_fd
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_implicit_tangent_boundary

    subroutine vector_residual_jvp(x, p, tp, jac_x, f_p_tp, context)
        real(dp), intent(in) :: x(:), p(:), tp(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        real(dp), intent(out) :: f_p_tp(size(x))
        class(*), intent(inout), optional :: context

        if (size(p) /= 2) error stop "vector residual expects two parameters"
        jac_x(1, 1) = 2.0_dp*x(1)
        jac_x(1, 2) = 1.0_dp
        jac_x(2, 1) = 1.0_dp
        jac_x(2, 2) = 2.0_dp*x(2)
        f_p_tp(1) = -tp(1)
        f_p_tp(2) = -tp(2)
    end subroutine vector_residual_jvp

    subroutine test_multiroot_implicit_adjoint_boundary(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: p(2), pp(2), pm(2), xstar(2), xp(2), xm(2)
        real(dp) :: u(2), jtu(2), jtu_fd(2), h
        type(fortnum_status_t) :: st
        integer :: i

        p = [2.0_dp, 2.0_dp]
        u = [1.3_dp, -0.4_dp]
        h = 1.0e-5_dp
        call solve_mr(p, xstar)
        call multiroot_implicit_vjp(vector_state_jacobian, &
            vector_parameter_vjp, xstar, p, u, jtu, st)

        do i = 1, 2
            pp = p
            pm = p
            pp(i) = pp(i) + h
            pm(i) = pm(i) - h
            call solve_mr(pp, xp)
            call solve_mr(pm, xm)
            jtu_fd(i) = dot_product(u, xp - xm)/(2.0_dp*h)
        end do

        if (.not. status_ok(st)) then
            write (error_unit, '(a)') &
                "FAIL [multiroot_implicit_adjoint] unexpected status error"
            nfail = nfail + 1
            return
        end if
        if (maxval(abs(jtu - jtu_fd)) > 1.0e-7_dp) then
            write (error_unit, '(a,2es24.16)') &
                "FAIL [multiroot_implicit_adjoint] implicit=", jtu
            write (error_unit, '(a,2es24.16)') &
                "FAIL [multiroot_implicit_adjoint] fd=", jtu_fd
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_implicit_adjoint_boundary

    subroutine vector_state_jacobian(x, p, jac_x, context)
        real(dp), intent(in) :: x(:), p(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        class(*), intent(inout), optional :: context

        if (size(p) /= 2) error stop "vector residual expects two parameters"
        jac_x(1, 1) = 2.0_dp*x(1)
        jac_x(1, 2) = 1.0_dp
        jac_x(2, 1) = 1.0_dp
        jac_x(2, 2) = 2.0_dp*x(2)
    end subroutine vector_state_jacobian

    subroutine vector_parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        if (size(x) /= 2) error stop "vector residual expects two states"
        f_p_t_u(1) = -u(1)
        f_p_t_u(2) = -u(2)
    end subroutine vector_parameter_vjp

    subroutine test_multiroot_preconditioner_hook(nfail)
        integer, intent(inout) :: nfail
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: jac_x(2, 2), f_p(2, 2), tp(2), u(2)
        real(dp) :: dx(2), jtu(2), dx_fd(2), xp(2), xm(2), pbase(2)
        real(dp) :: lhs, rhs
        type(diagonal_preconditioner_t) :: preconditioner
        type(fortnum_status_t) :: st

        pbase = [2.0_dp, 2.0_dp]
        jac_x = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        f_p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
        tp = [0.3_dp, -0.5_dp]
        u = [1.7_dp, -0.4_dp]
        preconditioner%inverse_diagonal = [0.5_dp, 0.5_dp]

        call multiroot_jvp(jac_x, f_p, tp, dx, st, diagonal_solve, preconditioner)
        call solve_mr(pbase + h*tp, xp)
        call solve_mr(pbase - h*tp, xm)
        dx_fd = (xp - xm)/(2.0_dp*h)
        if (maxval(abs(dx - dx_fd)) > 1.0e-7_dp) then
            write (error_unit, '(a)') "FAIL [mr_preconditioner] JVP vs FD"
            nfail = nfail + 1
        end if

        call multiroot_vjp(jac_x, f_p, u, jtu, st, diagonal_solve, preconditioner)
        lhs = dot_product(u, dx_fd)
        rhs = dot_product(tp, jtu)
        if (abs(lhs - rhs) > 1.0e-7_dp) then
            write (error_unit, '(a)') "FAIL [mr_preconditioner] VJP vs FD"
            nfail = nfail + 1
        end if
        if (preconditioner%calls /= 2) then
            write (error_unit, '(a,i0)') &
                "FAIL [mr_preconditioner] hook calls=", preconditioner%calls
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_preconditioner_hook

    subroutine test_multiroot_condition_diagnostic(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: jac_x(2, 2), f_p(2, 2), tp(2), dx(2), rcond
        real(dp) :: u(2), jtu(2), scalar_f_p(2), dxdp(2)
        type(fortnum_status_t) :: st

        f_p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
        tp = [0.3_dp, -0.5_dp]

        jac_x = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        call multiroot_jvp(jac_x, f_p, tp, dx, st, &
            reciprocal_condition=rcond)
        if (abs(rcond - 1.0_dp/3.0_dp) > 2.0e-15_dp) then
            write (error_unit, '(a,es24.16)') &
                "FAIL [mr_condition] expected 1/3, got ", rcond
            nfail = nfail + 1
        end if

        jac_x = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0e-8_dp], [2, 2])
        call multiroot_jvp(jac_x, f_p, tp, dx, st, &
            reciprocal_condition=rcond)
        if (abs(rcond - 1.0e-8_dp) > 1.0e-22_dp) then
            write (error_unit, '(a,es24.16)') &
                "FAIL [mr_condition] expected 1e-8, got ", rcond
            nfail = nfail + 1
        end if

        call multiroot_jvp(jac_x, f_p, tp, dx, st, &
            reciprocal_condition=rcond, minimum_reciprocal_condition=1.0e-6_dp)
        if (st%code /= FORTNUM_DOMAIN_ERROR .or. maxval(abs(dx)) /= 0.0_dp) then
            write (error_unit, '(a)') &
                "FAIL [mr_condition] unreliable JVP status"
            nfail = nfail + 1
        end if

        u = [0.4_dp, -0.2_dp]
        call multiroot_vjp(jac_x, f_p, u, jtu, st, &
            minimum_reciprocal_condition=1.0e-6_dp)
        if (st%code /= FORTNUM_DOMAIN_ERROR .or. maxval(abs(jtu)) /= 0.0_dp) then
            write (error_unit, '(a)') &
                "FAIL [mr_condition] unreliable VJP status"
            nfail = nfail + 1
        end if

        scalar_f_p = f_p(:, 1)
        call multiroot_grad(jac_x, scalar_f_p, dxdp, st, &
            minimum_reciprocal_condition=1.0e-6_dp)
        if (st%code /= FORTNUM_DOMAIN_ERROR .or. maxval(abs(dxdp)) /= 0.0_dp) then
            write (error_unit, '(a)') &
                "FAIL [mr_condition] unreliable gradient status"
            nfail = nfail + 1
        end if

        call multiroot_jvp(jac_x, f_p, tp, dx, st, &
            minimum_reciprocal_condition=1.0e-9_dp)
        if (.not. status_ok(st)) then
            write (error_unit, '(a)') &
                "FAIL [mr_condition] reliable threshold rejected"
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_condition_diagnostic

    subroutine diagonal_solve(a, b, x, info, context)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        class(*), intent(inout), optional :: context
        real(dp) :: pa(2, 2), pb(2), determinant
        integer :: i

        info = 1
        if (.not. present(context)) return
        select type (preconditioner => context)
        type is (diagonal_preconditioner_t)
            do i = 1, 2
                pa(i, :) = preconditioner%inverse_diagonal(i)*a(i, :)
                pb(i) = preconditioner%inverse_diagonal(i)*b(i)
            end do
            preconditioner%calls = preconditioner%calls + 1
        class default
            return
        end select

        determinant = pa(1, 1)*pa(2, 2) - pa(1, 2)*pa(2, 1)
        if (abs(determinant) <= epsilon(1.0_dp)) return
        x(1) = (pb(1)*pa(2, 2) - pa(1, 2)*pb(2))/determinant
        x(2) = (pa(1, 1)*pb(2) - pb(1)*pa(2, 1))/determinant
        info = 0
    end subroutine diagonal_solve

    ! Test 9: adjoint identity u.(M tp) = tp.(M^T u) for multiroot_jvp/vjp.
    subroutine test_multiroot_dot_product_id(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: jac_x(2, 2), f_p(2, 2), tp(2), u(2)
        real(dp) :: dx(2), jtu(2), lhs, rhs, e
        type(fortnum_status_t) :: st

        jac_x = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        f_p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
        tp = [0.3_dp, -0.5_dp]
        u = [1.7_dp, -0.4_dp]

        call multiroot_jvp(jac_x, f_p, tp, dx, st)
        call multiroot_vjp(jac_x, f_p, u, jtu, st)

        lhs = dot_product(u, dx)
        rhs = dot_product(tp, jtu)
        e = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)

        if (e > 1.0e-13_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [mr_dot_product_id] u.(Mv)=", lhs, &
                " v.(M^T u)=", rhs, " rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_dot_product_id

    ! Test 10: singular Jacobian guard returns FORTNUM_DOMAIN_ERROR.
    subroutine test_multiroot_singular(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: jac_x(2, 2), f_p(2, 2), f_pv(2), tp(2), u(2)
        real(dp) :: dx(2), jtu(2), dxdp(2)
        type(fortnum_status_t) :: st

        jac_x = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2])
        f_p = reshape([-1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp], [2, 2])
        f_pv = [-1.0_dp, 0.0_dp]
        tp = [1.0_dp, 0.0_dp]
        u = [1.0_dp, 0.0_dp]

        call multiroot_jvp(jac_x, f_p, tp, dx, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [mr_singular_jvp] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if

        call multiroot_vjp(jac_x, f_p, u, jtu, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [mr_singular_vjp] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if

        call multiroot_grad(jac_x, f_pv, dxdp, st)
        if (st%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, '(a)') &
                "FAIL [mr_singular_grad] expected FORTNUM_DOMAIN_ERROR"
            nfail = nfail + 1
        end if
    end subroutine test_multiroot_singular

    ! Internal: solve the n-dim system F(x, p_mr) = 0 via multiroot_hybrid.
    subroutine solve_mr(pv, xstar)
        real(dp), intent(in)  :: pv(2)
        real(dp), intent(out) :: xstar(2)
        type(fortnum_status_t) :: st
        p_mr = pv
        call multiroot_hybrid(fdf_mr, 2, [1.0_dp, 1.0_dp], xstar, st, &
            xtol=1.0e-14_dp, ftol=1.0e-14_dp)
    end subroutine solve_mr

    ! F(x) and analytic Jacobian for the n-dim system; p_mr via host association.
    subroutine fdf_mr(x, f, jac, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx
        f(1) = x(1)**2 + x(2) - p_mr(1)
        f(2) = x(1) + x(2)**2 - p_mr(2)
        jac(1, 1) = 2.0_dp * x(1); jac(1, 2) = 1.0_dp
        jac(2, 1) = 1.0_dp; jac(2, 2) = 2.0_dp * x(2)
    end subroutine fdf_mr

    ! Internal: solve x^2 = p_solve on [0, p_solve+1] via root_brent.
    ! p_solve is set by the caller via host association before calling this.
    subroutine solve_x2mp(p_val, xstar)
        real(dp), intent(in)  :: p_val
        real(dp), intent(out) :: xstar
        type(fortnum_status_t) :: st
        p_solve = p_val
        call root_brent(f_x2mp, 0.0_dp, p_solve + 1.0_dp, xstar, st, &
            ftol=1.0e-14_dp)
    end subroutine solve_x2mp

    ! f(x) = x^2 - p_solve, accessed via host association from the program.
    pure real(dp) function f_x2mp(x)
        real(dp), intent(in) :: x
        f_x2mp = x*x - p_solve
    end function f_x2mp

end program test_roots_ad
