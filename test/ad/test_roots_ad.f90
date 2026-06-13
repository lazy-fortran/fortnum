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
    use fortnum_roots,  only: root_brent, root_jvp, root_vjp, root_grad
    implicit none

    ! p_solve is shared between solve_x2mp and f_x2mp via host association.
    real(dp) :: p_solve

    integer :: nfail
    nfail = 0

    call test_grad_analytic(nfail)
    call test_grad_vs_fd(nfail)
    call test_jvp_scalar_p(nfail)
    call test_dot_product_id(nfail)
    call test_near_multiple_root(nfail)
    call test_vector_p(nfail)

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
        xstar  = sqrt(p)         ! 2.0
        f_x    = 2.0_dp * xstar  ! 4.0
        f_p    = -1.0_dp
        dxdp_exact = 1.0_dp / (2.0_dp * xstar)  ! 0.25

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
        dx_exact = (-f_p(1) * dp_vec(1)) / f_x  ! = 0.7 / (2 sqrt(2))

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
