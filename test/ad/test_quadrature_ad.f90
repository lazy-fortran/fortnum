program test_quadrature_ad
    ! Derivative tests for fortnum_quadrature (ad.md §4: transparent policy).
    ! Map under test: f -> I = sum_i w_i f_i   (linear, R^n -> R^1).
    ! Jacobian: J = w^T  (1 x n row vector; dI/df_i = w_i).
    !
    ! Tests:
    !   (1) gauss_legendre_grad returns the weight vector.
    !   (2) gauss_legendre_jvp matches central finite difference.
    !   (3) Dot-product identity  u.(J v) = v.(J^T u).
    !   (4) Known-function integration: sum w_i f(x_i) matches analytic integral.
    !   (5) gauss_gen_laguerre transparent policy: the same linear-map products
    !       apply to a generalized-Laguerre weight vector (ad.md transparent).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab, &
        gauss_legendre_jvp, gauss_legendre_vjp, gauss_legendre_grad, &
        gauss_gen_laguerre
    use fortnum_ad_test_utils, only: rel_err
    implicit none

    integer :: nfail
    nfail = 0

    call test_grad_equals_weights(nfail)
    call test_jvp_vs_fd(nfail)
    call test_dot_product_identity(nfail)
    call test_known_integral(nfail)
    call test_gen_laguerre_transparent(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! (1) gauss_legendre_grad must return the weight vector exactly.
    ! dI/df_i = w_i; grad is the u=1 special case of vjp.
    subroutine test_grad_equals_weights(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 7
        real(dp) :: x(n), w(n), grad(n)
        integer  :: i
        logical  :: ok

        call gauss_legendre(n, x, w)
        call gauss_legendre_grad(w, grad)

        ok = .true.
        do i = 1, n
            if (grad(i) /= w(i)) then
                write (error_unit, '(a,i0,a,es24.16,a,es24.16)') &
                    "FAIL [grad_equals_weights] i=", i, &
                    " grad=", grad(i), " w=", w(i)
                ok = .false.
            end if
        end do
        if (.not. ok) nfail = nfail + 1
    end subroutine test_grad_equals_weights

    ! (2) gauss_legendre_jvp vs central finite difference.
    ! Map: f_values -> I = sum_i w_i f_i.  FD perturbation is on f_values.
    ! Central FD: dI ~ (I(f+hv) - I(f-hv)) / (2h) = sum_i w_i v_i.
    subroutine test_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 6
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp), parameter :: tol = 1.0e-9_dp

        real(dp) :: x(n), w(n), f(n), v(n), jv(1)
        real(dp) :: fp(n), fm(n), I_p, I_m, fd_val, err
        integer  :: i

        call gauss_legendre(n, x, w)
        ! f_i = exp(x_i) at the Gauss nodes.
        do i = 1, n
            f(i) = exp(x(i))
        end do
        ! tangent: v_i = sin(i)  (generic, not aligned with w)
        do i = 1, n
            v(i) = sin(real(i, dp))
        end do

        call gauss_legendre_jvp(w, v, jv)

        ! FD: perturb f -> I(f +/- hv)
        fp = f + h*v
        fm = f - h*v
        I_p = dot_product(w, fp)
        I_m = dot_product(w, fm)
        fd_val = (I_p - I_m) / (2.0_dp*h)

        err = rel_err(jv(1), fd_val)
        if (err > tol) then
            write (error_unit, '(a,es12.4,a,es12.4,a,es12.4)') &
                "FAIL [jvp_vs_fd] jvp=", jv(1), " fd=", fd_val, &
                " rel_err=", err
            nfail = nfail + 1
        end if
    end subroutine test_jvp_vs_fd

    ! (3) Dot-product identity  u . (J v) = v . (J^T u).
    ! J = w^T; J v = dot(w,v) (scalar); J^T u = u(1)*w.
    ! Identity: u(1) * dot(w,v) == dot(v, u(1)*w) -- exact for any u,v,w.
    subroutine test_dot_product_identity(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 5
        real(dp), parameter :: tol = 4.0_dp*epsilon(1.0_dp)

        real(dp) :: x(n), w(n), v(n), u(1), jv(1), jtu(n)
        real(dp) :: lhs, rhs, e
        integer  :: i

        call gauss_legendre(n, x, w)
        do i = 1, n
            v(i) = real(i, dp) * 0.3_dp - 0.7_dp
        end do
        u(1) = 1.7_dp

        call gauss_legendre_jvp(w, v, jv)
        call gauss_legendre_vjp(w, u, jtu)

        lhs = u(1) * jv(1) ! u . (J v)
        rhs = dot_product(v, jtu) ! v . (J^T u)

        e = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)
        if (e > tol) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [dot_product_identity] lhs=", lhs, " rhs=", rhs, &
                " rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_dot_product_identity

    ! (4) Integrating a known function: sum w_i f(x_i) must match analytic value.
    ! Uses gauss_legendre_ab on [0, pi/2]; integral of sin(x) = 1.
    ! A 16-point rule integrates this to near machine precision.
    subroutine test_known_integral(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 16
        real(dp), parameter :: pi = 3.14159265358979324_dp
        real(dp), parameter :: a  = 0.0_dp, b = 0.5_dp*pi
        real(dp), parameter :: tol = 1.0e-13_dp

        real(dp) :: x(n), w(n), I_quad, err
        integer  :: i

        call gauss_legendre_ab(n, a, b, x, w)
        I_quad = 0.0_dp
        do i = 1, n
            I_quad = I_quad + w(i) * sin(x(i))
        end do

        ! Verify using the gradient: grad_i = w_i; I = dot(grad, f) = dot(w, f).
        ! This is identical to the direct sum -- the check is that the derivative
        ! interpretation (sum w_i f_i) matches the integral.
        err = abs(I_quad - 1.0_dp)
        if (err > tol) then
            write (error_unit, '(a,es24.16,a,es12.4)') &
                "FAIL [known_integral] I=", I_quad, " err=", err
            nfail = nfail + 1
        end if
    end subroutine test_known_integral

    ! (5) gauss_gen_laguerre is transparent w.r.t. integrand values, exactly as
    ! gauss_legendre: I = sum_i w_i f_i is linear with Jacobian w^T, regardless
    ! of how the rule (here generalized Laguerre, alpha = 5/2) produced w. The
    ! caller holds nodes/weights fixed (libneo transport calc_D_one_over_nu takes
    ! w, x as intent(in)), so the legendre jvp/vjp/grad products apply verbatim.
    ! Checks: grad equals the gen-Laguerre weights; jvp matches central FD;
    ! the dot-product adjoint identity holds on the gen-Laguerre weight vector.
    subroutine test_gen_laguerre_transparent(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 8
        real(dp), parameter :: alpha = 2.5_dp
        real(dp), parameter :: h = 1.0e-6_dp
        real(dp), parameter :: tol_fd = 1.0e-8_dp
        real(dp), parameter :: tol_id = 4.0_dp*epsilon(1.0_dp)

        real(dp) :: x(n), w(n), grad(n), v(n), u(1), jv(1), jtu(n)
        real(dp) :: f(n), fp(n), fm(n), fd_val, lhs, rhs, e
        integer  :: i
        logical  :: ok

        call gauss_gen_laguerre(n, alpha, x, w)

        ! grad must equal the weight vector.
        call gauss_legendre_grad(w, grad)
        ok = .true.
        do i = 1, n
            if (grad(i) /= w(i)) ok = .false.
        end do
        if (.not. ok) then
            write (error_unit, '(a)') "FAIL [gen_laguerre grad/=weights]"
            nfail = nfail + 1
        end if

        ! jvp vs central FD of I = sum_i w_i f_i on f_i = exp(-x_i/4).
        do i = 1, n
            f(i) = exp(-0.25_dp*x(i))
            v(i) = cos(real(i, dp))
        end do
        call gauss_legendre_jvp(w, v, jv)
        fp = f + h*v
        fm = f - h*v
        fd_val = (dot_product(w, fp) - dot_product(w, fm))/(2.0_dp*h)
        e = rel_err(jv(1), fd_val)
        if (e > tol_fd) then
            write (error_unit, '(a,es12.4,a,es12.4,a,es12.4)') &
                "FAIL [gen_laguerre jvp_vs_fd] jvp=", jv(1), " fd=", fd_val, &
                " rel_err=", e
            nfail = nfail + 1
        end if

        ! Dot-product adjoint identity u.(Jv) = v.(J^T u).
        u(1) = 0.9_dp
        call gauss_legendre_vjp(w, u, jtu)
        lhs = u(1)*jv(1)
        rhs = dot_product(v, jtu)
        e = abs(lhs - rhs)/max(abs(lhs), abs(rhs), 1.0_dp)
        if (e > tol_id) then
            write (error_unit, '(a,es12.4)') &
                "FAIL [gen_laguerre dot_identity] rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_gen_laguerre_transparent

end program test_quadrature_ad
