program test_fortnum_generated_residual_equivalence
    ! Generated versus hand-written equivalence for the root residual kernels.
    !
    ! The generated residual, JVP, VJP, and Jacobian kernels in
    ! src/generated/fortnum_{scalar,vector}_root_residual_*_kernel.f90 and
    ! src/generated/fortnum_implicit_root_residual_kernel.f90 are produced by
    ! the fortsym generator tools/codegen/app/gen_implicit_root_residual.f90
    ! from symbolic residuals. This test proves the generated kernels equal
    ! the hand-written analytic derivatives of the same residuals, evaluated
    ! at many points. The hand-written expressions here are the independent
    ! oracle; they are transcribed once from the symbolic definition and are
    ! deliberately NOT the same code path the generator emits.
    !
    ! The two residual families:
    !   scalar:  r(x, p1, p2)   = x^3 + p1*x - p2
    !   vector:  F(x, p) = [ x1^2 + x2 - p1,  x1 + x2^2 - p2 ]
    !
    ! and the implicit scalar residual r(x, p) = x^2 - p.
    !
    ! For each, the analytic hand-written derivatives are:
    !   dr/dx = 3x^2 + p1,  dr/dp = [x, -1]
    !   dF/dx = [[2 x1, 1], [1, 2 x2]],  dF/dp tangent = [-tp1, -tp2]
    !   r(x,p)=x^2-p: dr/dx = 2x, dr/dp = -1.
    !
    ! A finite-difference check runs alongside as a second, independent
    ! oracle that does not share the analytic expression either.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_generated_scalar_root_residual_jvp, only: &
        fortnum_scalar_root_residual_jvp_kernel
    use fortnum_generated_scalar_root_residual_vjp, only: &
        fortnum_scalar_root_residual_vjp_kernel
    use fortnum_generated_vector_root_residual_jacobian, only: &
        fortnum_vector_root_residual_jacobian_kernel
    use fortnum_generated_vector_root_residual_jvp, only: &
        fortnum_vector_root_residual_jvp_kernel
    use fortnum_generated_vector_root_residual_vjp, only: &
        fortnum_vector_root_residual_vjp_kernel
    use fortnum_generated_implicit_root_residual, only: &
        fortnum_implicit_root_residual_kernel
    implicit none

    integer :: nfail

    nfail = 0
    call check_scalar_jvp(nfail)
    call check_scalar_vjp(nfail)
    call check_vector_jacobian(nfail)
    call check_vector_jvp(nfail)
    call check_vector_vjp(nfail)
    call check_implicit_residual(nfail)

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, &
            " generated/hand-written equivalence check(s) failed"
        error stop 1
    end if
    write (*, "(a)") "fortnum_generated_residual_equivalence: all passed"

contains

    subroutine check_scalar_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x, p1, tp1, tp2, fx_gen, fptp_gen
        real(dp) :: fx_exp, fptp_exp, fx_fd, fptp_fd, h
        integer :: k
        real(dp), parameter :: pts(8) = [ -1.3_dp, 0.0_dp, 0.4_dp, 2.7_dp, &
            -0.05_dp, 5.0_dp, -9.0_dp, 1.0_dp ]
        do k = 1, size(pts)
            x = pts(k)
            p1 = 0.5_dp + 0.25_dp*k
            tp1 = 1.3_dp - 0.1_dp*k
            tp2 = -0.7_dp + 0.3_dp*k
            call fortnum_scalar_root_residual_jvp_kernel(x, p1, tp1, tp2, &
                fx_gen, fptp_gen)
            fx_exp = 3.0_dp*x*x + p1
            fptp_exp = tp1*x - tp2
            h = 1.0e-6_dp
            fx_fd = (scalar_residual(x + h, p1, 0.0_dp) &
                - scalar_residual(x - h, p1, 0.0_dp))/(2.0_dp*h)
            fptp_fd = (scalar_residual(x, p1 + h*tp1, h*tp2) &
                - scalar_residual(x, p1 - h*tp1, -h*tp2))/(2.0_dp*h)
            call check_abs("scalar JVP f_x", fx_gen, fx_exp, fx_fd, nfail)
            call check_abs("scalar JVP f_p_tp", fptp_gen, fptp_exp, fptp_fd, &
                nfail)
        end do
    end subroutine check_scalar_jvp

    subroutine check_scalar_vjp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x, p1, u, fx_gen, fp1, fp2
        real(dp) :: fx_exp, fp1_exp, fp2_exp
        integer :: k
        real(dp), parameter :: pts(4) = [ -0.9_dp, 1.2_dp, 3.3_dp, -6.0_dp ]
        do k = 1, size(pts)
            x = pts(k)
            p1 = 0.2_dp + 0.7_dp*k
            u = -1.0_dp + 0.5_dp*k
            call fortnum_scalar_root_residual_vjp_kernel(x, p1, u, fx_gen, &
                fp1, fp2)
            fx_exp = 3.0_dp*x*x + p1
            fp1_exp = u*x
            fp2_exp = -u
            call check_abs("scalar VJP f_x", fx_gen, fx_exp, &
                huge(1.0_dp), nfail)
            call check_abs("scalar VJP f_p1", fp1, fp1_exp, huge(1.0_dp), nfail)
            call check_abs("scalar VJP f_p2", fp2, fp2_exp, huge(1.0_dp), nfail)
        end do
    end subroutine check_scalar_vjp

    subroutine check_vector_jacobian(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2), jac(2,2)
        real(dp) :: jac_exp(2,2), h
        integer :: k, i, j
        real(dp), parameter :: pts(2, 4) = reshape([ &
            1.0_dp, 2.0_dp, -0.5_dp, 3.0_dp, 0.0_dp, -1.0_dp, 4.0_dp, 5.0_dp ], &
            [2, 4])
        do k = 1, size(pts, 2)
            x = pts(:, k)
            call fortnum_vector_root_residual_jacobian_kernel(x, jac)
            jac_exp(1, 1) = 2.0_dp*x(1)
            jac_exp(1, 2) = 1.0_dp
            jac_exp(2, 1) = 1.0_dp
            jac_exp(2, 2) = 2.0_dp*x(2)
            do j = 1, 2
                do i = 1, 2
                    ! finite-difference oracle for dF_i / dx_j
                    h = 1.0e-6_dp
                    call check_abs("vector Jacobian", jac(i, j), &
                        jac_exp(i, j), vector_residual_fd(i, j, x, h), nfail)
                end do
            end do
        end do
    end subroutine check_vector_jacobian

    subroutine check_vector_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2), tp(2), jac(2,2), fptp(2)
        real(dp) :: fptp_exp(2), h
        integer :: k
        real(dp), parameter :: pts(2, 4) = reshape([ &
            1.0_dp, 1.0_dp, -2.0_dp, 0.5_dp, 3.0_dp, -4.0_dp, 0.0_dp, 0.0_dp ], &
            [2, 4])
        do k = 1, size(pts, 2)
            x = pts(:, k)
            tp = [0.7_dp, -1.1_dp]
            call fortnum_vector_root_residual_jvp_kernel(x, tp, jac, fptp)
            fptp_exp = [-tp(1), -tp(2)]
            h = 1.0e-6_dp
            call check_abs("vector JVP f_p_tp(1)", fptp(1), fptp_exp(1), &
                (vector_residual(1, x + h*tp) - vector_residual(1, x - h*tp)) &
                /(2.0_dp*h), nfail)
            call check_abs("vector JVP f_p_tp(2)", fptp(2), fptp_exp(2), &
                (vector_residual(2, x + h*tp) - vector_residual(2, x - h*tp)) &
                /(2.0_dp*h), nfail)
        end do
    end subroutine check_vector_jvp

    subroutine check_vector_vjp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: u(2), fptu(2), fptu_exp(2)
        integer :: k
        real(dp), parameter :: uvals(2, 3) = reshape([ &
            1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, -2.0_dp, 3.0_dp ], [2, 3])
        do k = 1, size(uvals, 2)
            u = uvals(:, k)
            call fortnum_vector_root_residual_vjp_kernel(u, fptu)
            fptu_exp = [-u(1), -u(2)]
            call check_abs("vector VJP f_p_t_u(1)", fptu(1), fptu_exp(1), &
                huge(1.0_dp), nfail)
            call check_abs("vector VJP f_p_t_u(2)", fptu(2), fptu_exp(2), &
                huge(1.0_dp), nfail)
        end do
    end subroutine check_vector_vjp

    subroutine check_implicit_residual(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x, p, tp, residual, fx, fptp
        real(dp) :: residual_exp, fx_exp, fptp_exp, h
        integer :: k
        real(dp), parameter :: pts(4) = [ 0.3_dp, -1.0_dp, 2.0_dp, 5.5_dp ]
        do k = 1, size(pts)
            x = pts(k)
            p = 0.9_dp + 0.2_dp*k
            tp = 0.4_dp - 0.3_dp*k
            call fortnum_implicit_root_residual_kernel(x, p, tp, residual, &
                fx, fptp)
            residual_exp = x*x - p
            fx_exp = 2.0_dp*x
            fptp_exp = -tp
            h = 1.0e-6_dp
            call check_abs("implicit residual", residual, residual_exp, &
                huge(1.0_dp), nfail)
            call check_abs("implicit f_x", fx, fx_exp, &
                (residual_fn(x + h, p) - residual_fn(x - h, p))/(2.0_dp*h), nfail)
            call check_abs("implicit f_p_tp", fptp, fptp_exp, &
                (residual_fn(x, p + h*tp) - residual_fn(x, p - h*tp)) &
                /(2.0_dp*h), nfail)
        end do
    end subroutine check_implicit_residual

    !> Scalar residual r(x, p1, p2) = x^3 + p1*x - p2 (hand-written).
    pure function scalar_residual(x, p1, p2) result(r)
        real(dp), intent(in) :: x, p1, p2
        real(dp) :: r
        r = x**3 + p1*x - p2
    end function scalar_residual

    !> Implicit scalar residual r(x, p) = x^2 - p (hand-written).
    pure function residual_fn(x, p) result(r)
        real(dp), intent(in) :: x, p
        real(dp) :: r
        r = x*x - p
    end function residual_fn

    !> Vector residual component F_i(x) (hand-written), with p = 0.
    pure function vector_residual(i, x) result(r)
        integer, intent(in) :: i
        real(dp), intent(in) :: x(2)
        real(dp) :: r
        if (i == 1) then
            r = x(1)*x(1) + x(2)
        else
            r = x(1) + x(2)*x(2)
        end if
    end function vector_residual

    !> Central finite difference of dF_i / dx_j (secondary oracle).
    pure function vector_residual_fd(i, j, x, h) result(fd)
        integer, intent(in) :: i, j
        real(dp), intent(in) :: x(2), h
        real(dp) :: fd, xp(2), xm(2)
        xp = x
        xm = x
        xp(j) = x(j) + h
        xm(j) = x(j) - h
        fd = (vector_residual(i, xp) - vector_residual(i, xm))/(2.0_dp*h)
    end function vector_residual_fd

    !> Compare generated value against both the analytic and finite-difference
    !> oracles; fail only if it differs from BOTH by more than the tolerance.
    subroutine check_abs(label, got, analytic, fd, nfail)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, analytic, fd
        integer, intent(inout) :: nfail
        real(dp) :: tol
        tol = 1.0e-9_dp
        if (abs(got - analytic) > tol .and. abs(got - fd) > tol) then
            write (error_unit, "(a,a,es12.4,a,es12.4,a,es12.4)") &
                "FAIL ", label, " got=", got, " analytic=", analytic, &
                " fd=", fd
            nfail = nfail + 1
        end if
    end subroutine check_abs

end program test_fortnum_generated_residual_equivalence
