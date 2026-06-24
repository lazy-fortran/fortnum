program test_fortnum_multiroot_rc
    ! Behavioral tests for the reverse-communication multiroot stepper.
    !
    ! The RC stepper must reproduce the host multiroot_hybrid iterate path: the
    ! caller drives a NEED_FJ loop, evaluating the residual and analytic
    ! Jacobian inline, and the converged root must match multiroot_hybrid on the
    ! same fixture to tight tolerance.  A second case drives the RC loop on a
    ! 2x2 system and checks the residual at the returned root, and a third
    ! checks the singular-Jacobian failure path.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_multiroot_rc, only: multiroot_rc_t, multiroot_rc_init, &
                                    multiroot_step, MULTIROOT_NEED_FJ, MULTIROOT_DONE, &
                                    MULTIROOT_FAILED, MULTIROOT_RC_SINGULAR
    implicit none

    integer :: nfail
    nfail = 0

    call test_rc_vs_host_3x3(nfail)
    call test_rc_residual_2x2(nfail)
    call test_rc_singular(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! Residual and analytic Jacobian of the 3x3 test system, shared by the host
    ! callback and the inline RC evaluation so both solve the identical map:
    !   f1 = x1^2 + x2^2 + x3^2 - 3
    !   f2 = x1 + x2 - x3 - 1
    !   f3 = x1*x2*x3 - 1
    ! Root near (1,1,1).
    pure subroutine sys3(x, f, jac)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: f(3), jac(3, 3)
        f(1) = x(1)**2 + x(2)**2 + x(3)**2 - 3.0_dp
        f(2) = x(1) + x(2) - x(3) - 1.0_dp
        f(3) = x(1)*x(2)*x(3) - 1.0_dp
        jac(1, 1) = 2.0_dp*x(1); jac(1, 2) = 2.0_dp*x(2); jac(1, 3) = 2.0_dp*x(3)
        jac(2, 1) = 1.0_dp; jac(2, 2) = 1.0_dp; jac(2, 3) = -1.0_dp
        jac(3, 1) = x(2)*x(3); jac(3, 2) = x(1)*x(3); jac(3, 3) = x(1)*x(2)
    end subroutine sys3

    subroutine sys3_cb(x, f, jac, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx
        call sys3(x, f, jac)
    end subroutine sys3_cb

    subroutine test_rc_vs_host_3x3(nfail)
        integer, intent(inout) :: nfail
        type(multiroot_rc_t) :: st
        type(fortnum_status_t) :: status
        real(dp) :: x0(3), xhost(3), f(3), jac(3, 3)
        integer :: action, steps

        x0 = [1.5_dp, 0.8_dp, 1.3_dp]

        ! Host reference.
        call multiroot_hybrid(sys3_cb, 3, x0, xhost, status)
        if (status%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL host solve did not converge"
            return
        end if

        ! RC drive: evaluate F,J at st%x, step, repeat on NEED_FJ.
        call multiroot_rc_init(st, 3, x0)
        steps = 0
        do
            call sys3(st%x(1:3), f, jac)
            call multiroot_step(st, f, jac, action)
            steps = steps + 1
            if (action /= MULTIROOT_NEED_FJ) exit
            if (steps > 100000) exit
        end do

        if (action /= MULTIROOT_DONE) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL RC action=", action
            return
        end if

        ! RC root must match the host root.
        call report("rc_vs_host", maxval(abs(st%x(1:3) - xhost)), 1.0e-9_dp, nfail)
        ! Residual at the RC root must be tiny.
        call sys3(st%x(1:3), f, jac)
        call report("rc_residual", maxval(abs(f)), 1.0e-9_dp, nfail)
        write (*, "(a,es10.2,a,es10.2,a,i0)") &
            "rc_vs_host: |x_rc-x_host|_inf=", maxval(abs(st%x(1:3) - xhost)), &
            "  |F(x_rc)|_inf=", maxval(abs(f)), "  steps=", steps
    end subroutine test_rc_vs_host_3x3

    subroutine test_rc_residual_2x2(nfail)
        integer, intent(inout) :: nfail
        type(multiroot_rc_t) :: st
        real(dp) :: x0(2), f(2), jac(2, 2)
        integer :: action, steps
        ! f1 = x1^2 - x2 - 1 ; f2 = x1 + x2^2 - 3 ; root near (1.2, 0.5).
        x0 = [1.0_dp, 1.0_dp]
        call multiroot_rc_init(st, 2, x0)
        steps = 0
        do
            f(1) = st%x(1)**2 - st%x(2) - 1.0_dp
            f(2) = st%x(1) + st%x(2)**2 - 3.0_dp
            jac(1, 1) = 2.0_dp*st%x(1); jac(1, 2) = -1.0_dp
            jac(2, 1) = 1.0_dp; jac(2, 2) = 2.0_dp*st%x(2)
            call multiroot_step(st, f, jac, action)
            steps = steps + 1
            if (action /= MULTIROOT_NEED_FJ) exit
            if (steps > 100000) exit
        end do
        if (action /= MULTIROOT_DONE) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL RC 2x2 action=", action
            return
        end if
        f(1) = st%x(1)**2 - st%x(2) - 1.0_dp
        f(2) = st%x(1) + st%x(2)**2 - 3.0_dp
        call report("rc_2x2_residual", maxval(abs(f)), 1.0e-9_dp, nfail)
    end subroutine test_rc_residual_2x2

    subroutine test_rc_singular(nfail)
        integer, intent(inout) :: nfail
        type(multiroot_rc_t) :: st
        real(dp) :: x0(2), f(2), jac(2, 2)
        integer :: action
        ! A singular Jacobian (rank 1) at the start must fail with SINGULAR.
        x0 = [0.5_dp, 0.5_dp]
        call multiroot_rc_init(st, 2, x0)
        f = [1.0_dp, 1.0_dp]
        jac = reshape([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], [2, 2])  ! col2 = 2*col1
        call multiroot_step(st, f, jac, action)
        if (action /= MULTIROOT_FAILED) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL RC singular action=", action
            return
        end if
        if (st%fail_code /= MULTIROOT_RC_SINGULAR) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL RC singular fail_code=", st%fail_code
        end if
    end subroutine test_rc_singular

    subroutine report(label, err, atol, nfail)
        character(*), intent(in) :: label
        real(dp), intent(in) :: err, atol
        integer, intent(inout) :: nfail
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es12.4,a,es12.4)") "FAIL [", label, &
                "] err=", err, " atol=", atol
        end if
    end subroutine report

end program test_fortnum_multiroot_rc
