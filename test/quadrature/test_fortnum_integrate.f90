program test_fortnum_integrate
    ! Behavioral tests for the QAG/QAGS adaptive integrator. Analytic integrals
    ! check value and reported error; status tests check the domain/convergence
    ! split; trace tests check the frozen-subdivision hooks #40 needs.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_CONVERGENCE_ERROR
    use fortnum_integrate, only: integrate, integrate_qag, integrate_qags, &
        integrate_workspace_t, integrate_epstab_t, &
        integrate_result_t
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer :: nfail

    nfail = 0

    call test_polynomial()
    call test_exp_wrapper()
    call test_oscillatory()
    call test_lorentzian_qags()
    call test_ctx_threading()
    call test_trace_consistency()
    call test_value_invariant_to_workspace_reuse()
    call test_reject_reversed_interval()
    call test_reject_bad_key()
    call test_convergence_failure_low_limit()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_integrate: all behavioral tests passed"

contains

    subroutine check(cond, name)
        logical,      intent(in) :: cond
        character(*), intent(in) :: name
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: "//name
        end if
    end subroutine check

    subroutine check_close(got, ref, tol, name)
        real(dp),     intent(in) :: got, ref, tol
        character(*), intent(in) :: name
        if (.not. (abs(got - ref) <= tol)) then
            nfail = nfail + 1
            write (error_unit, "(a,3(a,es20.12e3))") "FAIL: "//name, &
                "  got=", got, "  ref=", ref, "  diff=", abs(got - ref)
        end if
    end subroutine check_close

    ! 3x^2 + 2x + 1 on [0,2] = 14, GK integrates a cubic exactly on one panel.
    function f_poly(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = 3.0_dp*x**2 + 2.0_dp*x + 1.0_dp
    end function f_poly

    function f_exp(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = exp(x)
    end function f_exp

    function f_osc(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = sin(10.0_dp*x)
    end function f_osc

    ! Peaked Lorentzian 1/(1 + ((x-c)/w)^2) with c, w threaded via ctx.
    function f_lorentz(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        real(dp) :: c, w
        c = 0.5_dp
        w = 1.0e-3_dp
        if (present(ctx)) then
            select type (ctx)
                type is (real(dp))
                w = ctx
            end select
        end if
        fx = 1.0_dp/(1.0_dp + ((x - c)/w)**2)
    end function f_lorentz

    subroutine test_polynomial()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qag(f_poly, 0.0_dp, 2.0_dp, 0.0_dp, 1.0e-10_dp, ws, &
            res, st)
        call check(status_ok(st), "polynomial status ok")
        call check_close(res%value, 14.0_dp, 1.0e-10_dp, "polynomial value")
        call check(res%abserr <= 1.0e-10_dp, "polynomial abserr small")
        call check(res%nsub >= 1, "polynomial nsub >= 1")
        call check(res%neval >= 21, "polynomial neval >= 21")
    end subroutine test_polynomial

    subroutine test_exp_wrapper()
        real(dp) :: value
        type(fortnum_status_t) :: st
        call integrate(f_exp, 0.0_dp, 1.0_dp, value, st)
        call check(status_ok(st), "exp wrapper status ok")
        call check_close(value, exp(1.0_dp) - 1.0_dp, 1.0e-8_dp, "exp value")
    end subroutine test_exp_wrapper

    ! sin(10x) on [0,pi] = (1 - cos(10 pi))/10 = 0.
    subroutine test_oscillatory()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: ref
        ! True value is 0, so a pure relative tolerance is unattainable; an
        ! absolute tolerance is the correct request for a zero-mean integrand.
        ref = (1.0_dp - cos(10.0_dp*PI))/10.0_dp
        call integrate_qag(f_osc, 0.0_dp, PI, 1.0e-10_dp, 1.0e-10_dp, ws, res, &
            st, key=31)
        call check(status_ok(st), "oscillatory status ok")
        call check_close(res%value, ref, 1.0e-9_dp, "oscillatory value")
        call check(res%key == 31, "oscillatory records key 31")
    end subroutine test_oscillatory

    ! Peaked Lorentzian on [0,1]: w*(atan((1-c)/w) + atan(c/w)).
    subroutine test_lorentzian_qags()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: w, c, ref
        w = 1.0e-3_dp
        c = 0.5_dp
        ref = w*(atan((1.0_dp - c)/w) + atan(c/w))
        call integrate_qags(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-8_dp, ws, &
            eps, res, st, ctx=w)
        call check(status_ok(st), "lorentzian qags status ok")
        call check_close(res%value, ref, 1.0e-6_dp, "lorentzian value")
    end subroutine test_lorentzian_qags

    ! Same integrand, two ctx widths give the two analytic values: the ctx
    ! channel actually reaches the integrand, no hidden state.
    subroutine test_ctx_threading()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: w, c, ref
        c = 0.5_dp
        w = 0.1_dp
        ref = w*(atan((1.0_dp - c)/w) + atan(c/w))
        call integrate_qags(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-9_dp, ws, &
            eps, res, st, ctx=w)
        call check_close(res%value, ref, 1.0e-8_dp, "ctx width 0.1 value")
    end subroutine test_ctx_threading

    ! Trace is the frozen subdivision: ordered left-to-right, contiguous, and
    ! its sub_r sum matches the reported value (the trace_rule schedule, #40).
    subroutine test_trace_consistency()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        integer  :: i
        logical  :: ordered, contiguous_trace
        real(dp) :: s
        call integrate_qag(f_osc, 0.0_dp, PI, 1.0e-10_dp, 1.0e-10_dp, ws, &
            res, st)
        call check(res%nsub >= 1, "trace nsub >= 1")
        ordered = .true.
        contiguous_trace = .true.
        do i = 1, res%nsub
            if (res%sub_b(i) < res%sub_a(i)) ordered = .false.
            if (i > 1) then
                if (abs(res%sub_a(i) - res%sub_b(i-1)) > 1.0e-13_dp) &
                    contiguous_trace = .false.
            end if
        end do
        call check(ordered, "trace subintervals oriented a<=b")
        call check(contiguous_trace, "trace subintervals contiguous")
        s = 0.0_dp
        do i = 1, res%nsub
            s = s + res%sub_r(i)
        end do
        call check_close(s, res%value, 1.0e-12_dp, "trace sub_r sums to value")
        call check(.not. res%extrapolated, "qag trace not extrapolated")
    end subroutine test_trace_consistency

    ! Reusing one workspace across calls must not change the result.
    subroutine test_value_invariant_to_workspace_reuse()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: v1, v2
        call integrate_qag(f_osc, 0.0_dp, PI, 1.0e-10_dp, 1.0e-10_dp, ws, &
            res, st)
        v1 = res%value
        call integrate_qag(f_poly, 0.0_dp, 2.0_dp, 0.0_dp, 1.0e-10_dp, ws, &
            res, st)
        call integrate_qag(f_osc, 0.0_dp, PI, 1.0e-10_dp, 1.0e-10_dp, ws, &
            res, st)
        v2 = res%value
        call check_close(v1, v2, 1.0e-14_dp, "workspace reuse invariant")
    end subroutine test_value_invariant_to_workspace_reuse

    subroutine test_reject_reversed_interval()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qag(f_poly, 2.0_dp, 0.0_dp, 0.0_dp, 1.0e-8_dp, ws, &
            res, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, &
            "reversed interval is domain error")
    end subroutine test_reject_reversed_interval

    subroutine test_reject_bad_key()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qag(f_poly, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-8_dp, ws, &
            res, st, key=7)
        call check(st%code == FORTNUM_DOMAIN_ERROR, "bad key is domain error")
    end subroutine test_reject_bad_key

    ! limit=1 on an integrand the single panel cannot resolve: convergence
    ! failure, not a domain error, and the value is still the best estimate.
    subroutine test_convergence_failure_low_limit()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qag(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, ws, &
            res, st, limit=1)
        call check(st%code == FORTNUM_CONVERGENCE_ERROR, &
            "limit=1 hard integrand is convergence error")
        call check(res%nsub >= 1, "failed run still records a trace")
    end subroutine test_convergence_failure_low_limit

end program test_fortnum_integrate
