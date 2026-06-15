program test_fortnum_integrate_singular
    ! Behavioral tests for the singular and infinite-interval variants:
    ! QAGS endpoint singularities, QAGP interior break points, QAGIU half-line
    ! and full-line transforms. Analytic integrals check value and the frozen
    ! trace; status tests check the domain/break/convergence split that a #40
    ! derivative product reads before differentiating at the subdivision.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok, &
                              FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
                              FORTNUM_CONVERGENCE_ERROR
    use fortnum_integrate, only: integrate_qags, integrate_qagp, &
                                 integrate_qagiu, integrate_workspace_t, &
                                 integrate_epstab_t, integrate_result_t
    implicit none

    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
    integer :: nfail

    nfail = 0

    call test_qags_inv_sqrt()
    call test_qags_log()
    call test_qagp_interior_singularity()
    call test_qagp_two_kinks()
    call test_qagp_empty_reduces_to_qags()
    call test_qagp_break_on_endpoint_is_domain_error()
    call test_qagiu_exp_decay()
    call test_qagiu_lorentz_halfline()
    call test_qagiu_lower_halfline()
    call test_qagiu_doubly_infinite()
    call test_qagiu_bad_inf_is_domain_error()
    call test_trace_frozen_qagp()
    call test_qags_extrapolation_matches_quadpack()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_integrate_singular: all behavioral tests passed"

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

    function f_inv_sqrt(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = 1.0_dp/sqrt(x)
    end function f_inv_sqrt

    function f_log(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = log(x)
    end function f_log

    ! Interior algebraic singularity at 0.5: |x - 0.5|^(-1/2), integral on
    ! [0,1] is 2*sqrt(0.5) + 2*sqrt(0.5) = 2*sqrt(2).
    function f_interior_sing(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = abs(x - 0.5_dp)**(-0.5_dp)
    end function f_interior_sing

    ! Two kinks at 0.3 and 0.7: |x-0.3| + |x-0.7| on [0,1].
    function f_two_kinks(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = abs(x - 0.3_dp) + abs(x - 0.7_dp)
    end function f_two_kinks

    function f_exp_decay(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = exp(-x)
    end function f_exp_decay

    function f_exp_grow(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = exp(x)
    end function f_exp_grow

    function f_lorentz0(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = 1.0_dp/(1.0_dp + x*x)
    end function f_lorentz0

    function f_gauss(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = exp(-x*x)
    end function f_gauss

    ! integral 1/sqrt(x) on [0,1] = 2.
    subroutine test_qags_inv_sqrt()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qags(f_inv_sqrt, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-8_dp, ws, &
                            eps, res, st)
        call check(status_ok(st), "qags 1/sqrt status ok")
        call check_close(res%value, 2.0_dp, 1.0e-7_dp, "qags 1/sqrt value")
    end subroutine test_qags_inv_sqrt

    ! integral log(x) on [0,1] = -1.
    subroutine test_qags_log()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qags(f_log, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-8_dp, ws, &
                            eps, res, st)
        call check(status_ok(st), "qags log status ok")
        call check_close(res%value, -1.0_dp, 1.0e-7_dp, "qags log value")
    end subroutine test_qags_log

    subroutine test_qagp_interior_singularity()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(1)
        pts = [0.5_dp]
        call integrate_qagp(f_interior_sing, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                            1.0e-8_dp, ws, eps, res, st)
        call check(status_ok(st), "qagp interior status ok")
        call check_close(res%value, 2.0_dp*sqrt(2.0_dp), 1.0e-6_dp, &
                         "qagp interior singularity value")
    end subroutine test_qagp_interior_singularity

    ! |x-0.3|+|x-0.7| on [0,1] = 0.58 (analytic piecewise integral).
    subroutine test_qagp_two_kinks()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(2)
        pts = [0.3_dp, 0.7_dp]
        call integrate_qagp(f_two_kinks, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                            1.0e-9_dp, ws, eps, res, st)
        call check(status_ok(st), "qagp two kinks status ok")
        call check_close(res%value, 0.58_dp, 1.0e-8_dp, "qagp two kinks value")
        ! Both break points seeded a panel boundary in the frozen trace.
        call check(res%nsub >= 3, "qagp two kinks nsub >= seeded panels")
    end subroutine test_qagp_two_kinks

    ! Empty break list: QAGP must equal QAGS on the same endpoint singularity.
    subroutine test_qagp_empty_reduces_to_qags()
        type(integrate_workspace_t) :: wsp, wss
        type(integrate_epstab_t)    :: epsp, epss
        type(integrate_result_t)    :: resp, ress
        type(fortnum_status_t)      :: stp, sts
        real(dp) :: pts(0)
        call integrate_qagp(f_inv_sqrt, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                            1.0e-8_dp, wsp, epsp, resp, stp)
        call integrate_qags(f_inv_sqrt, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-8_dp, &
                            wss, epss, ress, sts)
        call check(status_ok(stp), "qagp empty status ok")
        call check_close(resp%value, ress%value, 1.0e-12_dp, &
                         "qagp empty equals qags")
    end subroutine test_qagp_empty_reduces_to_qags

    subroutine test_qagp_break_on_endpoint_is_domain_error()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(1)
        pts = [0.0_dp]   ! coincides with a: break structure ambiguous
        call integrate_qagp(f_inv_sqrt, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                            1.0e-8_dp, ws, eps, res, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, &
                   "qagp break on endpoint is domain error")
    end subroutine test_qagp_break_on_endpoint_is_domain_error

    ! integral exp(-x) on [0,inf) = 1.
    subroutine test_qagiu_exp_decay()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qagiu(f_exp_decay, 0.0_dp, 1, 0.0_dp, 1.0e-8_dp, ws, &
                             eps, res, st)
        call check(status_ok(st), "qagiu exp decay status ok")
        call check_close(res%value, 1.0_dp, 1.0e-7_dp, "qagiu exp decay value")
    end subroutine test_qagiu_exp_decay

    ! integral 1/(1+x^2) on [0,inf) = pi/2.
    subroutine test_qagiu_lorentz_halfline()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qagiu(f_lorentz0, 0.0_dp, 1, 0.0_dp, 1.0e-8_dp, ws, &
                             eps, res, st)
        call check(status_ok(st), "qagiu lorentz status ok")
        call check_close(res%value, 0.5_dp*PI, 1.0e-7_dp, "qagiu lorentz value")
    end subroutine test_qagiu_lorentz_halfline

    ! integral exp(x) on (-inf,0] = 1.
    subroutine test_qagiu_lower_halfline()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qagiu(f_exp_grow, 0.0_dp, -1, 0.0_dp, 1.0e-8_dp, ws, &
                             eps, res, st)
        call check(status_ok(st), "qagiu lower halfline status ok")
        call check_close(res%value, 1.0_dp, 1.0e-7_dp, "qagiu lower value")
    end subroutine test_qagiu_lower_halfline

    ! integral exp(-x^2) on (-inf,inf) = sqrt(pi).
    subroutine test_qagiu_doubly_infinite()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qagiu(f_gauss, 0.0_dp, 2, 0.0_dp, 1.0e-8_dp, ws, &
                             eps, res, st)
        call check(status_ok(st), "qagiu doubly infinite status ok")
        call check_close(res%value, sqrt(PI), 1.0e-7_dp, &
                         "qagiu doubly infinite value")
    end subroutine test_qagiu_doubly_infinite

    subroutine test_qagiu_bad_inf_is_domain_error()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qagiu(f_exp_decay, 0.0_dp, 3, 0.0_dp, 1.0e-8_dp, ws, &
                             eps, res, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, &
                   "qagiu bad inf is domain error")
    end subroutine test_qagiu_bad_inf_is_domain_error

    ! Endpoint-singular QAGS at a tight tolerance must reach the QUADPACK
    ! dqagse accuracy with the QUADPACK subdivision economy: the Wynn-epsilon
    ! extrapolation drives both int log(x) and int 1/sqrt(x) on [0,1] to machine
    ! precision in 6 subintervals (231 GK21 evals), matching Netlib dqagse. The
    ! pre-extrapolation driver drifted to ~1e-12 in 34-66 subintervals at the
    ! same tolerance; this asserts both the value and the panel economy so a
    ! regression to the slow path (the consumer 600 s timeout) is caught.
    subroutine test_qags_extrapolation_matches_quadpack()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        call integrate_qags(f_log, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, ws, &
                            eps, res, st)
        call check(status_ok(st), "qags log tight status ok")
        call check_close(res%value, -1.0_dp, 1.0e-13_dp, "qags log tight value")
        call check(res%nsub <= 8, "qags log tight subdivision economy")
        call integrate_qags(f_inv_sqrt, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-10_dp, ws, &
                            eps, res, st)
        call check(status_ok(st), "qags 1/sqrt tight status ok")
        call check_close(res%value, 2.0_dp, 1.0e-11_dp, "qags 1/sqrt tight value")
        call check(res%nsub <= 8, "qags 1/sqrt tight subdivision economy")
    end subroutine test_qags_extrapolation_matches_quadpack

    ! The QAGP trace is the frozen subdivision: oriented, contiguous, and its
    ! left endpoint set contains the seeded break point (hook #40 re-walks).
    subroutine test_trace_frozen_qagp()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(1)
        integer  :: i
        logical  :: ordered, contiguous_trace, has_break
        pts = [0.5_dp]
        call integrate_qagp(f_interior_sing, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                            1.0e-8_dp, ws, eps, res, st)
        call check(res%nsub >= 2, "qagp trace nsub >= 2")
        ordered = .true.
        contiguous_trace = .true.
        has_break = .false.
        do i = 1, res%nsub
            if (res%sub_b(i) < res%sub_a(i)) ordered = .false.
            if (i > 1) then
                if (abs(res%sub_a(i) - res%sub_b(i-1)) > 1.0e-13_dp) &
                    contiguous_trace = .false.
            end if
            if (abs(res%sub_a(i) - 0.5_dp) <= 1.0e-13_dp) has_break = .true.
        end do
        call check(ordered, "qagp trace oriented a<=b")
        call check(contiguous_trace, "qagp trace contiguous")
        call check(has_break, "qagp trace keeps the seeded break boundary")
    end subroutine test_trace_frozen_qagp

end program test_fortnum_integrate_singular
