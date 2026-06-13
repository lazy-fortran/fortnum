program test_integrate_ad
    ! Trace-rule forward product for adaptive integration (issue #40).
    !
    ! Smooth case: I(p) = integral_0^1 exp(p x) dx = (exp(p) - 1)/p, with the
    ! known closed form dI/dp = integral_0^1 x exp(p x) dx
    !                         = (exp(p)(p - 1) + 1)/p^2.
    ! The tangent integrand is df/dp = x exp(p x). integrate_qag_jvp re-walks the
    ! frozen subdivision the primal accepted and quadratures the tangent on it.
    ! We check the trace JVP against (a) the analytic dI/dp and (b) the central
    ! finite difference of the primal I(p) computed on the SAME subdivision (the
    ! base and perturbed runs land on the same accepted panels), so the FD probes
    ! the frozen-trace derivative the JVP claims.
    !
    ! Non-smooth case: a parameter-placed inverse-square-root spike
    ! g(x,p) = 1/sqrt(|x - p|) on [0,1] with the singularity interior. The
    ! adaptive primal cannot reach tolerance and records a non-OK status; the
    ! frozen subdivision is then not a valid linearization point. The JVP must
    ! propagate that non-smooth/failure status and not return a product.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ad_test_utils, only: rel_err, check_smoothness, ad_status_t, &
        AD_SMOOTH, AD_NONSMOOTH
    use fortnum_integrate, only: integrate_integrand_t, integrate_workspace_t, &
        integrate_result_t, integrate_qag, integrate_qag_jvp
    implicit none

    ! Parameter for the smooth exponential case, threaded through ctx.
    type :: pbox_t
        real(dp) :: p = 0.0_dp
    end type pbox_t

    integer :: nfail
    nfail = 0

    call test_smooth_exp(nfail)
    call test_nonsmooth_spike(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! ---------------------------------------------------------- smooth case

    subroutine test_smooth_exp(nfail)
        integer, intent(inout) :: nfail
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st, jst
        type(pbox_t) :: box
        real(dp) :: p, di_jvp, di_exact, di_fd, e_exact, e_fd, h
        real(dp) :: ip, im
        integer  :: nsub_base

        p = 12.0_dp
        box%p = p

        ! Primal run freezes the accepted subdivision into res. A steep
        ! exponential plus the low-order key=15 rule forces several bisections,
        ! so the frozen trace has more than one panel and the JVP exercises the
        ! multi-panel trace-rule sum, not a single panel.
        call integrate_qag(g_exp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, ws, res, &
                           st, key=15, ctx=box)
        if (st%code /= FORTNUM_OK) then
            write (error_unit, '(a,i0)') &
                "FAIL [smooth_exp] primal status code=", st%code
            nfail = nfail + 1
            return
        end if
        nsub_base = res%nsub

        ! Trace JVP: quadrature of df/dp = x exp(p x) on the frozen subdivision.
        call integrate_qag_jvp(dg_exp_dp, res, di_jvp, jst, ctx=box)
        if (jst%code /= FORTNUM_OK) then
            write (error_unit, '(a,i0)') &
                "FAIL [smooth_exp] jvp status code=", jst%code
            nfail = nfail + 1
            return
        end if

        ! Analytic dI/dp = (exp(p)(p-1) + 1)/p^2.
        di_exact = (exp(p)*(p - 1.0_dp) + 1.0_dp)/(p*p)

        ! Central FD of the primal integral w.r.t. p on the same subdivision.
        ! The step is small enough that the perturbed runs reuse the same
        ! accepted panels (verified below via the trace count).
        h = 1.0e-6_dp
        ip = primal_integral(p + h, nfail, nsub_base)
        im = primal_integral(p - h, nfail, nsub_base)
        di_fd = (ip - im)/(2.0_dp*h)

        e_exact = rel_err(di_jvp, di_exact)
        e_fd    = rel_err(di_jvp, di_fd)

        write (*, '(a,es24.16)') "smooth_exp jvp dI/dp     = ", di_jvp
        write (*, '(a,es24.16)') "smooth_exp analytic dI/dp= ", di_exact
        write (*, '(a,es24.16)') "smooth_exp central-FD    = ", di_fd
        write (*, '(a,es12.4)')  "smooth_exp rel_err vs analytic = ", e_exact
        write (*, '(a,es12.4)')  "smooth_exp rel_err vs FD       = ", e_fd
        write (*, '(a,i0)')      "smooth_exp frozen nsub         = ", nsub_base

        if (e_exact > 1.0e-9_dp) then
            write (error_unit, '(a,es12.4)') &
                "FAIL [smooth_exp] jvp vs analytic rel_err=", e_exact
            nfail = nfail + 1
        end if
        if (e_fd > 1.0e-6_dp) then
            write (error_unit, '(a,es12.4)') &
                "FAIL [smooth_exp] jvp vs central-FD rel_err=", e_fd
            nfail = nfail + 1
        end if
    end subroutine test_smooth_exp

    ! Re-run the primal at a perturbed p and return the integral. Asserts the
    ! perturbation kept the SAME accepted subdivision (trace tag = nsub), so the
    ! FD difference probes the frozen-trace derivative, not a re-adaptation.
    function primal_integral(p, nfail, nsub_base) result(value)
        real(dp), intent(in)    :: p
        integer,  intent(inout) :: nfail
        integer,  intent(in)    :: nsub_base
        real(dp) :: value
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        type(pbox_t) :: box
        type(ad_status_t) :: smooth
        box%p = p
        call integrate_qag(g_exp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, ws, res, &
                           st, key=15, ctx=box)
        value = res%value
        ! Trace-tag check: same panel count => same accepted subdivision branch.
        smooth = check_smoothness(nsub_base, res%nsub, AD_SMOOTH)
        if (.not. smooth%ok) then
            write (error_unit, '(a,i0,a,i0)') &
                "FAIL [smooth_exp] FD perturbation changed subdivision: ", &
                nsub_base, " -> ", res%nsub
            nfail = nfail + 1
        end if
    end function primal_integral

    ! --------------------------------------------------------- non-smooth case

    subroutine test_nonsmooth_spike(nfail)
        integer, intent(inout) :: nfail
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st, jst
        type(pbox_t) :: box
        real(dp) :: di_jvp
        type(ad_status_t) :: smooth

        ! Interior inverse-square-root spike: the integrable singularity at p
        ! defeats plain QAG within the subdivision limit, so the primal records
        ! a non-OK status (the perturbation-moves-the-boundary case).
        box%p = 0.5_dp
        call integrate_qag(g_spike, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-12_dp, ws, &
                           res, st, limit=40, ctx=box)
        if (st%code == FORTNUM_OK) then
            write (error_unit, '(a)') &
                "FAIL [nonsmooth_spike] primal unexpectedly converged"
            nfail = nfail + 1
            return
        end if

        ! The JVP must report non-smoothness: it propagates the recorded
        ! non-OK status and returns no product.
        call integrate_qag_jvp(dg_spike_dp, res, di_jvp, jst, ctx=box)
        write (*, '(a,i0)') "nonsmooth_spike primal status code = ", st%code
        write (*, '(a,i0)') "nonsmooth_spike jvp    status code = ", jst%code

        ! Trace-tag verdict via the harness: tag the OK base point AD_SMOOTH and
        ! the observed derivative point by its status (non-OK => AD_NONSMOOTH),
        ! then assert the intended branch-change verdict.
        smooth = check_smoothness(AD_SMOOTH, &
                                  merge(AD_SMOOTH, AD_NONSMOOTH, &
                                        jst%code == FORTNUM_OK), &
                                  AD_NONSMOOTH)
        if (jst%code == FORTNUM_OK) then
            write (error_unit, '(a)') &
                "FAIL [nonsmooth_spike] jvp did not flag non-smoothness"
            nfail = nfail + 1
        end if
        if (.not. smooth%ok) then
            write (error_unit, '(a)') &
                "FAIL [nonsmooth_spike] smoothness verdict mismatch"
            nfail = nfail + 1
        end if
    end subroutine test_nonsmooth_spike

    ! ---------------------------------------------------------- integrands

    function g_exp(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        real(dp) :: p
        p = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (pbox_t)
                p = ctx%p
            end select
        end if
        fx = exp(p*x)
    end function g_exp

    ! df/dp of exp(p x) is x exp(p x).
    function dg_exp_dp(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        real(dp) :: p
        p = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (pbox_t)
                p = ctx%p
            end select
        end if
        fx = x*exp(p*x)
    end function dg_exp_dp

    function g_spike(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        real(dp) :: p, d
        p = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (pbox_t)
                p = ctx%p
            end select
        end if
        d = abs(x - p)
        if (d <= 0.0_dp) then
            fx = 0.0_dp
        else
            fx = 1.0_dp/sqrt(d)
        end if
    end function g_spike

    ! Formal tangent of the spike (sign of (x-p) times the p-derivative). Only
    ! reached if the primal converged, which it must not in this case.
    function dg_spike_dp(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        real(dp) :: p, d
        p = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
            type is (pbox_t)
                p = ctx%p
            end select
        end if
        d = abs(x - p)
        if (d <= 0.0_dp) then
            fx = 0.0_dp
        else
            fx = 0.5_dp*sign(1.0_dp, x - p)/(d*sqrt(d))
        end if
    end function dg_spike_dp

end program test_integrate_ad
