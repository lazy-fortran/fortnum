module fortnum_ode_cash_karp
    ! One Cash-Karp RK5(4) step: six stages, a fifth-order solution, and the
    ! embedded fourth-order solution that supplies the local error estimate.
    !
    ! Derivative policy: trace_rule (ad.md §1, §4).
    !   The stepper is a fixed linear map from stage values to the step update
    !   once the node/weight tableau is frozen. The adaptive schedule lives in
    !   the integrator (fortnum_ode); a trace_rule derivative differentiates
    !   this kernel at the step sizes the primal chose. Keeping the kernel
    !   separable and stateless is what makes that frozen-schedule rule cheap.
    !   Active arguments: y, the stage derivatives. Inactive: h (the step is a
    !   frozen control on the recorded mesh).
    !
    ! Tableau: Cash and Karp, "A Variable Order Runge-Kutta Method for Initial
    !   Value Problems with Rapidly Varying Right-Hand Sides", ACM TOMS 16
    !   (1990) 201-222. Nodes c, sub-diagonal matrix a, fifth-order weights b5,
    !   fourth-order weights b4. yerr = y5 - y4 is the local error estimate.
    !
    ! No module-level state. The caller owns every array; the step writes only
    ! its output arguments and the six stage slots it is handed.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    ! The tableau is generated from the verified Butcher data, not declared
    ! here. See tools/codegen/app/gen_rk54_cpu_tableau.f90.
    use fortnum_rk54_ck_tableau, only: &
        ck_c2, ck_c3, ck_c4, ck_c5, ck_c6, &
        ck_a21, ck_a31, ck_a32, ck_a41, ck_a42, ck_a43, &
        ck_a51, ck_a52, ck_a53, ck_a54, &
        ck_a61, ck_a62, ck_a63, ck_a64, ck_a65, &
        ck_b1, ck_b3, ck_b4, ck_b6, &
        ck_e1, ck_e3, ck_e4, ck_e5, ck_e6
    implicit none
    private

    public :: cash_karp_step

contains

    ! Advance one Cash-Karp step of size h from (t, y).
    !
    ! rhs is the user RHS (matching ode_rhs_t); ctx is its optional context.
    ! k1..k6 are caller-owned stage-derivative slots (length neq); ytmp is a
    ! caller-owned scratch state. y5 receives the fifth-order solution, yerr
    ! the local error estimate y5 - y4. nfev is incremented by six.
    !
    ! k1 may carry the first-stage derivative from a prior evaluation: when
    ! have_k1 is .true. the routine trusts k1 and skips its evaluation, saving
    ! one RHS call across a rejected-then-retried step. Cash-Karp is not FSAL,
    ! so the integrator only sets have_k1 when k1 already holds f(t, y).
    subroutine cash_karp_step(rhs, t, y, h, have_k1, &
            k1, k2, k3, k4, k5, k6, ytmp, y5, yerr, nfev, ctx)
        interface
            subroutine rhs(t, y, dydt, ctx)
                import :: dp
                real(dp), intent(in)  :: t
                real(dp), intent(in)  :: y(:)
                real(dp), intent(out) :: dydt(:)
                class(*), intent(in), optional :: ctx
            end subroutine rhs
        end interface
        real(dp), intent(in)    :: t
        real(dp), intent(in)    :: y(:)
        real(dp), intent(in)    :: h
        logical,  intent(in)    :: have_k1
        real(dp), intent(inout) :: k1(:)
        real(dp), intent(out)   :: k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), intent(out)   :: ytmp(:)
        real(dp), intent(out)   :: y5(:)
        real(dp), intent(out)   :: yerr(:)
        integer,  intent(inout) :: nfev
        class(*), intent(in), optional :: ctx

        if (.not. have_k1) then
            call rhs(t, y, k1, ctx)
            nfev = nfev + 1
        end if

        ytmp = y + h * (ck_a21 * k1)
        call rhs(t + ck_c2 * h, ytmp, k2, ctx)

        ytmp = y + h * (ck_a31 * k1 + ck_a32 * k2)
        call rhs(t + ck_c3 * h, ytmp, k3, ctx)

        ytmp = y + h * (ck_a41 * k1 + ck_a42 * k2 + ck_a43 * k3)
        call rhs(t + ck_c4 * h, ytmp, k4, ctx)

        ytmp = y + h * (ck_a51 * k1 + ck_a52 * k2 + ck_a53 * k3 + ck_a54 * k4)
        call rhs(t + ck_c5 * h, ytmp, k5, ctx)

        ytmp = y + h * (ck_a61 * k1 + ck_a62 * k2 + ck_a63 * k3 + ck_a64 * k4 + ck_a65 * k5)
        call rhs(t + ck_c6 * h, ytmp, k6, ctx)

        nfev = nfev + 5

        y5 = y + h * (ck_b1 * k1 + ck_b3 * k3 + ck_b4 * k4 + ck_b6 * k6)

        ! yerr = y5 - y4. The weights are b - bhat, differenced once in exact
        ! arithmetic by the generator, so the shared y term never appears and
        ! the subtraction is not restated here.
        yerr = h * (ck_e1 * k1 + ck_e3 * k3 + ck_e4 * k4 + ck_e5 * k5 &
            + ck_e6 * k6)
    end subroutine cash_karp_step

end module fortnum_ode_cash_karp
