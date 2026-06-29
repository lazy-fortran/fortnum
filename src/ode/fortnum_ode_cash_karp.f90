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
    implicit none
    private

    public :: cash_karp_step

    ! Stage nodes c2..c6 (c1 = 0 is implicit).
    real(dp), parameter :: c2 = 0.2_dp
    real(dp), parameter :: c3 = 0.3_dp
    real(dp), parameter :: c4 = 0.6_dp
    real(dp), parameter :: c5 = 1.0_dp
    real(dp), parameter :: c6 = 0.875_dp

    ! Coupling coefficients a(i,j), j < i.
    real(dp), parameter :: a21 = 0.2_dp
    real(dp), parameter :: a31 = 3.0_dp / 40.0_dp
    real(dp), parameter :: a32 = 9.0_dp / 40.0_dp
    real(dp), parameter :: a41 = 0.3_dp
    real(dp), parameter :: a42 = -0.9_dp
    real(dp), parameter :: a43 = 1.2_dp
    real(dp), parameter :: a51 = -11.0_dp / 54.0_dp
    real(dp), parameter :: a52 = 2.5_dp
    real(dp), parameter :: a53 = -70.0_dp / 27.0_dp
    real(dp), parameter :: a54 = 35.0_dp / 27.0_dp
    real(dp), parameter :: a61 = 1631.0_dp / 55296.0_dp
    real(dp), parameter :: a62 = 175.0_dp / 512.0_dp
    real(dp), parameter :: a63 = 575.0_dp / 13824.0_dp
    real(dp), parameter :: a64 = 44275.0_dp / 110592.0_dp
    real(dp), parameter :: a65 = 253.0_dp / 4096.0_dp

    ! Fifth-order solution weights (stage 2 has weight zero).
    real(dp), parameter :: b5_1 = 37.0_dp / 378.0_dp
    real(dp), parameter :: b5_3 = 250.0_dp / 621.0_dp
    real(dp), parameter :: b5_4 = 125.0_dp / 594.0_dp
    real(dp), parameter :: b5_6 = 512.0_dp / 1771.0_dp

    ! Embedded fourth-order solution weights.
    real(dp), parameter :: b4_1 = 2825.0_dp / 27648.0_dp
    real(dp), parameter :: b4_3 = 18575.0_dp / 48384.0_dp
    real(dp), parameter :: b4_4 = 13525.0_dp / 55296.0_dp
    real(dp), parameter :: b4_5 = 277.0_dp / 14336.0_dp
    real(dp), parameter :: b4_6 = 0.25_dp


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

        ytmp = y + h * (a21 * k1)
        call rhs(t + c2 * h, ytmp, k2, ctx)

        ytmp = y + h * (a31 * k1 + a32 * k2)
        call rhs(t + c3 * h, ytmp, k3, ctx)

        ytmp = y + h * (a41 * k1 + a42 * k2 + a43 * k3)
        call rhs(t + c4 * h, ytmp, k4, ctx)

        ytmp = y + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
        call rhs(t + c5 * h, ytmp, k5, ctx)

        ytmp = y + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
        call rhs(t + c6 * h, ytmp, k6, ctx)

        nfev = nfev + 5

        y5 = y + h * (b5_1 * k1 + b5_3 * k3 + b5_4 * k4 + b5_6 * k6)

        ! yerr = y5 - y4, formed from the weight differences to avoid cancelling
        ! the shared y term.
        yerr = h * ((b5_1 - b4_1) * k1 + (b5_3 - b4_3) * k3 &
            + (b5_4 - b4_4) * k4 - b4_5 * k5 + (b5_6 - b4_6) * k6)
    end subroutine cash_karp_step

end module fortnum_ode_cash_karp
