module fortnum_ode_dop853
    ! Adaptive Prince-Dormand RK8(7)-13M integrator (dop853 class). Twelve
    ! explicit stages give an eighth-order solution; two embedded error
    ! estimators (orders 5 and 3) drive the adaptive step. Replaces GSL
    ! gsl_odeiv_step_rk8pd for the very tight tolerances KiLCA needs.
    !
    ! Derivative policy: trace_rule (ad.md sec 1, 4; ode.md sec 1, 4).
    !   Same policy as the Cash-Karp integrator in fortnum_ode: the adaptive
    !   schedule is data-dependent, and a sensitivity differentiates the frozen
    !   accepted-step mesh (solution%t, solution%y, solution%h) with that
    !   schedule held fixed. The high-order forward/reverse products ride the
    !   recorded trace exactly as fortnum_ode's do; this module shares the
    !   ode_problem_t / ode_workspace_t / ode_solution_t carriers and the
    !   recorded mesh those products walk. Active: y0, ctx parameters. Inactive:
    !   rtol, atol, h0, hmin, hmax, max_steps.
    !
    ! Method and coefficients: Prince and Dormand, "High order embedded
    !   Runge-Kutta formulae", J. Comput. Appl. Math. 7 (1981) 67-75; Hairer,
    !   Norsett, Wanner, "Solving Ordinary Differential Equations I", 2nd ed.,
    !   II.5 (DOP853). The 13M node/coupling/weight constants and the order-5
    !   and order-3 error weights are the published RK8(7)13M values. PI step
    !   control and the starting-step estimate follow Hairer I, II.4.
    !
    ! No module-level state. The caller owns the workspace and the recorded
    ! trace; the step writes only its output arguments and stage slots.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_rhs_t
    implicit none
    private

    public :: dop853_step, ode_integrate_dop, ode_solve_dop

    ! PI controller for an order-8 method (Hairer I, II.4). The error estimate
    ! the step accepts is order 8 (err = err5^2 / sqrt(err5^2 + 0.01 err3^2)
    ! has the order-8 leading term), so the controller exponents use order 8.
    real(dp), parameter :: PI_ALPHA = 1.0_dp / 8.0_dp - 0.2_dp * 0.75_dp / 8.0_dp
    real(dp), parameter :: PI_BETA  = 0.2_dp * 0.75_dp / 8.0_dp
    real(dp), parameter :: SAFETY   = 0.9_dp
    real(dp), parameter :: FAC_MIN  = 0.333_dp
    real(dp), parameter :: FAC_MAX  = 6.0_dp
    real(dp), parameter :: TRACE_GROWTH = 2.0_dp

    ! Stage nodes c2..c12 (c1 = 0 is implicit).
    real(dp), parameter :: C2  = 0.05260015195876773_dp
    real(dp), parameter :: C3  = 0.0789002279381516_dp
    real(dp), parameter :: C4  = 0.1183503419072274_dp
    real(dp), parameter :: C5  = 0.2816496580927726_dp
    real(dp), parameter :: C6  = 0.3333333333333333_dp
    real(dp), parameter :: C7  = 0.25_dp
    real(dp), parameter :: C8  = 0.3076923076923077_dp
    real(dp), parameter :: C9  = 0.6512820512820513_dp
    real(dp), parameter :: C10 = 0.6_dp
    real(dp), parameter :: C11 = 0.8571428571428571_dp
    real(dp), parameter :: C12 = 1.0_dp

    ! Coupling coefficients a(i,j), j < i.
    real(dp), parameter :: A2_1 = 0.05260015195876773_dp
    real(dp), parameter :: A3_1 = 0.0197250569845379_dp
    real(dp), parameter :: A3_2 = 0.0591751709536137_dp
    real(dp), parameter :: A4_1 = 0.02958758547680685_dp
    real(dp), parameter :: A4_3 = 0.08876275643042054_dp
    real(dp), parameter :: A5_1 = 0.2413651341592667_dp
    real(dp), parameter :: A5_3 = -0.8845494793282861_dp
    real(dp), parameter :: A5_4 = 0.924834003261792_dp
    real(dp), parameter :: A6_1 = 0.037037037037037035_dp
    real(dp), parameter :: A6_4 = 0.17082860872947386_dp
    real(dp), parameter :: A6_5 = 0.12546768756682242_dp
    real(dp), parameter :: A7_1 = 0.037109375_dp
    real(dp), parameter :: A7_4 = 0.17025221101954405_dp
    real(dp), parameter :: A7_5 = 0.06021653898045596_dp
    real(dp), parameter :: A7_6 = -0.017578125_dp
    real(dp), parameter :: A8_1 = 0.03709200011850479_dp
    real(dp), parameter :: A8_4 = 0.17038392571223998_dp
    real(dp), parameter :: A8_5 = 0.10726203044637328_dp
    real(dp), parameter :: A8_6 = -0.015319437748624402_dp
    real(dp), parameter :: A8_7 = 0.008273789163814023_dp
    real(dp), parameter :: A9_1 = 0.6241109587160757_dp
    real(dp), parameter :: A9_4 = -3.3608926294469414_dp
    real(dp), parameter :: A9_5 = -0.868219346841726_dp
    real(dp), parameter :: A9_6 = 27.59209969944671_dp
    real(dp), parameter :: A9_7 = 20.154067550477894_dp
    real(dp), parameter :: A9_8 = -43.48988418106996_dp
    real(dp), parameter :: A10_1 = 0.47766253643826434_dp
    real(dp), parameter :: A10_4 = -2.4881146199716677_dp
    real(dp), parameter :: A10_5 = -0.590290826836843_dp
    real(dp), parameter :: A10_6 = 21.230051448181193_dp
    real(dp), parameter :: A10_7 = 15.279233632882423_dp
    real(dp), parameter :: A10_8 = -33.28821096898486_dp
    real(dp), parameter :: A10_9 = -0.020331201708508627_dp
    real(dp), parameter :: A11_1 = -0.9371424300859873_dp
    real(dp), parameter :: A11_4 = 5.186372428844064_dp
    real(dp), parameter :: A11_5 = 1.0914373489967295_dp
    real(dp), parameter :: A11_6 = -8.149787010746927_dp
    real(dp), parameter :: A11_7 = -18.52006565999696_dp
    real(dp), parameter :: A11_8 = 22.739487099350505_dp
    real(dp), parameter :: A11_9 = 2.4936055526796523_dp
    real(dp), parameter :: A11_10 = -3.0467644718982196_dp
    real(dp), parameter :: A12_1 = 2.273310147516538_dp
    real(dp), parameter :: A12_4 = -10.53449546673725_dp
    real(dp), parameter :: A12_5 = -2.0008720582248625_dp
    real(dp), parameter :: A12_6 = -17.9589318631188_dp
    real(dp), parameter :: A12_7 = 27.94888452941996_dp
    real(dp), parameter :: A12_8 = -2.8589982771350235_dp
    real(dp), parameter :: A12_9 = -8.87285693353063_dp
    real(dp), parameter :: A12_10 = 12.360567175794303_dp
    real(dp), parameter :: A12_11 = 0.6433927460157636_dp

    ! Eighth-order solution weights (stages 2..5 carry zero weight).
    real(dp), parameter :: B1  = 0.054293734116568765_dp
    real(dp), parameter :: B6  = 4.450312892752409_dp
    real(dp), parameter :: B7  = 1.8915178993145003_dp
    real(dp), parameter :: B8  = -5.801203960010585_dp
    real(dp), parameter :: B9  = 0.3111643669578199_dp
    real(dp), parameter :: B10 = -0.1521609496625161_dp
    real(dp), parameter :: B11 = 0.20136540080403034_dp
    real(dp), parameter :: B12 = 0.04471061572777259_dp

    ! Order-5 embedded error weights (err5 = sum E5_i k_i).
    real(dp), parameter :: E5_1  = 0.01312004499419488_dp
    real(dp), parameter :: E5_6  = -1.2251564463762044_dp
    real(dp), parameter :: E5_7  = -0.4957589496572502_dp
    real(dp), parameter :: E5_8  = 1.6643771824549864_dp
    real(dp), parameter :: E5_9  = -0.35032884874997366_dp
    real(dp), parameter :: E5_10 = 0.3341791187130175_dp
    real(dp), parameter :: E5_11 = 0.08192320648511571_dp
    real(dp), parameter :: E5_12 = -0.022355307863886294_dp

    ! Order-3 embedded error weights (err3 = sum E3_i k_i).
    real(dp), parameter :: E3_1  = -0.18980075407240762_dp
    real(dp), parameter :: E3_6  = 4.450312892752409_dp
    real(dp), parameter :: E3_7  = 1.8915178993145003_dp
    real(dp), parameter :: E3_8  = -5.801203960010585_dp
    real(dp), parameter :: E3_9  = -0.4226823213237919_dp
    real(dp), parameter :: E3_10 = -0.1521609496625161_dp
    real(dp), parameter :: E3_11 = 0.20136540080403034_dp
    real(dp), parameter :: E3_12 = 0.02265179219836082_dp

contains

    ! Advance one RK8(7)13M step of size h from (t, y). rhs matches ode_rhs_t.
    ! k1..k12 are caller-owned stage-derivative slots (length neq); ytmp is a
    ! caller-owned scratch state. y8 receives the eighth-order solution, err5
    ! and err3 the two embedded error vectors (already scaled by h). nfev is
    ! incremented by the stages evaluated.
    !
    ! have_k1 reuses an externally supplied first-stage derivative (the method
    ! is not FSAL, so the integrator only sets it when k1 already holds
    ! f(t, y)).
    subroutine dop853_step(rhs, t, y, h, have_k1, k1, k2, k3, k4, k5, k6, &
                           k7, k8, k9, k10, k11, k12, ytmp, y8, err5, err3, &
                           nfev, ctx)
        procedure(ode_rhs_t)            :: rhs
        real(dp), intent(in)            :: t
        real(dp), intent(in)            :: y(:)
        real(dp), intent(in)            :: h
        logical,  intent(in)            :: have_k1
        real(dp), intent(inout)         :: k1(:)
        real(dp), intent(out)           :: k2(:), k3(:), k4(:), k5(:), k6(:)
        real(dp), intent(out)           :: k7(:), k8(:), k9(:), k10(:)
        real(dp), intent(out)           :: k11(:), k12(:)
        real(dp), intent(out)           :: ytmp(:)
        real(dp), intent(out)           :: y8(:)
        real(dp), intent(out)           :: err5(:), err3(:)
        integer,  intent(inout)         :: nfev
        class(*), intent(in), optional  :: ctx

        if (.not. have_k1) then
            call rhs(t, y, k1, ctx)
            nfev = nfev + 1
        end if

        ytmp = y + h * (A2_1 * k1)
        call rhs(t + C2 * h, ytmp, k2, ctx)

        ytmp = y + h * (A3_1 * k1 + A3_2 * k2)
        call rhs(t + C3 * h, ytmp, k3, ctx)

        ytmp = y + h * (A4_1 * k1 + A4_3 * k3)
        call rhs(t + C4 * h, ytmp, k4, ctx)

        ytmp = y + h * (A5_1 * k1 + A5_3 * k3 + A5_4 * k4)
        call rhs(t + C5 * h, ytmp, k5, ctx)

        ytmp = y + h * (A6_1 * k1 + A6_4 * k4 + A6_5 * k5)
        call rhs(t + C6 * h, ytmp, k6, ctx)

        ytmp = y + h * (A7_1 * k1 + A7_4 * k4 + A7_5 * k5 + A7_6 * k6)
        call rhs(t + C7 * h, ytmp, k7, ctx)

        ytmp = y + h * (A8_1 * k1 + A8_4 * k4 + A8_5 * k5 + A8_6 * k6 &
                        + A8_7 * k7)
        call rhs(t + C8 * h, ytmp, k8, ctx)

        ytmp = y + h * (A9_1 * k1 + A9_4 * k4 + A9_5 * k5 + A9_6 * k6 &
                        + A9_7 * k7 + A9_8 * k8)
        call rhs(t + C9 * h, ytmp, k9, ctx)

        ytmp = y + h * (A10_1 * k1 + A10_4 * k4 + A10_5 * k5 + A10_6 * k6 &
                        + A10_7 * k7 + A10_8 * k8 + A10_9 * k9)
        call rhs(t + C10 * h, ytmp, k10, ctx)

        ytmp = y + h * (A11_1 * k1 + A11_4 * k4 + A11_5 * k5 + A11_6 * k6 &
                        + A11_7 * k7 + A11_8 * k8 + A11_9 * k9 + A11_10 * k10)
        call rhs(t + C11 * h, ytmp, k11, ctx)

        ytmp = y + h * (A12_1 * k1 + A12_4 * k4 + A12_5 * k5 + A12_6 * k6 &
                        + A12_7 * k7 + A12_8 * k8 + A12_9 * k9 + A12_10 * k10 &
                        + A12_11 * k11)
        call rhs(t + C12 * h, ytmp, k12, ctx)

        nfev = nfev + 11

        y8 = y + h * (B1 * k1 + B6 * k6 + B7 * k7 + B8 * k8 + B9 * k9 &
                      + B10 * k10 + B11 * k11 + B12 * k12)

        err5 = h * (E5_1 * k1 + E5_6 * k6 + E5_7 * k7 + E5_8 * k8 + E5_9 * k9 &
                    + E5_10 * k10 + E5_11 * k11 + E5_12 * k12)
        err3 = h * (E3_1 * k1 + E3_6 * k6 + E3_7 * k7 + E3_8 * k8 + E3_9 * k9 &
                    + E3_10 * k10 + E3_11 * k11 + E3_12 * k12)
    end subroutine dop853_step

    ! Integrate problem%rhs from t0 to t1 with adaptive RK8(7)13M. Records the
    ! accepted-step mesh into solution, matching ode_integrate's layout so the
    ! fortnum_ode trace_rule products can walk it. Events are not exposed.
    subroutine ode_integrate_dop(problem, workspace, solution, status)
        type(ode_problem_t),    intent(in)    :: problem
        type(ode_workspace_t),  intent(inout) :: workspace
        type(ode_solution_t),   intent(inout) :: solution
        type(fortnum_status_t), intent(out)   :: status

        integer  :: neq, nstep, cap
        real(dp) :: t, dir, span, h, hmax, hmin
        real(dp) :: err_norm, err_prev, fac, h_new
        logical  :: accepted, first_step, after_reject, final_step

        call status_set(status, FORTNUM_OK, "")
        solution%nsteps = 0
        solution%nrejected = 0
        solution%nfev = 0
        solution%event_found = .false.
        solution%t_event = 0.0_dp

        if (.not. validate_problem(problem, status)) then
            solution%status = status
            return
        end if

        neq = size(problem%y0)
        call ensure_workspace(workspace, neq)

        dir  = sign(1.0_dp, problem%t1 - problem%t0)
        span = abs(problem%t1 - problem%t0)
        if (problem%hmax > 0.0_dp) then
            hmax = min(problem%hmax, span)
        else
            hmax = span
        end if
        hmin = problem%hmin

        cap = 64
        call alloc_trace(solution, neq, cap)
        solution%t(1)   = problem%t0
        solution%y(:,1) = problem%y0
        nstep = 0

        if (span <= 0.0_dp) then
            call trim_trace(solution, 0, neq)
            solution%status = status
            return
        end if

        h = initial_step(problem, workspace, hmax)
        err_prev = 1.0_dp
        first_step = .true.
        after_reject = .false.
        t = problem%t0

        do
            if (nstep >= problem%max_steps) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ode_integrate_dop: exceeded max_steps")
                exit
            end if

            final_step = abs(problem%t1 - t) <= h
            if (final_step) h = abs(problem%t1 - t)
            if (h <= 0.0_dp) exit

            call dop853_step(problem%rhs, t, solution%y(:,nstep+1), dir * h, &
                .false., workspace%k1, workspace%k2, workspace%k3, &
                workspace%k4, workspace%k5, workspace%k6, workspace%k7, &
                workspace%k8, workspace%k9, workspace%k10, workspace%k11, &
                workspace%k12, workspace%ytmp, workspace%y8, &
                workspace%err5, workspace%err3, solution%nfev)

            err_norm = error_norm(solution%y(:,nstep+1), workspace%y8, &
                                  workspace%err5, workspace%err3, &
                                  problem%rtol, problem%atol)
            accepted = err_norm <= 1.0_dp

            if (accepted) then
                if (final_step) then
                    t = problem%t1
                else
                    t = t + dir * h
                end if
                nstep = nstep + 1
                if (nstep + 1 > cap) then
                    cap = max(cap + 1, nint(cap * TRACE_GROWTH))
                    call grow_trace(solution, neq, cap)
                end if
                solution%t(nstep+1) = t
                solution%y(:,nstep+1) = workspace%y8
                solution%h(nstep) = dir * h
                solution%err(nstep) = err_norm
                if (final_step) exit
            else
                solution%nrejected = solution%nrejected + 1
            end if

            fac = control_factor(err_norm, err_prev, first_step .or. after_reject)
            h_new = h * fac
            if (h_new > hmax) h_new = hmax

            if (accepted) then
                err_prev = max(err_norm, 1.0e-10_dp)
                first_step = .false.
                after_reject = .false.
            else
                after_reject = .true.
            end if

            if (hmin > 0.0_dp .and. h_new < hmin .and. .not. accepted) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ode_integrate_dop: step forced below hmin")
                exit
            end if
            h = h_new
        end do

        solution%nsteps = nstep
        call trim_trace(solution, nstep, neq)
        solution%status = status
    end subroutine ode_integrate_dop

    ! Flat call: build a problem, integrate with RK8(7)13M, hand back the trace.
    subroutine ode_solve_dop(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
        procedure(ode_rhs_t)               :: rhs
        real(dp),               intent(in) :: t0, t1
        real(dp),               intent(in) :: y0(:)
        real(dp), allocatable, intent(out) :: t_out(:)
        real(dp), allocatable, intent(out) :: y_out(:,:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional     :: rtol, atol

        type(ode_problem_t)   :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t)  :: solution
        integer :: npts

        problem%rhs => rhs
        problem%t0 = t0
        problem%t1 = t1
        problem%y0 = y0
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol

        call ode_integrate_dop(problem, workspace, solution, status)

        npts = solution%nsteps + 1
        if (allocated(solution%t)) then
            t_out = solution%t(1:npts)
            y_out = solution%y(:, 1:npts)
        else
            allocate(t_out(0))
            allocate(y_out(size(y0), 0))
        end if
    end subroutine ode_solve_dop

    ! --- internals ---

    logical function validate_problem(problem, status) result(ok)
        type(ode_problem_t),    intent(in)  :: problem
        type(fortnum_status_t), intent(out) :: status
        ok = .false.
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: rhs not associated")
            return
        end if
        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: y0 not allocated")
            return
        end if
        if (size(problem%y0) < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: neq < 1")
            return
        end if
        if (problem%rtol <= 0.0_dp .or. problem%atol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: tolerances must be positive")
            return
        end if
        if (problem%max_steps < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ode_integrate_dop: max_steps < 1")
            return
        end if
        ok = .true.
    end function validate_problem

    ! Allocate the twelve stage slots and temporaries when the size changes.
    ! Reuses ode_workspace_t: k1..k6 + ytmp are shared with Cash-Karp; the
    ! higher stages live in the workspace's spare slots.
    subroutine ensure_workspace(workspace, neq)
        type(ode_workspace_t), intent(inout) :: workspace
        integer,               intent(in)    :: neq
        if (workspace%neq == neq .and. allocated(workspace%k7)) return
        workspace%neq = neq
        call realloc(workspace%k1, neq)
        call realloc(workspace%k2, neq)
        call realloc(workspace%k3, neq)
        call realloc(workspace%k4, neq)
        call realloc(workspace%k5, neq)
        call realloc(workspace%k6, neq)
        call realloc(workspace%k7, neq)
        call realloc(workspace%k8, neq)
        call realloc(workspace%k9, neq)
        call realloc(workspace%k10, neq)
        call realloc(workspace%k11, neq)
        call realloc(workspace%k12, neq)
        call realloc(workspace%ytmp, neq)
        call realloc(workspace%y8, neq)
        call realloc(workspace%err5, neq)
        call realloc(workspace%err3, neq)
    end subroutine ensure_workspace

    subroutine realloc(a, n)
        real(dp), allocatable, intent(inout) :: a(:)
        integer,               intent(in)    :: n
        if (allocated(a)) then
            if (size(a) == n) return
            deallocate(a)
        end if
        allocate(a(n))
    end subroutine realloc

    ! Hairer-Norsett-Wanner DOP853 error norm: combine the order-5 and order-3
    ! estimates so the accepted error has the order-8 leading term. Both errN
    ! arrays already carry the h factor (see dop853_step).
    real(dp) function error_norm(y, y8, err5, err3, rtol, atol) result(en)
        real(dp), intent(in) :: y(:), y8(:), err5(:), err3(:)
        real(dp), intent(in) :: rtol, atol
        real(dp) :: sc, e5sq, e3sq, denom
        integer  :: i, n
        n = size(y)
        e5sq = 0.0_dp
        e3sq = 0.0_dp
        do i = 1, n
            sc = atol + rtol * max(abs(y(i)), abs(y8(i)))
            e5sq = e5sq + (err5(i) / sc)**2
            e3sq = e3sq + (err3(i) / sc)**2
        end do
        if (e5sq <= 0.0_dp .and. e3sq <= 0.0_dp) then
            en = 0.0_dp
            return
        end if
        denom = e5sq + 0.01_dp * e3sq
        en = sqrt(e5sq / (denom * real(n, dp))) * sqrt(e5sq)
    end function error_norm

    ! Step-growth factor. PI form on accepted history; pure I form when
    ! restarting (Hairer I, II.4). Clamped to [FAC_MIN, FAC_MAX].
    real(dp) function control_factor(err_norm, err_prev, restart) result(fac)
        real(dp), intent(in) :: err_norm, err_prev
        logical,  intent(in) :: restart
        real(dp) :: e
        e = max(err_norm, 1.0e-10_dp)
        if (restart) then
            fac = SAFETY * e**(-1.0_dp / 8.0_dp)
        else
            fac = SAFETY * e**(-PI_ALPHA) * err_prev**(PI_BETA)
        end if
        fac = max(FAC_MIN, min(FAC_MAX, fac))
    end function control_factor

    ! Starting step from the scaled RHS magnitude at t0 (Hairer I, II.4),
    ! tuned for an order-8 method via the 1/8 power. Returns a positive
    ! magnitude bounded by hmax.
    real(dp) function initial_step(problem, workspace, hmax) result(h)
        type(ode_problem_t),   intent(in)    :: problem
        type(ode_workspace_t), intent(inout) :: workspace
        real(dp),              intent(in)    :: hmax
        real(dp) :: d0, d1, sc
        integer  :: i, n
        if (problem%h0 > 0.0_dp) then
            h = min(problem%h0, hmax)
            return
        end if
        n = size(problem%y0)
        call problem%rhs(problem%t0, problem%y0, workspace%k1)
        d0 = 0.0_dp
        d1 = 0.0_dp
        do i = 1, n
            sc = problem%atol + problem%rtol * abs(problem%y0(i))
            d0 = d0 + (problem%y0(i) / sc)**2
            d1 = d1 + (workspace%k1(i) / sc)**2
        end do
        d0 = sqrt(d0 / real(n, dp))
        d1 = sqrt(d1 / real(n, dp))
        if (d1 <= 1.0e-10_dp) then
            h = 1.0e-6_dp
        else
            h = 0.01_dp * d0 / d1
        end if
        h = max(h, 1.0e-10_dp)
        h = min(h, hmax)
        if (h <= 0.0_dp) h = min(1.0e-4_dp, hmax)
    end function initial_step

    subroutine alloc_trace(solution, neq, cap)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: neq, cap
        if (allocated(solution%t)) deallocate(solution%t)
        if (allocated(solution%y)) deallocate(solution%y)
        if (allocated(solution%h)) deallocate(solution%h)
        if (allocated(solution%err)) deallocate(solution%err)
        allocate(solution%t(cap))
        allocate(solution%y(neq, cap))
        allocate(solution%h(cap))
        allocate(solution%err(cap))
    end subroutine alloc_trace

    subroutine grow_trace(solution, neq, cap)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: neq, cap
        real(dp), allocatable :: tt(:), yy(:,:), hh(:), ee(:)
        integer :: old
        old = size(solution%t)
        allocate(tt(cap))
        allocate(yy(neq, cap))
        allocate(hh(cap))
        allocate(ee(cap))
        tt(1:old) = solution%t
        yy(:,1:old) = solution%y
        hh(1:old) = solution%h
        ee(1:old) = solution%err
        call move_alloc(tt, solution%t)
        call move_alloc(yy, solution%y)
        call move_alloc(hh, solution%h)
        call move_alloc(ee, solution%err)
    end subroutine grow_trace

    subroutine trim_trace(solution, nstep, neq)
        type(ode_solution_t), intent(inout) :: solution
        integer,              intent(in)    :: nstep, neq
        real(dp), allocatable :: tt(:), yy(:,:), hh(:), ee(:)
        integer :: npts
        npts = nstep + 1
        allocate(tt(npts))
        allocate(yy(neq, npts))
        tt = solution%t(1:npts)
        yy = solution%y(:, 1:npts)
        call move_alloc(tt, solution%t)
        call move_alloc(yy, solution%y)
        if (nstep > 0) then
            allocate(hh(nstep))
            allocate(ee(nstep))
            hh = solution%h(1:nstep)
            ee = solution%err(1:nstep)
            call move_alloc(hh, solution%h)
            call move_alloc(ee, solution%err)
        else
            if (allocated(solution%h)) deallocate(solution%h)
            if (allocated(solution%err)) deallocate(solution%err)
            allocate(solution%h(0))
            allocate(solution%err(0))
        end if
    end subroutine trim_trace

end module fortnum_ode_dop853
