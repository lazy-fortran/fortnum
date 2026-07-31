module fortnum_ode_extrapolation
    ! Gragg-Bulirsch-Stoer extrapolation: an arbitrary-order explicit integrator
    ! for y' = F(t, y).
    !
    ! Why this rather than a transcribed high-order tableau. The benchmark this
    ! feeds needs a *very* high order explicit method to sit opposite the
    ! symplectic schemes. The obvious candidates are Verner 9(8) and Feagin
    ! RK10/12/14, but those are several hundred coefficients each, and a single
    ! mistyped digit yields an integrator that still runs, still converges, and
    ! is quietly of the wrong order -- the exact failure mode this project keeps
    ! tripping over. Extrapolation reaches any even order from first principles
    ! with no table at all: every constant here is derived in code from the step
    ! sequence. It is also not a substitute chosen for convenience -- Bulirsch-
    ! Stoer is a standard high-accuracy integrator in celestial mechanics, which
    ! is the literature this benchmark is answering.
    !
    ! Method. Gragg's modified midpoint rule over n substeps of H,
    !
    !     z_0 = y,  z_1 = z_0 + h F(t, z_0),  h = H/n
    !     z_{m+1} = z_{m-1} + 2h F(t + m h, z_m)
    !     S(n)    = (z_n + z_{n-1} + h F(t + H, z_n)) / 2
    !
    ! has an asymptotic error expansion in even powers of h. Aitken-Neville
    ! extrapolation of S over a sequence n_1 < n_2 < ... therefore gains two
    ! orders per column: column j has order 2j.
    !
    ! Derivative policy: trace_rule (ad.md sections 1 and 4), as for the other
    ! steppers here.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
                              FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t, ode_problem_t, ode_solution_t

    implicit none
    private

    public :: ode_integrate_gbs, ode_solve_gbs
    public :: gbs_modified_midpoint, gbs_step_sequence

    ! Column count. Column j carries order 2j, so 8 columns is order 16 -- the
    ! same ballpark as Feagin RK14, which is what this stands in for.
    integer, parameter, public :: GBS_MAX_COLS = 8

    real(dp), parameter :: SAFETY = 0.8_dp
    real(dp), parameter :: FAC_MIN = 0.1_dp
    real(dp), parameter :: FAC_MAX = 5.0_dp

contains

    ! Deuflhard's sequence 2, 4, 6, 8, ... Even entries are what make the error
    ! expansion contain only even powers of h.
    pure subroutine gbs_step_sequence(nseq)
        integer, intent(out) :: nseq(:)
        integer :: k
        do k = 1, size(nseq)
            nseq(k) = 2*k
        end do
    end subroutine gbs_step_sequence

    ! Gragg's modified midpoint rule: n substeps across H, returning the
    ! smoothed endpoint. nfev is incremented by n + 1.
    subroutine gbs_modified_midpoint(rhs, t, y, big_h, n, zprev, zcur, ftmp, &
                                     yout, nfev, ctx)
        procedure(ode_rhs_t) :: rhs
        real(dp), intent(in) :: t, big_h
        real(dp), intent(in) :: y(:)
        integer,  intent(in) :: n
        real(dp), intent(inout) :: zprev(:), zcur(:), ftmp(:)
        real(dp), intent(out) :: yout(:)
        integer, intent(inout) :: nfev
        class(*), intent(in), optional :: ctx

        real(dp) :: h, swap
        integer  :: m, k, neq

        neq = size(y)
        h = big_h/real(n, dp)

        zprev(1:neq) = y(1:neq)
        call rhs(t, y, ftmp, ctx)
        nfev = nfev + 1
        zcur(1:neq) = zprev(1:neq) + h*ftmp(1:neq)

        do m = 1, n - 1
            call rhs(t + real(m, dp)*h, zcur(1:neq), ftmp, ctx)
            nfev = nfev + 1
            do k = 1, neq
                swap = zprev(k) + 2.0_dp*h*ftmp(k)
                zprev(k) = zcur(k)
                zcur(k) = swap
            end do
        end do

        call rhs(t + big_h, zcur(1:neq), ftmp, ctx)
        nfev = nfev + 1
        yout(1:neq) = 0.5_dp*(zcur(1:neq) + zprev(1:neq) + h*ftmp(1:neq))
    end subroutine gbs_modified_midpoint

    ! One extrapolation step. Adds columns until the tolerance is met or the
    ! table is exhausted; ncols reports how many were needed, which is the
    ! effective order (2*ncols) actually used.
    subroutine gbs_step(rhs, t, y, big_h, rtol, atol, nseq, table, zprev, zcur, &
                        ftmp, ymid, yout, err, ncols, nfev, ctx)
        procedure(ode_rhs_t) :: rhs
        real(dp), intent(in) :: t, big_h, rtol, atol
        real(dp), intent(in) :: y(:)
        integer,  intent(in) :: nseq(:)
        real(dp), intent(inout) :: table(:, :)   ! (GBS_MAX_COLS, neq)
        real(dp), intent(inout) :: zprev(:), zcur(:), ftmp(:), ymid(:)
        real(dp), intent(out) :: yout(:)
        real(dp), intent(out) :: err
        integer,  intent(out) :: ncols
        integer,  intent(inout) :: nfev
        class(*), intent(in), optional :: ctx

        integer  :: k, j, m, neq
        real(dp) :: ratio, denom, sc, e, aux, cnew

        neq = size(y)
        err = huge(1.0_dp)
        ncols = 0

        do k = 1, GBS_MAX_COLS
            ! ymid is a distinct array from ftmp: passing the same actual
            ! argument for the scratch derivative and the result would alias.
            call gbs_modified_midpoint(rhs, t, y, big_h, nseq(k), zprev, zcur, &
                                       ftmp, ymid, nfev, ctx)

            ! Aitken-Neville in h^2, one component at a time so the single
            ! carried T_old(j-1) can be a scalar:
            !     T_new(j) = T_new(j-1) + (T_new(j-1) - T_old(j-1)) / R,
            !     R = (n_k / n_{k-j+1})^2 - 1.
            ! Every constant is derived from the step sequence; nothing tabulated.
            do m = 1, neq
                aux = table(1, m)          ! T_old(1), before it is overwritten
                table(1, m) = ymid(m)      ! T_new(1)
                do j = 2, k
                    ratio = real(nseq(k), dp)/real(nseq(k - j + 1), dp)
                    denom = ratio*ratio - 1.0_dp
                    cnew = table(j - 1, m) + (table(j - 1, m) - aux)/denom
                    aux = table(j, m)      ! T_old(j), needed by the next column
                    table(j, m) = cnew
                end do
            end do

            ncols = k
            if (k >= 2) then
                ! Error estimate: difference between the two highest columns.
                err = 0.0_dp
                do m = 1, neq
                    sc = atol + rtol*max(abs(y(m)), abs(table(k, m)))
                    if (sc <= 0.0_dp) sc = 1.0_dp
                    e = abs(table(k, m) - table(k - 1, m))/sc
                    if (e > err) err = e
                end do
                if (err <= 1.0_dp) exit
            end if
        end do

        yout(1:neq) = table(ncols, 1:neq)
    end subroutine gbs_step

    subroutine ode_integrate_gbs(problem, solution, status)
        type(ode_problem_t), intent(in) :: problem
        type(ode_solution_t), intent(inout) :: solution
        type(fortnum_status_t), intent(out) :: status

        real(dp), allocatable :: table(:, :), zprev(:), zcur(:), ftmp(:), ymid(:)
        real(dp), allocatable :: yout(:), yv(:)
        real(dp), allocatable :: t_rec(:), y_rec(:, :), h_rec(:), e_rec(:)
        integer, allocatable :: nseq(:)
        real(dp) :: t, h, hdir, err, fac, span, expo
        integer :: neq, nsteps, nrej, nfev, ncols, cap, rec_cap

        call status_set(status, FORTNUM_OK, "")

        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_gbs: y0 not allocated")
            return
        end if
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_gbs: rhs not associated")
            return
        end if
        if (problem%rtol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_gbs: rtol must be positive")
            return
        end if

        neq = size(problem%y0)
        span = problem%t1 - problem%t0
        if (span == 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_gbs: t1 must differ from t0")
            return
        end if
        hdir = sign(1.0_dp, span)

        allocate (nseq(GBS_MAX_COLS))
        call gbs_step_sequence(nseq)
        allocate (table(GBS_MAX_COLS, neq))
        allocate (zprev(neq), zcur(neq), ftmp(neq), ymid(neq), yout(neq), yv(neq))
        table = 0.0_dp

        ! Trace capacity grows by doubling rather than being allocated at
        ! max_steps up front. Callers that drive this one macro-step at a time
        ! (SIMPLE's orbit loop does) would otherwise pay a max_steps-sized
        ! allocation per call, which would dominate the very wall-clock the
        ! integrator benchmark is trying to measure.
        cap = max(problem%max_steps, 1)
        rec_cap = min(cap + 1, 256)
        allocate (t_rec(rec_cap), y_rec(neq, rec_cap), h_rec(rec_cap), e_rec(rec_cap))

        t = problem%t0
        yv = problem%y0
        h = problem%h0
        if (h == 0.0_dp) h = 1.0e-2_dp*abs(span)
        h = hdir*min(abs(h), abs(span))
        if (problem%hmax > 0.0_dp) h = hdir*min(abs(h), problem%hmax)

        nsteps = 0
        nrej = 0
        nfev = 0
        t_rec(1) = t
        y_rec(:, 1) = yv
        h_rec(1) = 0.0_dp
        e_rec(1) = 0.0_dp

        do while ((problem%t1 - t)*hdir > 0.0_dp)
            if (nsteps >= cap) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                                "ode_integrate_gbs: step limit reached")
                exit
            end if
            if (abs(h) > abs(problem%t1 - t)) h = problem%t1 - t

            call gbs_step(problem%rhs, t, yv, h, problem%rtol, problem%atol, &
                          nseq, table, zprev, zcur, ftmp, ymid, yout, err, ncols, &
                          nfev)

            if (err <= 1.0_dp .or. abs(h) <= abs(problem%hmin)) then
                t = t + h
                yv = yout
                nsteps = nsteps + 1
                if (nsteps + 1 > rec_cap) then
                    call grow_trace(t_rec, y_rec, h_rec, e_rec, rec_cap, cap + 1)
                end if
                t_rec(nsteps + 1) = t
                y_rec(:, nsteps + 1) = yv
                h_rec(nsteps + 1) = h
                e_rec(nsteps + 1) = err
            else
                nrej = nrej + 1
            end if

            ! The order actually used sets the exponent: column ncols has order
            ! 2*ncols, so the local error is O(h^(2*ncols+1)).
            expo = 1.0_dp/real(2*max(ncols, 1) + 1, dp)
            if (err > 0.0_dp) then
                fac = SAFETY*(1.0_dp/err)**expo
            else
                fac = FAC_MAX
            end if
            fac = max(FAC_MIN, min(FAC_MAX, fac))
            h = h*fac
            if (problem%hmax > 0.0_dp) h = hdir*min(abs(h), problem%hmax)
            if (problem%hmin > 0.0_dp) h = hdir*max(abs(h), problem%hmin)
        end do

        solution%nsteps = nsteps
        solution%nrejected = nrej
        solution%nfev = nfev
        if (allocated(solution%t)) deallocate (solution%t)
        if (allocated(solution%y)) deallocate (solution%y)
        if (allocated(solution%h)) deallocate (solution%h)
        if (allocated(solution%err)) deallocate (solution%err)
        allocate (solution%t(nsteps + 1), solution%y(neq, nsteps + 1))
        allocate (solution%h(nsteps + 1), solution%err(nsteps + 1))
        solution%t = t_rec(1:nsteps + 1)
        solution%y = y_rec(:, 1:nsteps + 1)
        solution%h = h_rec(1:nsteps + 1)
        solution%err = e_rec(1:nsteps + 1)
        solution%status = status
    end subroutine ode_integrate_gbs

    subroutine ode_solve_gbs(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
        procedure(ode_rhs_t) :: rhs
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:)
        real(dp), allocatable, intent(out) :: t_out(:)
        real(dp), allocatable, intent(out) :: y_out(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: rtol, atol

        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution

        problem%rhs => rhs
        problem%t0 = t0
        problem%t1 = t1
        problem%y0 = y0
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol

        call ode_integrate_gbs(problem, solution, status)

        if (allocated(solution%t)) then
            t_out = solution%t
            y_out = solution%y
        else
            allocate (t_out(1), y_out(size(y0), 1))
            t_out(1) = t0
            y_out(:, 1) = y0
        end if
    end subroutine ode_solve_gbs


    ! Double the recorded-trace capacity, never exceeding hard_cap.
    subroutine grow_trace(t_rec, y_rec, h_rec, e_rec, rec_cap, hard_cap)
        real(dp), allocatable, intent(inout) :: t_rec(:), y_rec(:, :)
        real(dp), allocatable, intent(inout) :: h_rec(:), e_rec(:)
        integer, intent(inout) :: rec_cap
        integer, intent(in) :: hard_cap

        real(dp), allocatable :: tt(:), yy(:, :), hh(:), ee(:)
        integer :: newcap, neq

        newcap = min(max(2*rec_cap, rec_cap + 1), hard_cap)
        if (newcap <= rec_cap) return
        neq = size(y_rec, 1)

        allocate (tt(newcap), yy(neq, newcap), hh(newcap), ee(newcap))
        tt(1:rec_cap) = t_rec(1:rec_cap)
        yy(:, 1:rec_cap) = y_rec(:, 1:rec_cap)
        hh(1:rec_cap) = h_rec(1:rec_cap)
        ee(1:rec_cap) = e_rec(1:rec_cap)
        call move_alloc(tt, t_rec)
        call move_alloc(yy, y_rec)
        call move_alloc(hh, h_rec)
        call move_alloc(ee, e_rec)
        rec_cap = newcap
    end subroutine grow_trace

end module fortnum_ode_extrapolation
