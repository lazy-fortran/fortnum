module fortnum_integrate
    ! Globally adaptive Gauss-Kronrod integration (QUADPACK QAG/QAGS pattern).
    !
    ! Driver per ADR docs/design/integrate.md. The per-panel GK estimate comes
    ! from fortnum_integrate_gk%gk_apply; this module owns only the global
    ! bisection bookkeeping, the Wynn epsilon extrapolation, the caller-owned
    ! work stack and epsilon table, and the recorded subdivision trace. The GK
    ! rule itself is not re-derived here.
    !
    ! Derivative policy: trace_rule (ad.md sections 1, 4). The subdivision is
    ! data-dependent; once the primal freezes it into integrate_result_t the
    ! integral is a fixed weighted sum of GK panel values, so a #40 derivative
    ! product differentiates at that frozen subdivision. The trace columns
    ! (sub_a/sub_b/sub_r/sub_e plus key and extrapolated) are the hooks #40
    ! re-walks; nothing here moves a boundary or re-adapts.
    !
    ! No module-level state and no global procedure pointer: the caller's
    ! integrand and its optional ctx are threaded to the ctx-free kernel through
    ! a host-associated wrapper that lives on the call stack for the duration of
    ! one driver call (ADR section 2). This lets a later derivative product call
    ! the integrand without a hidden side channel.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_integrate_gk, only: gk_apply
    use fortnum_status, only: fortnum_status_t, status_set, &
                              FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
                              FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: integrate_integrand_t
    public :: integrate_workspace_t, integrate_epstab_t, integrate_result_t
    public :: integrate_qag, integrate_qags, integrate
    public :: integrate_qagp, integrate_qagiu
    public :: integrate_qag_jvp

    abstract interface
        function integrate_integrand_t(x, ctx) result(fx)
            import :: dp
            real(dp), intent(in) :: x
            class(*), intent(in), optional :: ctx
            real(dp) :: fx
        end function integrate_integrand_t
    end interface

    ! Bridge interface: the kernel's gk_apply takes a ctx-free integrand. Each
    ! entry point passes a host-associated wrapper around (f, ctx) that matches
    ! this; the driver's procedure dummy is typed against it. No module pointer.
    abstract interface
        function panel_kernel_t(x) result(fx)
            import :: dp
            real(dp), intent(in) :: x
            real(dp) :: fx
        end function panel_kernel_t
    end interface

    ! QUADPACK alist/blist/rlist/elist workspace, one entry per subinterval.
    type :: integrate_workspace_t
        integer               :: limit = 0
        integer               :: last  = 0
        real(dp), allocatable :: a(:)
        real(dp), allocatable :: b(:)
        real(dp), allocatable :: r(:)
        real(dp), allocatable :: e(:)
        integer,  allocatable :: iord(:)
        integer,  allocatable :: level(:)
        integer,  allocatable :: ndin(:)
    end type integrate_workspace_t

    ! Wynn epsilon table (QUADPACK qelg workspace) for the QAGS path.
    type :: integrate_epstab_t
        integer  :: n = 0
        real(dp) :: tab(52) = 0.0_dp
        real(dp) :: result = 0.0_dp
        real(dp) :: abserr = huge(0.0_dp)
        real(dp) :: res3la(3) = 0.0_dp
        integer  :: nres = 0
    end type integrate_epstab_t

    ! Value, error, status, and the frozen accepted subdivision (the trace_rule
    ! schedule of ADR section 3.3).
    type :: integrate_result_t
        real(dp)               :: value  = 0.0_dp
        real(dp)               :: abserr = 0.0_dp
        integer                :: neval  = 0
        integer                :: nsub   = 0
        real(dp), allocatable  :: sub_a(:)
        real(dp), allocatable  :: sub_b(:)
        real(dp), allocatable  :: sub_r(:)
        real(dp), allocatable  :: sub_e(:)
        integer                :: key   = 21
        logical                :: extrapolated = .false.
        type(fortnum_status_t) :: status
    end type integrate_result_t

    real(dp), parameter :: epmach = epsilon(1.0_dp)
    real(dp), parameter :: uflow  = tiny(1.0_dp)
    real(dp), parameter :: oflow  = huge(1.0_dp)

    integer, parameter :: DEFAULT_LIMIT = 500

contains

    ! ------------------------------------------------------------------
    ! QAG: globally adaptive, selectable GK rule, no extrapolation.
    ! ------------------------------------------------------------------
    subroutine integrate_qag(f, a, b, epsabs, epsrel, workspace, result, &
                             status, key, limit, ctx)
        procedure(integrate_integrand_t)        :: f
        real(dp),                  intent(in)    :: a, b, epsabs, epsrel
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status
        integer,  intent(in), optional :: key, limit
        class(*), intent(in), optional :: ctx

        integer :: key_loc, limit_loc

        key_loc   = 21
        limit_loc = DEFAULT_LIMIT
        if (present(key))   key_loc   = key
        if (present(limit)) limit_loc = limit

        result%extrapolated = .false.
        result%key          = key_loc

        if (.not. valid_finite(a, b, epsabs, epsrel, limit_loc, key_loc, &
                               status)) then
            result%status = status
            return
        end if

        call ensure_workspace(workspace, limit_loc)
        call qag_core(.false.)
        result%status = status

    contains

        ! The driver core sits inside integrate_qag so panel_f sees f and ctx
        ! by host association without a module pointer. use_eps is false here.
        subroutine qag_core(use_eps)
            logical, intent(in) :: use_eps
            type(integrate_epstab_t) :: dummy_eps
            call driver(panel_f, a, b, epsabs, epsrel, key_loc, limit_loc, &
                        use_eps, workspace, dummy_eps, result, status)
        end subroutine qag_core

        function panel_f(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            fx = f(x, ctx)
        end function panel_f

    end subroutine integrate_qag

    ! ------------------------------------------------------------------
    ! QAGS: adaptive bisection plus Wynn epsilon extrapolation. Fixed GK21.
    ! ------------------------------------------------------------------
    subroutine integrate_qags(f, a, b, epsabs, epsrel, workspace, epstab, &
                              result, status, limit, ctx)
        procedure(integrate_integrand_t)        :: f
        real(dp),                  intent(in)    :: a, b, epsabs, epsrel
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status
        integer,  intent(in), optional :: limit
        class(*), intent(in), optional :: ctx

        integer :: limit_loc

        limit_loc = DEFAULT_LIMIT
        if (present(limit)) limit_loc = limit

        result%key = 21

        if (.not. valid_finite(a, b, epsabs, epsrel, limit_loc, 21, status)) then
            result%status = status
            return
        end if

        call ensure_workspace(workspace, limit_loc)
        call reset_epstab(epstab)
        call qags_core()
        result%status = status

    contains

        subroutine qags_core()
            call driver(panel_f, a, b, epsabs, epsrel, 21, limit_loc, &
                        .true., workspace, epstab, result, status)
        end subroutine qags_core

        function panel_f(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            fx = f(x, ctx)
        end function panel_f

    end subroutine integrate_qags

    ! ------------------------------------------------------------------
    ! QAGP: QAGS seeded with user break points so each known interior
    ! singularity or kink starts on a panel boundary. Break points coinciding
    ! with a or b are ambiguous (the break structure is not well defined for a
    ! frozen-subdivision derivative): reported as FORTNUM_DOMAIN_ERROR.
    ! ------------------------------------------------------------------
    subroutine integrate_qagp(f, a, b, points, epsabs, epsrel, workspace, &
                              epstab, result, status, limit, ctx)
        procedure(integrate_integrand_t)        :: f
        real(dp),                  intent(in)    :: a, b, epsabs, epsrel
        real(dp),                  intent(in)    :: points(:)
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status
        integer,  intent(in), optional :: limit
        class(*), intent(in), optional :: ctx

        real(dp) :: brk(size(points))
        real(dp) :: bnds(size(points) + 2)
        integer  :: limit_loc, nbrk, npan, i, j
        real(dp) :: bk, gap
        logical  :: dup

        limit_loc = DEFAULT_LIMIT
        if (present(limit)) limit_loc = limit

        result%key = 21

        if (.not. valid_finite(a, b, epsabs, epsrel, limit_loc, 21, status)) then
            result%status = status
            return
        end if

        ! Keep interior break points only; reject any that collapse onto an
        ! endpoint (ADR section 4.2 ambiguous break structure).
        gap = (b - a)*epmach*100.0_dp
        nbrk = 0
        do i = 1, size(points)
            bk = points(i)
            if (bk <= a + gap .or. bk >= b - gap) then
                if (bk > a - gap .and. bk < b + gap) then
                    call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        "integrate: break point coincides with an endpoint")
                    result%status = status
                    return
                end if
                cycle
            end if
            dup = .false.
            do j = 1, nbrk
                if (abs(brk(j) - bk) <= gap) dup = .true.
            end do
            if (.not. dup) then
                nbrk = nbrk + 1
                brk(nbrk) = bk
            end if
        end do

        ! Ascending sort of the kept break points (insertion; nbrk is small).
        do i = 2, nbrk
            bk = brk(i)
            j = i - 1
            do while (j >= 1)
                if (brk(j) <= bk) exit
                brk(j + 1) = brk(j)
                j = j - 1
            end do
            brk(j + 1) = bk
        end do

        bnds(1) = a
        do i = 1, nbrk
            bnds(i + 1) = brk(i)
        end do
        bnds(nbrk + 2) = b
        npan = nbrk + 1

        if (limit_loc < npan) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: limit below the seeded panel count")
            result%status = status
            return
        end if

        call ensure_workspace(workspace, limit_loc)
        call reset_epstab(epstab)
        do i = 1, npan
            workspace%a(i) = bnds(i)
            workspace%b(i) = bnds(i + 1)
            workspace%ndin(i) = 1   ! break-seeded panel (ADR section 3.1)
        end do
        call qagp_core()
        result%status = status

    contains

        subroutine qagp_core()
            call driver_seeded(panel_f, 21, epsabs, epsrel, limit_loc, npan, &
                               workspace, epstab, result, status)
        end subroutine qagp_core

        function panel_f(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            fx = f(x, ctx)
        end function panel_f

    end subroutine integrate_qagp

    ! ------------------------------------------------------------------
    ! QAGIU: semi-infinite or doubly infinite interval mapped to (0,1] by the
    ! QUADPACK dqagi transform, then driven by the QAGS machinery. The Jacobian
    ! 1/t**2 is folded into the host-associated panel wrapper; the endpoint
    ! singularity it puts at t = 0 is what the epsilon path is built for.
    !   inf = +1: [bound, +inf),  x = bound + (1-t)/t
    !   inf = -1: (-inf, bound],  x = bound - (1-t)/t
    !   inf = +2: (-inf, +inf), split at bound, both transforms, summed.
    ! ------------------------------------------------------------------
    subroutine integrate_qagiu(f, bound, inf, epsabs, epsrel, workspace, &
                               epstab, result, status, limit, ctx)
        procedure(integrate_integrand_t)        :: f
        real(dp),                  intent(in)    :: bound, epsabs, epsrel
        integer,                   intent(in)    :: inf
        type(integrate_workspace_t), intent(inout) :: workspace
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status
        integer,  intent(in), optional :: limit
        class(*), intent(in), optional :: ctx

        integer :: limit_loc, sgn

        limit_loc = DEFAULT_LIMIT
        if (present(limit)) limit_loc = limit

        result%key = 21

        ! The transformed interval is the finite (0,1]; reuse the finite
        ! validator on it. inf and a non-finite bound are the extra rejects.
        if (inf /= -1 .and. inf /= 1 .and. inf /= 2) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: inf must be one of -1, +1, +2")
            result%status = status
            return
        end if
        if (.not. (abs(bound) <= oflow)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: bound must be finite")
            result%status = status
            return
        end if
        if (.not. valid_finite(0.0_dp, 1.0_dp, epsabs, epsrel, limit_loc, 21, &
                               status)) then
            result%status = status
            return
        end if

        call ensure_workspace(workspace, limit_loc)
        call reset_epstab(epstab)
        sgn = inf
        call qagiu_core()
        result%status = status

    contains

        ! inf = +/-1 is one QAGS run on (0,1] with the sign-selected transform.
        ! inf = +2 is the sum of the +inf and -inf halves: two independent runs
        ! sharing the same transformed interval, each with its own sign, summed.
        subroutine qagiu_core()
            type(integrate_result_t) :: res2
            type(integrate_epstab_t) :: eps2
            type(integrate_workspace_t) :: ws2
            if (inf == 2) then
                ! Two independent semi-infinite integrals, summed. Each half is
                ! a full QAGS run on (0,1] with its own transform.
                call ensure_workspace(ws2, limit_loc)
                call reset_epstab(eps2)
                sgn = 1
                call driver(panel_f, 0.0_dp, 1.0_dp, epsabs, epsrel, 21, &
                            limit_loc, .true., workspace, epstab, result, &
                            status)
                sgn = -1
                call driver(panel_f, 0.0_dp, 1.0_dp, epsabs, epsrel, 21, &
                            limit_loc, .true., ws2, eps2, res2, status)
                result%value  = result%value + res2%value
                result%abserr = result%abserr + res2%abserr
                result%neval  = result%neval + res2%neval
            else
                call driver(panel_f, 0.0_dp, 1.0_dp, epsabs, epsrel, 21, &
                            limit_loc, .true., workspace, epstab, result, &
                            status)
            end if
        end subroutine qagiu_core

        ! Transformed integrand on t in (0,1]: x = bound + sgn*(1-t)/t,
        ! dx = dt/t**2. The Jacobian rides inside the wrapper, so the kernel
        ! never sees the transform (no module state).
        function panel_f(t) result(fx)
            real(dp), intent(in) :: t
            real(dp) :: fx, x, jac
            if (t <= 0.0_dp) then
                fx = 0.0_dp
                return
            end if
            x   = bound + real(sgn, dp)*(1.0_dp - t)/t
            jac = 1.0_dp/(t*t)
            fx  = f(x, ctx)*jac
        end function panel_f

    end subroutine integrate_qagiu

    ! ------------------------------------------------------------------
    ! Flat wrapper: finite interval, owns its own workspace and result.
    ! ------------------------------------------------------------------
    subroutine integrate(f, a, b, value, status, epsabs, epsrel, key, ctx)
        procedure(integrate_integrand_t)     :: f
        real(dp),               intent(in)    :: a, b
        real(dp),               intent(out)   :: value
        type(fortnum_status_t), intent(out)   :: status
        real(dp), intent(in), optional :: epsabs, epsrel
        integer,  intent(in), optional :: key
        class(*), intent(in), optional :: ctx

        type(integrate_workspace_t) :: workspace
        type(integrate_result_t)    :: result
        real(dp) :: epsabs_loc, epsrel_loc

        epsabs_loc = 0.0_dp
        epsrel_loc = 1.0e-8_dp
        if (present(epsabs)) epsabs_loc = epsabs
        if (present(epsrel)) epsrel_loc = epsrel

        call integrate_qag(f, a, b, epsabs_loc, epsrel_loc, workspace, result, &
                           status, key=key, ctx=ctx)
        value = result%value
    end subroutine integrate

    ! ------------------------------------------------------------------
    ! integrate_qag_jvp: forward-mode product for the trace_rule policy
    ! (ad.md sections 1, 4). The primal froze an accepted subdivision into
    ! `result`; the integral is then the fixed weighted sum of GK panel values
    !   I(p) = sum_panels  integral_{sub_a..sub_b}  f(x, p) dx.
    ! Differentiating at that frozen subdivision and pushing d/dp inside each
    ! panel integral gives
    !   dI/dp = sum_panels  integral_{sub_a..sub_b}  (df/dp)(x) dx,
    ! i.e. re-evaluate the SAME per-panel GK rule (same key, same frozen nodes)
    ! on the tangent integrand df/dp. dfdp is the directional derivative of the
    ! integrand along the parameter direction; for a vector parameter the caller
    ! supplies the contraction (df/dp).v as dfdp and gets the scalar dI/dp.v.
    !
    ! Active argument: dfdp (the integrand tangent). Inactive: the frozen trace
    ! (sub_a/sub_b/key/nsub) and the status side channel (ad.md section 3).
    !
    ! Status: the derivative is meaningful only on the trace the primal accepted.
    ! If the recorded primal status is not FORTNUM_OK (e.g. the bad-integrand
    ! code-3 path, where a perturbation would move or collapse a panel boundary),
    ! the frozen subdivision is not a valid linearization point and this reports
    ! the same non-smooth/failure status without computing a product. This reuses
    ! the existing status path (set_driver_status); no new code class.
    !
    ! VJP/grad: I(p) is scalar, so the reverse product is the same scalar
    ! sensitivity. integrate_qag_grad would equal integrate_qag_jvp with the
    ! identity tangent on the active parameter, and a VJP with a scalar seed u
    ! is u*dI/dp. Both are one scalar quadrature over the frozen trace, so a
    ! separate reverse entry adds no information over the forward product for a
    ! scalar output; the forward routine is the single derivative surface. An
    ! HVP needs the integrand's second parameter derivative d2f/dp2 on the frozen
    ! trace (same quadrature, second-order tangent); it is deferred until a
    ! second-order integrand-tangent contract exists (ad.md section 6, additive).
    ! ------------------------------------------------------------------
    subroutine integrate_qag_jvp(dfdp, result, di_dp, status, ctx)
        procedure(integrate_integrand_t)        :: dfdp
        type(integrate_result_t),  intent(in)    :: result
        real(dp),                  intent(out)   :: di_dp
        type(fortnum_status_t),    intent(out)   :: status
        class(*), intent(in), optional :: ctx

        real(dp) :: panel_r, panel_e, panel_resabs, panel_resasc
        integer  :: i, key_loc

        di_dp = 0.0_dp

        ! A perturbation that would change the accepted subdivision shows up as
        ! a non-OK recorded primal status; the frozen-trace derivative is then
        ! not well-posed. Propagate that verdict and stop (reuse status path).
        if (result%status%code /= FORTNUM_OK) then
            status = result%status
            return
        end if
        if (.not. allocated(result%sub_a) .or. result%nsub < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate_qag_jvp: empty frozen subdivision trace")
            return
        end if

        call status_set(status, FORTNUM_OK, "")
        key_loc = result%key

        ! Re-walk the frozen panels and apply the same GK rule to the tangent.
        ! gk_apply%result is the per-panel integral estimate; summing it over
        ! the recorded subdivision is the frozen-trace derivative of I.
        do i = 1, result%nsub
            call gk_apply(tangent_f, key_loc, result%sub_a(i), &
                          result%sub_b(i), panel_r, panel_e, &
                          panel_resabs, panel_resasc)
            di_dp = di_dp + panel_r
        end do

    contains

        ! Host-associated tangent wrapper: bridges the caller's df/dp and its
        ! optional ctx to the kernel's ctx-free integrand. Lives on this call
        ! stack; no module state, matching the primal panel_f pattern.
        function tangent_f(x) result(fx)
            real(dp), intent(in) :: x
            real(dp) :: fx
            fx = dfdp(x, ctx)
        end function tangent_f

    end subroutine integrate_qag_jvp

    ! ==================================================================
    ! Shared globally adaptive driver. panel_f is the ctx-free integrand
    ! bridged by the caller's host-associated wrapper; never a module pointer.
    ! ==================================================================
    subroutine driver(panel_f, a, b, epsabs, epsrel, key, limit, use_eps, &
                      ws, epstab, result, status)
        procedure(panel_kernel_t)               :: panel_f
        real(dp),                  intent(in)    :: a, b, epsabs, epsrel
        integer,                   intent(in)    :: key, limit
        logical,                   intent(in)    :: use_eps
        type(integrate_workspace_t), intent(inout) :: ws
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status

        real(dp) :: r0, e0, resabs0, resasc0
        real(dp) :: area, errsum, errbnd
        integer  :: neval_panel

        call status_set(status, FORTNUM_OK, "")
        neval_panel = panel_neval(key)

        ! First panel over the whole interval seeds the work stack.
        call gk_apply(panel_f, key, a, b, r0, e0, resabs0, resasc0)
        ws%last = 1
        ws%a(1) = a
        ws%b(1) = b
        ws%r(1) = r0
        ws%e(1) = e0
        ws%iord(1) = 1
        ws%level(1) = 0
        ws%ndin(1) = 0
        result%neval = neval_panel

        area   = r0
        errsum = e0
        errbnd = max(epsabs, epsrel*abs(area))

        ! Single-panel acceptance, with QUADPACK's resasc guard against a
        ! deceptively small estimate on an oscillatory smooth panel.
        if (e0 <= 100.0_dp*epmach*resabs0 .and. e0 > errbnd) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: roundoff dominates the first panel")
        else if (e0 <= errbnd .and. e0 /= resasc0) then
            call finalize(ws, result, area, errsum, .false., 0.0_dp)
            return
        else if (limit == 1) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: limit reached before tolerance")
        end if
        if (status%code /= FORTNUM_OK) then
            call finalize(ws, result, area, errsum, .false., 0.0_dp)
            return
        end if

        call adapt_loop(panel_f, key, epsabs, epsrel, limit, use_eps, &
                        neval_panel, area, errsum, ws, epstab, result, status)
    end subroutine driver

    ! ==================================================================
    ! Same globally adaptive driver, but the work stack is pre-seeded with a
    ! supplied subdivision (QAGP break points). The caller fills ws%a/ws%b for
    ! npan panels and marks the break-seeded ones in ws%ndin; this evaluates the
    ! GK pair on each seed panel, then runs the shared adaptive loop. use_eps is
    ! always true here (QAGP rides the QAGS extrapolation machinery).
    ! ==================================================================
    subroutine driver_seeded(panel_f, key, epsabs, epsrel, limit, npan, &
                             ws, epstab, result, status)
        procedure(panel_kernel_t)               :: panel_f
        integer,                   intent(in)    :: key, limit, npan
        real(dp),                  intent(in)    :: epsabs, epsrel
        type(integrate_workspace_t), intent(inout) :: ws
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status

        real(dp) :: r0, e0, resabs0, resasc0
        real(dp) :: area, errsum, errbnd
        integer  :: i, neval_panel

        call status_set(status, FORTNUM_OK, "")
        neval_panel = panel_neval(key)
        area   = 0.0_dp
        errsum = 0.0_dp
        do i = 1, npan
            call gk_apply(panel_f, key, ws%a(i), ws%b(i), r0, e0, &
                          resabs0, resasc0)
            ws%r(i) = r0
            ws%e(i) = e0
            ws%level(i) = 0
            area   = area + r0
            errsum = errsum + e0
        end do
        ws%last = npan
        result%neval = npan*neval_panel
        call reorder(ws)

        errbnd = max(epsabs, epsrel*abs(area))
        if (errsum <= errbnd) then
            call finalize(ws, result, area, errsum, .false., 0.0_dp)
            return
        else if (limit <= npan) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: limit reached before tolerance")
            call finalize(ws, result, area, errsum, .false., 0.0_dp)
            return
        end if

        call adapt_loop(panel_f, key, epsabs, epsrel, limit, .true., &
                        neval_panel, area, errsum, ws, epstab, result, status)
    end subroutine driver_seeded

    ! ==================================================================
    ! Shared bisection loop. Both seed paths hand off here after the work stack
    ! holds at least one panel and the single-panel acceptance test (if any) has
    ! passed. Owns the worst-panel pop, the GK re-evaluation on the two halves,
    ! the roundoff/limit detectors, the optional epsilon extrapolation, and the
    ! final value choice between direct sum and extrapolated value.
    ! ==================================================================
    subroutine adapt_loop(panel_f, key, epsabs, epsrel, limit, use_eps, &
                          neval_panel, area, errsum, ws, epstab, result, status)
        procedure(panel_kernel_t)               :: panel_f
        integer,                   intent(in)    :: key, limit, neval_panel
        real(dp),                  intent(in)    :: epsabs, epsrel
        logical,                   intent(in)    :: use_eps
        real(dp),                  intent(inout) :: area, errsum
        type(integrate_workspace_t), intent(inout) :: ws
        type(integrate_epstab_t),  intent(inout) :: epstab
        type(integrate_result_t),  intent(inout) :: result
        type(fortnum_status_t),    intent(out)   :: status

        real(dp) :: errbnd, errmax
        real(dp) :: a1, b1, a2, b2
        real(dp) :: r1, e1, resabs1, resasc1
        real(dp) :: r2, e2, resabs2, resasc2
        real(dp) :: r12, e12
        integer  :: maxerr, ncode, iroff1, iroff2, iroff3
        logical  :: extrap_better

        call status_set(status, FORTNUM_OK, "")
        ncode = 0
        iroff1 = 0
        iroff2 = 0
        iroff3 = 0
        extrap_better = .false.
        errbnd = max(epsabs, epsrel*abs(area))

        call eps_init(epstab, use_eps, area)

        do while (ws%last < limit)
            ! Pop the worst-error subinterval (iord(1)) and bisect it.
            maxerr = ws%iord(1)
            errmax = ws%e(maxerr)
            a1 = ws%a(maxerr)
            b1 = 0.5_dp*(ws%a(maxerr) + ws%b(maxerr))
            a2 = b1
            b2 = ws%b(maxerr)

            call gk_apply(panel_f, key, a1, b1, r1, e1, resabs1, resasc1)
            call gk_apply(panel_f, key, a2, b2, r2, e2, resabs2, resasc2)
            result%neval = result%neval + 2*neval_panel

            r12 = r1 + r2
            e12 = e1 + e2

            ! A bisection so deep that a GK node lands on the singularity
            ! itself (or its panel sum overflows) is the QUADPACK code-3
            ! "bad integrand behaviour" case: keep the last finite area and
            ! stop before it corrupts the running totals (ADR section 4.2).
            if (.not. is_finite(r12) .or. .not. is_finite(e12)) then
                ncode = 3
                exit
            end if

            errsum = errsum + e12 - errmax
            area   = area + r12 - ws%r(maxerr)

            call roundoff_count(ws%r(maxerr), r12, errmax, e12, ws%last, &
                                resasc1, e1, resasc2, e2, iroff1, iroff2, &
                                iroff3, use_eps)

            ! Replace maxerr with the left half, append the right half.
            ws%level(ws%last + 1) = ws%level(maxerr) + 1
            ws%level(maxerr)      = ws%level(maxerr) + 1
            ws%a(maxerr) = a1
            ws%b(maxerr) = b1
            ws%r(maxerr) = r1
            ws%e(maxerr) = e1
            ws%last = ws%last + 1
            ws%a(ws%last) = a2
            ws%b(ws%last) = b2
            ws%r(ws%last) = r2
            ws%e(ws%last) = e2
            ws%ndin(ws%last) = 0
            call reorder(ws)

            errbnd = max(epsabs, epsrel*abs(area))
            if (errsum <= errbnd) exit

            ! Convergence-failure detectors (QUADPACK iroffN, code 1, code 3).
            if (iroff1 + iroff2 >= 10 .or. iroff3 >= 20) ncode = 2
            if (iroff2 >= 5) ncode = 4
            if (ws%last == limit) ncode = 1
            if (max(abs(a1), abs(b2)) <= (1.0_dp + 100.0_dp*epmach)* &
                (abs(a2) + 1000.0_dp*uflow)) ncode = 3
            if (ncode /= 0) exit

            if (use_eps) call eps_step(epstab, ws, area)
        end do

        ! On the way out, choose between the direct panel sum and the best
        ! extrapolation when the direct sum did not meet the bound.
        if (use_eps .and. epstab%nres >= 1 .and. errsum > errbnd) then
            if (epstab%abserr < errsum) then
                area = epstab%result
                errsum = epstab%abserr
                extrap_better = .true.
            end if
        end if

        ! A panel that overflowed at a break-point singularity (is_finite guard,
        ! ncode 3) is not a non-smooth verdict if the epsilon table already met
        ! the tolerance: the bisection merely probed past the useful depth, and
        ! the extrapolated value is the converged answer. Clear that code.
        errbnd = max(epsabs, epsrel*abs(area))
        if (ncode == 3 .and. extrap_better .and. errsum <= errbnd) ncode = 0

        call set_driver_status(status, ncode, errsum, errbnd)
        call finalize(ws, result, area, errsum, extrap_better, epstab%result)
    end subroutine adapt_loop

    ! ------------------------------------------------------------------
    ! Helpers.
    ! ------------------------------------------------------------------

    ! True for an ordinary finite value; false for NaN or +/-Inf. Used to stop
    ! a singular bisection before a blown-up panel corrupts the running totals.
    pure logical function is_finite(x)
        real(dp), intent(in) :: x
        is_finite = (x == x) .and. (abs(x) <= oflow)
    end function is_finite

    pure integer function panel_neval(key)
        integer, intent(in) :: key
        select case (key)
        case (15);  panel_neval = 15
        case (31);  panel_neval = 31
        case (61);  panel_neval = 61
        case default; panel_neval = 21
        end select
    end function panel_neval

    ! Reject invalid input before any panel runs (ADR section 4.1).
    logical function valid_finite(a, b, epsabs, epsrel, limit, key, status)
        real(dp),               intent(in)  :: a, b, epsabs, epsrel
        integer,                intent(in)  :: limit, key
        type(fortnum_status_t), intent(out) :: status
        real(dp), parameter :: eps_floor = max(50.0_dp*epmach, 0.5e-28_dp)
        valid_finite = .false.
        if (b <= a) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: require a < b")
            return
        end if
        if (limit < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: limit must be >= 1")
            return
        end if
        if (key /= 15 .and. key /= 21 .and. key /= 31 .and. key /= 61) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: key must be one of 15, 21, 31, 61")
            return
        end if
        if (epsabs <= 0.0_dp .and. epsrel < eps_floor) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: epsabs/epsrel below the kernel floor")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
        valid_finite = .true.
    end function valid_finite

    ! Allocate the work-stack columns once when the capacity changes; reuse
    ! otherwise. No allocation in the bisection loop (ADR section 3.1).
    subroutine ensure_workspace(ws, limit)
        type(integrate_workspace_t), intent(inout) :: ws
        integer,                     intent(in)    :: limit
        if (ws%limit /= limit .or. .not. allocated(ws%a)) then
            if (allocated(ws%a))     deallocate (ws%a)
            if (allocated(ws%b))     deallocate (ws%b)
            if (allocated(ws%r))     deallocate (ws%r)
            if (allocated(ws%e))     deallocate (ws%e)
            if (allocated(ws%iord))  deallocate (ws%iord)
            if (allocated(ws%level)) deallocate (ws%level)
            if (allocated(ws%ndin))  deallocate (ws%ndin)
            allocate (ws%a(limit), ws%b(limit), ws%r(limit), ws%e(limit))
            allocate (ws%iord(limit), ws%level(limit), ws%ndin(limit))
            ws%limit = limit
        end if
        ws%last = 0
        ws%ndin = 0
    end subroutine ensure_workspace

    subroutine reset_epstab(epstab)
        type(integrate_epstab_t), intent(inout) :: epstab
        epstab%n      = 0
        epstab%tab    = 0.0_dp
        epstab%result = 0.0_dp
        epstab%abserr = oflow
        epstab%res3la = 0.0_dp
        epstab%nres   = 0
    end subroutine reset_epstab

    ! Maintain iord so iord(1) is the largest-error subinterval. The adaptive
    ! loop pops only iord(1), so the full descending sort the old code built was
    ! dead work beyond locating the maximum. Find that maximum in one O(last)
    ! pass with the same strict-greater, lowest-index tie-break the selection
    ! sort used, so iord(1) (hence the pop sequence and every downstream value)
    ! is identical; fill the rest with the identity to keep the array defined.
    subroutine reorder(ws)
        type(integrate_workspace_t), intent(inout) :: ws
        integer  :: i, k
        real(dp) :: emax
        do i = 1, ws%last
            ws%iord(i) = i
        end do
        if (ws%last < 2) return
        k    = 1
        emax = ws%e(1)
        do i = 2, ws%last
            if (ws%e(i) > emax) then
                emax = ws%e(i)
                k = i
            end if
        end do
        if (k /= 1) then
            ws%iord(1) = k
            ws%iord(k) = 1
        end if
    end subroutine reorder

    ! QUADPACK roundoff detectors. iroff1/iroff2 catch a bisection that fails
    ! to reduce the error; iroff3 is the extrapolation-phase guard.
    subroutine roundoff_count(rold, rnew, eold, enew, last, resasc1, e1, &
                              resasc2, e2, iroff1, iroff2, iroff3, use_eps)
        real(dp), intent(in)    :: rold, rnew, eold, enew, resasc1, e1
        real(dp), intent(in)    :: resasc2, e2
        integer,  intent(in)    :: last
        integer,  intent(inout) :: iroff1, iroff2, iroff3
        logical,  intent(in)    :: use_eps
        if (resasc1 /= e1 .and. resasc2 /= e2) then
            if (abs(rold - rnew) <= 1.0e-5_dp*abs(rnew) .and. &
                enew >= 0.99_dp*eold) then
                if (use_eps) then
                    iroff3 = iroff3 + 1
                else
                    iroff1 = iroff1 + 1
                end if
            end if
            if (last > 10 .and. enew > eold) iroff2 = iroff2 + 1
        end if
    end subroutine roundoff_count

    subroutine set_driver_status(status, ncode, errsum, errbnd)
        type(fortnum_status_t), intent(out) :: status
        integer,                intent(in)  :: ncode
        real(dp),               intent(in)  :: errsum, errbnd
        select case (ncode)
        case (0)
            if (errsum <= errbnd) then
                call status_set(status, FORTNUM_OK, "")
            else
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                                "integrate: tolerance not reached")
            end if
        case (1)
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: max subdivisions (limit) reached")
        case (2)
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: roundoff prevents reaching tolerance")
        case (3)
            ! Non-smooth: subdivision collapsed; the frozen-subdivision
            ! derivative is not well-posed (ADR section 4.2).
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "integrate: bad integrand behaviour, "// &
                            "non-integrable on this interval")
        case (4)
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: slow convergence, extrapolation stalled")
        case default
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "integrate: convergence failure")
        end select
    end subroutine set_driver_status

    ! Copy the work stack into the result trace in interval (left-to-right)
    ! order; this is the frozen trace_rule subdivision (ADR section 3.3).
    subroutine finalize(ws, result, area, errsum, extrapolated, eps_value)
        type(integrate_workspace_t), intent(in)    :: ws
        type(integrate_result_t),    intent(inout) :: result
        real(dp),                    intent(in)    :: area, errsum, eps_value
        logical,                     intent(in)    :: extrapolated
        integer :: i, j, n
        integer, allocatable :: ord(:)
        n = ws%last
        call ensure_trace(result, n)
        ! Sort indices by left endpoint for a left-to-right trace.
        allocate (ord(n))
        do i = 1, n
            ord(i) = i
        end do
        call sort_by_left(ws, ord, n)
        do i = 1, n
            j = ord(i)
            result%sub_a(i) = ws%a(j)
            result%sub_b(i) = ws%b(j)
            result%sub_r(i) = ws%r(j)
            result%sub_e(i) = ws%e(j)
        end do
        result%nsub   = n
        result%abserr = errsum
        result%extrapolated = extrapolated
        if (extrapolated) then
            result%value = eps_value
        else
            ! Sum fresh to avoid the running-update cancellation.
            result%value = sum(result%sub_r(1:n))
        end if
    end subroutine finalize

    subroutine sort_by_left(ws, ord, n)
        type(integrate_workspace_t), intent(in)    :: ws
        integer,                     intent(inout) :: ord(:)
        integer,                     intent(in)    :: n
        integer  :: i, j, k
        real(dp) :: amin
        do i = 1, n - 1
            k = i
            amin = ws%a(ord(i))
            do j = i + 1, n
                if (ws%a(ord(j)) < amin) then
                    amin = ws%a(ord(j))
                    k = j
                end if
            end do
            if (k /= i) then
                j = ord(i)
                ord(i) = ord(k)
                ord(k) = j
            end if
        end do
    end subroutine sort_by_left

    subroutine ensure_trace(result, n)
        type(integrate_result_t), intent(inout) :: result
        integer,                  intent(in)    :: n
        logical :: grow
        grow = .not. allocated(result%sub_a)
        if (.not. grow) grow = (size(result%sub_a) < n)
        if (grow) then
            if (allocated(result%sub_a)) deallocate (result%sub_a)
            if (allocated(result%sub_b)) deallocate (result%sub_b)
            if (allocated(result%sub_r)) deallocate (result%sub_r)
            if (allocated(result%sub_e)) deallocate (result%sub_e)
            allocate (result%sub_a(n), result%sub_b(n))
            allocate (result%sub_r(n), result%sub_e(n))
        end if
    end subroutine ensure_trace

    ! ---- Wynn epsilon extrapolation (QUADPACK qelg). ----

    subroutine eps_init(epstab, use_eps, area)
        type(integrate_epstab_t), intent(inout) :: epstab
        logical,                  intent(in)    :: use_eps
        real(dp),                 intent(in)    :: area
        if (.not. use_eps) return
        epstab%n = 1
        epstab%tab(1) = area
        epstab%abserr = oflow
        epstab%nres = 0
    end subroutine eps_init

    ! Append the current total area as a new sequence element and extrapolate.
    subroutine eps_step(epstab, ws, area)
        type(integrate_epstab_t),    intent(inout) :: epstab
        type(integrate_workspace_t), intent(in)    :: ws
        real(dp),                    intent(in)    :: area
        if (ws%last < 3) return
        epstab%n = epstab%n + 1
        if (epstab%n > 50) then
            ! Shift the column to keep qelg's width bound (QUADPACK).
            epstab%tab(1:49) = epstab%tab(3:51)
            epstab%n = 49
        end if
        epstab%tab(epstab%n) = area
        if (epstab%n >= 3) call qelg(epstab)
    end subroutine eps_step

    ! Wynn epsilon algorithm on epstab%tab(1:n). Produces epstab%result and
    ! epstab%abserr, with the QUADPACK three-term roundoff guard via res3la.
    subroutine qelg(epstab)
        type(integrate_epstab_t), intent(inout) :: epstab
        real(dp) :: tab(52)
        real(dp) :: e1, e2, e3, res, err
        real(dp) :: err2, err3, ss, epsinf, tol2, tol3
        integer  :: n, newelm, i, ib, ib2, k1

        n = epstab%n
        tab = epstab%tab
        epstab%nres = epstab%nres + 1
        if (n < 3) then
            epstab%result = tab(n)
            epstab%abserr = oflow
            return
        end if
        k1 = n
        res = tab(n)
        err = oflow
        newelm = (n - 1)/2
        do i = 1, newelm
            k1 = k1 - 1
            e1 = tab(k1 + 1)
            e2 = res
            e3 = tab(k1)
            err2 = abs(e2 - e1)
            tol2 = max(abs(e2), abs(e1))*epmach
            err3 = abs(e1 - e3)
            tol3 = max(abs(e1), abs(e3))*epmach
            if (err2 <= tol2 .and. err3 <= tol3) then
                res = e2
                err = err2 + err3
                exit
            end if
            ! Three-point Wynn lozenge step.
            ss = 1.0_dp/(e2 - e1) - 1.0_dp/(e1 - e3)
            epsinf = abs(ss*e1)
            if (epsinf <= 1.0e-4_dp) then
                res = e2
                err = err2 + err3
                exit
            end if
            res = e1 + 1.0_dp/ss
            tab(k1) = res
            err = err2 + abs(res - e2) + err3
            if (err <= epstab%abserr) then
                epstab%abserr = err
                epstab%result = res
            end if
        end do
        if (epstab%abserr == oflow) then
            epstab%result = res
            epstab%abserr = err
        end if
        ! Three-term roundoff floor (QUADPACK res3la).
        if (epstab%nres >= 4) then
            epstab%abserr = max(epstab%abserr, &
                5.0_dp*epmach*max(abs(epstab%result), abs(res)))
        end if
        epstab%res3la(1) = epstab%res3la(2)
        epstab%res3la(2) = epstab%res3la(3)
        epstab%res3la(3) = epstab%result
        ! Compress the worked column back into the table for the next call.
        do i = 1, newelm + 1
            ib = 2*i - 1
            ib2 = ib + 2
            if (ib2 <= n) tab(ib) = tab(ib2)
        end do
        if (mod(n, 2) == 0) then
            do i = 1, newelm
                k1 = 2*i
                if (k1 + 2 <= n + 2) tab(k1) = tab(k1 + 2)
            end do
        end if
        epstab%tab = tab
    end subroutine qelg

end module fortnum_integrate
