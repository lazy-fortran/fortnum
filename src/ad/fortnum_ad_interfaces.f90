module fortnum_ad_interfaces
    ! Backend-opaque derivative-product interfaces for optimizer-facing code
    ! (issue #41).
    !
    ! A downstream code consumes a fortnum derivative product through one of
    ! these abstract interfaces without knowing how the product was produced.
    ! Whether the Jacobian-vector product came from an analytic recurrence, an
    ! implicit-function rule, a frozen adaptive trace, an Enzyme-generated
    ! pass, or a finite-difference fallback is reported through the status, not
    ! through the call shape. The optimizer wires a kernel that matches one of
    ! these signatures and reads the backend and quality fields to decide how
    ! much to trust the result.
    !
    ! context is class(*): the kernel carries its own configuration and
    ! workspace through it (a derived type with profiles, grids, fixed traces).
    ! No module-level or global pointers; the optimizer owns the context and
    ! threads it through every call, so the same kernel code is reentrant and
    ! safe to call from parallel optimizer evaluations.
    use fortnum_kinds,  only: dp
    use fortnum_status, only: fortnum_status_t
    implicit none
    private

    ! ------------------------------------------------------------ backends
    !
    ! Which machinery produced the derivative product. Mirrors the policy
    ! classes of docs/design/ad.md §1 plus the finite-difference fallback that
    ! a kernel may use when no exact product exists. The optimizer treats the
    ! product as opaque; the backend tag is advisory metadata for logging,
    ! trust thresholds, and step-acceptance heuristics.
    integer, parameter, public :: FORTNUM_AD_BACKEND_NONE        = 0
    integer, parameter, public :: FORTNUM_AD_BACKEND_ANALYTIC    = 1
    integer, parameter, public :: FORTNUM_AD_BACKEND_IMPLICIT    = 2
    integer, parameter, public :: FORTNUM_AD_BACKEND_TRACE       = 3
    integer, parameter, public :: FORTNUM_AD_BACKEND_GENERATED   = 4
    integer, parameter, public :: FORTNUM_AD_BACKEND_FINITE_DIFF = 5

    ! ------------------------------------------------------------ quality
    !
    ! How good the product is, independent of the backend. EXACT means correct
    ! to rounding (analytic, implicit, transparent-on-trace). APPROXIMATE means
    ! a controlled truncation error (finite difference, a frozen trace differing
    ! from the true adaptive schedule at the perturbed point). NONSMOOTH flags a
    ! point where the derivative is not defined (a branch or event boundary,
    ! ad.md §3); the optimizer must not trust the product there.
    integer, parameter, public :: FORTNUM_AD_QUALITY_UNKNOWN     = 0
    integer, parameter, public :: FORTNUM_AD_QUALITY_EXACT       = 1
    integer, parameter, public :: FORTNUM_AD_QUALITY_APPROXIMATE = 2
    integer, parameter, public :: FORTNUM_AD_QUALITY_NONSMOOTH   = 3

    ! Derivative status carrier. Extends the side-channel role of
    ! fortnum_status_t (ad.md §3: status is inactive, never a differentiable
    ! output) with the backend and quality tags this layer needs. The embedded
    ! fortnum_status_t still reports the primal error code and message, so a
    ! caller that only checks status_ok keeps working unchanged.
    type, public :: fortnum_ad_status_t
        type(fortnum_status_t) :: status
        integer :: backend = FORTNUM_AD_BACKEND_NONE
        integer :: quality = FORTNUM_AD_QUALITY_UNKNOWN
    end type fortnum_ad_status_t

    public :: ad_status_set
    public :: ad_status_ok

    ! ----------------------------------------------------- kernel interfaces
    !
    ! Backend-opaque derivative products. n is the flat active-vector length
    ! (see fortnum_active_vector). All real arrays are contiguous explicit-shape
    ! to match the layout the Enzyme path is tested against first (ad.md §3).
    ! context is class(*), supplied by the optimizer and threaded unchanged.
    abstract interface

        ! Primal: y = f(x). y may be scalar (m = 1) or a residual vector.
        subroutine value_fn(n, x, y, context, status)
            import :: dp, fortnum_ad_status_t
            integer,                  intent(in)    :: n
            real(dp),                 intent(in)    :: x(n)
            real(dp),                 intent(out)   :: y(:)
            class(*),                 intent(inout) :: context
            type(fortnum_ad_status_t), intent(out)  :: status
        end subroutine value_fn

        ! Forward product: y = f(x), y_dot = J(x) x_dot.
        subroutine jvp_fn(n, x, x_dot, y, y_dot, context, status)
            import :: dp, fortnum_ad_status_t
            integer,                  intent(in)    :: n
            real(dp),                 intent(in)    :: x(n)
            real(dp),                 intent(in)    :: x_dot(n)
            real(dp),                 intent(out)   :: y(:)
            real(dp),                 intent(out)   :: y_dot(:)
            class(*),                 intent(inout) :: context
            type(fortnum_ad_status_t), intent(out)  :: status
        end subroutine jvp_fn

        ! Reverse product: x_bar = J(x)^T y_bar.
        subroutine vjp_fn(n, x, y_bar, x_bar, context, status)
            import :: dp, fortnum_ad_status_t
            integer,                  intent(in)    :: n
            real(dp),                 intent(in)    :: x(n)
            real(dp),                 intent(in)    :: y_bar(:)
            real(dp),                 intent(out)   :: x_bar(n)
            class(*),                 intent(inout) :: context
            type(fortnum_ad_status_t), intent(out)  :: status
        end subroutine vjp_fn

        ! Scalar objective and its gradient: f = f(x), g = grad f(x).
        subroutine grad_fn(n, x, f, g, context, status)
            import :: dp, fortnum_ad_status_t
            integer,                  intent(in)    :: n
            real(dp),                 intent(in)    :: x(n)
            real(dp),                 intent(out)   :: f
            real(dp),                 intent(out)   :: g(n)
            class(*),                 intent(inout) :: context
            type(fortnum_ad_status_t), intent(out)  :: status
        end subroutine grad_fn

        ! Hessian-vector product of a scalar objective: f = f(x),
        ! hv = (grad^2 f(x)) v.
        subroutine hvp_fn(n, x, v, f, hv, context, status)
            import :: dp, fortnum_ad_status_t
            integer,                  intent(in)    :: n
            real(dp),                 intent(in)    :: x(n)
            real(dp),                 intent(in)    :: v(n)
            real(dp),                 intent(out)   :: f
            real(dp),                 intent(out)   :: hv(n)
            class(*),                 intent(inout) :: context
            type(fortnum_ad_status_t), intent(out)  :: status
        end subroutine hvp_fn

    end interface

contains

    ! Sets the embedded primal status plus the backend and quality tags in one
    ! call, the common path for a kernel reporting a successful product.
    pure subroutine ad_status_set(s, code, msg, backend, quality)
        use fortnum_status, only: status_set
        type(fortnum_ad_status_t), intent(out) :: s
        integer,      intent(in) :: code
        character(*), intent(in) :: msg
        integer,      intent(in) :: backend
        integer,      intent(in) :: quality
        call status_set(s%status, code, msg)
        s%backend = backend
        s%quality = quality
    end subroutine ad_status_set

    ! True iff the primal status is OK and the quality is not the unusable
    ! NONSMOOTH verdict. An optimizer gates a step on this before trusting g.
    pure logical function ad_status_ok(s)
        use fortnum_status, only: status_ok
        type(fortnum_ad_status_t), intent(in) :: s
        ad_status_ok = status_ok(s%status) .and. &
            (s%quality /= FORTNUM_AD_QUALITY_NONSMOOTH)
    end function ad_status_ok

end module fortnum_ad_interfaces
