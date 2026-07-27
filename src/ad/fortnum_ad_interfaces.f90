module fortnum_ad_interfaces
    ! Backend-opaque derivative-product interfaces for optimizer-facing code
    ! (issue #41).
    !
    ! A downstream code consumes a fortnum derivative product through one of
    ! these abstract interfaces without knowing how the product was produced.
    ! Whether the Jacobian-vector product came from an analytical recurrence,
    ! an implicit-function rule, a frozen adaptive trace, an autodiff-generated
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
    integer, parameter, public :: FORTNUM_AD_BACKEND_NONE = 0
    integer, parameter, public :: FORTNUM_AD_BACKEND_ANALYTICAL = 1
    integer, parameter, public :: FORTNUM_AD_BACKEND_AUTODIFF = 2
    integer, parameter, public :: FORTNUM_AD_BACKEND_HYBRID = 3
    integer, parameter, public :: FORTNUM_AD_BACKEND_FINITE_DIFFERENCE_REFERENCE = 4

    ! Compatibility aliases. New public code uses autodiff, analytical, and
    ! hybrid; these names remain so existing downstream sources keep building.
    integer, parameter, public :: FORTNUM_AD_BACKEND_ANALYTIC = &
        FORTNUM_AD_BACKEND_ANALYTICAL
    integer, parameter, public :: FORTNUM_AD_BACKEND_IMPLICIT = &
        FORTNUM_AD_BACKEND_ANALYTICAL
    integer, parameter, public :: FORTNUM_AD_BACKEND_TRACE = &
        FORTNUM_AD_BACKEND_ANALYTICAL
    integer, parameter, public :: FORTNUM_AD_BACKEND_GENERATED = &
        FORTNUM_AD_BACKEND_AUTODIFF
    integer, parameter, public :: FORTNUM_AD_BACKEND_FINITE_DIFF = &
        FORTNUM_AD_BACKEND_FINITE_DIFFERENCE_REFERENCE

    ! ------------------------------------------------------------ quality
    !
    ! How good the product is, independent of the backend. EXACT means correct
    ! to rounding (analytical, implicit, or frozen-trace). APPROXIMATE means
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
    public :: ad_status_merge

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

    ! Merge provenance and quality when derivative products cross an operator
    ! boundary. Mixed mechanisms are hybrid. The first primal failure is
    ! retained, while derivative quality follows an explicit worst-case order:
    ! nonsmooth, unknown, approximate, exact.
    pure subroutine ad_status_merge(left, right, merged)
        use fortnum_status, only: FORTNUM_OK
        type(fortnum_ad_status_t), intent(in) :: left, right
        type(fortnum_ad_status_t), intent(out) :: merged

        if (left%status%code /= FORTNUM_OK) then
            merged%status = left%status
        else
            merged%status = right%status
        end if
        merged%backend = merge_backend(left%backend, right%backend)
        if (quality_rank(left%quality) >= quality_rank(right%quality)) then
            merged%quality = left%quality
        else
            merged%quality = right%quality
        end if
    end subroutine ad_status_merge

    pure integer function merge_backend(left, right)
        integer, intent(in) :: left, right
        if (left == FORTNUM_AD_BACKEND_NONE) then
            merge_backend = right
        else if (right == FORTNUM_AD_BACKEND_NONE .or. left == right) then
            merge_backend = left
        else
            merge_backend = FORTNUM_AD_BACKEND_HYBRID
        end if
    end function merge_backend

    pure integer function quality_rank(quality)
        integer, intent(in) :: quality
        select case (quality)
        case (FORTNUM_AD_QUALITY_NONSMOOTH)
            quality_rank = 3
        case (FORTNUM_AD_QUALITY_UNKNOWN)
            quality_rank = 2
        case (FORTNUM_AD_QUALITY_APPROXIMATE)
            quality_rank = 1
        case default
            quality_rank = 0
        end select
    end function quality_rank

end module fortnum_ad_interfaces
