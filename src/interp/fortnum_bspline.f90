module fortnum_bspline
    ! B-spline basis functions and derivatives over a breakpoint set, matching
    ! the GSL/FGSL convention used by NEO-2 (gsl_bspline_routines_mod.f90,
    ! collop_bspline.f90): "order" k is the spline order (k = degree + 1), the
    ! caller supplies nbreak breakpoints, and the augmented knot vector repeats
    ! the two end breakpoints with multiplicity k so the basis is clamped /
    ! interpolatory at the boundary. The number of basis functions (coefficients)
    ! is ncoef = nbreak + k - 2.
    !
    ! Reference algorithm (clean-room, from published sources, no GSL code):
    !   Cox-de Boor recursion for the basis B_{i,k}(x) and the all-orders
    !   derivative table; de Boor, "A Practical Guide to Splines" (2001),
    !   chapter X; Piegl & Tiller, "The NURBS Book", algorithms A2.2/A2.3
    !   (basis funs and derivative basis funs); DLMF 1.B-splines.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4):
    !   Default class: transparent inside a fixed knot span. Within one span the
    !   spline value sum_i c_i B_{i,k}(x) is a fixed polynomial in x; the active
    !   variable is x. The knot index (which span x falls in) is selected by
    !   grid_search and is primal_only: crossing a knot is a non-smooth event and
    !   the caller must hold the span fixed (see check_smoothness in the AD test).
    !   The derivative weights are produced by the analytic recurrence in
    !   bspline_eval_deriv; the JVP/VJP entry points below expose them.
    !   Active: x. Inactive: order k, the knot/breakpoint array, nderiv.
    !
    ! Out of fortnum scope: the Taylor extrapolation NEO-2 applies above the last
    ! breakpoint (collop_bspline_taylor) lives in the consumer, not here.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private

    public :: bspline_workspace_t
    public :: bspline_init
    public :: bspline_set_knots
    public :: bspline_eval_basis
    public :: bspline_eval_deriv
    public :: bspline_span_index
    public :: bspline_eval_jvp   ! d/dx [sum_i c_i B_i(x)] . vx   (x active)
    public :: bspline_eval_vjp   ! (d/dx [sum_i c_i B_i(x)])^T . u (x active)

    ! Caller-owned B-spline state. order, nbreak, ncoef and the augmented knot
    ! vector are filled by bspline_init / bspline_set_knots; no module-level
    ! mutable state exists, so two workspaces never race.
    type :: bspline_workspace_t
        integer               :: order  = 0   ! spline order k (degree + 1)
        integer               :: nbreak = 0   ! number of breakpoints
        integer               :: ncoef  = 0   ! number of basis funcs = nbreak+k-2
        real(dp), allocatable :: knots(:)     ! augmented knot vector, length nbreak+2(k-1)
        logical               :: knots_set = .false.
    end type bspline_workspace_t

contains

    ! Allocate and dimension a workspace for the given order and breakpoint
    ! count. order must be >= 2 (NEO-2 uses 3-5); nbreak must be >= 2. The knot
    ! vector is allocated here but only filled by bspline_set_knots.
    subroutine bspline_init(ws, order, nbreak, status)
        type(bspline_workspace_t), intent(out) :: ws
        integer,                   intent(in)  :: order
        integer,                   intent(in)  :: nbreak
        type(fortnum_status_t),    intent(out) :: status

        integer :: nknots

        call status_set(status, FORTNUM_OK, "")
        if (order < 2) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: order must be >= 2")
            return
        end if
        if (nbreak < 2) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: nbreak must be >= 2")
            return
        end if

        ws%order  = order
        ws%nbreak = nbreak
        ws%ncoef  = nbreak + order - 2
        nknots    = nbreak + 2*(order - 1)
        if (allocated(ws%knots)) deallocate (ws%knots)
        allocate (ws%knots(nknots))
        ws%knots     = 0.0_dp
        ws%knots_set = .false.
    end subroutine bspline_init

    ! Build the augmented (clamped) knot vector from nbreak breakpoints. The
    ! first and last breakpoint are repeated with multiplicity order; the
    ! interior breakpoints appear once. Breakpoints must be strictly increasing.
    subroutine bspline_set_knots(ws, breakpts, status)
        type(bspline_workspace_t), intent(inout) :: ws
        real(dp),                  intent(in)     :: breakpts(:)
        type(fortnum_status_t),    intent(out)    :: status

        integer :: k, m, i, pos

        call status_set(status, FORTNUM_OK, "")
        if (ws%order < 2 .or. ws%nbreak < 2) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: init before set_knots")
            return
        end if
        if (size(breakpts) /= ws%nbreak) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: breakpts size /= nbreak")
            return
        end if
        do i = 2, ws%nbreak
            if (breakpts(i) <= breakpts(i-1)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "bspline: breakpoints must be strictly increasing")
                return
            end if
        end do

        k = ws%order
        m = ws%nbreak
        ! Left clamp: k copies of breakpts(1).
        do i = 1, k
            ws%knots(i) = breakpts(1)
        end do
        ! Interior breakpoints 2..m-1, single multiplicity.
        pos = k
        do i = 2, m - 1
            pos = pos + 1
            ws%knots(pos) = breakpts(i)
        end do
        ! Right clamp: k copies of breakpts(m).
        do i = 1, k
            pos = pos + 1
            ws%knots(pos) = breakpts(m)
        end do
        ws%knots_set = .true.
    end subroutine bspline_set_knots

    ! Index l of the knot span [knots(l), knots(l+1)) that contains x, clamped to
    ! the valid range [k, nknots-k] so the basis is fully supported. Matches the
    ! FindSpan convention (Piegl & Tiller A2.1) for a clamped knot vector. This
    ! lookup is primal_only: the returned index is constant within a span.
    pure function bspline_span_index(ws, x) result(span)
        type(bspline_workspace_t), intent(in) :: ws
        real(dp),                  intent(in) :: x
        integer :: span

        integer :: k, nknots, lo, hi, mid

        k      = ws%order
        nknots = size(ws%knots)
        ! Valid coefficient-supporting spans are k..(nknots-k); the last
        ! breakpoint sits at knots(nknots-k+1), so clamp x at the top span.
        if (x >= ws%knots(nknots - k + 1)) then
            span = nknots - k
            return
        end if
        if (x <= ws%knots(k)) then
            span = k
            return
        end if
        lo = k
        hi = nknots - k + 1
        do
            if (hi - lo <= 1) exit
            mid = (lo + hi)/2
            if (x < ws%knots(mid)) then
                hi = mid
            else
                lo = mid
            end if
        end do
        span = lo
    end function bspline_span_index

    ! Evaluate all ncoef basis functions B_{i,k}(x) at x, returning values(ncoef).
    ! Only the k functions whose support contains the active span are nonzero;
    ! the rest are zero. Cox-de Boor recursion (Piegl & Tiller A2.2).
    subroutine bspline_eval_basis(ws, x, values, status)
        type(bspline_workspace_t), intent(in)  :: ws
        real(dp),                  intent(in)  :: x
        real(dp),                  intent(out) :: values(:)
        type(fortnum_status_t),    intent(out) :: status

        integer  :: k, span, j, r
        real(dp) :: nb(ws%order)
        real(dp) :: left(ws%order), right(ws%order)
        real(dp) :: saved, temp

        call status_set(status, FORTNUM_OK, "")
        if (.not. ws%knots_set) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: knots not set")
            return
        end if
        if (size(values) /= ws%ncoef) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: values size /= ncoef")
            return
        end if

        k    = ws%order
        span = bspline_span_index(ws, x)

        ! Nonzero basis functions on span: B_{span-k+1..span} in 1-based knot
        ! indexing. Compute the k values nb(1..k) for degree k-1.
        nb(1) = 1.0_dp
        do j = 1, k - 1
            left(j)  = x - ws%knots(span + 1 - j)
            right(j) = ws%knots(span + j) - x
            saved = 0.0_dp
            do r = 1, j
                temp  = nb(r)/(right(r) + left(j - r + 1))
                nb(r) = saved + right(r)*temp
                saved = left(j - r + 1)*temp
            end do
            nb(j + 1) = saved
        end do

        values = 0.0_dp
        ! nb(j) is the basis function with global coefficient index span-k+j.
        do j = 1, k
            values(span - k + j) = nb(j)
        end do
    end subroutine bspline_eval_basis

    ! Evaluate basis functions and their derivatives up to order nderiv at x.
    ! dvalues(0:nderiv, ncoef): row 0 is the value, row d is the d-th derivative.
    ! Algorithm A2.3 (Piegl & Tiller): the derivative basis-function table.
    subroutine bspline_eval_deriv(ws, x, nderiv, dvalues, status)
        type(bspline_workspace_t), intent(in)  :: ws
        real(dp),                  intent(in)  :: x
        integer,                   intent(in)  :: nderiv
        real(dp),                  intent(out) :: dvalues(0:, :)
        type(fortnum_status_t),    intent(out) :: status

        integer  :: k, span, j, r, s1, s2, rk, pk, j1, j2, d, col
        real(dp) :: ndu(ws%order, ws%order)
        real(dp) :: a(2, ws%order)
        real(dp) :: ders(0:nderiv, ws%order)
        real(dp) :: left(ws%order), right(ws%order)
        real(dp) :: saved, temp, dd, fac
        integer  :: p

        call status_set(status, FORTNUM_OK, "")
        if (.not. ws%knots_set) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: knots not set")
            return
        end if
        if (nderiv < 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: nderiv must be >= 0")
            return
        end if
        if (size(dvalues, 1) /= nderiv + 1 .or. size(dvalues, 2) /= ws%ncoef) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline: dvalues shape /= (0:nderiv, ncoef)")
            return
        end if

        k    = ws%order
        p    = k - 1            ! polynomial degree
        span = bspline_span_index(ws, x)

        ! ndu(1,1) holds the triangular table of basis values and knot diffs.
        ndu(1, 1) = 1.0_dp
        do j = 1, p
            left(j)  = x - ws%knots(span + 1 - j)
            right(j) = ws%knots(span + j) - x
            saved = 0.0_dp
            do r = 1, j
                ! lower triangle: knot differences
                ndu(j + 1, r) = right(r) + left(j - r + 1)
                temp = ndu(r, j)/ndu(j + 1, r)
                ! upper triangle: basis values
                ndu(r, j + 1) = saved + right(r)*temp
                saved = left(j - r + 1)*temp
            end do
            ndu(j + 1, j + 1) = saved
        end do

        ! Function values (degree p) in the last column.
        do j = 1, p + 1
            ders(0, j) = ndu(j, p + 1)
        end do

        ! Derivatives via the A2.3 recurrence over the columns of ndu.
        do r = 1, p + 1
            s1 = 1
            s2 = 2
            a(1, 1) = 1.0_dp
            do d = 1, nderiv
                dd = 0.0_dp
                rk = r - 1 - d
                pk = p - d
                if (r - 1 >= d) then
                    a(s2, 1) = a(s1, 1)/ndu(pk + 2, rk + 1)
                    dd = a(s2, 1)*ndu(rk + 1, pk + 1)
                end if
                if (rk >= -1) then
                    j1 = 1
                else
                    j1 = -rk
                end if
                if (r - 1 <= pk + 1) then
                    j2 = d - 1
                else
                    j2 = p - (r - 1)
                end if
                do j = j1, j2
                    a(s2, j + 1) = (a(s1, j + 1) - a(s1, j))/ndu(pk + 2, rk + j + 1)
                    dd = dd + a(s2, j + 1)*ndu(rk + j + 1, pk + 1)
                end do
                if (r - 1 <= pk) then
                    a(s2, d + 1) = -a(s1, d)/ndu(pk + 2, r)
                    dd = dd + a(s2, d + 1)*ndu(r, pk + 1)
                end if
                ders(d, r) = dd
                j = s1
                s1 = s2
                s2 = j
            end do
        end do

        ! Multiply by the factorial factor p!/(p-d)!.
        fac = real(p, dp)
        do d = 1, nderiv
            do r = 1, p + 1
                ders(d, r) = ders(d, r)*fac
            end do
            fac = fac*real(p - d, dp)
        end do

        ! Scatter the k nonzero functions into the global coefficient columns.
        dvalues = 0.0_dp
        do d = 0, nderiv
            do j = 1, k
                col = span - k + j
                dvalues(d, col) = ders(d, j)
            end do
        end do
    end subroutine bspline_eval_deriv

    ! JVP of the spline value s(x) = sum_i c_i B_{i,k}(x) w.r.t. x. Policy:
    ! transparent within a fixed span (ad.md §4). Active: x. Inactive: coef,
    ! knots, order. Valid only inside the span x sits in; crossing a knot is a
    ! non-smooth event the caller must guard (bspline_span_index is primal_only).
    !
    ! jv = (d s / d x) * vx = [sum_i c_i B'_{i,k}(x)] * vx
    subroutine bspline_eval_jvp(ws, x, coef, vx, jv, status)
        type(bspline_workspace_t), intent(in)  :: ws
        real(dp),                  intent(in)  :: x
        real(dp),                  intent(in)  :: coef(:)   ! spline coefficients (inactive)
        real(dp),                  intent(in)  :: vx        ! tangent for x
        real(dp),                  intent(out) :: jv
        type(fortnum_status_t),    intent(out) :: status

        real(dp) :: dvals(0:1, ws%ncoef)

        jv = 0.0_dp
        if (size(coef) /= ws%ncoef) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: coef size /= ncoef")
            return
        end if
        call bspline_eval_deriv(ws, x, 1, dvals, status)
        if (status%code /= FORTNUM_OK) return
        jv = dot_product(coef, dvals(1, :))*vx
    end subroutine bspline_eval_jvp

    ! VJP of the scalar spline value s(x) w.r.t. x. Scalar output collapses the
    ! VJP to the gradient times the output cotangent u.
    !
    ! jtu = u * (d s / d x) = u * [sum_i c_i B'_{i,k}(x)]
    subroutine bspline_eval_vjp(ws, x, coef, u, jtu, status)
        type(bspline_workspace_t), intent(in)  :: ws
        real(dp),                  intent(in)  :: x
        real(dp),                  intent(in)  :: coef(:)
        real(dp),                  intent(in)  :: u         ! output cotangent
        real(dp),                  intent(out) :: jtu       ! input cotangent for x
        type(fortnum_status_t),    intent(out) :: status

        real(dp) :: dvals(0:1, ws%ncoef)

        jtu = 0.0_dp
        if (size(coef) /= ws%ncoef) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "bspline: coef size /= ncoef")
            return
        end if
        call bspline_eval_deriv(ws, x, 1, dvals, status)
        if (status%code /= FORTNUM_OK) return
        jtu = u*dot_product(coef, dvals(1, :))
    end subroutine bspline_eval_vjp

end module fortnum_bspline
