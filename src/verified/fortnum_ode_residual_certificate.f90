!> A posteriori trajectory certificate for an ODE candidate y' = f(t, y).
!>
!> The candidate is binary64 node data: times t_k, values y_k and slopes
!> f_k. On each cell [t_k, t_{k+1}] of exact length h these define the exact
!> cubic Hermite reconstruction p(s), s in [0, 1], with p(0) = y_k,
!> p(1) = y_{k+1}, p'(0) = h f_k, p'(1) = h f_{k+1}. The reconstruction is
!> continuous with continuous first derivative across nodes, so it is
!> absolutely continuous and its defect
!>
!>   r(s) = p'(s)/h - f(t_k + h s, p(s))
!>
!> drives the error e = y - p through e' = f(y) - f(p) - r. The checker
!> never uses an exact solution and never trusts the candidate integrator.
!>
!>   - `hermite_cubic_enclosure` encloses the power coefficients of p in
!>     outward interval arithmetic from the binary64 node data.
!>   - `linear_residual_bound`: for f = A y + b(s) with a caller-supplied
!>     interval action of A and optional cubic forcing b, r is a cubic
!>     vector polynomial. Its interval Bernstein control vectors R_j give
!>     sup_s ||r(s)||_2 <= max_j ||R_j||_2 (convex-hull property) for the
!>     whole closed cell, not a sample. The bound is O(h^3) per cell for
!>     smooth solutions.
!>   - `nonlinear_residual_bound`: for general f, the cell is split into
!>     `nsub` equal parts; on each part the Bernstein hull of p and p' gives
!>     boxes X and D, the caller encloses f(T, X), and R = D/h - f(T, X).
!>     This bound is only first order, O(L h |y'| / nsub) with L a Lipschitz
!>     constant, because the box evaluation loses the cancellation between
!>     p'/h and f(p). It is valid but loose.
!>   - `residual_radius_step` accumulates the Euclidean error radius with a
!>     caller-supplied stability model: an upper bound mu on the logarithmic
!>     2-norm (one-sided Lipschitz constant) of f on a convex set containing
!>     the exact and reconstructed trajectories. mu <= 0 (unitary,
!>     skew-adjoint or contractive generators) uses amplification 1:
!>     e_{k+1} = e_k + h R_k. mu > 0 uses the Gronwall/Duhamel bound
!>     e_{k+1} = exp(mu h) e_k + R_k (exp(mu h) - 1)/mu. Both bounds are
!>     nondecreasing in the time inside the cell, so e_{k+1} bounds the error
!>     throughout the closed cell, not only at its right endpoint.
!>   - `certify_linear_trajectory` and `certify_nonlinear_trajectory` run
!>     the three steps over a stored trajectory with caller-owned scratch.
!>
!> Complex states use the real/imaginary split (re(:), im(:)) of length 2n;
!> its Euclidean norm is the complex Hilbert norm. The per-cell routines do
!> not allocate and keep no module state, so they are reentrant. Arbitrary
!> finite candidate data are accepted: a corrupted trace yields a large
!> radius rather than a small false one. Nonfinite data, invalid meshes,
!> callback failures and nonfinite bounds fail closed.
module fortnum_ode_residual_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_positive_inf
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*), operator(/), exp, mag, hull
    use fortnum_rounding, only: add_up, mul_up, sqrt_up
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, status_set
    implicit none
    private

    public :: residual_action_if, residual_rhs_box_if
    public :: hermite_cubic_enclosure, linear_residual_bound
    public :: nonlinear_residual_bound, residual_radius_step
    public :: certify_linear_trajectory, certify_nonlinear_trajectory

    abstract interface
        !> Enclose A x for every x in the box x (and every admissible A).
        subroutine residual_action_if(x, ax, ok, ctx)
            import :: interval_t
            type(interval_t), intent(in) :: x(:)
            type(interval_t), intent(out) :: ax(:)
            logical, intent(out) :: ok
            class(*), intent(in), optional :: ctx
        end subroutine residual_action_if

        !> Enclose f(t, x) for every t in the box t and x in the box x.
        subroutine residual_rhs_box_if(t, x, fx, ok, ctx)
            import :: interval_t
            type(interval_t), intent(in) :: t, x(:)
            type(interval_t), intent(out) :: fx(:)
            logical, intent(out) :: ok
            class(*), intent(in), optional :: ctx
        end subroutine residual_rhs_box_if
    end interface

contains

    !> Power coefficients c(:, 0:3) of the exact cubic Hermite p(s) with
    !> p(0) = y0, p(1) = y1, p'(0) = h f0, p'(1) = h f1, for every exact cell
    !> length in the enclosure h.
    pure subroutine hermite_cubic_enclosure(y0, y1, f0, f1, h, c)
        real(dp), intent(in) :: y0(:), y1(:), f0(:), f1(:)
        type(interval_t), intent(in) :: h
        type(interval_t), intent(out) :: c(:, 0:)
        type(interval_t) :: a, b, fa, fb, dy
        integer :: i

        do i = 1, size(y0)
            a = interval(y0(i))
            b = interval(y1(i))
            fa = interval(f0(i))
            fb = interval(f1(i))
            dy = b - a
            c(i, 0) = a
            c(i, 1) = h*fa
            c(i, 2) = 3.0_dp*dy - h*(2.0_dp*fa + fb)
            c(i, 3) = h*(fa + fb) - 2.0_dp*dy
        end do
    end subroutine hermite_cubic_enclosure

    !> Bound sup_{s in [0,1]} ||p'(s)/h - A p(s) - b(s)||_2 over the whole
    !> cell for the Hermite enclosure c and the linear action `action`.
    !> `forcing(:, 0:3)`, if present, encloses the power coefficients of a
    !> cubic forcing b(s) (a constant offset uses column 0 only). The action
    !> is called four times, on the coefficient vectors c(:, k); w(:, 0:3)
    !> is caller scratch. On failure bound = +inf and ok = .false.
    subroutine linear_residual_bound(action, c, h, w, bound, ok, ctx, forcing)
        procedure(residual_action_if) :: action
        type(interval_t), intent(in) :: c(:, 0:), h
        type(interval_t), intent(inout) :: w(:, 0:)
        real(dp), intent(out) :: bound
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx
        type(interval_t), intent(in), optional :: forcing(:, 0:)
        type(interval_t) :: ih, r(0:3), b(0:3)
        real(dp) :: sq(0:3), m
        integer :: i, k
        logical :: action_ok

        ok = .false.
        bound = ieee_value(0.0_dp, ieee_positive_inf)
        if (.not. (h%lo > 0.0_dp)) return
        ih = 1.0_dp/h
        do k = 0, 3
            call action(c(:, k), w(:, k), action_ok, ctx)
            if (.not. action_ok) return
        end do
        sq = 0.0_dp
        do i = 1, size(c, 1)
            r(0) = c(i, 1)*ih - w(i, 0)
            r(1) = 2.0_dp*c(i, 2)*ih - w(i, 1)
            r(2) = 3.0_dp*c(i, 3)*ih - w(i, 2)
            r(3) = interval(0.0_dp) - w(i, 3)
            if (present(forcing)) then
                do k = 0, 3
                    r(k) = r(k) - forcing(i, k)
                end do
            end if
            call power_to_bernstein3(r, b)
            do k = 0, 3
                m = mag(b(k))
                sq(k) = add_up(sq(k), mul_up(m, m))
            end do
        end do
        bound = sqrt_up(maxval(sq))
        ok = ieee_is_finite(bound)
        if (.not. ok) bound = ieee_value(0.0_dp, ieee_positive_inf)
    end subroutine linear_residual_bound

    !> Bound sup_s ||p'(s)/h - f(t0 + h s, p(s))||_2 for general f by
    !> splitting [0, 1] into nsub equal parts. On each part the Bernstein
    !> hulls of p and p' give boxes X and D and `rhs_box` encloses f on
    !> (T, X). First order and loose: see the module header.
    !> w(:, 0:2) is caller scratch.
    subroutine nonlinear_residual_bound(rhs_box, t0, c, h, nsub, w, bound, &
                                        ok, ctx)
        procedure(residual_rhs_box_if) :: rhs_box
        real(dp), intent(in) :: t0
        type(interval_t), intent(in) :: c(:, 0:), h
        integer, intent(in) :: nsub
        type(interval_t), intent(inout) :: w(:, 0:)
        real(dp), intent(out) :: bound
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx
        type(interval_t) :: ih, a, d, e, tbox, q(0:3), b(0:3), dq(0:2)
        real(dp) :: sq, m, worst
        integer :: i, j
        logical :: rhs_ok

        ok = .false.
        bound = ieee_value(0.0_dp, ieee_positive_inf)
        if (nsub < 1) return
        if (.not. (h%lo > 0.0_dp)) return
        ih = 1.0_dp/h
        d = interval(1.0_dp)/real(nsub, dp)
        worst = 0.0_dp
        do j = 1, nsub
            a = interval(real(j - 1, dp))/real(nsub, dp)
            e = a + d
            tbox = interval(t0) + h*interval(a%lo, e%hi)
            do i = 1, size(c, 1)
                call shifted_cubic(c(i, :), a, d, q, dq)
                call power_to_bernstein3(q, b)
                w(i, 0) = hull(hull(b(0), b(1)), hull(b(2), b(3)))
                w(i, 2) = hull(hull(dq(0), dq(0) + 0.5_dp*dq(1)), &
                               dq(0) + dq(1) + dq(2))
            end do
            call rhs_box(tbox, w(:, 0), w(:, 1), rhs_ok, ctx)
            if (.not. rhs_ok) return
            sq = 0.0_dp
            do i = 1, size(c, 1)
                m = mag(w(i, 2)*ih - w(i, 1))
                sq = add_up(sq, mul_up(m, m))
            end do
            worst = max(worst, sq)
        end do
        bound = sqrt_up(worst)
        ok = ieee_is_finite(bound)
        if (.not. ok) bound = ieee_value(0.0_dp, ieee_positive_inf)
    end subroutine nonlinear_residual_bound

    !> Error radius after one cell: e1 bounds ||y - p|| throughout the closed
    !> cell given e0 at its left node, residual bound rbound and an interval
    !> mu whose upper endpoint bounds the logarithmic 2-norm of f. mu%hi <= 0
    !> uses amplification 1; mu%hi > 0 the Gronwall/Duhamel factor
    !> exp(mu h) with phi(h) = (exp(mu h) - 1)/mu <= h min(exp(mu h),
    !> outward (exp(z) - 1)/z), z = mu h, both nondecreasing in z.
    pure subroutine residual_radius_step(mu, e0, h, rbound, e1, ok)
        type(interval_t), intent(in) :: mu, h
        real(dp), intent(in) :: e0, rbound
        real(dp), intent(out) :: e1
        logical, intent(out) :: ok
        type(interval_t) :: z, g, phi
        real(dp) :: growth, factor

        ok = .false.
        e1 = ieee_value(0.0_dp, ieee_positive_inf)
        if (.not. ieee_is_finite(mu%hi)) return
        if (.not. ieee_is_finite(h%hi) .or. .not. (h%lo > 0.0_dp)) return
        if (.not. ieee_is_finite(e0) .or. .not. ieee_is_finite(rbound)) return
        if (e0 < 0.0_dp .or. rbound < 0.0_dp) return
        if (mu%hi <= 0.0_dp) then
            e1 = add_up(e0, mul_up(h%hi, rbound))
        else
            z = interval(mul_up(mu%hi, h%hi))
            g = exp(z)
            growth = g%hi
            phi = (g - 1.0_dp)/z
            factor = growth
            if (ieee_is_finite(phi%hi)) factor = min(factor, phi%hi)
            e1 = add_up(mul_up(growth, e0), mul_up(mul_up(h%hi, factor), rbound))
        end if
        ok = ieee_is_finite(e1)
        if (.not. ok) e1 = ieee_value(0.0_dp, ieee_positive_inf)
    end subroutine residual_radius_step

    !> Certify a stored trace for f = A y. t(0:m) are strictly increasing
    !> exact binary64 nodes, y(:, 0:m) values, f(:, 0:m) slopes (normally
    !> the candidate's A y_k, but any finite data are admitted),
    !> initial_radius >= ||y(t_0) - y_0||_2. On success radius(0) =
    !> initial_radius and radius(k) bounds ||y(t) - p(t)||_2 for all t in
    !> [t_{k-1}, t_k]; residual(k) is the cell residual bound. c and w are
    !> caller scratch of shape (n, 0:3). Failure leaves +inf radii.
    subroutine certify_linear_trajectory(action, t, y, f, mu, initial_radius, &
                                         c, w, residual, radius, status, ctx)
        procedure(residual_action_if) :: action
        real(dp), intent(in) :: t(0:), y(:, 0:), f(:, 0:), initial_radius
        type(interval_t), intent(in) :: mu
        type(interval_t), intent(inout) :: c(:, 0:), w(:, 0:)
        real(dp), intent(out) :: residual(:), radius(0:)
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx
        type(interval_t) :: h
        integer :: k
        logical :: ok

        call validate_trace(t, y, f, mu, initial_radius, c, w, 4, residual, &
                            radius, status)
        if (status%code /= FORTNUM_OK) return
        radius(0) = initial_radius
        do k = 1, size(t) - 1
            h = interval(t(k)) - interval(t(k - 1))
            call hermite_cubic_enclosure(y(:, k - 1), y(:, k), f(:, k - 1), &
                                         f(:, k), h, c)
            call linear_residual_bound(action, c, h, w, residual(k), ok, ctx)
            if (.not. ok) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "linear residual bound failed or is not finite")
                return
            end if
            call residual_radius_step(mu, radius(k - 1), h, residual(k), &
                                      radius(k), ok)
            if (.not. ok) then
                radius(k) = ieee_value(0.0_dp, ieee_positive_inf)
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "error radius is not finite")
                return
            end if
        end do
    end subroutine certify_linear_trajectory

    !> As `certify_linear_trajectory` for a general right-hand side enclosed
    !> by `rhs_box`, with nsub >= 1 subcells per cell and scratch
    !> c(n, 0:3), w(n, 0:2). The residual bound is first order (module
    !> header).
    subroutine certify_nonlinear_trajectory(rhs_box, t, y, f, nsub, mu, &
                                            initial_radius, c, w, residual, &
                                            radius, status, ctx)
        procedure(residual_rhs_box_if) :: rhs_box
        real(dp), intent(in) :: t(0:), y(:, 0:), f(:, 0:), initial_radius
        integer, intent(in) :: nsub
        type(interval_t), intent(in) :: mu
        type(interval_t), intent(inout) :: c(:, 0:), w(:, 0:)
        real(dp), intent(out) :: residual(:), radius(0:)
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx
        type(interval_t) :: h
        integer :: k
        logical :: ok

        call validate_trace(t, y, f, mu, initial_radius, c, w, 3, residual, &
                            radius, status)
        if (status%code /= FORTNUM_OK) return
        if (nsub < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "nsub must be >= 1")
            return
        end if
        radius(0) = initial_radius
        do k = 1, size(t) - 1
            h = interval(t(k)) - interval(t(k - 1))
            call hermite_cubic_enclosure(y(:, k - 1), y(:, k), f(:, k - 1), &
                                         f(:, k), h, c)
            call nonlinear_residual_bound(rhs_box, t(k - 1), c, h, nsub, w, &
                                          residual(k), ok, ctx)
            if (.not. ok) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "nonlinear residual bound failed or is not finite")
                return
            end if
            call residual_radius_step(mu, radius(k - 1), h, residual(k), &
                                      radius(k), ok)
            if (.not. ok) then
                radius(k) = ieee_value(0.0_dp, ieee_positive_inf)
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "error radius is not finite")
                return
            end if
        end do
    end subroutine certify_nonlinear_trajectory

    subroutine validate_trace(t, y, f, mu, initial_radius, c, w, wcols, &
                              residual, radius, status)
        real(dp), intent(in) :: t(0:), y(:, 0:), f(:, 0:), initial_radius
        type(interval_t), intent(in) :: mu, c(:, 0:), w(:, 0:)
        integer, intent(in) :: wcols
        real(dp), intent(out) :: residual(:), radius(0:)
        type(fortnum_status_t), intent(out) :: status
        integer :: k, m, n

        residual = ieee_value(0.0_dp, ieee_positive_inf)
        radius = ieee_value(0.0_dp, ieee_positive_inf)
        m = size(t) - 1
        n = size(y, 1)
        status%code = FORTNUM_OK
        status%msg = ""
        if (m < 1 .or. n < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "empty trace")
        else if (size(y, 2) /= m + 1 .or. size(f, 1) /= n &
                 .or. size(f, 2) /= m + 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "trace shape mismatch")
        else if (size(c, 1) /= n .or. size(c, 2) < 4 .or. size(w, 1) /= n &
                 .or. size(w, 2) < wcols) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "scratch shape mismatch")
        else if (size(residual) /= m .or. size(radius) /= m + 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "output shape mismatch")
        else if (.not. all(ieee_is_finite(t))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "nonfinite time node")
        else if (.not. all(ieee_is_finite(y)) .or. &
                 .not. all(ieee_is_finite(f))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "nonfinite candidate")
        else if (.not. ieee_is_finite(initial_radius)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "invalid initial radius")
        else if (initial_radius < 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "invalid initial radius")
        else if (.not. ieee_is_finite(mu%hi)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, "nonfinite stability mu")
        end if
        if (status%code /= FORTNUM_OK) return
        do k = 1, m
            if (.not. (t(k) > t(k - 1))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "time nodes are not strictly increasing")
                return
            end if
        end do
    end subroutine validate_trace

    !> Coefficients in u of p(a + d u) (q) and p'(a + d u) (dq) for every
    !> a in the box a and d in the box d.
    pure subroutine shifted_cubic(c, a, d, q, dq)
        type(interval_t), intent(in) :: c(0:), a, d
        type(interval_t), intent(out) :: q(0:3), dq(0:2)
        type(interval_t) :: p0, p1, p2

        p0 = ((c(3)*a + c(2))*a + c(1))*a + c(0)
        p1 = (3.0_dp*c(3)*a + 2.0_dp*c(2))*a + c(1)
        p2 = 3.0_dp*c(3)*a + c(2)
        q(0) = p0
        q(1) = p1*d
        q(2) = p2*(d*d)
        q(3) = c(3)*(d*d*d)
        dq(0) = p1
        dq(1) = 2.0_dp*p2*d
        dq(2) = 3.0_dp*c(3)*(d*d)
    end subroutine shifted_cubic

    !> Bernstein coefficients of the cubic with power coefficients r.
    pure subroutine power_to_bernstein3(r, b)
        type(interval_t), intent(in) :: r(0:3)
        type(interval_t), intent(out) :: b(0:3)

        b(0) = r(0)
        b(1) = r(0) + r(1)/3.0_dp
        b(2) = r(0) + (2.0_dp*r(1) + r(2))/3.0_dp
        b(3) = r(0) + r(1) + r(2) + r(3)
    end subroutine power_to_bernstein3

end module fortnum_ode_residual_certificate
