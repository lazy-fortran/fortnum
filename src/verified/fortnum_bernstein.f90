!> Rigorous evaluation of polynomials given in Bernstein form,
!>   p(s) = sum_{j=0}^n coef(j) B_{j,n}(s),  B_{j,n}(s) = C(n,j) s^j (1-s)^(n-j),
!> on real and complex boxes, by local Taylor re-expansion at a box centre
!> (Corliss and Rihm 1996; the load-time grid cache and complex evaluation
!> generalize `cheb_field` from `gc-loss-certificate`, which is the module's
!> only current user).
!>
!> The Bernstein basis is well conditioned on [0, 1] (every B_{j,n} lies in
!> [0, 1] and they sum to 1), unlike a centred monomial expansion, whose
!> coefficients blow up away from the fit centre and wrap catastrophically in
!> interval or complex-dual arithmetic. Two independent rigorous bounds are
!> provided:
!>
!>   - `bernstein_hull`: for a real sub-box [a, b] subset [0, 1], the exact
!>     convex-hull property of Bezier curves bounds the curve (and its
!>     derivative) by the control points of two de Casteljau subdivisions,
!>     computed in interval arithmetic so the bound is outward rounded.
!>   - `bernstein_taylor_coeffs`: the local Taylor coefficients of p at a
!>     real centre c, from the forward-difference identity
!>       a_k = C(n,k) sum_{j=0}^{n-k} (Delta^k coef)_j B_{j,n-k}(c),
!>     rigorous because the forward differences are exact (real coefficients)
!>     and the basis values `bernstein_basis_table` returns are interval
!>     enclosures. Composing a_k with a (real or complex) box radius around c
!>     gives a Taylor-model enclosure of p on that box; `bernstein_grid_t`
!>     caches a_k at a fixed grid of centres so repeated evaluation (an
!>     ODE right-hand side sampled at many boxes) reuses the O(n^2)
!>     coefficient sum instead of recomputing it per box.
!>
!> `bernstein_lsq_fit` is an ordinary (non-rigorous) least-squares helper
!> that returns Bernstein-form coefficients from sampled data, for building
!> the polynomial models the rigorous routines above then enclose.
module fortnum_bernstein
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, cinterval_t, interval, cinterval, &
        operator(+), operator(-), operator(*), operator(/), mid
    use fortnum_linalg, only: dense_solve
    implicit none
    private

    public :: bernstein_basis_real, bernstein_diff_table, bernstein_basis_table
    public :: bernstein_taylor_coeffs, bernstein_hull, bernstein_lsq_fit
    public :: bernstein_eval_real_box, bernstein_eval_complex_box
    public :: bernstein_grid_t, bernstein_grid_build, bernstein_grid_index
    public :: bernstein_grid_eval_real, bernstein_grid_eval_complex

    !> Cache of local Taylor coefficients of one degree-n Bernstein
    !> polynomial at `ngrid` equispaced centres in [0, 1], built once and
    !> reused by every box evaluation (`bernstein_grid_eval_real/complex`).
    !> Any centre in [0, 1] gives an exact re-expansion of the same
    !> polynomial identity, so the grid trades enclosure width (distance
    !> from a box centre to the nearest cached grid point) for load-time
    !> cost; it never affects correctness.
    type :: bernstein_grid_t
        integer :: n = 0, ngrid = 0
        real(dp), allocatable :: centre(:)          ! (ngrid)
        type(interval_t), allocatable :: a(:, :)    ! (0:n, ngrid)
    end type bernstein_grid_t

contains

    !> Bernstein basis values and s-derivatives at a real point s, from
    !> shared powers of s and (1-s).
    pure subroutine bernstein_basis_real(deg, s, t, dt)
        integer, intent(in) :: deg
        real(dp), intent(in) :: s
        real(dp), intent(out) :: t(0:deg), dt(0:deg)
        real(dp) :: sp(0:deg), op(0:deg), bino, term1, term2
        integer :: j

        sp(0) = 1.0_dp
        op(0) = 1.0_dp
        do j = 1, deg
            sp(j) = sp(j - 1)*s
            op(j) = op(j - 1)*(1.0_dp - s)
        end do
        bino = 1.0_dp
        do j = 0, deg
            if (j > 0) bino = bino*real(deg - j + 1, dp)/real(j, dp)
            t(j) = bino*sp(j)*op(deg - j)
            if (j == 0) then
                term1 = 0.0_dp
            else
                term1 = real(j, dp)*sp(j - 1)*op(deg - j)
            end if
            if (j == deg) then
                term2 = 0.0_dp
            else
                term2 = real(deg - j, dp)*sp(j)*op(deg - j - 1)
            end if
            dt(j) = bino*(term1 - term2)
        end do
    end subroutine bernstein_basis_real

    !> Forward-difference triangle fd(k, j) = Delta^k coef(j), valid for
    !> j = 0..n-k. Exact (finitely many real subtractions, no enclosure
    !> needed on this side).
    pure subroutine bernstein_diff_table(n, coef, fd)
        integer, intent(in) :: n
        real(dp), intent(in) :: coef(0:n)
        real(dp), intent(out) :: fd(0:n, 0:n)
        integer :: k, j

        fd = 0.0_dp
        fd(0, :) = coef
        do k = 1, n
            do j = 0, n - k
                fd(k, j) = fd(k - 1, j + 1) - fd(k - 1, j)
            end do
        end do
    end subroutine bernstein_diff_table

    !> Bernstein basis values B_{j,m}(c) at a single real point c, for every
    !> degree m = 0..deg, from the triangular recursion
    !> B_{j,m} = (1-c) B_{j,m-1} + c B_{j-1,m-1}, in interval arithmetic (c a
    !> point interval): stable because the weights c, 1-c both lie in [0, 1]
    !> whenever c does.
    pure subroutine bernstein_basis_table(deg, c, btab)
        integer, intent(in) :: deg
        real(dp), intent(in) :: c
        type(interval_t), intent(out) :: btab(0:deg, 0:deg)
        type(interval_t) :: cc, omc
        integer :: m, j

        cc = interval(c)
        omc = interval(1.0_dp) - cc
        btab(0, 0) = interval(1.0_dp)
        do m = 1, deg
            btab(0, m) = omc*btab(0, m - 1)
            do j = 1, m - 1
                btab(j, m) = omc*btab(j, m - 1) + cc*btab(j - 1, m - 1)
            end do
            btab(m, m) = cc*btab(m - 1, m - 1)
        end do
    end subroutine bernstein_basis_table

    !> Local Taylor coefficients a_k = p^{(k)}(c)/k! of the degree-n
    !> Bernstein polynomial at the point c of `bernstein_basis_table`'s
    !> btab, rigorously enclosed since btab is an interval and fd is exact.
    pure subroutine bernstein_taylor_coeffs(n, fd, btab, a)
        integer, intent(in) :: n
        real(dp), intent(in) :: fd(0:n, 0:n)
        type(interval_t), intent(in) :: btab(0:n, 0:n)
        type(interval_t), intent(out) :: a(0:n)
        real(dp) :: bino
        integer :: k, j, nk

        bino = 1.0_dp
        do k = 0, n
            if (k > 0) bino = bino*real(n - k + 1, dp)/real(k, dp)
            nk = n - k
            a(k) = interval(0.0_dp)
            do j = 0, nk
                a(k) = a(k) + btab(j, nk)*fd(k, j)
            end do
            a(k) = a(k)*bino
        end do
    end subroutine bernstein_taylor_coeffs

    !> De Casteljau subdivision of a degree-n Bezier curve (interval control
    !> points p) at interval parameter t: control points of the two pieces
    !> on the sub-ranges t divides [0, 1] into, computed in interval
    !> arithmetic so an inexact t (from a real division) still yields a
    !> rigorous enclosure of the exact split.
    pure subroutine bezier_split_i(n, p, t, left, right)
        integer, intent(in) :: n
        type(interval_t), intent(in) :: p(0:n), t
        type(interval_t), intent(out) :: left(0:n), right(0:n)
        type(interval_t) :: q(0:n, 0:n), one_minus_t
        integer :: i, k

        one_minus_t = interval(1.0_dp) - t
        q(0:n, 0) = p
        do k = 1, n
            do i = 0, n - k
                q(i, k) = one_minus_t*q(i, k - 1) + t*q(i + 1, k - 1)
            end do
        end do
        do i = 0, n
            left(i) = q(0, i)
            right(i) = q(i, n - i)
        end do
    end subroutine bezier_split_i

    !> Control-point enclosure of the degree-n Bezier curve "coef" restricted
    !> to [a, b] subset [0, 1] and reparametrized to [0, 1]: two exact-in-
    !> exact-arithmetic de Casteljau subdivisions, both carried out in
    !> interval arithmetic. Their convex hull (`bernstein_hull`) bounds the
    !> curve on [a, b] and narrows properly as the box shrinks; it never
    !> inflates from cancellation between near-equal large Bernstein terms,
    !> at any width, unlike summing the shared basis directly over an
    !> s-interval would.
    pure subroutine bezier_restrict_i(n, coef, a, b, q)
        integer, intent(in) :: n
        real(dp), intent(in) :: coef(0:n), a, b
        type(interval_t), intent(out) :: q(0:n)
        type(interval_t) :: p(0:n), left1(0:n), right1(0:n)
        type(interval_t) :: left2(0:n), right2(0:n), t2
        integer :: i

        do i = 0, n
            p(i) = interval(coef(i))
        end do
        call bezier_split_i(n, p, interval(a), left1, right1)
        if (a >= 1.0_dp) then
            q = right1
            return
        end if
        t2 = (interval(b) - interval(a))/(interval(1.0_dp) - interval(a))
        call bezier_split_i(n, right1, t2, left2, right2)
        q = left2
    end subroutine bezier_restrict_i

    !> Rigorous value and s-derivative bound of the degree-n Bernstein sum
    !> "coef" for s in [a, b] subset [0, 1]: convex hull of the restricted
    !> control points, and of the restricted curve's own degree-(n-1)
    !> derivative control points n(q_{i+1}-q_i)/(b-a). A degenerate box
    !> (a == b) returns a zero-width derivative bound; callers needing a
    !> point derivative there should evaluate the analytic derivative
    !> directly instead.
    pure subroutine bernstein_hull(n, coef, a, b, vbound, dbound)
        integer, intent(in) :: n
        real(dp), intent(in) :: coef(0:n), a, b
        type(interval_t), intent(out) :: vbound, dbound
        type(interval_t) :: q(0:n), d(0:n - 1), width
        real(dp) :: vlo, vhi, dlo, dhi
        integer :: i

        call bezier_restrict_i(n, coef, a, b, q)
        vlo = q(0)%lo
        vhi = q(0)%hi
        do i = 1, n
            vlo = min(vlo, q(i)%lo)
            vhi = max(vhi, q(i)%hi)
        end do
        vbound = interval(vlo, vhi)

        width = interval(b) - interval(a)
        if (b <= a) then
            dbound = interval(0.0_dp)
            return
        end if
        do i = 0, n - 1
            d(i) = interval(real(n, dp))*(q(i + 1) - q(i))/width
        end do
        dlo = d(0)%lo
        dhi = d(0)%hi
        do i = 1, n - 1
            dlo = min(dlo, d(i)%lo)
            dhi = max(dhi, d(i)%hi)
        end do
        dbound = interval(dlo, dhi)
    end subroutine bernstein_hull

    !> Value of the degree-n Bernstein sum on a real box, by local Taylor
    !> re-expansion at the real centre c: a(k) are the rigorous local Taylor
    !> coefficients, box_s is an enclosure of s, and the Horner sum in
    !> (box_s - c) is evaluated in interval arithmetic.
    pure subroutine bernstein_eval_real_box(n, coef, c, box_s, value)
        integer, intent(in) :: n
        real(dp), intent(in) :: coef(0:n), c
        type(interval_t), intent(in) :: box_s
        type(interval_t), intent(out) :: value
        real(dp) :: fd(0:n, 0:n)
        type(interval_t) :: btab(0:n, 0:n), a(0:n), d
        integer :: k

        call bernstein_diff_table(n, coef, fd)
        call bernstein_basis_table(n, c, btab)
        call bernstein_taylor_coeffs(n, fd, btab, a)
        d = box_s - interval(c)
        value = a(n)
        do k = n - 1, 0, -1
            value = value*d + a(k)
        end do
    end subroutine bernstein_eval_real_box

    !> Value of the degree-n Bernstein sum on a complex box, by the same
    !> local Taylor re-expansion (the coefficients a(k) are real intervals;
    !> the polynomial extends analytically to complex s).
    pure subroutine bernstein_eval_complex_box(n, coef, c, box_z, value)
        integer, intent(in) :: n
        real(dp), intent(in) :: coef(0:n), c
        type(cinterval_t), intent(in) :: box_z
        type(cinterval_t), intent(out) :: value
        real(dp) :: fd(0:n, 0:n)
        type(interval_t) :: btab(0:n, 0:n), a(0:n)
        type(cinterval_t) :: d
        integer :: k

        call bernstein_diff_table(n, coef, fd)
        call bernstein_basis_table(n, c, btab)
        call bernstein_taylor_coeffs(n, fd, btab, a)
        d = box_z - cinterval(c, 0.0_dp)
        value = cinterval(a(n), interval(0.0_dp))
        do k = n - 1, 0, -1
            value = value*d + cinterval(a(k), interval(0.0_dp))
        end do
    end subroutine bernstein_eval_complex_box

    !> Build a grid of ngrid equispaced centres in [0, 1] and the local
    !> Taylor coefficients of "coef" at each one, once.
    subroutine bernstein_grid_build(n, coef, ngrid, grid)
        integer, intent(in) :: n, ngrid
        real(dp), intent(in) :: coef(0:n)
        type(bernstein_grid_t), intent(out) :: grid
        real(dp) :: fd(0:n, 0:n)
        type(interval_t) :: btab(0:n, 0:n)
        integer :: g

        grid%n = n
        grid%ngrid = ngrid
        allocate (grid%centre(ngrid), grid%a(0:n, ngrid))
        call bernstein_diff_table(n, coef, fd)
        do g = 1, ngrid
            grid%centre(g) = real(g - 1, dp)/real(max(ngrid - 1, 1), dp)
            call bernstein_basis_table(n, grid%centre(g), btab)
            call bernstein_taylor_coeffs(n, fd, btab, grid%a(:, g))
        end do
    end subroutine bernstein_grid_build

    !> Nearest grid index to a real centre c in [0, 1]. Any index in
    !> [1, ngrid] gives a valid re-expansion; the nearest one only minimizes
    !> the enclosure width.
    pure function bernstein_grid_index(grid, c) result(ig)
        type(bernstein_grid_t), intent(in) :: grid
        real(dp), intent(in) :: c
        integer :: ig

        ig = nint(c*real(grid%ngrid - 1, dp)) + 1
        ig = max(1, min(grid%ngrid, ig))
    end function bernstein_grid_index

    !> Evaluate on a real box using the cached coefficients nearest its
    !> midpoint, without recomputing the O(n^2) coefficient sum.
    pure subroutine bernstein_grid_eval_real(grid, box_s, value)
        type(bernstein_grid_t), intent(in) :: grid
        type(interval_t), intent(in) :: box_s
        type(interval_t), intent(out) :: value
        integer :: ig, k
        type(interval_t) :: d

        ig = bernstein_grid_index(grid, mid(box_s))
        d = box_s - interval(grid%centre(ig))
        value = grid%a(grid%n, ig)
        do k = grid%n - 1, 0, -1
            value = value*d + grid%a(k, ig)
        end do
    end subroutine bernstein_grid_eval_real

    !> Evaluate on a complex box using the cached coefficients nearest the
    !> real part of its midpoint.
    pure subroutine bernstein_grid_eval_complex(grid, box_z, value)
        type(bernstein_grid_t), intent(in) :: grid
        type(cinterval_t), intent(in) :: box_z
        type(cinterval_t), intent(out) :: value
        integer :: ig, k
        type(cinterval_t) :: d

        ig = bernstein_grid_index(grid, mid(box_z%re))
        d = box_z - cinterval(grid%centre(ig), 0.0_dp)
        value = cinterval(grid%a(grid%n, ig), interval(0.0_dp))
        do k = grid%n - 1, 0, -1
            value = value*d + cinterval(grid%a(k, ig), interval(0.0_dp))
        end do
    end subroutine bernstein_grid_eval_complex

    !> Least-squares Bernstein-form fit: coefficients coef(0:deg) minimizing
    !> sum_i (p(s_i) - y(i))^2 for p(s) = sum_j coef(j) B_{j,deg}(s), by the
    !> normal equations of the Bernstein design matrix. Not itself a
    !> rigorous routine (an ordinary floating-point least-squares solve);
    !> the fitted coefficients are the input the routines above enclose.
    !> rms is the residual root-mean-square at the sample points.
    subroutine bernstein_lsq_fit(deg, s, y, coef, rms, info)
        integer, intent(in) :: deg
        real(dp), intent(in) :: s(:), y(:)
        real(dp), intent(out) :: coef(0:deg)
        real(dp), intent(out) :: rms
        integer, intent(out) :: info
        real(dp), allocatable :: design(:, :), ata(:, :), aty(:), x(:)
        real(dp) :: t(0:deg), dt(0:deg), resid
        integer :: m, i, j

        m = size(s)
        info = -1
        if (m < deg + 1 .or. size(y) /= m) return
        allocate (design(m, 0:deg), ata(0:deg, 0:deg), aty(0:deg), x(0:deg))
        do i = 1, m
            call bernstein_basis_real(deg, s(i), t, dt)
            design(i, :) = t
        end do
        do j = 0, deg
            do i = 0, deg
                ata(i, j) = dot_product(design(:, i), design(:, j))
            end do
            aty(j) = dot_product(design(:, j), y)
        end do
        call dense_solve(ata, aty, x, info)
        if (info /= 0) return
        coef = x
        rms = 0.0_dp
        do i = 1, m
            call bernstein_basis_real(deg, s(i), t, dt)
            resid = dot_product(t, coef) - y(i)
            rms = rms + resid*resid
        end do
        rms = sqrt(rms/real(m, dp))
    end subroutine bernstein_lsq_fit

end module fortnum_bernstein
