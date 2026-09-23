!> Validated integration of ODEs y' = f(y) with an abstract right-hand side.
!>
!> Generalizes the fixed-dimension `lohner7`/`lohner8`/`flow_enclosure`
!> Lohner-QR integrators from `gc-loss-certificate` (Lohner 1988; Moore
!> 1966) to an arbitrary dimension `n` and an abstract right-hand side
!> `ode_rhs_t`, so a single implementation serves any autonomous system
!> instead of one hand-written copy per state dimension. The right-hand
!> side supplies two things:
!>
!>   - `eval_box(x, f, df, ok)`: enclosures of f and its Jacobian on a real
!>     box, used by the a priori (Picard) enclosure, the Gronwall bound, and
!>     the first-order Lohner mean-value step.
!>   - `eval_taylor(k, y, fk)`: the k-th order-by-order Taylor coefficient of
!>     f(y(t)) given the coefficients 0..k-1 already computed, as `idual_t`
!>     so the box Jacobian of the coefficient with respect to the seeded
!>     initial condition comes along for free (`fortnum_idual`,
!>     `fortnum_taylor_series`). Used by the high-order Taylor predictor.
!>
!> A Lohner step represents the flow image of an initial box as a moving
!> parallelepiped, x(t) in x_c + A (C - c) + B e, where C is the fixed
!> initial box, c its centre, A the accumulated linear part, and B e a
!> re-orthonormalized error frame. Each step: (1) computes an a priori
!> (Picard) enclosure of the box over the step, which bounds the mean-value
!> Jacobian; (2) advances x_c by a (possibly high-order) Taylor predictor;
!> (3) forms the mean-value step Jacobian Q = I + h A_apriori (1 + eps) and
!> propagates A, B through it; (4) re-orthonormalizes B (dynamically
!> reordered Gram-Schmidt) to control frame growth and encloses its inverse
!> rigorously (`fortnum_verified_linalg`'s general matrix-ball inverse,
!> applied to the real frame with zero input error) to fold the new error
!> into e without inflating A.
module fortnum_validated_ode
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*), operator(/), operator(**), exp, mid, &
        width, mag, hull, contains_zero
    use fortnum_idual, only: idual_t, idual_var, idual_const, operator(+), &
        operator(-), operator(*), operator(/)
    use fortnum_verified_linalg, only: mball_t, verified_inverse, &
        interval_matmul, interval_matvec
    implicit none
    private

    public :: ode_rhs_t
    public :: picard_apriori, gronwall_bound, lohner_step_jacobian
    public :: eps_ball_matrix
    public :: lohner_qr_frame, lohner_inverse_enclosure
    public :: taylor_lohner_predictor
    public :: lohner_state_t, lohner_state_init, lohner_step, lohner_integrate
    public :: event_fn_if, event_crossing_newton
    public :: section_fn_if, section_crossing
    public :: stop_unresolved, stop_crossed, stop_avoided, stopping_enclosure

    !> Abstract autonomous right-hand side y' = f(y).
    type, abstract :: ode_rhs_t
    contains
        procedure(ode_box_if), deferred :: eval_box
        procedure(ode_taylor_if), deferred :: eval_taylor
    end type ode_rhs_t

    abstract interface
        !> Enclosures of f(x) and its Jacobian df/dx on the box x.
        subroutine ode_box_if(this, x, f, df, ok)
            import :: ode_rhs_t, interval_t
            class(ode_rhs_t), intent(in) :: this
            type(interval_t), intent(in) :: x(:)
            type(interval_t), intent(out) :: f(:), df(:, :)
            logical, intent(out) :: ok
        end subroutine ode_box_if

        !> k-th order Taylor coefficient fk of f(y(t)) given the solution's
        !> own coefficients y(0:k-1, :) (idual_t so a seeded initial
        !> condition's box Jacobian propagates through the composition).
        subroutine ode_taylor_if(this, k, y, fk)
            import :: ode_rhs_t, idual_t
            class(ode_rhs_t), intent(in) :: this
            integer, intent(in) :: k
            type(idual_t), intent(in) :: y(0:, :)
            type(idual_t), intent(out) :: fk(:)
        end subroutine ode_taylor_if
    end interface

    !> Scalar event function along time: enclosures of g(t) and dg/dt on t.
    abstract interface
        subroutine event_fn_if(t, g, dgdt)
            import :: interval_t
            type(interval_t), intent(in) :: t
            type(interval_t), intent(out) :: g, dgdt
        end subroutine event_fn_if
    end interface

    !> Scalar section function g(y) and its gradient dg/dy, enclosed on a
    !> real box y: the crossing condition is g(y) = 0. Generalizes
    !> `lohner_time`'s xi = 0 coordinate section (there g is simply a state
    !> component, dg/dy a fixed unit vector) to an arbitrary differentiable
    !> function of the state.
    abstract interface
        subroutine section_fn_if(y, g, dgdy)
            import :: interval_t
            type(interval_t), intent(in) :: y(:)
            type(interval_t), intent(out) :: g
            type(interval_t), intent(out) :: dgdy(:)
        end subroutine section_fn_if
    end interface

    !> `stopping_enclosure` classification of m_T = max_{0<=t<=T} g(y(t))
    !> for a box of initial conditions: `stop_crossed` when a rigorous
    !> lower bound on m_T is already > 0 (some trajectory certainly hits
    !> g > 0 by or before the decisive step); `stop_avoided` when the run
    !> completes to T with a rigorous upper bound on m_T still < 0 (every
    !> trajectory from the box stays g < 0 throughout); `stop_unresolved`
    !> otherwise (the bounds straddle 0, or an a priori box or width-bound
    !> failure cut the run short before T).
    integer, parameter :: stop_unresolved = 0, stop_crossed = 1, &
        stop_avoided = 2

    !> Lohner parallelepiped state: x(t) in xc + A (dlt) + B e, dlt fixed at
    !> the initial cell minus its centre.
    type :: lohner_state_t
        integer :: n = 0
        real(dp) :: t = 0.0_dp
        real(dp), allocatable :: xc(:), a(:, :), b(:, :)
        type(interval_t), allocatable :: dlt(:), e(:), bx(:)
    end type lohner_state_t

contains

    !> Picard a priori box: an enclosure y of the flow image over [0, h]
    !> starting from `box`, found by inflating the trial box
    !> y1 = box + [0, h] f(y0) until it is contained in y0 (Moore 1966).
    !> f and df are the right-hand side enclosures at the returned y.
    subroutine picard_apriori(rhs, box, h, y, f, df, ok)
        class(ode_rhs_t), intent(in) :: rhs
        type(interval_t), intent(in) :: box(:)
        real(dp), intent(in) :: h
        type(interval_t), intent(out) :: y(:), f(:), df(:, :)
        logical, intent(out) :: ok
        integer :: n, it, i
        type(interval_t) :: th
        type(interval_t) :: y0(size(box)), y1(size(box)), f0(size(box))
        type(interval_t) :: df0(size(box), size(box))
        real(dp) :: w

        n = size(box)
        th = interval(0.0_dp, h)
        call rhs%eval_box(box, f0, df0, ok)
        if (.not. ok) return
        y1 = box + th*f0
        do it = 1, 20
            do i = 1, n
                w = 0.1_dp*(width(y1(i)) - width(box(i))) + 0.05_dp*h*mag(f0(i)) &
                    + 1.0e-13_dp*mag(y1(i)) + 1.0e-300_dp
                y0(i) = interval(y1(i)%lo - w, y1(i)%hi + w)
            end do
            call rhs%eval_box(y0, f0, df0, ok)
            if (.not. ok) return
            y1 = box + th*f0
            if (all(y1%lo >= y0%lo) .and. all(y1%hi <= y0%hi)) then
                y = y1
                call rhs%eval_box(y, f, df, ok)
                return
            end if
        end do
        ok = .false.
    end subroutine picard_apriori

    !> Gronwall variational bound: an interval enclosing an upper bound of
    !> the operator norm ||dPhi_h/dx0|| given a Jacobian enclosure df_box on
    !> an a priori box over the step, from the logarithmic-norm inequality
    !> ||dPhi_h/dx0|| <= exp(L h),  L = max_i sum_j |df_box(i,j)|.
    pure function gronwall_bound(df_box, h) result(bound)
        type(interval_t), intent(in) :: df_box(:, :)
        real(dp), intent(in) :: h
        type(interval_t) :: bound
        real(dp) :: l
        integer :: i

        l = 0.0_dp
        do i = 1, size(df_box, 1)
            l = max(l, sum(mag(df_box(i, :))))
        end do
        bound = exp(interval(0.0_dp, l*h))
    end function gronwall_bound

    !> Isotropic matrix ball M = I + [-eps, eps] enclosing the mean-value
    !> propagator D Phi_t for t in [0, h], eps a rigorous bound on
    !> exp(h ||a||_inf) - 1 by the linear majorant x + x^2 for
    !> x = h ||a||_inf (1 + delta) <= 0.5, following `lohner7`'s
    !> `step_jacobian`. Shared by the first-order step Jacobian and the
    !> order-q variational Jacobian's remainder term.
    pure function eps_ball_matrix(a, h) result(m)
        type(interval_t), intent(in) :: a(:, :)
        real(dp), intent(in) :: h
        type(interval_t) :: m(size(a, 1), size(a, 1))
        real(dp) :: norm, x, eps
        integer :: i, j, n

        n = size(a, 1)
        norm = 0.0_dp
        do i = 1, n
            norm = max(norm, sum(mag(a(i, :))))
        end do
        x = h*norm*(1.0_dp + 1.0e-12_dp)
        if (x > 0.5_dp) then
            eps = huge(1.0_dp)
        else
            eps = (x + x*x)*(1.0_dp + 1.0e-12_dp)
        end if
        do j = 1, n
            do i = 1, n
                m(i, j) = interval(-eps, eps)
            end do
            m(j, j) = interval(1.0_dp - eps, 1.0_dp + eps)
        end do
    end function eps_ball_matrix

    !> Mean-value one-step Jacobian enclosure Q = I + h A (1 + eps) for a
    !> Jacobian enclosure `a` on the a priori box, with eps a rigorous bound
    !> on the local (Lagrange-remainder-free) linearization error, following
    !> `lohner7`'s `step_jacobian`. First-order-in-h fallback; see
    !> `taylor_lohner_predictor`'s `q_jac` output for the order-q enclosure.
    pure function lohner_step_jacobian(a, h) result(q)
        type(interval_t), intent(in) :: a(:, :)
        real(dp), intent(in) :: h
        type(interval_t) :: q(size(a, 1), size(a, 1))
        type(interval_t) :: m(size(a, 1), size(a, 1)), am(size(a, 1), size(a, 1))
        integer :: i, j, n

        n = size(a, 1)
        m = eps_ball_matrix(a, h)
        am = interval_matmul(a, m)
        do j = 1, n
            do i = 1, n
                q(i, j) = h*am(i, j)
            end do
            q(j, j) = q(j, j) + 1.0_dp
        end do
    end function lohner_step_jacobian

    !> Dynamically reordered Gram-Schmidt frame of the real matrix m,
    !> columns weighted by norm2(column) * width(e(column)) so directions
    !> that carry more current uncertainty are orthonormalized first
    !> (Lohner 1988): controls long-term frame growth better than a fixed
    !> column order.
    subroutine lohner_qr_frame(m, e, b)
        real(dp), intent(in) :: m(:, :)
        type(interval_t), intent(in) :: e(:)
        real(dp), intent(out) :: b(size(m, 1), size(m, 1))
        real(dp) :: w(size(m, 1)), c(size(m, 1), size(m, 1)), nrm
        integer :: ord(size(m, 1)), j, k, n

        n = size(m, 1)
        do j = 1, n
            w(j) = norm2(m(:, j))*max(width(e(j)), tiny(1.0_dp))
        end do
        do j = 1, n
            ord(j) = maxloc(w, 1)
            w(ord(j)) = -1.0_dp
        end do
        c = m(:, ord)
        do j = 1, n
            do k = 1, j - 1
                c(:, j) = c(:, j) - dot_product(c(:, k), c(:, j))*c(:, k)
            end do
            nrm = norm2(c(:, j))
            if (nrm < 1.0e-12_dp) then
                c(:, j) = 0.0_dp
                c(mod(j - 1, n) + 1, j) = 1.0_dp
                do k = 1, j - 1
                    c(:, j) = c(:, j) - dot_product(c(:, k), c(:, j))*c(:, k)
                end do
                nrm = norm2(c(:, j))
            end if
            c(:, j) = c(:, j)/nrm
        end do
        b = c
    end subroutine lohner_qr_frame

    !> Rigorous enclosure of b^-1 for a real matrix b, from
    !> `fortnum_verified_linalg`'s general complex matrix-ball inverse
    !> (b embedded with zero imaginary part, zero input error): entry (i, j)
    !> is bounded by center(i,j) +/- rho, valid because |M_ij| <= ||M||_2
    !> for every matrix M in the ball.
    subroutine lohner_inverse_enclosure(b, binv, ok)
        real(dp), intent(in) :: b(:, :)
        type(interval_t), intent(out) :: binv(size(b, 1), size(b, 1))
        logical, intent(out) :: ok
        complex(dp) :: zb(size(b, 1), size(b, 1))
        type(mball_t) :: g
        real(dp) :: re
        integer :: i, j, n

        n = size(b, 1)
        zb = cmplx(b, 0.0_dp, dp)
        call verified_inverse(zb, 0.0_dp, 0.0_dp, g, ok)
        if (.not. ok) return
        do j = 1, n
            do i = 1, n
                re = real(g%c(i, j), dp)
                binv(i, j) = interval(re - g%rho, re + g%rho)
            end do
        end do
    end subroutine lohner_inverse_enclosure

    !> Order-q Taylor predictor of the point flow starting at xc, plus a
    !> rigorous Lagrange-style remainder from the order-(q+1) coefficient
    !> evaluated on the a priori box (`picard_apriori`'s output), and the
    !> point Jacobian d(xpred)/d(xc) carried by the idual_t coefficients'
    !> gradients. Generalizes the order-2 predictor hand-derived in
    !> `lohner7::flow7` (h f(x) + h^2/2 f'(x) f(x)) to any order q >= 1.
    !>
    !> With the optional `ay` (a Jacobian bound of f on `apriori_box`, as
    !> returned by `picard_apriori`), also returns `q_jac`, the order-q
    !> variational Jacobian enclosure of `lohner_taylor`:
    !>   Q = sum_{k=0}^{q} D y_k(Y) h^k + D y_{q+1}(Y) M h^{q+1},
    !> D y_k(Y) the box Jacobian of the k-th Taylor coefficient (the idual_t
    !> gradient carried by the box coefficients `y_bx` below, seeded on
    !> `apriori_box`), and M = `eps_ball_matrix(ay, h)` the Gronwall
    !> enclosure of D Phi_t for t in [0, h]. This is exact to order q in h
    !> (only the order-(q+1) remainder term uses the cruder Gronwall
    !> enclosure), unlike `lohner_step_jacobian`'s Q = I + h A (1 + eps),
    !> which is first-order in h regardless of q.
    !> Optional `ybx_out`: the order-(q+1) box Taylor coefficient values
    !> (the `%v` part only, no idual gradient), i.e. the array
    !> `y_bx(:, :)%v` computed below. Shared with `lohner_step`'s local
    !> step-time polynomial consumers (`section_crossing`,
    !> `stopping_enclosure`): the Horner form from `ybx_out(0:q, :)` plus
    !> the Lagrange remainder from `ybx_out(q+1, :)` is a rigorous
    !> enclosure of the state at any local step-time t in [0, h] for every
    !> point of the a priori box (the `y_bx` recursion starts from the
    !> whole box, not just its centre), following `lohner_time`'s
    !> `local_state` generalized from a point-plus-remainder enclosure to a
    !> full box-valid one.
    subroutine taylor_lohner_predictor(rhs, xc, apriori_box, h, q, xpred, &
        remainder, jac_pred, ok, ay, q_jac, ybx_out)
        class(ode_rhs_t), intent(in) :: rhs
        real(dp), intent(in) :: xc(:)
        type(interval_t), intent(in) :: apriori_box(:)
        real(dp), intent(in) :: h
        integer, intent(in) :: q
        real(dp), intent(out) :: xpred(:), jac_pred(:, :)
        type(interval_t), intent(out) :: remainder(:)
        logical, intent(out) :: ok
        type(interval_t), intent(in), optional :: ay(:, :)
        type(interval_t), intent(out), optional :: q_jac(:, :)
        type(interval_t), intent(out), optional :: ybx_out(0:, :)
        integer :: n, i, j, k
        type(idual_t) :: y_pt(0:q + 1, size(xc)), fk_pt(size(xc))
        type(idual_t) :: y_bx(0:q + 1, size(xc)), fk_bx(size(xc))
        real(dp) :: hp
        type(interval_t) :: jk(size(xc), size(xc)), m(size(xc), size(xc))
        type(interval_t) :: hq1

        n = size(xc)
        ok = .true.
        do i = 1, n
            y_pt(0, i) = idual_var(interval(xc(i)), i, n)
            ! Seeded as an independent variable per box component (not
            ! idual_const) so d(j) carries the box Jacobian D y_k(Y) needed
            ! by the optional `q_jac` output below; the remainder term only
            ! reads %v, unaffected by carrying the gradient too.
            y_bx(0, i) = idual_var(apriori_box(i), i, n)
        end do
        do k = 1, q + 1
            call rhs%eval_taylor(k - 1, y_pt(0:k - 1, :), fk_pt)
            call rhs%eval_taylor(k - 1, y_bx(0:k - 1, :), fk_bx)
            do i = 1, n
                y_pt(k, i) = fk_pt(i)/idual_const(interval(real(k, dp)), n)
                y_bx(k, i) = fk_bx(i)/idual_const(interval(real(k, dp)), 0)
            end do
        end do
        if (present(ybx_out)) then
            do k = 0, q + 1
                ybx_out(k, :) = y_bx(k, :)%v
            end do
        end if

        xpred = xc
        jac_pred = 0.0_dp
        do j = 1, n
            jac_pred(j, j) = 1.0_dp
        end do
        hp = 1.0_dp
        do k = 1, q
            hp = hp*h
            do i = 1, n
                xpred(i) = xpred(i) + hp*mid(y_pt(k, i)%v)
                do j = 1, n
                    jac_pred(i, j) = jac_pred(i, j) + hp*mid(y_pt(k, i)%d(j))
                end do
            end do
        end do
        do i = 1, n
            remainder(i) = interval(h**(q + 1))*y_bx(q + 1, i)%v
        end do

        if (present(q_jac)) then
            do j = 1, n
                do i = 1, n
                    q_jac(i, j) = y_bx(0, i)%d(j)
                end do
            end do
            hp = 1.0_dp
            do k = 1, q
                hp = hp*h
                do j = 1, n
                    do i = 1, n
                        q_jac(i, j) = q_jac(i, j) + interval(hp)*y_bx(k, i)%d(j)
                    end do
                end do
            end do
            if (present(ay)) then
                do j = 1, n
                    do i = 1, n
                        jk(i, j) = y_bx(q + 1, i)%d(j)
                    end do
                end do
                m = eps_ball_matrix(ay, h)
                jk = interval_matmul(jk, m)
                hq1 = interval(h**(q + 1))
                do j = 1, n
                    do i = 1, n
                        q_jac(i, j) = q_jac(i, j) + hq1*jk(i, j)
                    end do
                end do
            end if
        end if
    end subroutine taylor_lohner_predictor

    !> Initialize a Lohner integration from an initial box cell: centre at
    !> the box midpoint, identity linear part, identity frame, zero error.
    subroutine lohner_state_init(state, cell)
        type(lohner_state_t), intent(out) :: state
        type(interval_t), intent(in) :: cell(:)
        integer :: n, i

        n = size(cell)
        state%n = n
        state%t = 0.0_dp
        allocate (state%xc(n), state%a(n, n), state%b(n, n))
        allocate (state%dlt(n), state%e(n), state%bx(n))
        do i = 1, n
            state%xc(i) = mid(cell(i))
        end do
        state%a = 0.0_dp
        state%b = 0.0_dp
        do i = 1, n
            state%a(i, i) = 1.0_dp
            state%b(i, i) = 1.0_dp
        end do
        state%dlt = cell - interval(state%xc)
        state%e = interval(0.0_dp)
        state%bx = cell
    end subroutine lohner_state_init

    !> One Lohner-QR step of size h and predictor order q, following
    !> `lohner7::flow7`'s per-step update generalized to n dimensions and
    !> order q: a priori box, order-q centre prediction with remainder,
    !> mean-value Jacobian propagation, frame re-orthonormalization, and
    !> rigorous inverse-enclosed error-frame update.
    !> Optional `apriori_y`/`apriori_f`: this step's a priori box (`y` in
    !> `picard_apriori`, valid for every t in the step, not just its
    !> endpoints) and the right-hand side there. Optional `ybx_out`:
    !> pass-through of `taylor_lohner_predictor`'s local step-time box
    !> Taylor coefficients (see its header). Consumed by `section_crossing`
    !> and `stopping_enclosure` to build a priori-box section bounds and
    !> locate crossings within the step without a second a priori/Taylor
    !> computation.
    subroutine lohner_step(state, rhs, h, q, ok, variational, apriori_y, &
        apriori_f, ybx_out)
        type(lohner_state_t), intent(inout) :: state
        class(ode_rhs_t), intent(in) :: rhs
        real(dp), intent(in) :: h
        integer, intent(in) :: q
        logical, intent(out) :: ok
        logical, intent(in), optional :: variational
        type(interval_t), intent(out), optional :: apriori_y(:), apriori_f(:)
        type(interval_t), intent(out), optional :: ybx_out(0:, :)
        integer :: n
        type(interval_t) :: y(state%n), fy(state%n), ay(state%n, state%n)
        type(interval_t) :: xb(state%n), qmat(state%n, state%n)
        type(interval_t) :: qa(state%n, state%n), qb(state%n, state%n)
        type(interval_t) :: v(state%n), mbmat(state%n, state%n)
        type(interval_t) :: binv(state%n, state%n), lin(state%n)
        type(interval_t) :: remainder(state%n)
        real(dp) :: xpred(state%n), jac_pred(state%n, state%n)
        real(dp) :: am_new(state%n, state%n), bm_new(state%n, state%n)
        real(dp) :: x_new(state%n)
        logical :: use_var

        n = state%n
        use_var = .true.
        if (present(variational)) use_var = variational
        call picard_apriori(rhs, state%bx, h, y, fy, ay, ok)
        if (.not. ok) return
        if (present(apriori_y)) apriori_y = y
        if (present(apriori_f)) apriori_f = fy
        if (use_var) then
            call taylor_lohner_predictor(rhs, state%xc, y, h, q, xpred, &
                remainder, jac_pred, ok, ay=ay, q_jac=qmat, &
                ybx_out=ybx_out)
        else
            call taylor_lohner_predictor(rhs, state%xc, y, h, q, xpred, &
                remainder, jac_pred, ok, ybx_out=ybx_out)
        end if
        if (.not. ok) return
        xb = interval(xpred) + remainder

        if (.not. use_var) qmat = lohner_step_jacobian(ay, h)
        qa = interval_matmul(qmat, interval(state%a))
        qb = interval_matmul(qmat, interval(state%b))
        am_new = mid(qa)
        call lohner_qr_frame(mid(qb), state%e, bm_new)
        call lohner_inverse_enclosure(bm_new, binv, ok)
        if (.not. ok) return

        x_new = mid(xb)
        v = (xb - interval(x_new)) &
            + interval_matvec(qa - interval(am_new), state%dlt)
        mbmat = interval_matmul(binv, qb)
        state%e = interval_matvec(binv, v) + interval_matvec(mbmat, state%e)
        lin = interval(x_new) + interval_matvec(interval(am_new), state%dlt) &
            + interval_matvec(interval(bm_new), state%e)
        state%bx = hull(min_hull(y, lin), interval(x_new))
        state%xc = x_new
        state%a = am_new
        state%b = bm_new
        state%t = state%t + h
    end subroutine lohner_step

    !> Tighter-of-two-enclosures box, elementwise: valid because both `y`
    !> (the a priori bound) and `lin` (the propagated affine bound) contain
    !> the true flow image, so their intersection does too, whenever they
    !> overlap (they always do here since both contain the same new centre).
    pure elemental function min_hull(y, lin) result(r)
        type(interval_t), intent(in) :: y, lin
        type(interval_t) :: r

        r = interval(max(y%lo, lin%lo), min(y%hi, lin%hi))
    end function min_hull

    !> Integrate from an initial box `cell` over [0, tend] with step at most
    !> hmax and predictor order q, returning a rigorous enclosure of the
    !> image of every point of cell at time tend.
    subroutine lohner_integrate(rhs, cell, tend, hmax, q, out, ok, variational)
        class(ode_rhs_t), intent(in) :: rhs
        type(interval_t), intent(in) :: cell(:)
        real(dp), intent(in) :: tend, hmax
        integer, intent(in) :: q
        type(interval_t), intent(out) :: out(size(cell))
        logical, intent(out) :: ok
        logical, intent(in), optional :: variational
        type(lohner_state_t) :: state
        integer :: nstep, i
        real(dp) :: h

        call lohner_state_init(state, cell)
        nstep = max(1, ceiling(tend/hmax))
        h = tend/real(nstep, dp)
        ok = .true.
        do i = 1, nstep
            call lohner_step(state, rhs, h, q, ok, variational)
            if (.not. ok) return
        end do
        out = interval(state%xc) + interval_matvec(interval(state%a), state%dlt) &
            + interval_matvec(interval(state%b), state%e)
    end subroutine lohner_integrate

    !> Interval-Newton bracket refinement of a root of a scalar event
    !> function g on [t0, t1], given a sign change at the endpoints. Each
    !> iteration replaces the bracket t by
    !>   N(t) = mid(t) - g(mid(t))/g'(t)  intersected with t,
    !> which contains every root of g in t (interval Newton, Moore 1966);
    !> falls back to bisection when g'(t) contains zero. Returns ok = .false.
    !> if no sign change is found at the initial endpoints.
    subroutine event_crossing_newton(event, t0, t1, troot, ok)
        procedure(event_fn_if) :: event
        real(dp), intent(in) :: t0, t1
        type(interval_t), intent(out) :: troot
        logical, intent(out) :: ok
        type(interval_t) :: t, g0, g1, dgdt, gm, cand
        real(dp) :: tm
        integer :: it

        call event(interval(t0), g0, dgdt)
        call event(interval(t1), g1, dgdt)
        ok = (g0%lo <= 0.0_dp .and. g1%hi >= 0.0_dp) .or. &
            (g0%hi >= 0.0_dp .and. g1%lo <= 0.0_dp)
        if (.not. ok) return

        t = interval(t0, t1)
        do it = 1, 80
            tm = mid(t)
            call event(interval(tm), gm, dgdt)
            call event(t, g0, dgdt)
            if (dgdt%lo <= 0.0_dp .and. dgdt%hi >= 0.0_dp) then
                if (gm%lo <= 0.0_dp .and. gm%hi >= 0.0_dp) then
                    troot = interval(tm)
                    ok = .true.
                    return
                end if
                if (gm%lo > 0.0_dp .eqv. g0%lo > 0.0_dp) then
                    t = interval(tm, t%hi)
                else
                    t = interval(t%lo, tm)
                end if
                cycle
            end if
            cand = interval(tm) - gm/dgdt
            cand = interval(max(cand%lo, t%lo), min(cand%hi, t%hi))
            if (cand%lo > cand%hi) then
                ok = .false.
                return
            end if
            if (width(cand) < 1.0e-13_dp*max(1.0_dp, abs(tm))) then
                troot = cand
                ok = .true.
                return
            end if
            t = cand
        end do
        troot = t
        ok = .true.
    end subroutine event_crossing_newton

    !> Rigorous dg/dt = dg/dy . f(y) by the chain rule, given a box y, the
    !> right-hand side value f on y (as returned by `picard_apriori` or
    !> `lohner_step`'s `apriori_f`), and the section's own gradient dg/dy.
    !> Valid on any box, not just the a priori one: called on both the whole
    !> step's a priori box (transversality) and on narrowing local-time
    !> brackets (`locate_section_crossing`'s Newton derivative).
    subroutine section_rate(rhs, section, y, dgdt, ok)
        class(ode_rhs_t), intent(in) :: rhs
        procedure(section_fn_if) :: section
        type(interval_t), intent(in) :: y(:)
        type(interval_t), intent(out) :: dgdt
        logical, intent(out) :: ok
        type(interval_t) :: g, dgdy(size(y)), f(size(y)), df(size(y), size(y))
        integer :: i

        call rhs%eval_box(y, f, df, ok)
        if (.not. ok) return
        call section(y, g, dgdy)
        dgdt = dgdy(1)*f(1)
        do i = 2, size(y)
            dgdt = dgdt + dgdy(i)*f(i)
        end do
    end subroutine section_rate

    !> Rigorous state enclosure at local step-time t in [0, h], valid for
    !> every point of the a priori box the step's `ybx_out` was seeded from
    !> (not just the step's centre trajectory): Horner evaluation of the
    !> order-0..q box Taylor coefficients plus the Lagrange remainder from
    !> the box coefficient at order q + 1. Unlike `taylor_lohner_predictor`'s
    !> `xpred` (which uses the tight point coefficients `y_pt` for orders
    !> 0..q and the box coefficient only for the remainder, to predict the
    !> centre trajectory as precisely as possible), this generalizes
    !> `lohner_time`'s `local_state` from a single-trajectory enclosure to a
    !> genuine whole-box one, needed to certify a crossing for every point
    !> of a `section_crossing` cell that is not itself a near-point box.
    pure function section_local_state(ybx, q, t) result(s)
        type(interval_t), intent(in) :: ybx(0:, :)
        integer, intent(in) :: q
        type(interval_t), intent(in) :: t
        type(interval_t) :: s(size(ybx, 2))
        integer :: k

        s = ybx(q, :)
        do k = q - 1, 0, -1
            s = ybx(k, :) + t*s
        end do
        s = s + ybx(q + 1, :)*t**(q + 1)
    end function section_local_state

    !> Interval-Newton bracket refinement of the root t of g(state(t)) = 0
    !> on [0, h], given a sign change and transversality already certified
    !> by the caller over the whole step (`section_crossing`). Each
    !> iteration narrows t by the interval Newton map t <- mid(t) -
    !> g(state(mid(t))) / (dg/dt on state(t)), following `lohner_time`'s
    !> `locate_crossing` generalized from a coordinate section to an
    !> abstract `section_fn_if` and the derivative recomputed by the chain
    !> rule (`section_rate`) instead of a fixed state component. Fails
    !> (`ok = .false.`) if the derivative bracket ever contains zero or the
    !> bracket becomes empty, exactly as `lohner_time` does, rather than
    !> falling back to bisection: transversality was already certified over
    !> the whole step, so this should not happen away from numerical noise.
    subroutine locate_section_crossing(rhs, section, ybx, q, h, tau, &
        ztau, ok)
        class(ode_rhs_t), intent(in) :: rhs
        procedure(section_fn_if) :: section
        type(interval_t), intent(in) :: ybx(0:, :)
        integer, intent(in) :: q
        real(dp), intent(in) :: h
        type(interval_t), intent(out) :: tau, ztau(size(ybx, 2))
        logical, intent(out) :: ok
        type(interval_t) :: t, tnew, sm(size(ybx, 2)), st(size(ybx, 2))
        type(interval_t) :: gm, dgdy(size(ybx, 2)), dgdt
        real(dp) :: m
        integer :: it
        logical :: ok2

        t = interval(0.0_dp, h)
        ok = .true.
        do it = 1, 80
            m = mid(t)
            sm = section_local_state(ybx, q, interval(m, m))
            st = section_local_state(ybx, q, t)
            call section_rate(rhs, section, st, dgdt, ok2)
            if (.not. ok2 .or. contains_zero(dgdt)) then
                ok = .false.
                return
            end if
            call section(sm, gm, dgdy)
            tnew = interval(m, m) - gm/dgdt
            tnew = interval(max(tnew%lo, t%lo), min(tnew%hi, t%hi))
            if (tnew%lo > tnew%hi) then
                ok = .false.
                return
            end if
            if (width(tnew) >= width(t) .and. width(t) < 1.0e-15_dp*h) then
                t = tnew
                exit
            end if
            t = tnew
        end do
        tau = t
        ztau = section_local_state(ybx, q, tau)
    end subroutine locate_section_crossing

    !> Validated first-return map to a section g(y) = 0 with a prescribed
    !> crossing direction, generalizing `lohner_time::bounce_return_v`'s
    !> adaptive-step/state-machine/interval-Newton crossing localization
    !> from a fixed 8D complex-time state with a coordinate section (xi = 0)
    !> to an abstract `ode_rhs_t` of any dimension and an abstract scalar
    !> `section_fn_if`. From an initial box `cell0`, integrates forward with
    !> adaptive step (shrinking by half on a priori-box failure and retrying
    !> from the same running state, since `lohner_step` only mutates it on
    !> success; growing back by 1.2x towards h0 on success -- following
    !> `bounce_return_v`). `direction > 0` looks for g increasing through
    !> zero (dg/dt > 0 required); `direction < 0` for g decreasing through
    !> zero (dg/dt < 0). At each step the rigorous node enclosure
    !> `state%bx` gives g there (`section`); the FIRST step whose node
    !> value g_end is certainly on the prescribed "after" side (g_end%lo > 0
    !> for direction > 0, g_end%hi < 0 for direction < 0) while the
    !> preceding node g_prev was not yet certainly on that side (g_prev%lo
    !> <= 0, respectively g_prev%hi >= 0 -- true whether g_prev was
    !> certainly "before" or still straddling zero) is the one checked for
    !> transversality, over the WHOLE step's a priori box (`section_rate` on
    !> `lohner_step`'s `apriori_y`/`apriori_f`, a rigorous superset of the
    !> trajectory for every t in the step): only if dg/dt there is bounded
    !> away from 0 in the prescribed direction is the crossing located, by
    !> `locate_section_crossing` on that step's local Taylor polynomial
    !> (valid over the step regardless of exactly where within it the sign
    !> actually changes). No earlier crossing is missed, because every
    !> prior step's node pair failed this test, and every accepted step's a
    !> priori box is a valid enclosure of the trajectory throughout the
    !> step.
    !>
    !> Returns `tau` (enclosure of the crossing time from t = 0) and `ztau`
    !> (enclosure of the state there): `locate_section_crossing` builds its
    !> local polynomial from `ybx_out`, the box Taylor coefficients seeded
    !> on the whole a priori box for the detected step, not just the step's
    !> centre trajectory, so (tau, ztau) is a rigorous enclosure of the
    !> crossing event certified at that step. This differs from
    !> `lohner_time`'s `local_state`, which mixed a tight centre polynomial
    !> with only the top-order box coefficient (sufficient there because
    !> `bounce_return_v`'s `cell0` is effectively a single physical orbit
    !> carried with interval bookkeeping, not a set of visibly distinct
    !> initial conditions).
    !>
    !> Caution for a genuinely wide `cell0`: the node test fires at the
    !> FIRST step whose end is certainly past the section for the WHOLE
    !> running box, so if different points of `cell0` cross several steps
    !> apart, a point that crossed earlier is not guaranteed to lie within
    !> the returned (tau, ztau) -- it can cross before `tau` starts. The
    !> guarantee is exact for a point cell0 and degrades gracefully as
    !> width grows; callers needing every point's individual crossing
    !> bounded should keep `cell0` narrow enough (relative to h0 and dg/dt)
    !> that the whole box crosses within one step, or should split `cell0`
    !> and call `section_crossing` per sub-box. `nsteps`/`nshrink`
    !> optionally report step and halving counts for diagnostics.
    subroutine section_crossing(rhs, section, cell0, h0, q, direction, &
        nmax, tau, ztau, ok, nsteps, nshrink)
        class(ode_rhs_t), intent(in) :: rhs
        procedure(section_fn_if) :: section
        type(interval_t), intent(in) :: cell0(:)
        real(dp), intent(in) :: h0
        integer, intent(in) :: q, direction, nmax
        type(interval_t), intent(out) :: tau, ztau(size(cell0))
        logical, intent(out) :: ok
        integer, intent(out), optional :: nsteps, nshrink
        type(lohner_state_t) :: state
        type(interval_t) :: y(size(cell0)), fy(size(cell0))
        type(interval_t) :: ybx(0:q + 1, size(cell0))
        type(interval_t) :: g_prev, g_end, dgdy_dummy(size(cell0)), dgdt
        real(dp) :: h, hused, tcur
        integer :: j, nshrink_
        logical :: stepok, crossing

        call lohner_state_init(state, cell0)
        call section(state%bx, g_prev, dgdy_dummy)
        h = h0
        tcur = 0.0_dp
        nshrink_ = 0
        ok = .true.
        do j = 1, nmax
            do
                call lohner_step(state, rhs, h, q, stepok, &
                    apriori_y=y, apriori_f=fy, ybx_out=ybx)
                if (stepok) exit
                h = 0.5_dp*h
                nshrink_ = nshrink_ + 1
                if (h < h0*1.0e-8_dp) then
                    ok = .false.
                    if (present(nsteps)) nsteps = j - 1
                    if (present(nshrink)) nshrink = nshrink_
                    return
                end if
            end do
            hused = h
            h = min(h0, h*1.2_dp)
            call section(state%bx, g_end, dgdy_dummy)
            if (direction > 0) then
                crossing = g_prev%lo <= 0.0_dp .and. g_end%lo > 0.0_dp
            else
                crossing = g_prev%hi >= 0.0_dp .and. g_end%hi < 0.0_dp
            end if
            if (crossing) then
                call section_rate(rhs, section, y, dgdt, ok)
                if (ok .and. ((direction > 0 .and. dgdt%lo > 0.0_dp) .or. &
                    (direction < 0 .and. dgdt%hi < 0.0_dp))) then
                    call locate_section_crossing(rhs, section, ybx, q, &
                        hused, tau, ztau, ok)
                    if (present(nsteps)) nsteps = j
                    if (present(nshrink)) nshrink = nshrink_
                    if (.not. ok) return
                    tau = interval(tcur, tcur) + tau
                    return
                end if
            end if
            g_prev = g_end
            tcur = tcur + hused
        end do
        ok = .false.
        if (present(nsteps)) nsteps = nmax
        if (present(nshrink)) nshrink = nshrink_
    end subroutine section_crossing

    !> Enclosure of m_T = max_{0<=t<=T} g(y(t)) for the flow image of a box
    !> of initial conditions `cell`, generalizing `flow_enclosure`'s
    !> `enclose_cell` from a fixed 4D guiding-centre state and a
    !> `gc_system_t`-specific wall function to an abstract `ode_rhs_t` and
    !> `section_fn_if`. Steps at fixed size `h` (following `enclose_cell`,
    !> no adaptivity): each step's a priori box (`lohner_step`'s
    !> `apriori_y`, a rigorous superset of the trajectory for every t in
    !> the step) gives an upper bound on g there via `section`, accumulated
    !> into the running upper bound of m_T; the step-end refined box
    !> `state%bx` (a tighter enclosure, valid exactly at the grid time)
    !> gives a lower bound. `status` becomes `stop_crossed` as soon as the
    !> accumulated lower bound exceeds 0 (returned immediately, without
    !> completing to `tmax`); otherwise, if the run completes to `tmax`
    !> with the accumulated upper bound still below 0, `stop_avoided`;
    !> otherwise `stop_unresolved` (includes the case where `lohner_step`
    !> fails before `tmax`, `ok = .false.`, still returning the partial
    !> bounds and step count accumulated so far).
    subroutine stopping_enclosure(rhs, section, cell, h, tmax, q, status, &
        m, ok, nsteps)
        class(ode_rhs_t), intent(in) :: rhs
        procedure(section_fn_if) :: section
        type(interval_t), intent(in) :: cell(:)
        real(dp), intent(in) :: h, tmax
        integer, intent(in) :: q
        integer, intent(out) :: status
        type(interval_t), intent(out) :: m
        logical, intent(out) :: ok
        integer, intent(out), optional :: nsteps
        type(lohner_state_t) :: state
        type(interval_t) :: y(size(cell))
        type(interval_t) :: g, dgdy(size(cell))
        real(dp) :: t, hs
        integer :: j, nmax
        logical :: stepok

        call lohner_state_init(state, cell)
        call section(state%bx, g, dgdy)
        m = interval(g%lo, -huge(1.0_dp))
        status = stop_unresolved
        ok = .true.
        if (present(nsteps)) nsteps = 0
        if (g%lo > 0.0_dp) then
            status = stop_crossed
            return
        end if
        t = 0.0_dp
        nmax = ceiling(tmax/h) + 1
        do j = 1, nmax
            hs = min(h, tmax - t)
            if (hs <= 0.0_dp) exit
            call lohner_step(state, rhs, hs, q, stepok, apriori_y=y)
            if (.not. stepok) then
                ok = .false.
                if (present(nsteps)) nsteps = j - 1
                return
            end if
            call section(y, g, dgdy)
            m%hi = max(m%hi, g%hi)
            call section(state%bx, g, dgdy)
            m%lo = max(m%lo, g%lo)
            t = t + hs
            if (present(nsteps)) nsteps = j
            if (m%lo > 0.0_dp) then
                status = stop_crossed
                return
            end if
        end do
        if (m%hi < 0.0_dp) status = stop_avoided
    end subroutine stopping_enclosure

end module fortnum_validated_ode
