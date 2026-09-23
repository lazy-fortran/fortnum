!> Fixture right-hand sides for test_fortnum_validated_ode: type-bound
!> procedures must be module procedures in Fortran, so the ode_rhs_t
!> extensions live here and the program below just uses them.
module test_fortnum_validated_ode_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*), operator(/), sin, cos
    use fortnum_idual, only: idual_t, operator(+), operator(-), operator(*), &
        operator(/), sin, cos
    use fortnum_validated_ode, only: ode_rhs_t
    implicit none

    !> Rotation: y1' = -y2, y2' = y1; exact flow is rotation by t.
    type, extends(ode_rhs_t) :: rotation_rhs_t
    contains
        procedure :: eval_box => rotation_box
        procedure :: eval_taylor => rotation_taylor
    end type rotation_rhs_t

    !> Shear: y1' = y2, y2' = 0; exact flow is y1 + t y2, y2.
    type, extends(ode_rhs_t) :: shear_rhs_t
    contains
        procedure :: eval_box => shear_box
        procedure :: eval_taylor => shear_taylor
    end type shear_rhs_t

    !> y' = y; exact flow is exp(t) y0.
    type, extends(ode_rhs_t) :: exp_rhs_t
    contains
        procedure :: eval_box => exp_box
        procedure :: eval_taylor => exp_taylor
    end type exp_rhs_t

    !> Logistic: y' = y - y^2 (nonlinear, scalar). Exact flow
    !> Phi_t(y0) = y0 e^t / (1 - y0 + y0 e^t), a smooth strictly convex map
    !> on (0, 1) so the variational Jacobian dPhi_t/dy0 is not constant --
    !> an oracle a linear system (rotation, shear) cannot exercise.
    type, extends(ode_rhs_t) :: logistic_rhs_t
    contains
        procedure :: eval_box => logistic_box
        procedure :: eval_taylor => logistic_taylor
    end type logistic_rhs_t

    !> Nonlinear pendulum: y1' = y2, y2' = -sin(y1) (y1 the angle, y2 the
    !> angular velocity). Oracle for `section_crossing`: released from rest
    !> at angle theta0 (a turning point), energy conservation
    !> 0.5 y2^2 - cos(y1) = -cos(theta0) gives the exact velocity
    !> -sqrt(2(1 - cos theta0)) at the first return to y1 = 0, and the
    !> substitution sin(y1/2) = sin(theta0/2) sin(psi) gives the exact
    !> return time as the complete elliptic integral of the first kind
    !> K(sin(theta0/2)) (`elliptic_k_agm`, independent AGM evaluation).
    type, extends(ode_rhs_t) :: pendulum_rhs_t
    contains
        procedure :: eval_box => pendulum_box
        procedure :: eval_taylor => pendulum_taylor
    end type pendulum_rhs_t

contains

    subroutine rotation_box(this, x, f, df, ok)
        class(rotation_rhs_t), intent(in) :: this
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: f(:), df(:, :)
        logical, intent(out) :: ok
        f(1) = -x(2)
        f(2) = x(1)
        df(1, 1) = interval(0.0_dp)
        df(1, 2) = interval(-1.0_dp)
        df(2, 1) = interval(1.0_dp)
        df(2, 2) = interval(0.0_dp)
        ok = .true.
    end subroutine rotation_box

    subroutine rotation_taylor(this, k, y, fk)
        class(rotation_rhs_t), intent(in) :: this
        integer, intent(in) :: k
        type(idual_t), intent(in) :: y(0:, :)
        type(idual_t), intent(out) :: fk(:)
        fk(1) = -y(k, 2)
        fk(2) = y(k, 1)
    end subroutine rotation_taylor

    pure function exact_rotation(x0, t) result(x)
        real(qp), intent(in) :: x0(2), t
        real(qp) :: x(2)
        x(1) = cos(t)*x0(1) - sin(t)*x0(2)
        x(2) = sin(t)*x0(1) + cos(t)*x0(2)
    end function exact_rotation

    subroutine shear_box(this, x, f, df, ok)
        class(shear_rhs_t), intent(in) :: this
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: f(:), df(:, :)
        logical, intent(out) :: ok
        f(1) = x(2)
        f(2) = interval(0.0_dp)
        df(1, 1) = interval(0.0_dp)
        df(1, 2) = interval(1.0_dp)
        df(2, 1) = interval(0.0_dp)
        df(2, 2) = interval(0.0_dp)
        ok = .true.
    end subroutine shear_box

    subroutine shear_taylor(this, k, y, fk)
        class(shear_rhs_t), intent(in) :: this
        integer, intent(in) :: k
        type(idual_t), intent(in) :: y(0:, :)
        type(idual_t), intent(out) :: fk(:)
        fk(1) = y(k, 2)
        fk(2) = 0.0_dp*y(k, 1)
    end subroutine shear_taylor

    pure function exact_shear(x0, t) result(x)
        real(qp), intent(in) :: x0(2), t
        real(qp) :: x(2)
        x(1) = x0(1) + t*x0(2)
        x(2) = x0(2)
    end function exact_shear

    subroutine exp_box(this, x, f, df, ok)
        class(exp_rhs_t), intent(in) :: this
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: f(:), df(:, :)
        logical, intent(out) :: ok
        f(1) = x(1)
        df(1, 1) = interval(1.0_dp)
        ok = .true.
    end subroutine exp_box

    subroutine exp_taylor(this, k, y, fk)
        class(exp_rhs_t), intent(in) :: this
        integer, intent(in) :: k
        type(idual_t), intent(in) :: y(0:, :)
        type(idual_t), intent(out) :: fk(:)
        fk(1) = y(k, 1)
    end subroutine exp_taylor

    subroutine logistic_box(this, x, f, df, ok)
        class(logistic_rhs_t), intent(in) :: this
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: f(:), df(:, :)
        logical, intent(out) :: ok
        f(1) = x(1) - x(1)*x(1)
        df(1, 1) = interval(1.0_dp) - interval(2.0_dp)*x(1)
        ok = .true.
    end subroutine logistic_box

    !> f_k = a_k - sum_{i=0}^{k} a_i a_{k-i}, the order-k Taylor coefficient
    !> of y(t) - y(t)^2 from the Cauchy product of y with itself.
    subroutine logistic_taylor(this, k, y, fk)
        class(logistic_rhs_t), intent(in) :: this
        integer, intent(in) :: k
        type(idual_t), intent(in) :: y(0:, :)
        type(idual_t), intent(out) :: fk(:)
        type(idual_t) :: conv
        integer :: i
        conv = y(0, 1)*y(k, 1)
        do i = 1, k
            conv = conv + y(i, 1)*y(k - i, 1)
        end do
        fk(1) = y(k, 1) - conv
    end subroutine logistic_taylor

    pure function exact_logistic(x0, t) result(x)
        real(qp), intent(in) :: x0, t
        real(qp) :: x, et
        et = exp(t)
        x = x0*et/(1.0_qp - x0 + x0*et)
    end function exact_logistic

    !> dPhi_t/dy0 = e^t / (1 - y0 + y0 e^t)^2, by direct differentiation of
    !> exact_logistic.
    pure function exact_logistic_jac(x0, t) result(dj)
        real(qp), intent(in) :: x0, t
        real(qp) :: dj, et, den
        et = exp(t)
        den = 1.0_qp - x0 + x0*et
        dj = et/(den*den)
    end function exact_logistic_jac

    subroutine pendulum_box(this, x, f, df, ok)
        class(pendulum_rhs_t), intent(in) :: this
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: f(:), df(:, :)
        logical, intent(out) :: ok
        f(1) = x(2)
        f(2) = -sin(x(1))
        df(1, 1) = interval(0.0_dp)
        df(1, 2) = interval(1.0_dp)
        df(2, 1) = -cos(x(1))
        df(2, 2) = interval(0.0_dp)
        ok = .true.
    end subroutine pendulum_box

    !> Order-k Taylor coefficient of y2' = -sin(y1(t)) by the standard
    !> Faa-di-Bruno recursion for sin/cos of a Taylor series: with
    !> s = sin(y1(t)), c = cos(y1(t)), s_0 = sin(y1_0), c_0 = cos(y1_0), and
    !> for m >= 1, m s_m = sum_{j=0}^{m-1} (m-j) y1_(m-j) c_j,
    !> m c_m = -sum_{j=0}^{m-1} (m-j) y1_(m-j) s_j. Recomputed from scratch
    !> up to order k on every call (k stays small in practice), needing no
    !> persistent state across the interface's per-order calls.
    subroutine pendulum_taylor(this, k, y, fk)
        class(pendulum_rhs_t), intent(in) :: this
        integer, intent(in) :: k
        type(idual_t), intent(in) :: y(0:, :)
        type(idual_t), intent(out) :: fk(:)
        type(idual_t) :: s(0:k), c(0:k), acc_s, acc_c
        integer :: m, j
        s(0) = sin(y(0, 1))
        c(0) = cos(y(0, 1))
        do m = 1, k
            acc_s = real(m, dp)*y(m, 1)*c(0)
            acc_c = real(m, dp)*y(m, 1)*s(0)
            do j = 1, m - 1
                acc_s = acc_s + real(m - j, dp)*y(m - j, 1)*c(j)
                acc_c = acc_c + real(m - j, dp)*y(m - j, 1)*s(j)
            end do
            s(m) = acc_s/real(m, dp)
            c(m) = -acc_c/real(m, dp)
        end do
        fk(1) = y(k, 2)
        fk(2) = -s(k)
    end subroutine pendulum_taylor

    !> Complete elliptic integral of the first kind K(k) by the
    !> arithmetic-geometric mean (Gauss, quadratic convergence): an
    !> independent oracle sharing no code with the Taylor recursion above.
    pure function elliptic_k_agm(k) result(kk)
        real(qp), intent(in) :: k
        real(qp) :: kk, a, b, an
        integer :: it
        a = 1.0_qp
        b = sqrt(1.0_qp - k*k)
        do it = 1, 80
            an = 0.5_qp*(a + b)
            b = sqrt(a*b)
            a = an
            if (abs(a - b) < 1.0e-32_qp) exit
        end do
        kk = acos(-1.0_qp)/(2.0_qp*a)
    end function elliptic_k_agm

    !> Plain real cos/sin, disambiguated from the interval_t/idual_t `cos`,
    !> `sin` generics imported above (which have no real specific and would
    !> otherwise shadow the intrinsics for a plain-real argument).
    pure elemental function qcos(x) result(y)
        real(qp), intent(in) :: x
        real(qp) :: y
        intrinsic :: cos
        y = cos(x)
    end function qcos

    pure elemental function qsin(x) result(y)
        real(qp), intent(in) :: x
        real(qp) :: y
        intrinsic :: sin
        y = sin(x)
    end function qsin

    !> Section g(y) = y2 for the rotation system: crossing the y1-axis.
    subroutine rotation_section_y2(y, g, dgdy)
        type(interval_t), intent(in) :: y(:)
        type(interval_t), intent(out) :: g
        type(interval_t), intent(out) :: dgdy(:)
        g = y(2)
        dgdy(1) = interval(0.0_dp)
        dgdy(2) = interval(1.0_dp)
    end subroutine rotation_section_y2

    !> Section g(y) = y1 for the pendulum: crossing the angle-zero line.
    subroutine pendulum_section_y1(y, g, dgdy)
        type(interval_t), intent(in) :: y(:)
        type(interval_t), intent(out) :: g
        type(interval_t), intent(out) :: dgdy(:)
        g = y(1)
        dgdy(1) = interval(1.0_dp)
        dgdy(2) = interval(0.0_dp)
    end subroutine pendulum_section_y1

end module test_fortnum_validated_ode_fixtures

!> Validated ODE integration, checked against closed-form flows: the
!> rotation system y1' = -y2, y2' = y1 (exact flow: rotation by t) and the
!> shear system y1' = y2, y2' = 0 (exact flow: y1 + t y2, y2), both
!> integrated by `lohner_integrate` and checked by sampling many points of
!> the initial box, applying the exact closed-form flow, and requiring
!> containment in the returned enclosure. `taylor_lohner_predictor` is
!> checked separately against exp(h) on y' = y. `event_crossing_newton` is
!> checked against the exact root pi/2 of cos(t), the harmonic-oscillator
!> crossing time.
program test_fortnum_validated_ode
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, interval, width, sin, cos, &
        operator(-)
    use fortnum_validated_ode, only: ode_rhs_t, lohner_integrate, &
        taylor_lohner_predictor, event_crossing_newton, picard_apriori, &
        lohner_step_jacobian, section_crossing
    use test_fortnum_validated_ode_fixtures, only: rotation_rhs_t, &
        shear_rhs_t, exp_rhs_t, logistic_rhs_t, pendulum_rhs_t, &
        exact_rotation, exact_shear, exact_logistic, exact_logistic_jac, &
        rotation_section_y2, pendulum_section_y1, elliptic_k_agm, qcos, qsin
    implicit none

    integer :: nfail, seed_size
    integer, allocatable :: seed(:)
    type(rotation_rhs_t) :: rot
    type(shear_rhs_t) :: shr
    type(exp_rhs_t) :: expo
    type(logistic_rhs_t) :: logi
    type(pendulum_rhs_t) :: pend

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 5
    call random_seed(put=seed)

    call check_flow(rot, exact_rotation, "rotation", nfail)
    call check_flow(shr, exact_shear, "shear", nfail)
    call check_taylor_predictor(expo, nfail)
    call check_variational_jacobian(logi, nfail)
    call check_event(nfail)
    call check_section_crossing_harmonic(rot, nfail)
    call check_section_crossing_pendulum(pend, nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " validated_ode test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_validated_ode: all tests passed"

contains

    subroutine check_flow(rhs, exact, name, nfail)
        class(ode_rhs_t), intent(in) :: rhs
        interface
            pure function exact(x0, t) result(x)
                import qp
                real(qp), intent(in) :: x0(2), t
                real(qp) :: x(2)
            end function exact
        end interface
        character(*), intent(in) :: name
        integer, intent(inout) :: nfail
        type(interval_t) :: cell(2), out(2)
        real(dp) :: r(3)
        real(qp) :: x0(2), xt(2)
        logical :: ok
        integer :: i, j, nbad

        nbad = 0
        cell = [interval(0.9_dp, 1.1_dp), interval(-0.1_dp, 0.1_dp)]
        call lohner_integrate(rhs, cell, 0.6_dp, 0.05_dp, 3, out, ok)
        if (.not. ok) then
            nbad = nbad + 1
        else
            do i = 1, 400
                call random_number(r)
                x0(1) = real(cell(1)%lo, qp) + real(cell(1)%hi - cell(1)%lo, qp) &
                    *real(r(1), qp)
                x0(2) = real(cell(2)%lo, qp) + real(cell(2)%hi - cell(2)%lo, qp) &
                    *real(r(2), qp)
                xt = exact(x0, 0.6_qp)
                do j = 1, 2
                    if (.not. (real(out(j)%lo, qp) <= xt(j) .and. &
                        xt(j) <= real(out(j)%hi, qp))) nbad = nbad + 1
                end do
            end do
        end if
        call require(nbad == 0, name//": lohner_integrate encloses the exact "// &
            "closed-form flow of 400 sampled initial points", nfail)
    end subroutine check_flow

    subroutine check_taylor_predictor(rhs, nfail)
        class(ode_rhs_t), intent(in) :: rhs
        integer, intent(inout) :: nfail
        type(interval_t) :: box(1), remainder(1)
        real(dp) :: xc(1), xpred(1), jac(1, 1), h
        real(qp) :: ref
        logical :: ok

        xc(1) = 1.0_dp
        h = 0.2_dp
        box(1) = interval(0.9_dp, 1.1_dp)
        call taylor_lohner_predictor(rhs, xc, box, h, 6, xpred, remainder, jac, ok)
        ref = exp(real(h, qp))
        call require(ok .and. &
            real(xpred(1), qp) + real(remainder(1)%lo, qp) <= ref .and. &
            ref <= real(xpred(1), qp) + real(remainder(1)%hi, qp), &
            "order-6 Taylor-Lohner predictor of y'=y encloses exp(h)", nfail)
        call require(abs(jac(1, 1) - real(exp(real(h, qp)), dp)) < 1.0e-4_dp, &
            "predictor Jacobian matches d(exp(h) y0)/dy0 = exp(h)", nfail)
    end subroutine check_taylor_predictor

    !> Order-q variational Jacobian enclosure (`taylor_lohner_predictor`'s
    !> `q_jac`) against the closed-form logistic-map Jacobian: containment
    !> at several orders, width shrinking as h shrinks with q fixed, and a
    !> tighter enclosure than the first-order fallback `lohner_step_jacobian`
    !> (Q = I + h A (1 + eps)) in the moderately nonlinear regime where the
    !> resolved Taylor terms pay for the extra interval dependency they
    !> introduce. (At both very small h, where the box-Jacobian term A
    !> dominates both enclosures equally, and very large h, where the
    !> Taylor series itself is not yet resolving the nonlinearity, the two
    !> methods are not comparable in general; the tested step size is
    !> chosen where the order-q advantage is expected and verified.)
    subroutine check_variational_jacobian(rhs, nfail)
        class(ode_rhs_t), intent(in) :: rhs
        integer, intent(inout) :: nfail
        type(interval_t) :: box(1), y(1), fy(1), ay(1, 1), qjac(1, 1)
        type(interval_t) :: q1(1, 1), remainder(1)
        real(dp) :: xc(1), xpred(1), jac_pred(1, 1), h
        real(qp) :: refj
        logical :: ok
        integer :: q
        real(dp) :: w_prev, w_now

        xc(1) = 0.2_dp
        box(1) = interval(0.15_dp, 0.25_dp)

        do q = 2, 4
            call picard_apriori(rhs, box, 0.3_dp, y, fy, ay, ok)
            call require(ok, "logistic picard_apriori succeeds", nfail)
            call taylor_lohner_predictor(rhs, xc, y, 0.3_dp, q, xpred, &
                remainder, jac_pred, ok, ay=ay, q_jac=qjac)
            call require(ok, "order-q variational predictor succeeds", nfail)
            refj = exact_logistic_jac(0.2_qp, 0.3_qp)
            call require(real(qjac(1, 1)%lo, qp) <= refj .and. &
                refj <= real(qjac(1, 1)%hi, qp), &
                "order-q variational Jacobian encloses the exact logistic "// &
                "dPhi/dy0", nfail)
        end do

        ! Width shrinks with h at fixed order q = 3 (local variational
        ! error is O(h^(q+1)); the composed one-step bound should shrink
        ! well within a factor 0.5 per halving of h).
        w_prev = -1.0_dp
        do q = 1, 4
            h = 0.4_dp/real(2**q, dp)
            call picard_apriori(rhs, box, h, y, fy, ay, ok)
            call taylor_lohner_predictor(rhs, xc, y, h, 3, xpred, &
                remainder, jac_pred, ok, ay=ay, q_jac=qjac)
            w_now = width(qjac(1, 1))
            if (w_prev > 0.0_dp) then
                call require(w_now < 0.5_dp*w_prev, &
                    "order-q variational Jacobian width shrinks as h halves", &
                    nfail)
            end if
            w_prev = w_now
        end do

        ! At h = 0.15 the order-4 variational Jacobian is tighter than the
        ! first-order fallback (measured ratio about 0.89; margin 0.95
        ! leaves headroom against compiler/library rounding differences).
        h = 0.15_dp
        call picard_apriori(rhs, box, h, y, fy, ay, ok)
        call taylor_lohner_predictor(rhs, xc, y, h, 4, xpred, remainder, &
            jac_pred, ok, ay=ay, q_jac=qjac)
        q1 = lohner_step_jacobian(ay, h)
        call require(width(qjac(1, 1)) < 0.95_dp*width(q1(1, 1)), &
            "order-4 variational Jacobian is tighter than the first-order "// &
            "I + h A (1 + eps) enclosure at h = 0.15", nfail)
    end subroutine check_variational_jacobian

    ! ---- harmonic-oscillator event: g(t) = cos(t), root at pi/2 ----

    subroutine cos_event(t, g, dgdt)
        type(interval_t), intent(in) :: t
        type(interval_t), intent(out) :: g, dgdt
        g = cos(t)
        dgdt = -sin(t)
    end subroutine cos_event

    subroutine check_event(nfail)
        integer, intent(inout) :: nfail
        type(interval_t) :: troot
        logical :: ok
        real(qp) :: piq_half

        call event_crossing_newton(cos_event, 1.0_dp, 2.0_dp, troot, ok)
        piq_half = acos(-1.0_qp)/2.0_qp
        call require(ok .and. real(troot%lo, qp) <= piq_half .and. &
            piq_half <= real(troot%hi, qp) .and. width(troot) < 1.0e-9_dp, &
            "interval-Newton event crossing brackets pi/2, cos(t)'s root", nfail)
    end subroutine check_event

    !> `section_crossing` on the harmonic oscillator (rotation system),
    !> against the exact closed form: for a box of points near radius 1,
    !> angle phi0 (small spread in both a and b), the section g = y2 = 0 is
    !> next crossed (rotation preserves a sin(t + phi) with phi = atan2(b,
    !> a)) at t = pi - phi (direction < 0, downward) or t = -phi (direction
    !> > 0, upward, for phi0 < 0), exactly. Checked against 200 sampled
    !> points of the box, both directions.
    subroutine check_section_crossing_harmonic(rhs, nfail)
        class(ode_rhs_t), intent(in) :: rhs
        integer, intent(inout) :: nfail
        type(interval_t) :: cell(2), tau, ztau(2)
        real(dp) :: r(2)
        real(qp) :: a, b, phi, texact, piq
        real(qp) :: xt(2)
        logical :: ok
        integer :: i, nbad

        piq = acos(-1.0_qp)

        ! Downward crossing (direction < 0): box near angle 0.3.
        cell = [interval(real(qcos(0.3_qp), dp) - 0.005_dp, &
            real(qcos(0.3_qp), dp) + 0.005_dp), &
            interval(real(qsin(0.3_qp), dp) - 0.005_dp, &
            real(qsin(0.3_qp), dp) + 0.005_dp)]
        call section_crossing(rhs, rotation_section_y2, cell, 0.05_dp, 4, &
            -1, 300, tau, ztau, ok)
        nbad = 0
        if (.not. ok) then
            nbad = 1
        else
            do i = 1, 200
                call random_number(r)
                a = real(cell(1)%lo, qp) + real(cell(1)%hi - cell(1)%lo, qp) &
                    *real(r(1), qp)
                b = real(cell(2)%lo, qp) + real(cell(2)%hi - cell(2)%lo, qp) &
                    *real(r(2), qp)
                phi = atan2(b, a)
                texact = piq - phi
                if (.not. (real(tau%lo, qp) <= texact .and. &
                    texact <= real(tau%hi, qp))) nbad = nbad + 1
                xt = exact_rotation([a, b], texact)
                if (.not. (real(ztau(1)%lo, qp) <= xt(1) .and. &
                    xt(1) <= real(ztau(1)%hi, qp) .and. &
                    real(ztau(2)%lo, qp) <= xt(2) .and. &
                    xt(2) <= real(ztau(2)%hi, qp))) nbad = nbad + 1
            end do
        end if
        call require(nbad == 0, "section_crossing: downward y2 = 0 "// &
            "crossing of the rotation flow encloses 200 sampled exact "// &
            "crossing times and states", nfail)

        ! Upward crossing (direction > 0): box near angle -0.325 (off a
        ! step boundary so the exact crossing sits well inside the
        ! detected step's window, not at its edge).
        cell = [interval(real(qcos(-0.325_qp), dp) - 0.005_dp, &
            real(qcos(-0.325_qp), dp) + 0.005_dp), &
            interval(real(qsin(-0.325_qp), dp) - 0.005_dp, &
            real(qsin(-0.325_qp), dp) + 0.005_dp)]
        call section_crossing(rhs, rotation_section_y2, cell, 0.05_dp, 4, &
            1, 300, tau, ztau, ok)
        nbad = 0
        if (.not. ok) then
            nbad = 1
        else
            do i = 1, 200
                call random_number(r)
                a = real(cell(1)%lo, qp) + real(cell(1)%hi - cell(1)%lo, qp) &
                    *real(r(1), qp)
                b = real(cell(2)%lo, qp) + real(cell(2)%hi - cell(2)%lo, qp) &
                    *real(r(2), qp)
                phi = atan2(b, a)
                texact = -phi
                if (.not. (real(tau%lo, qp) <= texact .and. &
                    texact <= real(tau%hi, qp))) nbad = nbad + 1
                xt = exact_rotation([a, b], texact)
                if (.not. (real(ztau(1)%lo, qp) <= xt(1) .and. &
                    xt(1) <= real(ztau(1)%hi, qp) .and. &
                    real(ztau(2)%lo, qp) <= xt(2) .and. &
                    xt(2) <= real(ztau(2)%hi, qp))) nbad = nbad + 1
            end do
        end if
        call require(nbad == 0, "section_crossing: upward y2 = 0 "// &
            "crossing of the rotation flow encloses 200 sampled exact "// &
            "crossing times and states", nfail)
    end subroutine check_section_crossing_harmonic

    !> `section_crossing` on the nonlinear pendulum, against the return
    !> time and velocity at y1 = 0 from a turning point at theta0 = 0.8,
    !> checked with a real128 AGM elliptic-integral oracle
    !> (`elliptic_k_agm`) and exact energy conservation, independent of the
    !> Taylor-coefficient recursion the integrator uses.
    subroutine check_section_crossing_pendulum(rhs, nfail)
        class(ode_rhs_t), intent(in) :: rhs
        integer, intent(inout) :: nfail
        type(interval_t) :: cell(2), tau, ztau(2)
        real(qp) :: theta0, texact, vexact, k
        logical :: ok

        theta0 = 0.8_qp
        k = qsin(theta0/2.0_qp)
        texact = elliptic_k_agm(k)
        vexact = -sqrt(2.0_qp*(1.0_qp - qcos(theta0)))

        cell = [interval(real(theta0, dp) - 1.0e-7_dp, &
            real(theta0, dp) + 1.0e-7_dp), interval(-1.0e-7_dp, 1.0e-7_dp)]
        call section_crossing(rhs, pendulum_section_y1, cell, 0.05_dp, 5, &
            -1, 200, tau, ztau, ok)
        call require(ok .and. real(tau%lo, qp) <= texact .and. &
            texact <= real(tau%hi, qp), "section_crossing: pendulum "// &
            "quarter-period encloses the AGM elliptic-integral reference "// &
            "K(sin(theta0/2))", nfail)
        call require(ok .and. real(ztau(1)%lo, qp) <= 0.0_qp .and. &
            0.0_qp <= real(ztau(1)%hi, qp), &
            "section_crossing: pendulum crossing state has y1 = 0", nfail)
        call require(ok .and. real(ztau(2)%lo, qp) <= vexact .and. &
            vexact <= real(ztau(2)%hi, qp), "section_crossing: pendulum "// &
            "crossing velocity encloses the exact energy-conservation "// &
            "value", nfail)
    end subroutine check_section_crossing_pendulum

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail
        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_validated_ode
