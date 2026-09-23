!> Fixture right-hand sides for test_fortnum_validated_ode: type-bound
!> procedures must be module procedures in Fortran, so the ode_rhs_t
!> extensions live here and the program below just uses them.
module test_fortnum_validated_ode_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*), operator(/), sin, cos
    use fortnum_idual, only: idual_t, operator(+), operator(-), operator(*), &
        operator(/)
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
        taylor_lohner_predictor, event_crossing_newton
    use test_fortnum_validated_ode_fixtures, only: rotation_rhs_t, &
        shear_rhs_t, exp_rhs_t, exact_rotation, exact_shear
    implicit none

    integer :: nfail, seed_size
    integer, allocatable :: seed(:)
    type(rotation_rhs_t) :: rot
    type(shear_rhs_t) :: shr
    type(exp_rhs_t) :: expo

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 5
    call random_seed(put=seed)

    call check_flow(rot, exact_rotation, "rotation", nfail)
    call check_flow(shr, exact_shear, "shear", nfail)
    call check_taylor_predictor(expo, nfail)
    call check_event(nfail)

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
