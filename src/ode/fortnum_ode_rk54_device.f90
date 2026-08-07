module fortnum_ode_rk54_device
    ! Allocation-free, reverse-communication RK5(4) stepping for four-state
    ! GPU workloads. The application evaluates its RHS; fortnum owns the
    ! tableaux, error norm, and adaptive controller. No procedure pointer or
    ! polymorphic argument enters a device routine.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_rk54_ck_finish, only: fortnum_rk54_ck_finish
    use fortnum_generated_rk54_ck_stage2, only: fortnum_rk54_ck_stage2
    use fortnum_generated_rk54_ck_stage3, only: fortnum_rk54_ck_stage3
    use fortnum_generated_rk54_ck_stage4, only: fortnum_rk54_ck_stage4
    use fortnum_generated_rk54_ck_stage5, only: fortnum_rk54_ck_stage5
    use fortnum_generated_rk54_ck_stage6, only: fortnum_rk54_ck_stage6
    use fortnum_generated_rk54_dp_finish, only: fortnum_rk54_dp_finish
    use fortnum_generated_rk54_dp_stage2, only: fortnum_rk54_dp_stage2
    use fortnum_generated_rk54_dp_stage3, only: fortnum_rk54_dp_stage3
    use fortnum_generated_rk54_dp_stage4, only: fortnum_rk54_dp_stage4
    use fortnum_generated_rk54_dp_stage5, only: fortnum_rk54_dp_stage5
    use fortnum_generated_rk54_dp_stage6, only: fortnum_rk54_dp_stage6
    use fortnum_generated_rk54_dp_stage7, only: fortnum_rk54_dp_stage7
    implicit none
    private

    integer, parameter, public :: RK54_CASH_KARP = 1
    integer, parameter, public :: RK54_DORMAND_PRINCE = 2
    integer, parameter, public :: RK54_NEED_RHS = 1
    integer, parameter, public :: RK54_ACCEPTED = 2
    integer, parameter, public :: RK54_REJECTED = 3
    integer, parameter, public :: RK54_FAILED = 4

    real(dp), parameter :: SAFETY = 0.9_dp
    real(dp), parameter :: FAC_MIN = 0.2_dp
    real(dp), parameter :: FAC_MAX = 5.0_dp
    real(dp), parameter :: PI_ALPHA = 0.7_dp/5.0_dp
    real(dp), parameter :: PI_BETA = 0.4_dp/5.0_dp

    ! Keep device-local derived types free of component default initializers.
    ! NVHPC 26.5 emits duplicate hidden device globals when two kernels use a
    ! default-initialized type. Callers set every control field explicitly;
    ! rk54_initialize4 initializes every state field.
    type, public :: rk54_controls4_t
        integer :: method
        real(dp) :: rtol
        real(dp) :: atol(4)
        real(dp) :: hmin
        real(dp) :: hmax
    end type rk54_controls4_t

    type, public :: rk54_state4_t
        real(dp) :: t
        real(dp) :: y(4)
        real(dp) :: h
        real(dp) :: k(4, 7)
        real(dp) :: trial(4)
        real(dp) :: error(4)
        real(dp) :: previous_error
        real(dp) :: last_error
        integer :: stage
        integer :: nfev
        integer :: naccepted
        integer :: nrejected
        logical :: first_step
        logical :: after_reject
    end type rk54_state4_t

    public :: rk54_initialize4, rk54_request4, rk54_supply4
    public :: rk54_firm3d_error_norm4

contains

    !NVF$ INLINE
    pure subroutine rk54_initialize4(state, t, y, h)
        !$acc routine seq
        type(rk54_state4_t), intent(out) :: state
        real(dp), intent(in) :: t, y(4), h

        state%t = t
        state%y = y
        state%h = h
        state%k = 0.0_dp
        state%trial = 0.0_dp
        state%error = 0.0_dp
        state%previous_error = 1.0_dp
        state%last_error = 0.0_dp
        state%stage = 1
        state%nfev = 0
        state%naccepted = 0
        state%nrejected = 0
        state%first_step = .true.
        state%after_reject = .false.
    end subroutine rk54_initialize4

    !NVF$ INLINE
    pure subroutine rk54_request4(state, controls, t_eval, y_eval, request)
        !$acc routine seq
        type(rk54_state4_t), intent(in) :: state
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(out) :: t_eval, y_eval(4)
        integer, intent(out) :: request

        request = RK54_NEED_RHS
        select case (controls%method)
        case (RK54_CASH_KARP)
            call request_cash_karp(state, t_eval, y_eval, request)
        case (RK54_DORMAND_PRINCE)
            call request_dormand_prince(state, t_eval, y_eval, request)
        case default
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
        end select
    end subroutine rk54_request4

    !NVF$ INLINE
    pure subroutine rk54_supply4(state, controls, derivative, t_eval, y_eval, &
            request)
        !$acc routine seq
        type(rk54_state4_t), intent(inout) :: state
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(in) :: derivative(4)
        real(dp), intent(out) :: t_eval, y_eval(4)
        integer, intent(out) :: request
        integer :: last_stage

        if (state%stage < 1 .or. state%stage > 7) then
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
            return
        end if
        if (controls%rtol < 0.0_dp .or. any(controls%atol <= 0.0_dp)) then
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
            return
        end if

        state%k(:, state%stage) = derivative
        state%nfev = state%nfev + 1
        if (controls%method == RK54_CASH_KARP) then
            last_stage = 6
        else if (controls%method == RK54_DORMAND_PRINCE) then
            last_stage = 7
        else
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
            return
        end if

        if (state%stage < last_stage) then
            state%stage = state%stage + 1
            call rk54_request4(state, controls, t_eval, y_eval, request)
            return
        end if

        if (controls%method == RK54_CASH_KARP) then
            call finish_cash_karp(state)
            state%last_error = cash_karp_error_norm4(state, controls)
        else
            call finish_dormand_prince(state)
            state%last_error = rk54_firm3d_error_norm4(state%y, &
                state%k(:, 1), state%error, state%h, controls%rtol, &
                controls%atol)
        end if
        call finish_attempt(state, controls, t_eval, y_eval, request)
    end subroutine rk54_supply4

    !NVF$ INLINE
    pure subroutine request_cash_karp(state, t_eval, y_eval, request)
        !$acc routine seq
        type(rk54_state4_t), intent(in) :: state
        real(dp), intent(out) :: t_eval, y_eval(4)
        integer, intent(inout) :: request

        select case (state%stage)
        case (1)
            t_eval = state%t
            y_eval = state%y
        case (2)
            t_eval = state%t + state%h/5.0_dp
            call fortnum_rk54_ck_stage2(state%y, state%h, state%k(:, 1), y_eval)
        case (3)
            t_eval = state%t + 3.0_dp*state%h/10.0_dp
            call fortnum_rk54_ck_stage3(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), y_eval)
        case (4)
            t_eval = state%t + 3.0_dp*state%h/5.0_dp
            call fortnum_rk54_ck_stage4(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), y_eval)
        case (5)
            t_eval = state%t + state%h
            call fortnum_rk54_ck_stage5(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), state%k(:, 4), y_eval)
        case (6)
            t_eval = state%t + 7.0_dp*state%h/8.0_dp
            call fortnum_rk54_ck_stage6(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), state%k(:, 4), &
                state%k(:, 5), y_eval)
        case default
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
        end select
    end subroutine request_cash_karp

    !NVF$ INLINE
    pure subroutine request_dormand_prince(state, t_eval, y_eval, request)
        !$acc routine seq
        type(rk54_state4_t), intent(in) :: state
        real(dp), intent(out) :: t_eval, y_eval(4)
        integer, intent(inout) :: request

        select case (state%stage)
        case (1)
            t_eval = state%t
            y_eval = state%y
        case (2)
            t_eval = state%t + state%h/5.0_dp
            call fortnum_rk54_dp_stage2(state%y, state%h, state%k(:, 1), y_eval)
        case (3)
            t_eval = state%t + 3.0_dp*state%h/10.0_dp
            call fortnum_rk54_dp_stage3(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), y_eval)
        case (4)
            t_eval = state%t + 4.0_dp*state%h/5.0_dp
            call fortnum_rk54_dp_stage4(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), y_eval)
        case (5)
            t_eval = state%t + 8.0_dp*state%h/9.0_dp
            call fortnum_rk54_dp_stage5(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), state%k(:, 4), y_eval)
        case (6)
            t_eval = state%t + state%h
            call fortnum_rk54_dp_stage6(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), state%k(:, 4), &
                state%k(:, 5), y_eval)
        case (7)
            t_eval = state%t + state%h
            call fortnum_rk54_dp_stage7(state%y, state%h, state%k(:, 1), &
                state%k(:, 2), state%k(:, 3), state%k(:, 4), &
                state%k(:, 5), state%k(:, 6), y_eval)
        case default
            t_eval = state%t
            y_eval = state%y
            request = RK54_FAILED
        end select
    end subroutine request_dormand_prince

    !NVF$ INLINE
    pure subroutine finish_cash_karp(state)
        !$acc routine seq
        type(rk54_state4_t), intent(inout) :: state

        call fortnum_rk54_ck_finish(state%y, state%h, state%k(:, 1), &
            state%k(:, 3), state%k(:, 4), state%k(:, 5), state%k(:, 6), &
            state%trial, state%error)
    end subroutine finish_cash_karp

    !NVF$ INLINE
    pure subroutine finish_dormand_prince(state)
        !$acc routine seq
        type(rk54_state4_t), intent(inout) :: state

        call fortnum_rk54_dp_finish(state%y, state%h, state%k(:, 1), &
            state%k(:, 3), state%k(:, 4), state%k(:, 5), state%k(:, 6), &
            state%k(:, 7), state%trial, state%error)
    end subroutine finish_dormand_prince

    ! FIRM3D/CATAPULT's exact four-component maximum norm.
    !NVF$ INLINE
    pure function rk54_firm3d_error_norm4(y, k1, error, h, rtol, atol) result(norm)
        !$acc routine seq
        real(dp), intent(in) :: y(4), k1(4), error(4), h, rtol, atol(4)
        real(dp) :: norm, scale
        integer :: i

        norm = 0.0_dp
        do i = 1, 4
            scale = atol(i) + rtol*(abs(y(i)) + h*abs(k1(i)))
            norm = max(norm, abs(error(i))/scale)
        end do
    end function rk54_firm3d_error_norm4

    !NVF$ INLINE
    pure function cash_karp_error_norm4(state, controls) result(norm)
        !$acc routine seq
        type(rk54_state4_t), intent(in) :: state
        type(rk54_controls4_t), intent(in) :: controls
        real(dp) :: norm, scale, sumsq
        integer :: i

        sumsq = 0.0_dp
        do i = 1, 4
            scale = controls%atol(i) + controls%rtol* &
                max(abs(state%y(i)), abs(state%trial(i)))
            sumsq = sumsq + (state%error(i)/scale)**2
        end do
        norm = sqrt(sumsq/4.0_dp)
    end function cash_karp_error_norm4

    !NVF$ INLINE
    pure subroutine finish_attempt(state, controls, t_eval, y_eval, request)
        !$acc routine seq
        type(rk54_state4_t), intent(inout) :: state
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(out) :: t_eval, y_eval(4)
        integer, intent(out) :: request
        real(dp) :: factor, hnew, error_floor
        logical :: accepted, at_minimum

        error_floor = max(state%last_error, 1.0e-300_dp)
        if (controls%method == RK54_DORMAND_PRINCE) then
            factor = SAFETY*error_floor**(-1.0_dp/3.0_dp)
            factor = max(FAC_MIN, min(FAC_MAX, factor))
            if (state%last_error > 0.5_dp) then
                if (state%last_error < 1.0_dp) factor = 1.0_dp
            end if
        else if (state%first_step .or. state%after_reject) then
            factor = SAFETY*max(state%last_error, 1.0e-10_dp)**(-1.0_dp/5.0_dp)
            factor = max(FAC_MIN, min(FAC_MAX, factor))
        else
            factor = SAFETY*max(state%last_error, 1.0e-10_dp)**(-PI_ALPHA)* &
                state%previous_error**PI_BETA
            factor = max(FAC_MIN, min(FAC_MAX, factor))
        end if
        hnew = max(controls%hmin, min(controls%hmax, state%h*factor))
        at_minimum = controls%hmin > 0.0_dp
        if (at_minimum) at_minimum = state%h <= controls%hmin
        accepted = state%last_error <= 1.0_dp .or. at_minimum

        if (accepted) then
            state%t = state%t + state%h
            state%y = state%trial
            if (controls%method == RK54_DORMAND_PRINCE) &
                state%k(:, 1) = state%k(:, 7)
            state%h = hnew
            state%previous_error = max(state%last_error, 1.0e-10_dp)
            state%first_step = .false.
            state%after_reject = .false.
            state%naccepted = state%naccepted + 1
            request = RK54_ACCEPTED
        else
            state%h = hnew
            state%after_reject = .true.
            state%nrejected = state%nrejected + 1
            request = RK54_REJECTED
        end if
        if (controls%method == RK54_DORMAND_PRINCE) then
            ! The accepted step's seventh stage is f(t+h, y_new). After a
            ! rejection the original first stage is still f(t, y). Both are
            ! valid first stages for the next attempt.
            state%stage = 2
        else
            state%stage = 1
        end if
        t_eval = state%t
        y_eval = state%y
    end subroutine finish_attempt

end module fortnum_ode_rk54_device
