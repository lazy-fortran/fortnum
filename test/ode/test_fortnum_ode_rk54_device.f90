program test_fortnum_ode_rk54_device
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ode_rk54_device, only: rk54_controls4_t, rk54_state4_t, &
        rk54_initialize4, rk54_request4, rk54_supply4, &
        RK54_CASH_KARP, RK54_DORMAND_PRINCE, RK54_NEED_RHS, &
        RK54_DORMAND_PRINCE_TUNED, RK54_ACCEPTED, RK54_REJECTED, RK54_FAILED
    implicit none

    integer :: failures

    failures = 0
    call test_dormand_prince_external_fixture(RK54_DORMAND_PRINCE, &
        "DOPRI")
    call test_dormand_prince_external_fixture(RK54_DORMAND_PRINCE_TUNED, &
        "tuned DOPRI")
    call test_dormand_prince_fsal(RK54_DORMAND_PRINCE, "DOPRI")
    call test_dormand_prince_fsal(RK54_DORMAND_PRINCE_TUNED, "tuned DOPRI")
    call test_exact_exponential(RK54_CASH_KARP)
    call test_exact_exponential(RK54_DORMAND_PRINCE)
    call test_exact_exponential(RK54_DORMAND_PRINCE_TUNED)
    call test_device_exponential(RK54_CASH_KARP)
    call test_device_exponential(RK54_DORMAND_PRINCE)
    call test_device_exponential(RK54_DORMAND_PRINCE_TUNED)
    call test_tuned_controller()
    call test_rejection_and_minimum_step()
    if (failures /= 0) then
        write (error_unit, "(i0,a)") failures, " RK54 device test(s) failed"
        error stop 1
    end if

contains

    subroutine test_dormand_prince_external_fixture(method, label)
        ! Golden values were evaluated at 80 decimal digits from the published
        ! Dormand-Prince tableau, independently of the generated Fortran.
        real(dp), parameter :: expected_y(4) = [ &
            0.9900498337491683333333333333333333333333_dp, &
            2.040402680053546666666666666666666666667_dp, &
            -0.970445533548715_dp, &
            0.5012515638028975425211588541666666666667_dp]
        real(dp), parameter :: expected_error(4) = [ &
            8.11587949622882842e-14_dp, -5.13183995298760436e-12_dp, &
            -1.98803366824318315e-11_dp, -3.94476118437125968e-17_dp]
        real(dp), parameter :: expected_norm = &
            9.79326930169055741e-6_dp
        real(dp), parameter :: lambda(4) = [-1.0_dp, 2.0_dp, -3.0_dp, 0.25_dp]
        integer, intent(in) :: method
        character(*), intent(in) :: label
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4)
        integer :: request

        controls%method = method
        controls%rtol = 1.0e-6_dp
        controls%atol = [1.0e-6_dp, 1.0e-6_dp, 1.0e-6_dp, 13.0_dp]
        controls%hmin = 0.0_dp
        controls%hmax = 0.01_dp
        call rk54_initialize4(state, 0.0_dp, [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp], &
            0.01_dp)
        call rk54_request4(state, controls, t_eval, y_eval, request)
        do while (request == RK54_NEED_RHS)
            derivative = lambda*y_eval
            call rk54_supply4(state, controls, derivative, t_eval, y_eval, request)
        end do

        call check(request == RK54_ACCEPTED, label//" fixture accepts")
        call check(maxval(abs(state%y - expected_y)) < 3.0e-16_dp, &
            label//" fixture fifth-order state")
        if (maxval(abs(state%error - expected_error)) >= 3.0e-26_dp) then
            write (error_unit, "(a,4es24.16)") "actual error: ", state%error
        end if
        call check(maxval(abs(state%error - expected_error)) < 1.0e-18_dp, &
            label//" fixture embedded error")
        if (abs(state%last_error - expected_norm) >= 3.0e-19_dp) then
            write (error_unit, "(a,es24.16)") "actual norm: ", state%last_error
        end if
        call check(abs(state%last_error - expected_norm) < 5.0e-12_dp, &
            label//" FIRM3D maximum error norm")
        call check(state%nfev == 7, label//" fixture has seven RHS evaluations")
    end subroutine test_dormand_prince_external_fixture

    subroutine test_dormand_prince_fsal(method, label)
        ! Two fixed steps of the published DOPRI tableau require 7 + 6 RHS
        ! evaluations because stage 7 of the first step is stage 1 of the next.
        real(dp), parameter :: lambda(4) = [-1.0_dp, 2.0_dp, -3.0_dp, 0.25_dp]
        real(dp), parameter :: initial(4) = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp]
        integer, intent(in) :: method
        character(*), intent(in) :: label
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4)
        integer :: request, step

        controls%method = method
        controls%rtol = 1.0e-6_dp
        controls%atol = 1.0e-6_dp
        controls%hmin = 0.01_dp
        controls%hmax = 0.01_dp
        call rk54_initialize4(state, 0.0_dp, initial, 0.01_dp)
        do step = 1, 2
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                derivative = lambda*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, y_eval, &
                    request)
            end do
            call check(request == RK54_ACCEPTED, label//" FSAL fixture accepts")
        end do
        call check(state%nfev == 13, label//" FSAL uses 13 RHS evaluations")
        call check(maxval(abs(state%y - initial*exp(0.02_dp*lambda))) < 2.0e-12_dp, &
            label//" FSAL agrees with exact exponential")
    end subroutine test_dormand_prince_fsal

    subroutine test_exact_exponential(method)
        integer, intent(in) :: method
        real(dp), parameter :: lambda(4) = [-1.0_dp, 0.5_dp, 1.0_dp, -2.0_dp]
        real(dp), parameter :: initial(4) = [1.0_dp, -2.0_dp, 0.25_dp, 3.0_dp]
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4), remaining
        integer :: request, attempts

        controls%method = method
        controls%rtol = 1.0e-10_dp
        controls%atol = 1.0e-12_dp
        controls%hmin = 1.0e-14_dp
        controls%hmax = 0.2_dp
        call rk54_initialize4(state, 0.0_dp, initial, 0.02_dp)
        attempts = 0
        do while (state%t < 1.0_dp)
            remaining = 1.0_dp - state%t
            if (state%h > remaining) state%h = remaining
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                derivative = lambda*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, y_eval, &
                    request)
            end do
            call check(request == RK54_ACCEPTED .or. request == RK54_REJECTED, &
                "adaptive request terminates normally")
            attempts = attempts + 1
            if (attempts > 10000) exit
        end do
        call check(attempts <= 10000, "adaptive trace finishes")
        call check(maxval(abs(state%y - initial*exp(lambda))) < 2.0e-9_dp, &
            "adaptive trace agrees with exact exponential")
    end subroutine test_exact_exponential

    subroutine test_device_exponential(method)
        integer, parameter :: batch_size = 128
        integer, intent(in) :: method
        real(dp), parameter :: lambda(4) = [-1.0_dp, 0.5_dp, 1.0_dp, -2.0_dp]
        real(dp), parameter :: initial(4) = [1.0_dp, -2.0_dp, 0.25_dp, 3.0_dp]
        real(dp) :: result(4, batch_size), final_time(batch_size)
        integer :: status(batch_size), attempts(batch_size), particle
        integer :: request
        type(rk54_controls4_t) :: controls

        controls%method = method
        controls%rtol = 1.0e-10_dp
        controls%atol = 1.0e-12_dp
        controls%hmin = 1.0e-14_dp
        controls%hmax = 0.2_dp
        result = 0.0_dp
        final_time = 0.0_dp
        status = RK54_FAILED
        attempts = 0

        !$acc parallel loop copyout(result, final_time, status, attempts) &
        !$acc& firstprivate(controls)
        do particle = 1, batch_size
            call integrate_exponential_device(controls, lambda, initial, &
                result(:, particle), final_time(particle), request, &
                attempts(particle))
            status(particle) = request
        end do
        !$acc end parallel loop

        call check(all(status == RK54_ACCEPTED), &
            "device adaptive requests terminate normally")
        call check(all(attempts < 10000), "device adaptive traces finish")
        call check(maxval(abs(final_time - 1.0_dp)) < 2.0e-15_dp, &
            "device adaptive traces reach final time")
        call check(maxval(abs(result - spread(initial*exp(lambda), 2, &
            batch_size))) < 2.0e-9_dp, &
            "device adaptive traces agree with exact exponential")
    end subroutine test_device_exponential

    subroutine test_tuned_controller()
        ! The external fixture has a maximum error of 1.98803366824318315e-11.
        ! With zero relative tolerance and atol=1e-11, the FIRM3D norm is the
        ! known value below. The tuned first-step rule is 0.9*error**(-1/5),
        ! while the compatibility controller uses error**(-1/3).
        real(dp), parameter :: expected_norm = 1.98803366824318315_dp
        real(dp), parameter :: initial_h = 0.01_dp
        real(dp), parameter :: lambda(4) = [-1.0_dp, 2.0_dp, -3.0_dp, 0.25_dp]
        real(dp), parameter :: initial(4) = [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp]
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: tuned, original
        real(dp) :: t_eval, y_eval(4), derivative(4), expected_h
        integer :: request

        controls%rtol = 0.0_dp
        controls%atol = 1.0e-11_dp
        controls%hmin = 0.0_dp
        controls%hmax = 1.0_dp

        controls%method = RK54_DORMAND_PRINCE_TUNED
        call rk54_initialize4(tuned, 0.0_dp, initial, initial_h)
        call rk54_request4(tuned, controls, t_eval, y_eval, request)
        do while (request == RK54_NEED_RHS)
            derivative = lambda*y_eval
            call rk54_supply4(tuned, controls, derivative, t_eval, y_eval, &
                request)
        end do
        expected_h = initial_h*0.9_dp*expected_norm**(-1.0_dp/5.0_dp)
        call check(request == RK54_REJECTED, "tuned controller rejects")
        call check(abs(tuned%h - expected_h) < 2.0e-10_dp, &
            "tuned controller uses fifth-order exponent")
        call check(tuned%nfev == 7, "tuned controller evaluates seven stages")

        controls%method = RK54_DORMAND_PRINCE
        call rk54_initialize4(original, 0.0_dp, initial, initial_h)
        call rk54_request4(original, controls, t_eval, y_eval, request)
        do while (request == RK54_NEED_RHS)
            derivative = lambda*y_eval
            call rk54_supply4(original, controls, derivative, t_eval, y_eval, &
                request)
        end do
        expected_h = initial_h*0.9_dp*expected_norm**(-1.0_dp/3.0_dp)
        call check(request == RK54_REJECTED, "compatibility controller rejects")
        call check(abs(original%h - expected_h) < 2.0e-10_dp, &
            "compatibility controller keeps cubic exponent")
        call check(tuned%h > original%h, &
            "tuned controller takes the larger rejected-step retry")
    end subroutine test_tuned_controller

    subroutine integrate_exponential_device(controls, lambda, initial, result, &
            final_time, status, attempts)
        !$acc routine seq
        type(rk54_controls4_t), intent(in) :: controls
        real(dp), intent(in) :: lambda(4), initial(4)
        real(dp), intent(out) :: result(4), final_time
        integer, intent(out) :: status, attempts

        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4), remaining
        integer :: rhs_calls

        call rk54_initialize4(state, 0.0_dp, initial, 0.02_dp)
        status = RK54_ACCEPTED
        do attempts = 1, 10000
            if (state%t >= 1.0_dp) exit
            remaining = 1.0_dp - state%t
            if (state%h > remaining) state%h = remaining
            call rk54_request4(state, controls, t_eval, y_eval, status)
            rhs_calls = 0
            do while (status == RK54_NEED_RHS)
                derivative = lambda*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, &
                    y_eval, status)
                rhs_calls = rhs_calls + 1
                if (rhs_calls > 7) then
                    status = RK54_FAILED
                    exit
                end if
            end do
            if (status == RK54_FAILED) exit
        end do
        attempts = min(attempts, 10000)
        result = state%y
        final_time = state%t
    end subroutine integrate_exponential_device

    subroutine test_rejection_and_minimum_step()
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4)
        integer :: request, previous_nfev

        controls%method = RK54_DORMAND_PRINCE
        controls%rtol = 1.0e-15_dp
        controls%atol = 1.0e-15_dp
        controls%hmin = 0.02_dp
        controls%hmax = 0.1_dp
        call rk54_initialize4(state, 0.0_dp, &
            [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 0.1_dp)
        call rk54_request4(state, controls, t_eval, y_eval, request)
        do while (request == RK54_NEED_RHS)
            derivative = 100.0_dp*y_eval
            call rk54_supply4(state, controls, derivative, t_eval, y_eval, request)
        end do
        call check(request == RK54_REJECTED, "large-error step rejects")
        call check(abs(state%h - 0.02_dp) < epsilon(1.0_dp), &
            "rejection clamps to minimum step")

        do while (request == RK54_REJECTED)
            previous_nfev = state%nfev
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                derivative = 100.0_dp*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, y_eval, &
                    request)
            end do
            call check(state%nfev - previous_nfev == 6, &
                "DOPRI rejection reuses its unchanged first stage")
        end do
        call check(request == RK54_ACCEPTED, &
            "FIRM3D policy accepts at minimum step")
    end subroutine test_rejection_and_minimum_step

    subroutine check(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (condition) return
        failures = failures + 1
        write (error_unit, "(a)") "FAIL: "//label
    end subroutine check

end program test_fortnum_ode_rk54_device
