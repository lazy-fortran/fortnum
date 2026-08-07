program test_fortnum_ode_rk54_device
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ode_rk54_device, only: rk54_controls4_t, rk54_state4_t, &
        rk54_initialize4, rk54_request4, rk54_supply4, &
        RK54_CASH_KARP, RK54_DORMAND_PRINCE, RK54_NEED_RHS, &
        RK54_ACCEPTED, RK54_REJECTED
    implicit none

    integer :: failures

    failures = 0
    call test_dormand_prince_external_fixture()
    call test_exact_exponential(RK54_CASH_KARP)
    call test_exact_exponential(RK54_DORMAND_PRINCE)
    call test_rejection_and_minimum_step()
    if (failures /= 0) then
        write (error_unit, "(i0,a)") failures, " RK54 device test(s) failed"
        error stop 1
    end if

contains

    subroutine test_dormand_prince_external_fixture()
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
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4)
        integer :: request

        controls%method = RK54_DORMAND_PRINCE
        controls%rtol = 1.0e-6_dp
        controls%atol = [1.0e-6_dp, 1.0e-6_dp, 1.0e-6_dp, 13.0_dp]
        controls%hmax = 0.01_dp
        call rk54_initialize4(state, 0.0_dp, [1.0_dp, 2.0_dp, -1.0_dp, 0.5_dp], &
            0.01_dp)
        call rk54_request4(state, controls, t_eval, y_eval, request)
        do while (request == RK54_NEED_RHS)
            derivative = lambda*y_eval
            call rk54_supply4(state, controls, derivative, t_eval, y_eval, request)
        end do

        call check(request == RK54_ACCEPTED, "DOPRI fixture accepts")
        call check(maxval(abs(state%y - expected_y)) < 3.0e-16_dp, &
            "DOPRI fixture fifth-order state")
        if (maxval(abs(state%error - expected_error)) >= 3.0e-26_dp) then
            write (error_unit, "(a,4es24.16)") "actual error: ", state%error
        end if
        call check(maxval(abs(state%error - expected_error)) < 1.0e-18_dp, &
            "DOPRI fixture embedded error")
        if (abs(state%last_error - expected_norm) >= 3.0e-19_dp) then
            write (error_unit, "(a,es24.16)") "actual norm: ", state%last_error
        end if
        call check(abs(state%last_error - expected_norm) < 5.0e-12_dp, &
            "FIRM3D maximum error norm")
        call check(state%nfev == 7, "DOPRI fixture has seven RHS evaluations")
    end subroutine test_dormand_prince_external_fixture

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

    subroutine test_rejection_and_minimum_step()
        type(rk54_controls4_t) :: controls
        type(rk54_state4_t) :: state
        real(dp) :: t_eval, y_eval(4), derivative(4)
        integer :: request

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
            call rk54_request4(state, controls, t_eval, y_eval, request)
            do while (request == RK54_NEED_RHS)
                derivative = 100.0_dp*y_eval
                call rk54_supply4(state, controls, derivative, t_eval, y_eval, &
                    request)
            end do
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
