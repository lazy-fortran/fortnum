program test_concurrency_ode
    ! M5.1 concurrency: ode (per-thread workspace + solution).
    ! ode_integrate and ode_at carry caller-owned ode_workspace_t and write a
    ! caller-owned ode_solution_t / y_out. Thread safety rests on each thread
    ! owning its own workspace and solution; the rhs pointer is read-only.
    ! This test integrates a family of decay problems (varying rate via ctx)
    ! serially, then with one private workspace per thread, and asserts
    ! bit-for-bit equality of the solution at the final time.
    !
    ! PRIMAL concurrency only; derivative-product ode concurrency (sensitivity
    ! integration) lands in M6 (#40).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate
    use fortnum_ode_wrapper, only: ode_at
    implicit none

    integer, parameter :: nprob = 256
    integer :: nfail, j
    real(dp) :: rates(nprob)
    real(dp) :: ser_final(nprob), par_final(nprob)
    real(dp) :: ser_at(11, nprob), par_at(11, nprob)

    nfail = 0
    do j = 1, nprob
        rates(j) = 0.5_dp + 1.5_dp * real(j - 1, dp) / real(nprob, dp)
    end do

    call serial_integrate()
    call parallel_integrate()
    call check_exact1(ser_final, par_final, "ode_integrate final parallel == serial")

    call serial_at()
    call parallel_at()
    call check_exact2(ser_at, par_at, "ode_at samples parallel == serial")

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_ode: all tests passed"

contains

    ! Decay y' = -y.
    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
            dydt(1) = -y(1)
        end subroutine rhs_decay

        ! ode_integrate/ode_at call problem%rhs without forwarding a ctx, so the
        ! per-problem rate cannot ride through ctx here. Encode it instead as a
        ! per-problem final time: with rhs y'=-y, the state at t1 is exp(-t1), so
        ! varying t1 = 3*rate gives a distinct integration per j. The point of the
        ! test is per-thread workspace isolation, not the rate channel.
        subroutine build_problem(prob, rate)
            type(ode_problem_t), intent(out) :: prob
            real(dp),            intent(in)  :: rate
            prob%rhs  => rhs_decay
            prob%t0   = 0.0_dp
            prob%t1   = 3.0_dp * rate
            prob%y0   = [1.0_dp]
            prob%rtol = 1.0e-9_dp
            prob%atol = 1.0e-11_dp
        end subroutine build_problem

        subroutine serial_integrate()
            type(ode_problem_t)    :: prob
            type(ode_workspace_t)  :: ws
            type(ode_solution_t)   :: sol
            type(fortnum_status_t) :: st
            do j = 1, nprob
                call build_problem(prob, rates(j))
                call ode_integrate(prob, ws, sol, st)
                ser_final(j) = sol%y(1, sol%nsteps)
            end do
        end subroutine serial_integrate

        subroutine parallel_integrate()
            type(ode_problem_t)    :: prob
            type(ode_workspace_t)  :: ws
            type(ode_solution_t)   :: sol
            type(fortnum_status_t) :: st
            !$omp parallel do default(shared) private(prob, ws, sol, st, j) &
            !$omp   schedule(static)
            do j = 1, nprob
                call build_problem(prob, rates(j))
                call ode_integrate(prob, ws, sol, st)
                par_final(j) = sol%y(1, sol%nsteps)
            end do
            !$omp end parallel do
        end subroutine parallel_integrate

        subroutine serial_at()
            type(ode_problem_t)    :: prob
            type(ode_workspace_t)  :: ws
            type(fortnum_status_t) :: st
            real(dp), allocatable  :: y_out(:,:)
            real(dp) :: t_eval(11)
            integer  :: i
            do j = 1, nprob
                call build_problem(prob, rates(j))
                do i = 1, 11
                    t_eval(i) = prob%t1 * real(i - 1, dp) / 10.0_dp
                end do
                call ode_at(prob, t_eval, ws, y_out, st)
                ser_at(:, j) = y_out(1, :)
            end do
        end subroutine serial_at

        subroutine parallel_at()
            type(ode_problem_t)    :: prob
            type(ode_workspace_t)  :: ws
            type(fortnum_status_t) :: st
            real(dp), allocatable  :: y_out(:,:)
            real(dp) :: t_eval(11)
            integer  :: i
            !$omp parallel do default(shared) private(prob, ws, st, y_out, t_eval, i, j) &
            !$omp   schedule(static)
            do j = 1, nprob
                call build_problem(prob, rates(j))
                do i = 1, 11
                    t_eval(i) = prob%t1 * real(i - 1, dp) / 10.0_dp
                end do
                call ode_at(prob, t_eval, ws, y_out, st)
                par_at(:, j) = y_out(1, :)
            end do
            !$omp end parallel do
        end subroutine parallel_at

        subroutine check_exact1(ref, got, name)
            real(dp),     intent(in) :: ref(:), got(:)
            character(*), intent(in) :: name
            integer :: k
            do k = 1, size(ref)
                if (ref(k) /= got(k)) then
                    nfail = nfail + 1
                    write (error_unit, "(a,a,a,i0)") "FAIL: ", name, " at index ", k
                    return
                end if
            end do
        end subroutine check_exact1

        subroutine check_exact2(ref, got, name)
            real(dp),     intent(in) :: ref(:,:), got(:,:)
            character(*), intent(in) :: name
            integer :: a, b
            do b = 1, size(ref, 2)
                do a = 1, size(ref, 1)
                    if (ref(a, b) /= got(a, b)) then
                        nfail = nfail + 1
                        write (error_unit, "(a,a)") "FAIL: ", name
                        return
                    end if
                end do
            end do
        end subroutine check_exact2

    end program test_concurrency_ode
