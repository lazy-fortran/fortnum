program test_fortnum_ode_wrapper
    ! Behavioural tests for fortnum_ode_wrapper (ode_at).
    !   - decay y' = -y at fixed output points matches exp(-t)
    !   - growth y' = y at fixed output points matches exp(t)
    !   - harmonic oscillator conserves energy across all output points
    !   - single output point (n=1) returns y0 without integrating
    !   - bad input (t_eval not monotone) is rejected with FORTNUM_DOMAIN_ERROR
    !   - bad input (null rhs) is rejected with FORTNUM_DOMAIN_ERROR

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_workspace_t
    use fortnum_ode_wrapper, only: ode_at
    implicit none

    integer :: nfail
    real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp

    nfail = 0
    call check_decay(nfail)
    call check_growth(nfail)
    call check_harmonic_energy(nfail)
    call check_single_point(nfail)
    call check_nonmonotone(nfail)
    call check_null_rhs(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_fortnum_ode_wrapper: all tests passed"
    stop 0

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)
    end subroutine rhs_decay

    subroutine rhs_growth(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = y(1)
    end subroutine rhs_growth

    ! Harmonic oscillator y1'=y2, y2'=-y1; energy = 0.5*(y1^2+y2^2).
    subroutine rhs_osc(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) =  y(2)
        dydt(2) = -y(1)
    end subroutine rhs_osc

    subroutine check_decay(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(11), exact, err
        integer  :: i

        do i = 1, 11
            t_eval(i) = (i-1) * 0.5_dp   ! 0.0, 0.5, ..., 5.0
        end do

        problem%rhs  => rhs_decay
        problem%y0   = [1.0_dp]
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-11_dp

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "FAIL check_decay: code=", status%code
            nfail = nfail + 1
            return
        end if
        if (size(y_out, 2) /= 11) then
            write (error_unit, "(a,i0)") "FAIL check_decay: n_out=", size(y_out, 2)
            nfail = nfail + 1
            return
        end if
        do i = 1, 11
            exact = exp(-t_eval(i))
            err   = abs(y_out(1,i) - exact)
            if (err > 1.0e-7_dp) then
                write (error_unit, "(a,i0,a,es12.4,a,es12.4)") &
                    "FAIL check_decay i=", i, " err=", err, " t=", t_eval(i)
                nfail = nfail + 1
            end if
        end do
    end subroutine check_decay

    subroutine check_growth(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(7), exact, err
        integer  :: i

        do i = 1, 7
            t_eval(i) = (i-1) * 0.5_dp   ! 0.0, 0.5, ..., 3.0
        end do

        problem%rhs  => rhs_growth
        problem%y0   = [1.0_dp]
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-11_dp

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "FAIL check_growth: code=", status%code
            nfail = nfail + 1
            return
        end if
        do i = 1, 7
            exact = exp(t_eval(i))
            err   = abs(y_out(1,i) - exact)
            if (err > 1.0e-7_dp) then
                write (error_unit, "(a,i0,a,es12.4)") &
                    "FAIL check_growth i=", i, " err=", err
                nfail = nfail + 1
            end if
        end do
    end subroutine check_growth

    subroutine check_harmonic_energy(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(13), e0, ek, maxdrift
        integer  :: i

        do i = 1, 13
            t_eval(i) = (i-1) * (PI / 2.0_dp)   ! quarter-period steps
        end do

        problem%rhs  => rhs_osc
        problem%y0   = [1.0_dp, 0.0_dp]
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "FAIL check_harmonic_energy: code=", status%code
            nfail = nfail + 1
            return
        end if

        e0 = 0.5_dp * (y_out(1,1)**2 + y_out(2,1)**2)
        maxdrift = 0.0_dp
        do i = 1, 13
            ek = 0.5_dp * (y_out(1,i)**2 + y_out(2,i)**2)
            maxdrift = max(maxdrift, abs(ek - e0))
        end do
        if (maxdrift > 1.0e-6_dp) then
            write (error_unit, "(a,es12.4)") &
                "FAIL check_harmonic_energy: maxdrift=", maxdrift
            nfail = nfail + 1
        end if
    end subroutine check_harmonic_energy

    subroutine check_single_point(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(1)
        logical  :: ok

        t_eval(1) = 2.0_dp
        problem%rhs  => rhs_decay
        problem%y0   = [3.0_dp]

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "FAIL check_single_point: code=", status%code
            nfail = nfail + 1
            return
        end if
        ok = size(y_out, 2) == 1
        if (ok) ok = abs(y_out(1, 1) - 3.0_dp) <= 0.0_dp
        if (.not. ok) then
            write (error_unit, "(a)") "FAIL check_single_point: y0 not preserved"
            nfail = nfail + 1
        end if
    end subroutine check_single_point

    subroutine check_nonmonotone(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(3)

        ! t_eval is not monotone (goes forward then back).
        t_eval = [0.0_dp, 1.0_dp, 0.5_dp]
        problem%rhs => rhs_decay
        problem%y0  = [1.0_dp]

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_nonmonotone: expected DOMAIN_ERROR, got code=", status%code
            nfail = nfail + 1
        end if
    end subroutine check_nonmonotone

    subroutine check_null_rhs(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: workspace
        type(fortnum_status_t) :: status
        real(dp), allocatable  :: y_out(:,:)
        real(dp) :: t_eval(2)

        t_eval = [0.0_dp, 1.0_dp]
        problem%y0 = [1.0_dp]
        ! rhs stays null (default-initialized).

        call ode_at(problem, t_eval, workspace, y_out, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a,i0)") &
                "FAIL check_null_rhs: expected DOMAIN_ERROR, got code=", status%code
            nfail = nfail + 1
        end if
    end subroutine check_null_rhs

end program test_fortnum_ode_wrapper
