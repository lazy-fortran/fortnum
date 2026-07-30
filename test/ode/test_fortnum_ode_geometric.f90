program test_fortnum_ode_geometric
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ode_geometric, only: &
        geometric_strang_step, geometric_yoshida_fourth_step
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, status_set
    implicit none

    type(fortnum_status_t) :: status
    real(dp) :: basis_q(2), basis_p(2), determinant, energy, energy_error
    real(dp) :: initial(2), state(2)
    integer :: step

    basis_q = [1.0_dp, 0.0_dp]
    basis_p = [0.0_dp, 1.0_dp]
    call geometric_strang_step( &
        basis_q, 0.2_dp, drift, kick, status)
    call geometric_strang_step( &
        basis_p, 0.2_dp, drift, kick, status)
    determinant = basis_q(1)*basis_p(2) - basis_p(1)*basis_q(2)
    call check(abs(determinant - 1.0_dp) < 4.0e-15_dp, &
        "Strang split preserves canonical phase area")

    initial = [0.7_dp, -0.4_dp]
    state = initial
    call geometric_yoshida_fourth_step( &
        state, 0.17_dp, drift, kick, status)
    call geometric_yoshida_fourth_step( &
        state, -0.17_dp, drift, kick, status)
    call check(maxval(abs(state - initial)) < 2.0e-14_dp, &
        "Symmetric fourth-order composition is time reversible")

    state = initial
    energy = 0.5_dp*dot_product(state, state)
    energy_error = 0.0_dp
    do step = 1, 20000
        call geometric_strang_step(state, 0.03_dp, drift, kick, status)
        energy_error = max(energy_error, &
            abs(0.5_dp*dot_product(state, state) - energy))
    end do
    call check(status%code == FORTNUM_OK .and. energy_error < 8.0e-5_dp, &
        "Symplectic split has bounded long-time oscillator energy error")

contains

    subroutine drift(values, time_step, local_status, context)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(in) :: time_step
        type(fortnum_status_t), intent(out) :: local_status
        class(*), intent(in), optional :: context

        values(1) = values(1) + time_step*values(2)
        call status_set(local_status, FORTNUM_OK, "")
        associate(unused => context)
        end associate
    end subroutine drift

    subroutine kick(values, time_step, local_status, context)
        real(dp), intent(inout) :: values(:)
        real(dp), intent(in) :: time_step
        type(fortnum_status_t), intent(out) :: local_status
        class(*), intent(in), optional :: context

        values(2) = values(2) - time_step*values(1)
        call status_set(local_status, FORTNUM_OK, "")
        associate(unused => context)
        end associate
    end subroutine kick

    subroutine check(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (condition) return
        write(error_unit, "(a)") "FAIL: "//message
        stop 1
    end subroutine check

end program test_fortnum_ode_geometric
