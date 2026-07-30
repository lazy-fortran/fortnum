module fortnum_ode_geometric
    !! Symmetric composition of exactly solvable geometric subflows.
    !!
    !! If each supplied flow is symplectic (or Poisson), Strang composition
    !! and the Yoshida triple jump preserve that structure. This is the
    !! propagator pattern used by Hamiltonian-splitting plasma codes such as
    !! STRUPHY. Negative substeps in the fourth-order method make it suitable
    !! for reversible Hamiltonian flows, not dissipative operators.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: &
        fortnum_status_t, FORTNUM_DOMAIN_ERROR, FORTNUM_OK, status_set
    implicit none
    private

    public :: geometric_flow_t
    public :: geometric_strang_step
    public :: geometric_yoshida_fourth_step

    abstract interface
        subroutine geometric_flow_t(state, time_step, status, context)
            import :: dp, fortnum_status_t
            real(dp), intent(inout) :: state(:)
            real(dp), intent(in) :: time_step
            type(fortnum_status_t), intent(out) :: status
            class(*), intent(in), optional :: context
        end subroutine geometric_flow_t
    end interface

contains

    subroutine geometric_strang_step( &
            state, time_step, flow_a, flow_b, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: time_step
        procedure(geometric_flow_t) :: flow_a, flow_b
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context

        call status_set( &
            status, FORTNUM_DOMAIN_ERROR, &
            "Geometric Strang step requires a finite time step")
        if (.not. finite(time_step)) return
        call flow_a(state, 0.5_dp*time_step, status, context)
        if (status%code /= FORTNUM_OK) return
        call flow_b(state, time_step, status, context)
        if (status%code /= FORTNUM_OK) return
        call flow_a(state, 0.5_dp*time_step, status, context)
    end subroutine geometric_strang_step

    subroutine geometric_yoshida_fourth_step( &
            state, time_step, flow_a, flow_b, status, context)
        real(dp), intent(inout) :: state(:)
        real(dp), intent(in) :: time_step
        procedure(geometric_flow_t) :: flow_a, flow_b
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: context

        real(dp), parameter :: cube_root_two = 2.0_dp**(1.0_dp/3.0_dp)
        real(dp), parameter :: outer_weight = 1.0_dp/(2.0_dp - cube_root_two)
        real(dp), parameter :: middle_weight = &
            -cube_root_two/(2.0_dp - cube_root_two)

        call geometric_strang_step( &
            state, outer_weight*time_step, flow_a, flow_b, status, context)
        if (status%code /= FORTNUM_OK) return
        call geometric_strang_step( &
            state, middle_weight*time_step, flow_a, flow_b, status, context)
        if (status%code /= FORTNUM_OK) return
        call geometric_strang_step( &
            state, outer_weight*time_step, flow_a, flow_b, status, context)
    end subroutine geometric_yoshida_fourth_step

    pure logical function finite(value) result(is_finite)
        real(dp), intent(in) :: value

        is_finite = value == value .and. abs(value) <= huge(value)
    end function finite

end module fortnum_ode_geometric
