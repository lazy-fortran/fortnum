module fortnum_capi_rk8pd
    ! C ABI for the re-entrant adaptive rk8pd integrator (fortnum_ode_rk8pd).
    ! KiLCA's calc_back evolves its background ODE continuously across the radial
    ! grid; this surface mirrors a loop over an accepted-step evolve advance so the
    ! caller carries the adaptive step from one output abscissa to the next.
    !
    ! The evolve state is opaque to C: create returns a void* handle holding the
    ! re-entrant rk8pd_state_t plus the C rhs callback and its void* context;
    ! integrate_to advances the carried (t, y) to the next abscissa; destroy
    ! frees it. Handle lifetime is caller-managed; no module-level state, so
    ! distinct handles never race. The callback ABI matches fortnum_ode_rhs in
    ! fortnum.h (the same fortnum_capi c_ode_rhs signature).

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_funptr, c_ptr, &
        c_loc, c_f_pointer, c_f_procpointer, c_null_ptr, c_associated

    use fortnum_status, only: fortnum_status_t
    use fortnum_ode, only: ode_rhs_t
    use fortnum_ode_rk8pd, only: rk8pd_state_t, rk8pd_evolve_init, &
        rk8pd_evolve_apply

    implicit none
    private

    public :: fortnum_rk8pd_create, fortnum_rk8pd_integrate_to
    public :: fortnum_rk8pd_destroy

    ! C rhs callback ABI: dydt = f(t, y) with a leading abscissa, the equation
    ! count, and the opaque user context (matches fortnum_ode_rhs / c_ode_rhs).
    abstract interface
        subroutine c_ode_rhs(t, n, y, dydt, ctx) bind(c)
            import :: c_int, c_double, c_ptr
            real(c_double), value :: t
            integer(c_int), value :: n
            real(c_double), intent(in)  :: y(*)
            real(c_double), intent(out) :: dydt(*)
            type(c_ptr),    value :: ctx
        end subroutine c_ode_rhs
    end interface

    ! Heap object behind the void* handle: the re-entrant state, the equation
    ! count, the tolerances, and the C callback pair threaded into every step.
    type :: rk8pd_handle_t
        type(rk8pd_state_t) :: state
        integer             :: neq = 0
        real(dp)            :: eps_abs = 0.0_dp
        real(dp)            :: eps_rel = 0.0_dp
        integer             :: max_steps = 0
        type(c_funptr)      :: rhs
        type(c_ptr)         :: ctx
    end type rk8pd_handle_t

contains

    ! Allocate a re-entrant rk8pd evolve state for neq equations with first step
    ! magnitude h0 and the error-per-step tolerances (eps_abs, eps_rel) used by
    ! the y_new error-per-step control. rhs is the C callback, ctx its opaque user
    ! pointer, forwarded unchanged to every evaluation. The starting abscissa is
    ! carried by the caller through fortnum_rk8pd_integrate_to. Returns an opaque
    ! handle (c_null_ptr on a domain error); the handle owns the state until
    ! fortnum_rk8pd_destroy frees it.
    function fortnum_rk8pd_create(rhs, neq, h0, eps_abs, eps_rel, &
            max_steps, ctx) result(handle) bind(c, name="fortnum_rk8pd_create")
        type(c_funptr), value :: rhs
        integer(c_int), value :: neq, max_steps
        real(c_double), value :: h0, eps_abs, eps_rel
        type(c_ptr),    value :: ctx
        type(c_ptr)           :: handle
        type(rk8pd_handle_t), pointer :: h
        type(fortnum_status_t)        :: status
        allocate (h)
        call rk8pd_evolve_init(h%state, int(neq), real(h0, dp), status)
        if (status%code /= 0) then
            deallocate (h)
            handle = c_null_ptr
            return
        end if
        h%neq       = int(neq)
        h%eps_abs   = real(eps_abs, dp)
        h%eps_rel   = real(eps_rel, dp)
        h%max_steps = int(max_steps)
        h%rhs       = rhs
        h%ctx       = ctx
        handle = c_loc(h)
    end function fortnum_rk8pd_create

    ! Advance the carried solution from *t to t1, continuing the adaptive
    ! schedule (one an accepted-step evolve advance pass per call). On entry y[neq] holds
    ! the state at *t; on return *t == t1 and y holds the solution there, with
    ! the carried step ready for the next abscissa. Returns the fortnum status
    ! code (0 == FORTNUM_OK).
    function fortnum_rk8pd_integrate_to(handle, t, t1, y) result(code) &
            bind(c, name="fortnum_rk8pd_integrate_to")
        type(c_ptr),    value         :: handle
        real(c_double), intent(inout) :: t
        real(c_double), value         :: t1
        real(c_double), intent(inout) :: y(*)
        integer(c_int) :: code
        type(rk8pd_handle_t), pointer :: h
        type(fortnum_status_t)        :: status
        real(dp) :: tt
        real(dp), allocatable :: yy(:)
        integer  :: nfev
        if (.not. c_associated(handle)) then
            code = 1
            return
        end if
        call c_f_pointer(handle, h)
        allocate (yy(h%neq))
        yy = real(y(1:h%neq), dp)
        tt = real(t, dp)
        nfev = 0
        call rk8pd_evolve_apply(rhs_bridge, h%state, tt, real(t1, dp), yy, &
            h%eps_abs, h%eps_rel, h%max_steps, nfev, status)
        y(1:h%neq) = real(yy, c_double)
        t = real(tt, c_double)
        code = int(status%code, c_int)
    contains
        subroutine rhs_bridge(tb, yb, dydt, ctx)
            real(dp), intent(in)  :: tb
            real(dp), intent(in)  :: yb(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            procedure(c_ode_rhs), pointer :: cf
            real(c_double) :: yc(size(yb)), dc(size(yb))
            call c_f_procpointer(h%rhs, cf)
            yc = real(yb, c_double)
            call cf(real(tb, c_double), int(size(yb), c_int), yc, dc, h%ctx)
            dydt = real(dc, dp)
        end subroutine rhs_bridge
    end function fortnum_rk8pd_integrate_to

    ! Free the evolve state behind the handle.
    subroutine fortnum_rk8pd_destroy(handle) &
            bind(c, name="fortnum_rk8pd_destroy")
        type(c_ptr), value :: handle
        type(rk8pd_handle_t), pointer :: h
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, h)
        deallocate (h)
    end subroutine fortnum_rk8pd_destroy

end module fortnum_capi_rk8pd
