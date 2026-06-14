module fortnum_capi_bspline
    ! C ABI for the B-spline basis (fortnum_bspline), matching the GSL bspline
    ! workflow NEO-2 uses (collop_bspline.f90). The stateful workspace is opaque
    ! to C: create returns a void* handle, the eval/deriv calls take it back, and
    ! destroy frees it. Handle lifetime is caller-managed; no module-level state,
    ! so distinct handles never race.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_loc, &
        c_f_pointer, c_null_ptr, c_associated

    use fortnum_status, only: fortnum_status_t
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_basis, bspline_eval_deriv, &
        bspline_span_index

    implicit none
    private

    public :: fortnum_bspline_create, fortnum_bspline_set_knots
    public :: fortnum_bspline_ncoef, fortnum_bspline_eval_basis
    public :: fortnum_bspline_eval_deriv, fortnum_bspline_span_index
    public :: fortnum_bspline_destroy

contains

    ! Allocate a workspace for the given order and breakpoint count; returns an
    ! opaque handle (c_null_ptr on a domain error). The handle owns the Fortran
    ! workspace until fortnum_bspline_destroy frees it.
    function fortnum_bspline_create(order, nbreak) result(handle) &
            bind(c, name="fortnum_bspline_create")
        integer(c_int), value :: order, nbreak
        type(c_ptr)           :: handle
        type(bspline_workspace_t), pointer :: ws
        type(fortnum_status_t)             :: status
        allocate (ws)
        call bspline_init(ws, int(order), int(nbreak), status)
        if (status%code /= 0) then
            deallocate (ws)
            handle = c_null_ptr
            return
        end if
        handle = c_loc(ws)
    end function fortnum_bspline_create

    ! Build the clamped knot vector from nbreak strictly increasing breakpoints.
    function fortnum_bspline_set_knots(handle, nbreak, breakpts) result(code) &
            bind(c, name="fortnum_bspline_set_knots")
        type(c_ptr),    value      :: handle
        integer(c_int), value      :: nbreak
        real(c_double), intent(in) :: breakpts(nbreak)
        integer(c_int) :: code
        type(bspline_workspace_t), pointer :: ws
        type(fortnum_status_t)             :: status
        if (.not. c_associated(handle)) then
            code = 1
            return
        end if
        call c_f_pointer(handle, ws)
        call bspline_set_knots(ws, real(breakpts, dp), status)
        code = int(status%code, c_int)
    end function fortnum_bspline_set_knots

    ! Number of basis functions (coefficients) = nbreak + order - 2.
    function fortnum_bspline_ncoef(handle) result(ncoef) &
            bind(c, name="fortnum_bspline_ncoef")
        type(c_ptr), value :: handle
        integer(c_int)     :: ncoef
        type(bspline_workspace_t), pointer :: ws
        ncoef = 0
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ws)
        ncoef = int(ws%ncoef, c_int)
    end function fortnum_bspline_ncoef

    ! Evaluate all ncoef basis functions B_{i,k}(x) into values[ncoef].
    function fortnum_bspline_eval_basis(handle, x, ncoef, values) result(code) &
            bind(c, name="fortnum_bspline_eval_basis")
        type(c_ptr),    value       :: handle
        real(c_double), value       :: x
        integer(c_int), value       :: ncoef
        real(c_double), intent(out) :: values(ncoef)
        integer(c_int) :: code
        type(bspline_workspace_t), pointer :: ws
        type(fortnum_status_t)             :: status
        if (.not. c_associated(handle)) then
            code = 1
            return
        end if
        call c_f_pointer(handle, ws)
        call bspline_eval_basis(ws, real(x, dp), values, status)
        code = int(status%code, c_int)
    end function fortnum_bspline_eval_basis

    ! Evaluate basis functions and derivatives up to order nderiv. dvalues is
    ! row-major (nderiv+1) x ncoef: dvalues[d*ncoef + i] is the d-th derivative
    ! of basis function i (d = 0 is the value).
    function fortnum_bspline_eval_deriv(handle, x, nderiv, ncoef, dvalues) &
            result(code) bind(c, name="fortnum_bspline_eval_deriv")
        type(c_ptr),    value       :: handle
        real(c_double), value       :: x
        integer(c_int), value       :: nderiv, ncoef
        ! Declared (ncoef, 0:nderiv) so element (i,d) sits at the C row-major
        ! offset d*ncoef + (i-1), matching the documented dvalues[d*ncoef+i].
        real(c_double), intent(out) :: dvalues(ncoef, 0:nderiv)
        integer(c_int) :: code
        type(bspline_workspace_t), pointer :: ws
        type(fortnum_status_t)             :: status
        real(dp), allocatable              :: tmp(:, :)
        integer :: d, i
        if (.not. c_associated(handle)) then
            code = 1
            return
        end if
        call c_f_pointer(handle, ws)
        allocate (tmp(0:nderiv, ncoef))
        call bspline_eval_deriv(ws, real(x, dp), int(nderiv), tmp, status)
        do i = 1, ncoef
            do d = 0, nderiv
                dvalues(i, d) = tmp(d, i)
            end do
        end do
        code = int(status%code, c_int)
    end function fortnum_bspline_eval_deriv

    ! 1-based knot span index containing x (Piegl-Tiller FindSpan).
    function fortnum_bspline_span_index(handle, x) result(span) &
            bind(c, name="fortnum_bspline_span_index")
        type(c_ptr),    value :: handle
        real(c_double), value :: x
        integer(c_int)        :: span
        type(bspline_workspace_t), pointer :: ws
        span = 0
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ws)
        span = int(bspline_span_index(ws, real(x, dp)), c_int)
    end function fortnum_bspline_span_index

    ! Free the workspace behind the handle.
    subroutine fortnum_bspline_destroy(handle) &
            bind(c, name="fortnum_bspline_destroy")
        type(c_ptr), value :: handle
        type(bspline_workspace_t), pointer :: ws
        if (.not. c_associated(handle)) return
        call c_f_pointer(handle, ws)
        deallocate (ws)
    end subroutine fortnum_bspline_destroy

end module fortnum_capi_bspline
