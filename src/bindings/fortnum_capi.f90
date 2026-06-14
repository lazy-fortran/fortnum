module fortnum_capi
    ! C ABI for fortnum. bind(c, name="fortnum_*") wrappers let C and C++
    ! consumers (KAMEL KiLCA + QL-Balance, MEPHIT C sources) link fortnum in
    ! place of the external backend. Each wrapper forwards to a pure Fortran routine and
    ! translates the argument convention: C double <-> real(dp), C int <->
    ! integer, c_double_complex <-> complex(dp), assumed-size C arrays by
    ! reference, and the fortnum_status code as the integer return value.
    !
    ! Callback-driven routines (quadrature, root finding, multiroot, central
    ! difference, ODE) take a C function pointer plus an opaque void* context.
    ! A stack-local carrier type passes the (c_funptr, c_ptr) pair through the
    ! Fortran routine's class(*) ctx slot; a host-associated bridge calls back
    ! through c_f_procpointer. No module-level state, so two callers never race.
    !
    ! No __enzyme_* symbol is exposed here (ad.md sec 2): the C surface is the
    ! primal foo_* routines only.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_double_complex, &
        c_funptr, c_ptr, c_f_procpointer, c_null_ptr

    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR

    use fortnum_special_bessel, only: bessel_in, bessel_in_array, bessel_kn
    use fortnum_special_dawson, only: dawson
    use fortnum_special_gamma, only: gamma_lower, gamma_reg_p
    use fortnum_special_erf_cbind, only: fortnum_erf, fortnum_erfc
    use fortnum_special_complex_bessel, only: bessel_j_complex, &
        bessel_i_complex, bessel_i_complex_array, &
        bessel_k_complex, bessel_k_complex_array
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1, hyperg_1f1_a1

    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    use fortnum_levin, only: levin_u_accel
    use fortnum_integrate, only: integrate_qag, integrate_qags, &
        integrate_qagp, integrate_qagiu, &
        integrate_workspace_t, integrate_epstab_t, integrate_result_t

    use fortnum_roots, only: root_brent
    use fortnum_multiroot, only: multiroot_hybrid, deriv_central, argsort

    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t
    use fortnum_ode_dop853, only: ode_integrate_dop, ode_solve_dop

    implicit none
    private

    ! C callback ABIs the consumers pass.
    ! Declared pure so the bridge can satisfy fortnum_roots' pure root_fn_t
    ! contract; the math callbacks fortnum exercises are referentially
    ! transparent. Impure-where-impure-allowed (quadrature) still accepts it.
    abstract interface
        pure function c_scalar_fn(x, ctx) result(fx) bind(c)
            import :: c_double, c_ptr
            real(c_double), value :: x
            type(c_ptr),    value :: ctx
            real(c_double)        :: fx
        end function c_scalar_fn
    end interface

    abstract interface
        subroutine c_vector_fn(n, x, f, ctx) bind(c)
            import :: c_int, c_double, c_ptr
            integer(c_int), value :: n
            real(c_double), intent(in)  :: x(*)
            real(c_double), intent(out) :: f(*)
            type(c_ptr),    value :: ctx
        end subroutine c_vector_fn
    end interface

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

    ! Stack-local carrier threading a C callback + opaque ctx through a Fortran
    ! routine's class(*) ctx slot. Holds the C function pointer and the user
    ! void*; the bridge below re-derives the procedure and forwards the call.
    type :: scalar_bridge_t
        type(c_funptr) :: fn
        type(c_ptr)    :: ctx
    end type scalar_bridge_t

    type :: vector_bridge_t
        type(c_funptr) :: fn
        type(c_ptr)    :: ctx
        integer        :: n
    end type vector_bridge_t

    public :: fortnum_bessel_in, fortnum_bessel_in_array, fortnum_bessel_kn
    public :: fortnum_dawson, fortnum_gamma_lower, fortnum_gamma_reg_p
    public :: capi_erf, capi_erfc
    public :: fortnum_bessel_j_complex
    public :: fortnum_bessel_i_complex, fortnum_bessel_i_complex_array
    public :: fortnum_bessel_k_complex, fortnum_bessel_k_complex_array
    public :: fortnum_hyperg_1f1, fortnum_hyperg_1f1_a1
    public :: fortnum_integrate_qag, fortnum_integrate_qags
    public :: fortnum_integrate_qagp, fortnum_integrate_qagiu
    public :: fortnum_gauss_legendre, fortnum_gauss_legendre_ab
    public :: fortnum_levin_u_accel, fortnum_multiroot_hybrid
    public :: fortnum_root_brent, fortnum_deriv_central, fortnum_argsort
    public :: fortnum_ode_integrate_dop, fortnum_ode_solve_dop

contains

    ! ------------------------------------------------------------------ special

    function fortnum_bessel_in(n, x) result(value) bind(c, name="fortnum_bessel_in")
        integer(c_int), value :: n
        real(c_double), value :: x
        real(c_double)        :: value
        value = bessel_in(int(n), real(x, dp))
    end function fortnum_bessel_in

    subroutine fortnum_bessel_in_array(nmax, x, values) &
            bind(c, name="fortnum_bessel_in_array")
        integer(c_int), value       :: nmax
        real(c_double), value       :: x
        real(c_double), intent(out) :: values(0:nmax)
        call bessel_in_array(int(nmax), real(x, dp), values)
    end subroutine fortnum_bessel_in_array

    function fortnum_bessel_kn(n, x) result(value) bind(c, name="fortnum_bessel_kn")
        integer(c_int), value :: n
        real(c_double), value :: x
        real(c_double)        :: value
        value = bessel_kn(int(n), real(x, dp))
    end function fortnum_bessel_kn

    function fortnum_dawson(x) result(value) bind(c, name="fortnum_dawson")
        real(c_double), value :: x
        real(c_double)        :: value
        value = dawson(real(x, dp))
    end function fortnum_dawson

    function fortnum_gamma_lower(a, x) result(value) &
            bind(c, name="fortnum_gamma_lower")
        real(c_double), value :: a, x
        real(c_double)        :: value
        value = gamma_lower(real(a, dp), real(x, dp))
    end function fortnum_gamma_lower

    function fortnum_gamma_reg_p(a, x) result(value) &
            bind(c, name="fortnum_gamma_reg_p")
        real(c_double), value :: a, x
        real(c_double)        :: value
        value = gamma_reg_p(real(a, dp), real(x, dp))
    end function fortnum_gamma_reg_p

    function capi_erf(x) result(value) bind(c, name="fortnum_erf")
        real(c_double), value :: x
        real(c_double)        :: value
        value = fortnum_erf(real(x, dp))
    end function capi_erf

    function capi_erfc(x) result(value) bind(c, name="fortnum_erfc")
        real(c_double), value :: x
        real(c_double)        :: value
        value = fortnum_erfc(real(x, dp))
    end function capi_erfc

    ! --------------------------------------------------------- complex Bessel

    function fortnum_bessel_j_complex(order, z, result) result(code) &
            bind(c, name="fortnum_bessel_j_complex")
        integer(c_int),         value       :: order
        complex(c_double_complex), value    :: z
        complex(c_double_complex), intent(out) :: result
        integer(c_int)         :: code
        complex(dp)            :: r
        type(fortnum_status_t) :: status
        call bessel_j_complex(int(order), cmplx(z, kind=dp), r, status)
        result = r
        code = int(status%code, c_int)
    end function fortnum_bessel_j_complex

    function fortnum_bessel_i_complex(order, z, scaled, result) result(code) &
            bind(c, name="fortnum_bessel_i_complex")
        integer(c_int),         value       :: order
        complex(c_double_complex), value    :: z
        integer(c_int),         value       :: scaled
        complex(c_double_complex), intent(out) :: result
        integer(c_int)         :: code
        complex(dp)            :: r
        type(fortnum_status_t) :: status
        call bessel_i_complex(int(order), cmplx(z, kind=dp), scaled /= 0, r, &
            status)
        result = r
        code = int(status%code, c_int)
    end function fortnum_bessel_i_complex

    function fortnum_bessel_i_complex_array(order0, nseq, z, scaled, result) &
            result(code) bind(c, name="fortnum_bessel_i_complex_array")
        integer(c_int),         value       :: order0, nseq
        complex(c_double_complex), value    :: z
        integer(c_int),         value       :: scaled
        complex(c_double_complex), intent(out) :: result(nseq)
        integer(c_int)         :: code
        complex(dp)            :: seq(nseq)
        type(fortnum_status_t) :: status
        call bessel_i_complex_array(int(order0), int(nseq), cmplx(z, kind=dp), &
            scaled /= 0, seq, status)
        result = seq
        code = int(status%code, c_int)
    end function fortnum_bessel_i_complex_array

    function fortnum_bessel_k_complex(order, z, scaled, result) result(code) &
            bind(c, name="fortnum_bessel_k_complex")
        integer(c_int),         value       :: order
        complex(c_double_complex), value    :: z
        integer(c_int),         value       :: scaled
        complex(c_double_complex), intent(out) :: result
        integer(c_int)         :: code
        complex(dp)            :: r
        type(fortnum_status_t) :: status
        call bessel_k_complex(int(order), cmplx(z, kind=dp), scaled /= 0, r, &
            status)
        result = r
        code = int(status%code, c_int)
    end function fortnum_bessel_k_complex

    function fortnum_bessel_k_complex_array(order0, nseq, z, scaled, result) &
            result(code) bind(c, name="fortnum_bessel_k_complex_array")
        integer(c_int),         value       :: order0, nseq
        complex(c_double_complex), value    :: z
        integer(c_int),         value       :: scaled
        complex(c_double_complex), intent(out) :: result(nseq)
        integer(c_int)         :: code
        complex(dp)            :: seq(nseq)
        type(fortnum_status_t) :: status
        call bessel_k_complex_array(int(order0), int(nseq), cmplx(z, kind=dp), &
            scaled /= 0, seq, status)
        result = seq
        code = int(status%code, c_int)
    end function fortnum_bessel_k_complex_array

    ! ----------------------------------------------------------- 1F1 (Kummer)

    function fortnum_hyperg_1f1(a, b, z, result) result(code) &
            bind(c, name="fortnum_hyperg_1f1")
        complex(c_double_complex), value    :: a, b, z
        complex(c_double_complex), intent(out) :: result
        integer(c_int)         :: code
        complex(dp)            :: r
        type(fortnum_status_t) :: status
        call hyperg_1f1(cmplx(a, kind=dp), cmplx(b, kind=dp), &
            cmplx(z, kind=dp), r, status)
        result = r
        code = int(status%code, c_int)
    end function fortnum_hyperg_1f1

    function fortnum_hyperg_1f1_a1(b, z, result) result(code) &
            bind(c, name="fortnum_hyperg_1f1_a1")
        complex(c_double_complex), value    :: b, z
        complex(c_double_complex), intent(out) :: result
        integer(c_int)         :: code
        complex(dp)            :: r
        type(fortnum_status_t) :: status
        call hyperg_1f1_a1(cmplx(b, kind=dp), cmplx(z, kind=dp), r, status)
        result = r
        code = int(status%code, c_int)
    end function fortnum_hyperg_1f1_a1

    ! ---------------------------------------------------------- fixed quadrature

    subroutine fortnum_gauss_legendre(n, x, w) &
            bind(c, name="fortnum_gauss_legendre")
        integer(c_int), value       :: n
        real(c_double), intent(out) :: x(n), w(n)
        call gauss_legendre(int(n), x, w)
    end subroutine fortnum_gauss_legendre

    subroutine fortnum_gauss_legendre_ab(n, a, b, x, w) &
            bind(c, name="fortnum_gauss_legendre_ab")
        integer(c_int), value       :: n
        real(c_double), value       :: a, b
        real(c_double), intent(out) :: x(n), w(n)
        call gauss_legendre_ab(int(n), real(a, dp), real(b, dp), x, w)
    end subroutine fortnum_gauss_legendre_ab

    function fortnum_levin_u_accel(terms, n, sum_accel, abserr) result(code) &
            bind(c, name="fortnum_levin_u_accel")
        integer(c_int), value       :: n
        real(c_double), intent(in)  :: terms(n)
        real(c_double), intent(out) :: sum_accel, abserr
        integer(c_int)         :: code
        type(fortnum_status_t) :: status
        call levin_u_accel(real(terms, dp), int(n), sum_accel, abserr, status)
        code = int(status%code, c_int)
    end function fortnum_levin_u_accel

    ! ------------------------------------------------------- adaptive quadrature

    function fortnum_integrate_qag(f, a, b, epsabs, epsrel, key, value, abserr) &
            result(code) bind(c, name="fortnum_integrate_qag")
        type(c_funptr), value       :: f
        real(c_double), value       :: a, b, epsabs, epsrel
        integer(c_int), value       :: key
        real(c_double), intent(out) :: value, abserr
        integer(c_int) :: code
        code = run_finite_integral(0, f, a, b, epsabs, epsrel, key, value, abserr)
    end function fortnum_integrate_qag

    function fortnum_integrate_qags(f, a, b, epsabs, epsrel, value, abserr) &
            result(code) bind(c, name="fortnum_integrate_qags")
        type(c_funptr), value       :: f
        real(c_double), value       :: a, b, epsabs, epsrel
        real(c_double), intent(out) :: value, abserr
        integer(c_int) :: code
        code = run_finite_integral(1, f, a, b, epsabs, epsrel, 21_c_int, value, &
            abserr)
    end function fortnum_integrate_qags

    function fortnum_integrate_qagp(f, a, b, points, npts, epsabs, epsrel, &
            value, abserr) result(code) bind(c, name="fortnum_integrate_qagp")
        type(c_funptr), value       :: f
        real(c_double), value       :: a, b, epsabs, epsrel
        integer(c_int), value       :: npts
        real(c_double), intent(in)  :: points(npts)
        real(c_double), intent(out) :: value, abserr
        integer(c_int) :: code
        type(scalar_bridge_t), target :: br
        type(integrate_workspace_t)   :: ws
        type(integrate_epstab_t)      :: eps
        type(integrate_result_t)      :: res
        type(fortnum_status_t)        :: status
        br%fn  = f
        br%ctx = c_null_ptr
        call integrate_qagp(integrand, real(a, dp), real(b, dp), &
            real(points, dp), real(epsabs, dp), real(epsrel, dp), ws, eps, res, &
            status, ctx=br)
        value  = res%value
        abserr = res%abserr
        code   = int(status%code, c_int)
    end function fortnum_integrate_qagp

    function fortnum_integrate_qagiu(f, bound, inf, epsabs, epsrel, value, &
            abserr) result(code) bind(c, name="fortnum_integrate_qagiu")
        type(c_funptr), value       :: f
        real(c_double), value       :: bound, epsabs, epsrel
        integer(c_int), value       :: inf
        real(c_double), intent(out) :: value, abserr
        integer(c_int) :: code
        type(scalar_bridge_t), target :: br
        type(integrate_workspace_t)   :: ws
        type(integrate_epstab_t)      :: eps
        type(integrate_result_t)      :: res
        type(fortnum_status_t)        :: status
        br%fn  = f
        br%ctx = c_null_ptr
        call integrate_qagiu(integrand, real(bound, dp), int(inf), &
            real(epsabs, dp), real(epsrel, dp), ws, eps, res, status, ctx=br)
        value  = res%value
        abserr = res%abserr
        code   = int(status%code, c_int)
    end function fortnum_integrate_qagiu

    ! Shared QAG/QAGS finite-interval driver. mode 0 = QAG, 1 = QAGS.
    function run_finite_integral(mode, f, a, b, epsabs, epsrel, key, value, &
            abserr) result(code)
        integer,        intent(in)  :: mode
        type(c_funptr), intent(in)  :: f
        real(c_double), intent(in)  :: a, b, epsabs, epsrel
        integer(c_int), intent(in)  :: key
        real(c_double), intent(out) :: value, abserr
        integer(c_int) :: code
        type(scalar_bridge_t), target :: br
        type(integrate_workspace_t)   :: ws
        type(integrate_epstab_t)      :: eps
        type(integrate_result_t)      :: res
        type(fortnum_status_t)        :: status
        br%fn  = f
        br%ctx = c_null_ptr
        if (mode == 0) then
            call integrate_qag(integrand, real(a, dp), real(b, dp), &
                real(epsabs, dp), real(epsrel, dp), ws, res, status, &
                key=int(key), ctx=br)
        else
            call integrate_qags(integrand, real(a, dp), real(b, dp), &
                real(epsabs, dp), real(epsrel, dp), ws, eps, res, status, ctx=br)
        end if
        value  = res%value
        abserr = res%abserr
        code   = int(status%code, c_int)
    end function run_finite_integral

    ! Bridge for the integrate_integrand_t signature; recovers the C callback
    ! and its void* from the scalar_bridge_t threaded through ctx.
    function integrand(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        procedure(c_scalar_fn), pointer :: cf
        fx = 0.0_dp
        if (.not. present(ctx)) return
        select type (ctx)
        type is (scalar_bridge_t)
            call c_f_procpointer(ctx%fn, cf)
            fx = cf(real(x, c_double), ctx%ctx)
        end select
    end function integrand

    ! ---------------------------------------------------------------- roots

    function fortnum_root_brent(f, a, b, xtol, ftol, max_iter, root) &
            result(code) bind(c, name="fortnum_root_brent")
        type(c_funptr), value       :: f
        real(c_double), value       :: a, b, xtol, ftol
        integer(c_int), value       :: max_iter
        real(c_double), intent(out) :: root
        integer(c_int) :: code
        type(scalar_bridge_t), target :: br
        type(fortnum_status_t)        :: status
        real(dp)                      :: x
        procedure(c_scalar_fn), pointer :: cf
        br%fn  = f
        br%ctx = c_null_ptr
        ! Resolve the C procedure pointer here (c_f_procpointer is impure); the
        ! pure root_fn that root_brent demands then only calls through cf.
        call c_f_procpointer(br%fn, cf)
        call root_brent(root_fn, real(a, dp), real(b, dp), x, status, &
            xtol=real(xtol, dp), ftol=real(ftol, dp), max_iter=int(max_iter))
        root = x
        code = int(status%code, c_int)
    contains
        pure function root_fn(xv) result(yv)
            real(dp), intent(in) :: xv
            real(dp) :: yv
            yv = cf(real(xv, c_double), br%ctx)
        end function root_fn
    end function fortnum_root_brent

    function fortnum_deriv_central(f, x, h, result, abserr) result(code) &
            bind(c, name="fortnum_deriv_central")
        type(c_funptr), value       :: f
        real(c_double), value       :: x, h
        real(c_double), intent(out) :: result, abserr
        integer(c_int) :: code
        type(scalar_bridge_t), target :: br
        type(fortnum_status_t)        :: status
        real(dp)                      :: r, e
        br%fn  = f
        br%ctx = c_null_ptr
        call deriv_central(deriv_fn, real(x, dp), real(h, dp), r, e, status, &
            ctx=br)
        result = r
        abserr = e
        code   = int(status%code, c_int)
    contains
        function deriv_fn(xv, ctx) result(yv)
            real(dp), intent(in) :: xv
            class(*), intent(in), optional :: ctx
            real(dp) :: yv
            procedure(c_scalar_fn), pointer :: cf
            call c_f_procpointer(br%fn, cf)
            yv = cf(real(xv, c_double), br%ctx)
        end function deriv_fn
    end function fortnum_deriv_central

    subroutine fortnum_argsort(x, n, perm) bind(c, name="fortnum_argsort")
        integer(c_int), value       :: n
        real(c_double), intent(in)  :: x(n)
        integer(c_int), intent(out) :: perm(n)
        integer :: p(n)
        call argsort(real(x, dp), p)
        ! argsort returns 1-based indices; expose 0-based for C consumers.
        perm = int(p - 1, c_int)
    end subroutine fortnum_argsort

    ! ------------------------------------------------------------ multiroot

    function fortnum_multiroot_hybrid(fdf, n, x0, xtol, ftol, max_iter, x) &
            result(code) bind(c, name="fortnum_multiroot_hybrid")
        type(c_funptr), value       :: fdf
        integer(c_int), value       :: n
        real(c_double), intent(in)  :: x0(n)
        real(c_double), value       :: xtol, ftol
        integer(c_int), value       :: max_iter
        real(c_double), intent(out) :: x(n)
        integer(c_int) :: code
        type(vector_bridge_t), target :: br
        type(fortnum_status_t)        :: status
        real(dp)                      :: xout(n)
        br%fn  = fdf
        br%ctx = c_null_ptr
        br%n   = int(n)
        call multiroot_hybrid(system, int(n), real(x0, dp), xout, status, &
            xtol=real(xtol, dp), ftol=real(ftol, dp), max_iter=int(max_iter), &
            ctx=br)
        x    = xout
        code = int(status%code, c_int)
    contains
        ! Builds the analytic Jacobian by central differences of the C residual
        ! callback, matching the KiLCA "_hybrids" usage on a _hybrid signature.
        subroutine system(xv, fv, jac, ctx)
            real(dp), intent(in)  :: xv(:)
            real(dp), intent(out) :: fv(:)
            real(dp), intent(out) :: jac(:, :)
            class(*), intent(in), optional :: ctx
            real(dp), parameter :: cube_root_eps = 6.0554544523933429e-6_dp
            real(dp) :: xp(size(xv)), xm(size(xv)), fp(size(xv)), fm(size(xv)), hh
            integer :: j
            call eval_residual(xv, fv)
            do j = 1, size(xv)
                hh = cube_root_eps*max(abs(xv(j)), 1.0_dp)
                xp = xv; xm = xv
                xp(j) = xv(j) + hh
                xm(j) = xv(j) - hh
                call eval_residual(xp, fp)
                call eval_residual(xm, fm)
                jac(:, j) = (fp - fm)/(2.0_dp*hh)
            end do
        end subroutine system

        subroutine eval_residual(xv, fv)
            real(dp), intent(in)  :: xv(:)
            real(dp), intent(out) :: fv(:)
            procedure(c_vector_fn), pointer :: cf
            real(c_double) :: xc(size(xv)), fc(size(xv))
            call c_f_procpointer(br%fn, cf)
            xc = real(xv, c_double)
            call cf(int(size(xv), c_int), xc, fc, br%ctx)
            fv = real(fc, dp)
        end subroutine eval_residual
    end function fortnum_multiroot_hybrid

    ! ------------------------------------------------------------------- ode

    function fortnum_ode_integrate_dop(rhs, neq, t0, t1, y0, rtol, atol, &
            max_steps, npts_cap, t_out, y_out, npts) result(code) &
            bind(c, name="fortnum_ode_integrate_dop")
        type(c_funptr), value       :: rhs
        integer(c_int), value       :: neq, max_steps, npts_cap
        real(c_double), value       :: t0, t1, rtol, atol
        real(c_double), intent(in)  :: y0(neq)
        real(c_double), intent(out) :: t_out(npts_cap)
        real(c_double), intent(out) :: y_out(neq, npts_cap)
        integer(c_int), intent(out) :: npts
        integer(c_int) :: code
        code = run_ode_dop(rhs, neq, t0, t1, y0, rtol, atol, max_steps, &
            npts_cap, t_out, y_out, npts)
    end function fortnum_ode_integrate_dop

    ! ode_solve_dop is the allocatable-output flat call; the C ABI cannot return
    ! an allocatable, so this fixed-buffer entry caps the recorded mesh at
    ! npts_cap and reports the true point count in npts. Both fortnum entry
    ! points (ode_integrate_dop / ode_solve_dop) run the same integrator; the C
    ! surface exposes the buffer-bounded one.
    function fortnum_ode_solve_dop(rhs, neq, t0, t1, y0, rtol, atol, &
            npts_cap, t_out, y_out, npts) result(code) &
            bind(c, name="fortnum_ode_solve_dop")
        type(c_funptr), value       :: rhs
        integer(c_int), value       :: neq, npts_cap
        real(c_double), value       :: t0, t1, rtol, atol
        real(c_double), intent(in)  :: y0(neq)
        real(c_double), intent(out) :: t_out(npts_cap)
        real(c_double), intent(out) :: y_out(neq, npts_cap)
        integer(c_int), intent(out) :: npts
        integer(c_int) :: code
        code = run_ode_dop(rhs, neq, t0, t1, y0, rtol, atol, 100000_c_int, &
            npts_cap, t_out, y_out, npts)
    end function fortnum_ode_solve_dop

    ! Shared dop853 driver for both C ode entry points. Bridges the C rhs and
    ! copies the recorded mesh into the caller's fixed buffers (capped).
    function run_ode_dop(rhs, neq, t0, t1, y0, rtol, atol, max_steps, &
            npts_cap, t_out, y_out, npts) result(code)
        type(c_funptr), intent(in)  :: rhs
        integer(c_int), intent(in)  :: neq, max_steps, npts_cap
        real(c_double), intent(in)  :: t0, t1, rtol, atol
        real(c_double), intent(in)  :: y0(neq)
        real(c_double), intent(out) :: t_out(npts_cap)
        real(c_double), intent(out) :: y_out(neq, npts_cap)
        integer(c_int), intent(out) :: npts
        integer(c_int) :: code
        type(vector_bridge_t), target :: br
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: status
        integer :: nout
        br%fn  = rhs
        br%ctx = c_null_ptr
        br%n   = int(neq)
        problem%rhs => rhs_bridge
        problem%t0  = real(t0, dp)
        problem%t1  = real(t1, dp)
        problem%y0  = real(y0, dp)
        problem%rtol = real(rtol, dp)
        problem%atol = real(atol, dp)
        problem%max_steps = int(max_steps)
        call ode_integrate_dop(problem, ws, sol, status)
        nout = min(sol%nsteps + 1, int(npts_cap))
        npts = int(nout, c_int)
        if (nout >= 1 .and. allocated(sol%t)) then
            t_out(1:nout)    = sol%t(1:nout)
            y_out(:, 1:nout) = sol%y(:, 1:nout)
        end if
        code = int(status%code, c_int)
    contains
        subroutine rhs_bridge(t, y, dydt, ctx)
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            procedure(c_ode_rhs), pointer :: cf
            real(c_double) :: yc(size(y)), dc(size(y))
            call c_f_procpointer(br%fn, cf)
            yc = real(y, c_double)
            call cf(real(t, c_double), int(size(y), c_int), yc, dc, br%ctx)
            dydt = real(dc, dp)
        end subroutine rhs_bridge
    end function run_ode_dop

end module fortnum_capi
