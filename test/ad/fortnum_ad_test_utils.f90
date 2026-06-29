module fortnum_ad_test_utils
    ! Reusable derivative-test machinery for fortnum. Every module that adds a
    ! derivative product (issue #40) checks it here against an independent
    ! estimate: central finite difference, the dot-product (adjoint) identity,
    ! the complex-step derivative, and an optional Enzyme-vs-analytic compare.
    !
    ! The callables a module supplies are vector maps f: R^n -> R^m and their
    ! products: a JVP applies J(x), a VJP applies J(x)^T. The harness does not
    ! know how a product is produced (analytic_rule, transparent, ...); it only
    ! checks that the product the module returns matches the reference within a
    ! stated tolerance. See docs/design/ad.md for the policy classes and naming.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    implicit none
    private

    ! ------------------------------------------------------------ public API
    public :: vector_fn_i ! primal map  y = f(x)
    public :: jvp_fn_i ! forward product  Jv = J(x) v
    public :: vjp_fn_i ! reverse product  Jtu = J(x)^T u
    public :: scalar_complex_fn_i ! analytic scalar f for complex-step

    public :: fd_jvp ! central-FD estimate of J(x) v
    public :: fd_jvp_step ! safeguarded step size for a base point
    public :: check_jvp_vs_fd ! forward-product check (JVP vs central FD)
    public :: dot_product_identity ! mandatory adjoint identity u.(Jv)=v.(J^T u)
    public :: complex_step_deriv ! df/dx via Im f(x+ih)/h
    public :: check_complex_step ! complex-step vs analytic-derivative compare
    public :: rel_err ! scaled error helper used by the checks
    public :: ad_status_t ! smoothness/branch status carrier
    public :: AD_SMOOTH, AD_NONSMOOTH
    public :: check_smoothness ! same-trace passes, branch change reports
    public :: enzyme_available ! .true. only in FORTNUM_ENABLE_ENZYME builds
    public :: check_enzyme_vs_analytic ! inert unless built with Enzyme

    ! Smoothness verdicts for a same-trace vs branch-change derivative check.
    integer, parameter :: AD_SMOOTH    = 0
    integer, parameter :: AD_NONSMOOTH = 1

    ! Outcome of a smoothness/event check. ok is .true. when the observed
    ! verdict matched what the caller asserted; verdict carries which case it
    ! actually saw.
    type :: ad_status_t
        logical :: ok      = .true.
        integer :: verdict = AD_SMOOTH
    end type ad_status_t

    ! ------------------------------------------------------------ callables
    !
    ! All shapes are contiguous explicit-shape arrays, matching the simplest
    ! layout the Enzyme path is tested against first (ad.md section 3).
    abstract interface
        ! Primal vector map y = f(x), x in R^n, y in R^m.
        subroutine vector_fn_i(x, y)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(out) :: y(:)
        end subroutine vector_fn_i

        ! Forward product jv = J(x) v, with v, jv conforming to x, y.
        subroutine jvp_fn_i(x, v, jv)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(in)  :: v(:)
            real(dp), intent(out) :: jv(:)
        end subroutine jvp_fn_i

        ! Reverse product jtu = J(x)^T u, with u in R^m, jtu in R^n.
        subroutine vjp_fn_i(x, u, jtu)
            import :: dp
            real(dp), intent(in)  :: x(:)
            real(dp), intent(in)  :: u(:)
            real(dp), intent(out) :: jtu(:)
        end subroutine vjp_fn_i

        ! Analytic scalar map accepting a complex argument, used so the
        ! complex-step estimate is exact to machine precision.
        pure function scalar_complex_fn_i(z) result(w)
            import :: dp
            complex(dp), intent(in) :: z
            complex(dp)             :: w
        end function scalar_complex_fn_i
    end interface

contains

    ! ---------------------------------------------------------------- errors

    ! Error scaled by the magnitude of the reference. Falls back to absolute
    ! error when the reference is near zero so a tiny true value does not
    ! inflate the ratio.
    pure function rel_err(got, want) result(e)
        real(dp), intent(in) :: got, want
        real(dp)             :: e
        real(dp) :: scale
        scale = max(abs(want), 1.0_dp)
        e = abs(got - want) / scale
    end function rel_err

    ! ------------------------------------------------------- finite difference

    ! Safeguarded central-difference step for a directional derivative at x
    ! along v. The classic central-FD optimum scales as eps^(1/3); scale it by
    ! the base-point magnitude so the perturbation stays meaningful for large
    ! or small x, and floor it so it never collapses to zero.
    pure function fd_jvp_step(x, v) result(h)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: v(:)
        real(dp)             :: h
        real(dp), parameter :: cube_root_eps = 6.0554544523933429e-6_dp
        real(dp) :: xscale, vscale
        xscale = max(maxval(abs(x)), 1.0_dp)
        vscale = max(maxval(abs(v)), tiny(1.0_dp))
        h = cube_root_eps * xscale / vscale
    end function fd_jvp_step

    ! Central-FD estimate of the forward product J(x) v:
    !   J(x) v ~ ( f(x + h v) - f(x - h v) ) / (2 h).
    ! Independent of any module JVP; check_jvp_vs_fd compares the two.
    subroutine fd_jvp(f, x, v, jv, h_in)
        procedure(vector_fn_i)         :: f
        real(dp), intent(in)           :: x(:)
        real(dp), intent(in)           :: v(:)
        real(dp), intent(out)          :: jv(:)
        real(dp), intent(in), optional :: h_in
        real(dp) :: h
        real(dp), allocatable :: yp(:), ym(:)

        h = fd_jvp_step(x, v)
        if (present(h_in)) h = h_in

        allocate (yp(size(jv)), ym(size(jv)))
        call f(x + h*v, yp)
        call f(x - h*v, ym)
        jv = (yp - ym) / (2.0_dp*h)
    end subroutine fd_jvp

    ! Forward-product check: a module JVP routine against the central-FD
    ! estimate at x along v. Returns .true. when every component agrees within
    ! tol (scaled error); otherwise reports the worst component on error_unit.
    function check_jvp_vs_fd(label, f, jvp, x, v, tol, h_in) result(ok)
        character(*), intent(in)       :: label
        procedure(vector_fn_i)         :: f
        procedure(jvp_fn_i)            :: jvp
        real(dp), intent(in)           :: x(:)
        real(dp), intent(in)           :: v(:)
        real(dp), intent(in)           :: tol
        real(dp), intent(in), optional :: h_in
        logical                        :: ok
        real(dp), allocatable :: jv_ad(:), jv_fd(:)
        real(dp) :: worst
        integer  :: i, m

        m = size(x) ! square-Jacobian default; callers size jv to m
        allocate (jv_ad(m), jv_fd(m))
        call jvp(x, v, jv_ad)
        if (present(h_in)) then
            call fd_jvp(f, x, v, jv_fd, h_in=h_in)
        else
            call fd_jvp(f, x, v, jv_fd)
        end if

        worst = 0.0_dp
        do i = 1, m
            worst = max(worst, rel_err(jv_ad(i), jv_fd(i)))
        end do
        ok = (worst <= tol)
        if (.not. ok) then
            write (error_unit, '(a,a,a,es12.4,a,es12.4,a)') &
                "FAIL [", label, "] jvp-vs-fd worst rel_err=", worst, &
                " tol=", tol, ""
        end if
    end function check_jvp_vs_fd

    ! ------------------------------------------------- dot-product identity

    ! Mandatory adjoint identity. For any u, v the forward and reverse products
    ! must satisfy  u . (J v) == v . (J^T u)  exactly up to rounding, because
    ! both equal u^T J v. A mismatch means JVP and VJP are not transposes, the
    ! single most common reverse-mode bug. Returns .true. on agreement within
    ! tol (error scaled by the magnitude of the bilinear form).
    function dot_product_identity(label, jvp, vjp, x, u, v, tol) result(ok)
        character(*), intent(in) :: label
        procedure(jvp_fn_i)      :: jvp
        procedure(vjp_fn_i)      :: vjp
        real(dp), intent(in)     :: x(:)
        real(dp), intent(in)     :: u(:) ! in R^m (output space)
        real(dp), intent(in)     :: v(:) ! in R^n (input space)
        real(dp), intent(in)     :: tol
        logical                  :: ok
        real(dp), allocatable :: jv(:), jtu(:)
        real(dp) :: lhs, rhs, e

        allocate (jv(size(u)), jtu(size(v)))
        call jvp(x, v, jv)
        call vjp(x, u, jtu)
        lhs = dot_product(u, jv) ! u^T (J v)
        rhs = dot_product(v, jtu) ! v^T (J^T u)
        e = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)
        ok = (e <= tol)
        if (.not. ok) then
            write (error_unit, '(a,a,a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [", label, "] adjoint u.(Jv)=", lhs, &
                " v.(J^T u)=", rhs, " rel_err=", e
        end if
    end function dot_product_identity

    ! --------------------------------------------------------- complex step

    ! Complex-step derivative of a real analytic scalar function whose code
    ! accepts a complex argument:  f'(x) ~ Im( f(x + i h) ) / h, with no
    ! subtractive cancellation, so h can be tiny and the result is accurate to
    ! machine precision. The function must be analytic at x.
    pure function complex_step_deriv(f, x, h_in) result(d)
        procedure(scalar_complex_fn_i) :: f
        real(dp), intent(in)           :: x
        real(dp), intent(in), optional :: h_in
        real(dp)                       :: d
        real(dp) :: h
        h = 1.0e-200_dp
        if (present(h_in)) h = h_in
        d = aimag(f(cmplx(x, h, kind=dp))) / h
    end function complex_step_deriv

    ! Complex-step estimate against a caller-supplied analytic derivative
    ! value. Because complex-step has no cancellation, tol may be a few ulp.
    function check_complex_step(label, f, x, dwant, tol, h_in) result(ok)
        character(*), intent(in)       :: label
        procedure(scalar_complex_fn_i) :: f
        real(dp), intent(in)           :: x
        real(dp), intent(in)           :: dwant
        real(dp), intent(in)           :: tol
        real(dp), intent(in), optional :: h_in
        logical                        :: ok
        real(dp) :: dgot, e

        if (present(h_in)) then
            dgot = complex_step_deriv(f, x, h_in=h_in)
        else
            dgot = complex_step_deriv(f, x)
        end if
        e = rel_err(dgot, dwant)
        ok = (e <= tol)
        if (.not. ok) then
            write (error_unit, '(a,a,a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [", label, "] complex-step got=", dgot, &
                " want=", dwant, " rel_err=", e
        end if
    end function check_complex_step

    ! --------------------------------------------------- smoothness / events

    ! Branch/event status check. A derivative is valid only on the trace the
    ! primal took; crossing a branch or event boundary is a non-smooth point
    ! and the derivative the module reports there is meaningless. The caller
    ! passes the trace tag at the base point and at the perturbed point and
    ! asserts which case it intends:
    !   expect == AD_SMOOTH    -> the two tags must agree (same trace -> pass);
    !   expect == AD_NONSMOOTH -> the tags must differ (branch change -> the
    !                             status must report non-smoothness).
    ! The returned status carries ok and the verdict actually observed.
    pure function check_smoothness(trace_base, trace_pert, expect) result(s)
        integer, intent(in) :: trace_base
        integer, intent(in) :: trace_pert
        integer, intent(in) :: expect
        type(ad_status_t)   :: s
        if (trace_base == trace_pert) then
            s%verdict = AD_SMOOTH
        else
            s%verdict = AD_NONSMOOTH
        end if
        s%ok = (s%verdict == expect)
    end function check_smoothness

    ! ----------------------------------------------------- Enzyme compare

    ! True only when fortnum is built with the Enzyme pass available. The
    ! preprocessor define is set by the build (FORTNUM_ENABLE_ENZYME); without
    ! it the whole Enzyme path stays inert and tests skip cleanly.
    pure function enzyme_available() result(yes)
        logical :: yes
#ifdef FORTNUM_ENABLE_ENZYME
        yes = .true.
#else
        yes = .false.
#endif
    end function enzyme_available

    ! Enzyme-vs-analytic forward-product compare. Identical contract to
    ! check_jvp_vs_fd but the reference is an analytic JVP instead of finite
    ! difference, used to validate a generated (transparent-path) product. In a
    ! non-Enzyme build enzyme_available() is .false. and callers skip; the
    ! routine still compiles and, if called, reports the worst component.
    function check_enzyme_vs_analytic(label, enzyme_jvp, analytic_jvp, &
            x, v, tol) result(ok)
        character(*), intent(in) :: label
        procedure(jvp_fn_i)      :: enzyme_jvp
        procedure(jvp_fn_i)      :: analytic_jvp
        real(dp), intent(in)     :: x(:)
        real(dp), intent(in)     :: v(:)
        real(dp), intent(in)     :: tol
        logical                  :: ok
        real(dp), allocatable :: je(:), ja(:)
        real(dp) :: worst
        integer  :: i, n

        n = size(x)
        allocate (je(n), ja(n))
        call enzyme_jvp(x, v, je)
        call analytic_jvp(x, v, ja)
        worst = 0.0_dp
        do i = 1, n
            worst = max(worst, rel_err(je(i), ja(i)))
        end do
        ok = (worst <= tol)
        if (.not. ok) then
            write (error_unit, '(a,a,a,es12.4,a,es12.4)') &
                "FAIL [", label, "] enzyme-vs-analytic worst rel_err=", &
                worst, " tol=", tol
        end if
    end function check_enzyme_vs_analytic

end module fortnum_ad_test_utils
