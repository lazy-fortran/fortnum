module fortnum_homotopy
    ! Polynomial homotopy continuation for square systems F(x; c) = 0 in n
    ! complex unknowns. Each term has its own coefficient c(k), so the family
    ! "same supports, arbitrary coefficients" is parametrised by c, and paths
    ! follow straight coefficient segments c(t) = ca + t (cb - ca).
    !
    ! Every linear solve is equilibrated: rows by the largest term magnitude of
    ! their equation, columns by the current |x_j|. Corrections are measured in
    ! the relative max-norm max_j |dx_j| / |x_j|, which keeps solutions whose
    ! coordinates span many orders of magnitude trackable in double precision.
    !
    ! Predictor: a coordinatewise Pade approximant of the order-8 Taylor
    ! polynomial of the path at the accepted point. The Taylor coefficients come
    ! from truncated power-series evaluation of the terms and reuse the LU
    ! factors of the last corrector iteration, so a step costs about two
    ! factorizations. Corrector: Newton. A step is accepted only when the first
    ! relative correction is below max_first_correction, which keeps path
    ! jumping rare, and Newton converges or stagnates at roundoff. The step
    ! size aims at a first correction of step_target.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit, real128
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_kinds, only: dp
    use fortnum_multiroot, only: argsort
    use fortnum_rng, only: rng_t, rng_seed, rng_split, rng_normal, rng_uniform
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: poly_system_t, poly_system_init, poly_system_add_term
    public :: poly_system_finalize, poly_system_eval
    public :: homotopy_options_t, homotopy_stats_t
    public :: homotopy_track, homotopy_plant, homotopy_monodromy, homotopy_move

    ! Fixed-size term buffers keep evaluation free of allocation.
    integer, parameter :: MAX_TERM_FACTORS = 32
    integer, parameter :: TAYLOR_ORDER = 8
    real(dp), parameter :: TINY_SCALE = 1.0e-300_dp
    ! A non-contracting Newton correction is roundoff once it is this small.
    real(dp), parameter :: ROUNDOFF_ACCEPT = 1.0e-7_dp
    ! Extra Newton iterations once the residual is formed in quad precision.
    integer, parameter :: EXTENDED_ITERS = 4
    ! Rounds of re-tracking paths whose endpoints coincide.
    integer, parameter :: JUMP_ROUNDS = 3
    integer, parameter :: qp = real128
    real(dp), parameter :: STEP_SAFETY = 0.8_dp
    real(dp), parameter :: STEP_GROW_MAX = 4.0_dp
    real(dp), parameter :: STEP_SHRINK_MIN = 0.1_dp
    real(dp), parameter :: TWO_PI = 8.0_dp * atan(1.0_dp)

    type :: poly_system_t
        integer :: n = 0
        integer :: ncoef = 0
        integer :: nterm = 0
        integer :: max_degree = 0
        logical :: finalized = .false.
        ! Terms in insertion order while the system is built.
        integer, allocatable :: add_eq(:), add_coef(:), add_ptr(:), add_vars(:)
        ! Packed by equation: terms of equation i are eq_ptr(i):eq_ptr(i+1)-1,
        ! in insertion order.
        integer, allocatable :: eq_ptr(:), coef(:)
        ! Distinct variables with exponents, for values and Jacobians.
        integer, allocatable :: fac_ptr(:), fac_var(:), fac_exp(:)
        ! Variables with multiplicity, for truncated power series.
        integer, allocatable :: deg_ptr(:), deg_var(:)
        ! Partial power-series products of terms with two or more factors.
        integer, allocatable :: prod_ptr(:)
    end type poly_system_t

    type :: homotopy_options_t
        ! Corrector tolerance along the path (relative max-norm).
        real(dp) :: tol = 1.0e-10_dp
        ! Refinement tolerance at the end of a segment.
        real(dp) :: final_tol = 1.0e-14_dp
        ! A step is rejected when the first relative correction exceeds this.
        real(dp) :: max_first_correction = 1.0e-3_dp
        ! Step-size control aims at this first relative correction.
        real(dp) :: step_target = 1.0e-5_dp
        real(dp) :: h_init = 0.05_dp
        ! A path fails when its step falls below h_min times the distance in t
        ! to the nearer end of the segment.
        real(dp) :: h_min = 1.0e-13_dp
        integer :: max_steps = 20000
        integer :: newton_iters = 4
        integer :: final_newton_iters = 8
        integer :: max_loops = 50
        integer :: stagnant_loops = 3
        ! Expected number of solutions; 0 means unknown.
        integer :: target = 0
        integer :: seed = 1
        integer :: verbosity = 0
        ! Two points are one solution when they agree to this relative
        ! distance in every coordinate.
        real(dp) :: dup_tol = 1.0e-6_dp
        ! homotopy_move re-tracks a path that fails in double precision with
        ! every operation in quad precision.
        logical :: quad_fallback = .true.
    end type homotopy_options_t

    type :: homotopy_stats_t
        integer :: loops = 0
        integer(int64) :: tracks = 0
        integer(int64) :: failed = 0
        integer(int64) :: steps = 0
        integer(int64) :: factorizations = 0
        ! Monodromy endpoints of one loop that coincide. The loop permutes the
        ! fibre, so each one is a path that jumped to another branch.
        integer(int64) :: jumps = 0
    end type homotopy_stats_t

    type :: track_work_t
        complex(dp), allocatable :: c(:), dc(:)
        complex(dp), allocatable :: jac(:, :), lu(:, :)
        complex(dp), allocatable :: f(:), z(:), xp(:)
        complex(dp), allocatable :: ser(:, :), pp(:, :)
        ! rs holds reciprocal row scales after scale_and_factor.
        real(dp), allocatable :: rs(:), cs(:)
        integer, allocatable :: ipiv(:)
        integer(int64) :: nfact = 0
        ! Set once Newton on this path has needed the quad-precision residual;
        ! from then on every residual of the path is formed that way, so the
        ! step control sees prediction error rather than roundoff.
        logical :: ext = .false.
    end type track_work_t

contains

    subroutine poly_system_init(sys, n)
        type(poly_system_t), intent(out) :: sys
        integer, intent(in) :: n

        if (n < 1) error stop "poly_system_init: n must be positive"
        sys%n = n
        allocate (sys%add_eq(16), sys%add_coef(16), sys%add_ptr(17), &
            sys%add_vars(64))
        sys%add_ptr(1) = 1
    end subroutine poly_system_init

    ! Add the term c(coef_index) * prod(x(vars)) to equation eq. A repeated
    ! index raises the power; an empty vars gives a constant term.
    subroutine poly_system_add_term(sys, eq, vars, coef_index)
        type(poly_system_t), intent(inout) :: sys
        integer, intent(in) :: eq
        integer, intent(in) :: vars(:)
        integer, intent(in) :: coef_index

        integer :: k, nv, last

        if (sys%n < 1) error stop "poly_system_add_term: system not initialised"
        if (sys%finalized) error stop "poly_system_add_term: system finalised"
        if (eq < 1 .or. eq > sys%n) error stop "poly_system_add_term: bad eq"
        if (coef_index < 1) error stop "poly_system_add_term: bad coef_index"
        nv = size(vars)
        if (nv > 0) then
            if (minval(vars) < 1 .or. maxval(vars) > sys%n) &
                error stop "poly_system_add_term: variable index out of range"
        end if
        k = sys%nterm + 1
        if (k > size(sys%add_eq)) then
            call grow_int(sys%add_ptr, 2 * size(sys%add_eq) + 1)
            call grow_int(sys%add_coef, 2 * size(sys%add_eq))
            call grow_int(sys%add_eq, 2 * size(sys%add_eq))
        end if
        last = sys%add_ptr(k) + nv - 1
        if (last > size(sys%add_vars)) then
            call grow_int(sys%add_vars, max(2 * size(sys%add_vars), last))
        end if
        sys%add_eq(k) = eq
        sys%add_coef(k) = coef_index
        sys%add_vars(sys%add_ptr(k):last) = vars
        sys%add_ptr(k + 1) = last + 1
        sys%nterm = k
        sys%ncoef = max(sys%ncoef, coef_index)
    end subroutine poly_system_add_term

    ! Pack the terms by equation (stable, so the first term of an equation is
    ! the first one added) and split each monomial into distinct variables.
    subroutine poly_system_finalize(sys)
        type(poly_system_t), intent(inout) :: sys

        integer, allocatable :: pos(:), order(:), vtmp(:)
        integer :: i, k, q, nv, nt, p, lo, a, b, tmp

        if (sys%finalized) return
        if (sys%n < 1) error stop "poly_system_finalize: system not initialised"
        nt = sys%nterm
        allocate (sys%eq_ptr(sys%n + 1), pos(sys%n), order(max(nt, 1)))
        pos = 0
        do k = 1, nt
            pos(sys%add_eq(k)) = pos(sys%add_eq(k)) + 1
        end do
        if (any(pos == 0)) error stop "poly_system_finalize: equation without terms"
        sys%eq_ptr(1) = 1
        do i = 1, sys%n
            sys%eq_ptr(i + 1) = sys%eq_ptr(i) + pos(i)
        end do
        pos = sys%eq_ptr(1:sys%n)
        do k = 1, nt
            order(pos(sys%add_eq(k))) = k
            pos(sys%add_eq(k)) = pos(sys%add_eq(k)) + 1
        end do

        nv = 0
        do k = 1, nt
            nv = max(nv, sys%add_ptr(k + 1) - sys%add_ptr(k))
        end do
        allocate (vtmp(max(nv, 1)))
        p = sys%add_ptr(nt + 1) - 1
        allocate (sys%prod_ptr(nt + 1))
        allocate (sys%coef(nt), sys%deg_ptr(nt + 1), sys%deg_var(max(p, 1)), &
            sys%fac_ptr(nt + 1), sys%fac_var(max(p, 1)), sys%fac_exp(max(p, 1)))
        sys%deg_ptr(1) = 1
        sys%fac_ptr(1) = 1
        sys%prod_ptr(1) = 1
        sys%max_degree = 0
        do q = 1, nt
            k = order(q)
            sys%coef(q) = sys%add_coef(k)
            lo = sys%add_ptr(k)
            nv = sys%add_ptr(k + 1) - lo
            vtmp(1:nv) = sys%add_vars(lo:lo + nv - 1)
            do a = 2, nv
                tmp = vtmp(a)
                b = a - 1
                do while (b >= 1)
                    if (vtmp(b) <= tmp) exit
                    vtmp(b + 1) = vtmp(b)
                    b = b - 1
                end do
                vtmp(b + 1) = tmp
            end do
            sys%deg_var(sys%deg_ptr(q):sys%deg_ptr(q) + nv - 1) = vtmp(1:nv)
            sys%deg_ptr(q + 1) = sys%deg_ptr(q) + nv
            sys%prod_ptr(q + 1) = sys%prod_ptr(q) + max(nv - 1, 0)
            p = sys%fac_ptr(q) - 1
            do a = 1, nv
                if (a > 1) then
                    if (vtmp(a) == vtmp(a - 1)) then
                        sys%fac_exp(p) = sys%fac_exp(p) + 1
                        cycle
                    end if
                end if
                p = p + 1
                sys%fac_var(p) = vtmp(a)
                sys%fac_exp(p) = 1
            end do
            sys%fac_ptr(q + 1) = p + 1
            if (sys%fac_ptr(q + 1) - sys%fac_ptr(q) > MAX_TERM_FACTORS) &
                error stop "poly_system_finalize: too many distinct variables in a term"
            sys%max_degree = max(sys%max_degree, nv)
        end do
        deallocate (sys%add_eq, sys%add_coef, sys%add_ptr, sys%add_vars)
        sys%finalized = .true.
    end subroutine poly_system_finalize

    ! f = F(x; c), jac = dF/dx, ft = F(x; dc), rowscale(i) = largest
    ! |c_k m_k(x)| over the terms of equation i.
    pure subroutine poly_system_eval(sys, x, c, dc, f, jac, ft, rowscale)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x(:), c(:), dc(:)
        complex(dp), intent(out) :: f(:), jac(:, :), ft(:)
        real(dp), intent(out) :: rowscale(:)

        if (.not. sys%finalized) error stop "poly_system_eval: system not finalised"
        call eval_terms(sys, x, c, f, jac, rowscale, dc, ft)
    end subroutine poly_system_eval

    ! F(x; c) with every product and sum in quad precision, rounded to double
    ! once at the end.
    pure subroutine eval_residual_ext(sys, x, c, f)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x(:), c(:)
        complex(dp), intent(out) :: f(:)

        complex(qp) :: acc, mono
        integer :: i, k, l, lo, m

        do i = 1, sys%n
            acc = (0.0_qp, 0.0_qp)
            do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                mono = cmplx(c(sys%coef(k)), kind=qp)
                lo = sys%fac_ptr(k)
                m = sys%fac_ptr(k + 1) - lo
                do l = 1, m
                    mono = mono * cmplx(x(sys%fac_var(lo + l - 1)), kind=qp) &
                        **sys%fac_exp(lo + l - 1)
                end do
                acc = acc + mono
            end do
            f(i) = cmplx(acc, kind=dp)
        end do
    end subroutine eval_residual_ext

    pure subroutine eval_terms(sys, x, c, f, jac, rs, dc, ft)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x(:), c(:)
        complex(dp), intent(out) :: f(:), jac(:, :)
        real(dp), intent(out) :: rs(:)
        complex(dp), intent(in), optional :: dc(:)
        complex(dp), intent(out), optional :: ft(:)

        complex(dp) :: fv(MAX_TERM_FACTORS), pw(MAX_TERM_FACTORS)
        complex(dp) :: pre(MAX_TERM_FACTORS + 1)
        complex(dp) :: ck, mono, suf, cm, p
        integer :: i, j, k, l, lo, m, v, e
        logical :: want_ft

        want_ft = present(ft)
        do j = 1, sys%n
            do i = 1, sys%n
                jac(i, j) = (0.0_dp, 0.0_dp)
            end do
            f(j) = (0.0_dp, 0.0_dp)
            rs(j) = 0.0_dp
        end do
        if (want_ft) ft(1:sys%n) = (0.0_dp, 0.0_dp)
        do i = 1, sys%n
            do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                ck = c(sys%coef(k))
                lo = sys%fac_ptr(k)
                m = sys%fac_ptr(k + 1) - lo
                if (m == 0) then
                    mono = (1.0_dp, 0.0_dp)
                else if (m == 1) then
                    v = sys%fac_var(lo)
                    e = sys%fac_exp(lo)
                    if (e == 1) then
                        mono = x(v)
                        jac(i, v) = jac(i, v) + ck
                    else
                        p = x(v)**(e - 1)
                        mono = p * x(v)
                        jac(i, v) = jac(i, v) + ck * (real(e, dp) * p)
                    end if
                else
                    pre(1) = (1.0_dp, 0.0_dp)
                    do l = 1, m
                        v = sys%fac_var(lo + l - 1)
                        e = sys%fac_exp(lo + l - 1)
                        if (e == 1) then
                            pw(l) = (1.0_dp, 0.0_dp)
                            fv(l) = x(v)
                        else
                            pw(l) = x(v)**(e - 1)
                            fv(l) = pw(l) * x(v)
                        end if
                        pre(l + 1) = pre(l) * fv(l)
                    end do
                    mono = pre(m + 1)
                    suf = (1.0_dp, 0.0_dp)
                    do l = m, 1, -1
                        v = sys%fac_var(lo + l - 1)
                        e = sys%fac_exp(lo + l - 1)
                        jac(i, v) = jac(i, v) &
                            + ck * (real(e, dp) * pw(l)) * (pre(l) * suf)
                        suf = suf * fv(l)
                    end do
                end if
                cm = ck * mono
                f(i) = f(i) + cm
                rs(i) = max(rs(i), modulus(cm))
                if (want_ft) ft(i) = ft(i) + dc(sys%coef(k)) * mono
            end do
        end do
    end subroutine eval_terms

    ! |z| without the cost of hypot; falls back to it where squares would
    ! overflow or underflow.
    elemental real(dp) function modulus(z) result(a)
        complex(dp), intent(in) :: z

        real(dp) :: re, im, big

        re = real(z, dp)
        im = aimag(z)
        big = max(abs(re), abs(im))
        if (big < 1.0e150_dp .and. big > 1.0e-150_dp) then
            a = sqrt(re * re + im * im)
        else
            a = abs(z)
        end if
    end function modulus

    pure function term_value(sys, k, x) result(mono)
        type(poly_system_t), intent(in) :: sys
        integer, intent(in) :: k
        complex(dp), intent(in) :: x(:)
        complex(dp) :: mono

        integer :: l

        mono = (1.0_dp, 0.0_dp)
        do l = sys%fac_ptr(k), sys%fac_ptr(k + 1) - 1
            mono = mono * x(sys%fac_var(l))**sys%fac_exp(l)
        end do
    end function term_value

    subroutine work_init(w, sys)
        type(track_work_t), intent(out) :: w
        type(poly_system_t), intent(in) :: sys

        integer :: n

        n = sys%n
        allocate (w%c(sys%ncoef), w%dc(sys%ncoef), w%jac(n, n), w%lu(n, n), &
            w%f(n), w%z(n), w%xp(n), w%ser(0:TAYLOR_ORDER, n), &
            w%rs(n), w%cs(n), w%ipiv(n), &
            w%pp(0:TAYLOR_ORDER, max(1, sys%prod_ptr(sys%nterm + 1) - 1)))
        w%nfact = 0
    end subroutine work_init

    ! Equilibrate J at x with the row scales from evaluation and the column
    ! scales |x_j|, then factor. Leaves reciprocal row scales in w%rs.
    pure subroutine scale_and_factor(n, x, w, ok)
        integer, intent(in) :: n
        complex(dp), intent(in) :: x(:)
        type(track_work_t), intent(inout) :: w
        logical, intent(out) :: ok

        integer :: i, j

        do i = 1, n
            w%rs(i) = 1.0_dp / max(w%rs(i), TINY_SCALE)
            w%cs(i) = max(modulus(x(i)), TINY_SCALE)
        end do
        do j = 1, n
            do i = 1, n
                w%lu(i, j) = w%jac(i, j) * (w%cs(j) * w%rs(i))
            end do
        end do
        w%nfact = w%nfact + 1
        call lu_factor(n, w%lu, w%ipiv, ok)
    end subroutine scale_and_factor

    ! Right-looking LU with partial pivoting, column oriented, in place.
    pure subroutine lu_factor(n, a, ipiv, ok)
        integer, intent(in) :: n
        complex(dp), intent(inout) :: a(:, :)
        integer, intent(out) :: ipiv(:)
        logical, intent(out) :: ok

        integer :: i, j, k, p
        real(dp) :: best, v
        complex(dp) :: t, pinv

        ok = .false.
        do k = 1, n
            p = k
            best = abs(real(a(k, k), dp)) + abs(aimag(a(k, k)))
            do i = k + 1, n
                v = abs(real(a(i, k), dp)) + abs(aimag(a(i, k)))
                if (v > best) then
                    best = v
                    p = i
                end if
            end do
            ipiv(k) = p
            if (.not. (best > 0.0_dp)) return
            if (best > huge(1.0_dp)) return
            if (p /= k) then
                do j = 1, n
                    t = a(k, j)
                    a(k, j) = a(p, j)
                    a(p, j) = t
                end do
            end if
            pinv = (1.0_dp, 0.0_dp) / a(k, k)
            do i = k + 1, n
                a(i, k) = a(i, k) * pinv
            end do
            do j = k + 1, n
                t = a(k, j)
                do i = k + 1, n
                    a(i, j) = a(i, j) - t * a(i, k)
                end do
            end do
        end do
        ok = .true.
    end subroutine lu_factor

    pure subroutine lu_solve(n, a, ipiv, b)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a(:, :)
        integer, intent(in) :: ipiv(:)
        complex(dp), intent(inout) :: b(:)

        integer :: i, j
        complex(dp) :: t

        do j = 1, n
            if (ipiv(j) /= j) then
                t = b(j)
                b(j) = b(ipiv(j))
                b(ipiv(j)) = t
            end if
        end do
        do j = 1, n - 1
            t = b(j)
            do i = j + 1, n
                b(i) = b(i) - t * a(i, j)
            end do
        end do
        do j = n, 1, -1
            b(j) = b(j) / a(j, j)
            t = b(j)
            do i = 1, j - 1
                b(i) = b(i) - t * a(i, j)
            end do
        end do
    end subroutine lu_solve

    ! Newton at fixed c with equilibrated solves. first and last are the
    ! relative max-norms of the first and last corrections. The iteration
    ! fails if the first correction exceeds max_first; it converges below tol,
    ! when the predicted next correction is below tol, or when a correction
    ! stops contracting while the previous one was already at roundoff level.
    pure subroutine correct(sys, x, c, maxit, tol, max_first, w, ok, first, last)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(inout) :: x(:)
        complex(dp), intent(in) :: c(:)
        integer, intent(in) :: maxit
        real(dp), intent(in) :: tol, max_first
        type(track_work_t), intent(inout) :: w
        logical, intent(out) :: ok
        real(dp), intent(out) :: first, last

        real(dp) :: nrm, prev, total, a
        integer :: it, j, n, itmax
        logical :: fok, ext

        n = sys%n
        ok = .false.
        first = huge(1.0_dp)
        last = huge(1.0_dp)
        prev = huge(1.0_dp)
        ext = w%ext
        itmax = maxit
        it = 0
        do while (it < itmax)
            it = it + 1
            call eval_terms(sys, x, c, w%f, w%jac, w%rs)
            if (ext) call eval_residual_ext(sys, x, c, w%f)
            call scale_and_factor(n, x, w, fok)
            if (.not. fok) return
            do j = 1, n
                w%z(j) = -w%f(j) * w%rs(j)
            end do
            call lu_solve(n, w%lu, w%ipiv, w%z)
            nrm = 0.0_dp
            total = 0.0_dp
            do j = 1, n
                a = modulus(w%z(j))
                nrm = max(nrm, a)
                total = total + a
            end do
            if (.not. ieee_is_finite(total)) return
            do j = 1, n
                x(j) = x(j) + w%cs(j) * w%z(j)
            end do
            if (it == 1) then
                first = nrm
                if (first > max_first) return
            end if
            last = nrm
            if (nrm < tol) then
                ok = .true.
                return
            end if
            if (it > 1) then
                if (nrm > 0.5_dp * prev) then
                    ! Stalled above roundoff: the residual in double cannot
                    ! resolve an ill-conditioned point. Iterative refinement
                    ! with the residual in quad precision and the same double
                    ! Jacobian reaches double accuracy instead.
                    if (.not. ext .and. prev >= ROUNDOFF_ACCEPT) then
                        ext = .true.
                        w%ext = .true.
                        itmax = it + EXTENDED_ITERS
                        prev = nrm
                        cycle
                    end if
                    ok = (prev < ROUNDOFF_ACCEPT)
                    return
                end if
                if (nrm * (nrm / prev) < tol) then
                    ok = .true.
                    return
                end if
            end if
            prev = nrm
        end do
    end subroutine correct

    ! Taylor coefficients ser(0:K, :) of the path through (x, t), where w%c is
    ! c(t) and the factors in w belong to J at (a point within tol of) x.
    ! Order j solves J x_j = -[s^j] F(x_0 + ... + x_{j-1} s^{j-1}; c + s dc).
    ! w%pp keeps the partial products P_l = X_{v1} ... X_{vl} of every term,
    ! so each order costs O(j) per factor instead of a full series product:
    ! coefficient j-1 only lacks the terms with x_{j-1}, and coefficient j is
    ! formed with x_j = 0.
    pure subroutine taylor_series(sys, x, w)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x(:)
        type(track_work_t), intent(inout) :: w

        complex(dp) :: r, s, dl, dprev, p0prev, prev
        integer :: i, j, k, l, a, v, v1, lo, d, n, q0

        n = sys%n
        do v = 1, n
            w%ser(0, v) = x(v)
            do j = 1, TAYLOR_ORDER
                w%ser(j, v) = (0.0_dp, 0.0_dp)
            end do
        end do
        do j = 1, TAYLOR_ORDER
            do i = 1, n
                r = (0.0_dp, 0.0_dp)
                do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                    lo = sys%deg_ptr(k)
                    d = sys%deg_ptr(k + 1) - lo
                    if (d == 0) then
                        if (j == 1) r = r + w%dc(sys%coef(k))
                        cycle
                    end if
                    v1 = sys%deg_var(lo)
                    if (d == 1) then
                        r = r + w%dc(sys%coef(k)) * w%ser(j - 1, v1)
                        cycle
                    end if
                    q0 = sys%prod_ptr(k) - 2
                    if (j == 1) then
                        prev = w%ser(0, v1)
                        do l = 2, d
                            v = sys%deg_var(lo + l - 1)
                            prev = prev * w%ser(0, v)
                            w%pp(0, q0 + l) = prev
                            w%pp(1, q0 + l) = (0.0_dp, 0.0_dp)
                        end do
                    else
                        ! Complete coefficient j-1 now that x_{j-1} is known.
                        dprev = w%ser(j - 1, v1)
                        p0prev = w%ser(0, v1)
                        do l = 2, d
                            v = sys%deg_var(lo + l - 1)
                            dl = p0prev * w%ser(j - 1, v) + dprev * w%ser(0, v)
                            w%pp(j - 1, q0 + l) = w%pp(j - 1, q0 + l) + dl
                            dprev = dl
                            p0prev = w%pp(0, q0 + l)
                        end do
                        ! Coefficient j with x_j = 0.
                        v = sys%deg_var(lo + 1)
                        s = (0.0_dp, 0.0_dp)
                        do a = 1, j - 1
                            s = s + w%ser(a, v1) * w%ser(j - a, v)
                        end do
                        w%pp(j, q0 + 2) = s
                        do l = 3, d
                            v = sys%deg_var(lo + l - 1)
                            s = w%pp(j, q0 + l - 1) * w%ser(0, v)
                            do a = 1, j - 1
                                s = s + w%pp(a, q0 + l - 1) * w%ser(j - a, v)
                            end do
                            w%pp(j, q0 + l) = s
                        end do
                    end if
                    r = r + w%c(sys%coef(k)) * w%pp(j, q0 + d) &
                        + w%dc(sys%coef(k)) * w%pp(j - 1, q0 + d)
                end do
                w%z(i) = -r * w%rs(i)
            end do
            call lu_solve(n, w%lu, w%ipiv, w%z)
            do v = 1, n
                w%ser(j, v) = w%cs(v) * w%z(v)
            end do
        end do
    end subroutine taylor_series

    ! Coordinatewise Pade [K-1/1] approximant of the Taylor polynomial at h.
    ! Its pole -1/b follows a nearby singularity of the path; when that pole
    ! lies within 2h the plain Taylor polynomial is used instead.
    pure subroutine predict(n, w, h)
        integer, intent(in) :: n
        type(track_work_t), intent(inout) :: w
        real(dp), intent(in) :: h

        complex(dp) :: p, b
        integer :: j, v

        do v = 1, n
            if (modulus(w%ser(TAYLOR_ORDER, v)) * h &
                < 0.5_dp * modulus(w%ser(TAYLOR_ORDER - 1, v))) then
                b = -w%ser(TAYLOR_ORDER, v) / w%ser(TAYLOR_ORDER - 1, v)
                p = w%ser(TAYLOR_ORDER - 1, v) + b * w%ser(TAYLOR_ORDER - 2, v)
                do j = TAYLOR_ORDER - 2, 1, -1
                    p = p * h + (w%ser(j, v) + b * w%ser(j - 1, v))
                end do
                p = p * h + w%ser(0, v)
                w%xp(v) = p / (1.0_dp + b * h)
            else
                p = w%ser(TAYLOR_ORDER, v)
                do j = TAYLOR_ORDER - 1, 0, -1
                    p = p * h + w%ser(j, v)
                end do
                w%xp(v) = p
            end if
        end do
    end subroutine predict

    ! Step factor from the first corrector norm, taking the predictor error
    ! to scale like h**(TAYLOR_ORDER + 1).
    pure real(dp) function step_factor(first, target) result(fac)
        real(dp), intent(in) :: first, target

        if (first > 0.0_dp) then
            fac = STEP_SAFETY * (target / first)**(1.0_dp / (TAYLOR_ORDER + 1))
        else
            fac = STEP_GROW_MAX
        end if
    end function step_factor

    pure subroutine track_segment(sys, x, ca, cb, opts, w, ok, nsteps, acc, &
            t_reached, s_reached)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(inout) :: x(:)
        complex(dp), intent(in) :: ca(:), cb(:)
        type(homotopy_options_t), intent(in) :: opts
        type(track_work_t), intent(inout) :: w
        logical, intent(out) :: ok
        integer, intent(out) :: nsteps
        real(dp), intent(out) :: acc
        ! On failure x is the last accepted point and these give its position.
        real(dp), intent(out), optional :: t_reached, s_reached

        ! The position is held as t in the first half of the segment and as the
        ! remaining distance s = 1 - t in the second, and the coefficients are
        ! formed from the nearer end. Both ends are then resolved down to
        ! tiny(1.0_dp) without cancellation, which a start fibre at another
        ! coefficient scale and an endpoint whose coordinates are many orders
        ! of magnitude apart both need.
        real(dp) :: t, s, h, tn, sn, first, last
        logical :: cok, have_series
        integer :: k, n

        n = sys%n
        ok = .false.
        nsteps = 0
        acc = huge(1.0_dp)
        w%ext = .false.
        do k = 1, sys%ncoef
            w%dc(k) = cb(k) - ca(k)
            w%c(k) = ca(k)
        end do
        ! The start point's last factorization seeds the first predictor.
        call correct(sys, x, w%c, opts%newton_iters, opts%tol, huge(1.0_dp), &
            w, cok, first, last)
        if (.not. cok) return
        t = 0.0_dp
        s = 1.0_dp
        if (present(t_reached)) t_reached = t
        if (present(s_reached)) s_reached = s
        h = opts%h_init
        have_series = .false.
        do while (s > 0.0_dp)
            if (present(t_reached)) t_reached = t
            if (present(s_reached)) s_reached = s
            if (nsteps >= opts%max_steps) return
            if (h < opts%h_min * min(t, s)) return
            if (t < 0.5_dp) then
                if (.not. (t + h > t)) return
            else
                if (.not. (s - h < s)) return
            end if
            if (.not. have_series) then
                call taylor_series(sys, x, w)
                have_series = .true.
            end if
            nsteps = nsteps + 1
            if (h >= s) then
                h = s
                tn = 1.0_dp
                sn = 0.0_dp
                w%c(1:sys%ncoef) = cb(1:sys%ncoef)
            else if (t + h <= 0.5_dp) then
                tn = t + h
                sn = 1.0_dp - tn
                do k = 1, sys%ncoef
                    w%c(k) = ca(k) + tn * w%dc(k)
                end do
            else
                sn = s - h
                tn = 1.0_dp - sn
                do k = 1, sys%ncoef
                    w%c(k) = cb(k) - sn * w%dc(k)
                end do
            end if
            call predict(n, w, h)
            call correct(sys, w%xp, w%c, opts%newton_iters, opts%tol, &
                opts%max_first_correction, w, cok, first, last)
            if (cok) then
                x(1:n) = w%xp(1:n)
                t = tn
                s = sn
                have_series = .false.
                h = h * min(STEP_GROW_MAX, max(0.5_dp, &
                    step_factor(first, opts%step_target)))
            else if (first < huge(1.0_dp)) then
                h = h * min(0.5_dp, max(STEP_SHRINK_MIN, &
                    step_factor(first, opts%step_target)))
            else
                h = 0.5_dp * h
            end if
        end do
        ! The final refinement always forms the residual in quad precision: at
        ! an ill-conditioned point the residual in double can round to zero
        ! long before the point is accurate, which would end Newton early and
        ! report a false accuracy.
        w%ext = .true.
        call correct(sys, x, cb, opts%final_newton_iters, opts%final_tol, &
            huge(1.0_dp), w, cok, first, last)
        acc = last
        ok = .true.
    end subroutine track_segment

    ! Fallback tracker with every operation in quad precision, for paths on
    ! which double precision is exhausted (scaled condition numbers near
    ! 1/epsilon(1.0_dp)). It is the double tracker's method in its simplest
    ! form: an RK4 predictor on dx/dt = -J^{-1} dF/dt, a Newton corrector with
    ! the same row and column scaling, step halving and doubling, and the
    ! position held as t in the first half of the segment and 1 - t in the
    ! second. acc is the final relative Newton correction.
    subroutine track_segment_qp(sys, xd, ca, cb, opts, ok, nsteps, acc, t0, s0)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(inout) :: xd(:)
        complex(dp), intent(in) :: ca(:), cb(:)
        type(homotopy_options_t), intent(in) :: opts
        logical, intent(out) :: ok
        integer, intent(out) :: nsteps
        real(dp), intent(out) :: acc
        ! Start position on the segment (default: its start); xd lies on it.
        real(dp), intent(in), optional :: t0, s0

        ! The result is rounded to double, so quad Newton aims just below double
        ! accuracy; at scaled condition numbers near 1e14, quad precision itself
        ! cannot resolve corrections much below 1e-20.
        real(qp), parameter :: TOL = 1.0e-16_qp, FINAL_TOL = 1.0e-19_qp
        real(qp), parameter :: ACCEPT = 1.0e-13_qp, FINAL_ACCEPT = 1.0e-16_qp
        real(qp), parameter :: GROW_BELOW = 1.0e-6_qp, H_MIN = 1.0e-30_qp
        complex(qp), allocatable :: x(:), x0(:), xp(:), k1(:), k2(:), k3(:), k4(:)
        complex(qp), allocatable :: caq(:), cbq(:), dq(:), c(:)
        real(qp) :: t, s, h, tn, sn, first, last
        logical :: good
        integer :: n, nc

        n = sys%n
        nc = sys%ncoef
        allocate (x(n), x0(n), xp(n), k1(n), k2(n), k3(n), k4(n))
        allocate (caq(nc), cbq(nc), dq(nc), c(nc))
        caq = cmplx(ca(1:nc), kind=qp)
        cbq = cmplx(cb(1:nc), kind=qp)
        dq = cbq - caq
        x = cmplx(xd(1:n), kind=qp)
        ok = .false.
        nsteps = 0
        acc = huge(1.0_dp)
        t = 0.0_qp
        s = 1.0_qp
        if (present(t0) .and. present(s0)) then
            if (s0 < 0.5_dp) then
                s = real(s0, qp)
                t = 1.0_qp - s
            else
                t = real(t0, qp)
                s = 1.0_qp - t
            end if
        end if
        call coef_at(t, s, c)
        call newton_qp(sys, x, c, opts%final_newton_iters, TOL, huge(1.0_qp), &
            ACCEPT, good, first, last)
        if (.not. good) return
        h = min(real(opts%h_init, qp), 0.1_qp * max(min(t, s), H_MIN))
        if (t == 0.0_qp) h = real(opts%h_init, qp)
        do while (s > 0.0_qp)
            if (nsteps >= opts%max_steps) return
            if (h < H_MIN * min(t, s)) return
            nsteps = nsteps + 1
            if (h >= s) then
                h = s
                tn = 1.0_qp
                sn = 0.0_qp
            else if (t + h <= 0.5_qp) then
                tn = t + h
                sn = 1.0_qp - tn
            else
                sn = s - h
                tn = 1.0_qp - sn
            end if
            x0 = x
            call coef_at(t, s, c)
            call velocity_qp(sys, x0, c, dq, k1, good)
            if (good) then
                call coef_at(t + 0.5_qp * h, s - 0.5_qp * h, c)
                call velocity_qp(sys, x0 + 0.5_qp * h * k1, c, dq, k2, good)
            end if
            if (good) call velocity_qp(sys, x0 + 0.5_qp * h * k2, c, dq, k3, good)
            if (good) then
                call coef_at(tn, sn, c)
                call velocity_qp(sys, x0 + h * k3, c, dq, k4, good)
            end if
            if (good) then
                xp = x0 + h / 6.0_qp * (k1 + 2.0_qp * k2 + 2.0_qp * k3 + k4)
                call newton_qp(sys, xp, c, opts%newton_iters, TOL, &
                    real(opts%max_first_correction, qp), ACCEPT, good, first, last)
            end if
            if (good) then
                x = xp
                t = tn
                s = sn
                if (first < GROW_BELOW) h = 2.0_qp * h
            else
                x = x0
                h = 0.5_qp * h
            end if
        end do
        call newton_qp(sys, x, cbq, opts%final_newton_iters, FINAL_TOL, &
            huge(1.0_qp), FINAL_ACCEPT, good, first, last)
        xd(1:n) = cmplx(x, kind=dp)
        acc = real(last, dp)
        ok = .true.

    contains

        ! c(t) from the nearer end of the segment.
        pure subroutine coef_at(tt, ss, cc)
            real(qp), intent(in) :: tt, ss
            complex(qp), intent(out) :: cc(:)

            if (tt <= 0.5_qp) then
                cc = caq + tt * dq
            else
                cc = cbq - ss * dq
            end if
        end subroutine coef_at
    end subroutine track_segment_qp

    ! F(x; c), J = dF/dx, ft = F(x; dc) and the row scales, in quad precision.
    pure subroutine eval_qp(sys, x, c, dc, f, jac, ft, rs)
        type(poly_system_t), intent(in) :: sys
        complex(qp), intent(in) :: x(:), c(:), dc(:)
        complex(qp), intent(out) :: f(:), jac(:, :), ft(:)
        real(qp), intent(out) :: rs(:)

        complex(qp) :: mono, part, cm
        integer :: i, k, l, j, lo, m, v, e

        f = (0.0_qp, 0.0_qp)
        ft = (0.0_qp, 0.0_qp)
        jac = (0.0_qp, 0.0_qp)
        rs = 0.0_qp
        do i = 1, sys%n
            do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                lo = sys%fac_ptr(k)
                m = sys%fac_ptr(k + 1) - lo
                mono = (1.0_qp, 0.0_qp)
                do l = 1, m
                    mono = mono * x(sys%fac_var(lo + l - 1))**sys%fac_exp(lo + l - 1)
                end do
                cm = c(sys%coef(k)) * mono
                f(i) = f(i) + cm
                ft(i) = ft(i) + dc(sys%coef(k)) * mono
                rs(i) = max(rs(i), abs(cm))
                do l = 1, m
                    v = sys%fac_var(lo + l - 1)
                    e = sys%fac_exp(lo + l - 1)
                    part = c(sys%coef(k)) * real(e, qp) * x(v)**(e - 1)
                    do j = 1, m
                        if (j == l) cycle
                        part = part * x(sys%fac_var(lo + j - 1))**sys%fac_exp(lo + j - 1)
                    end do
                    jac(i, v) = jac(i, v) + part
                end do
            end do
        end do
    end subroutine eval_qp

    ! Solve (D_r J D_c) z = D_r b with row scales 1/rs and column scales |x|,
    ! by Gaussian elimination with partial pivoting; returns D_c z.
    pure subroutine solve_scaled_qp(jac, rs, x, b, dx, good)
        complex(qp), intent(in) :: jac(:, :), b(:), x(:)
        real(qp), intent(in) :: rs(:)
        complex(qp), intent(out) :: dx(:)
        logical, intent(out) :: good

        complex(qp) :: a(size(b), size(b)), z(size(b)), tmp, piv, fac
        real(qp) :: cs(size(b)), r(size(b)), best
        integer :: i, j, k, p, n

        n = size(b)
        good = .false.
        dx = (0.0_qp, 0.0_qp)
        r = 1.0_qp / max(rs, maxval(rs) * 1.0e-4000_qp, tiny(1.0_qp))
        cs = max(abs(x), tiny(1.0_qp))
        do j = 1, n
            do i = 1, n
                a(i, j) = jac(i, j) * cs(j) * r(i)
            end do
        end do
        z = b * r
        do k = 1, n
            p = k
            best = abs(a(k, k))
            do i = k + 1, n
                if (abs(a(i, k)) > best) then
                    best = abs(a(i, k))
                    p = i
                end if
            end do
            if (.not. (best > 0.0_qp)) return
            if (p /= k) then
                do j = 1, n
                    tmp = a(k, j)
                    a(k, j) = a(p, j)
                    a(p, j) = tmp
                end do
                tmp = z(k)
                z(k) = z(p)
                z(p) = tmp
            end if
            piv = a(k, k)
            do i = k + 1, n
                fac = a(i, k) / piv
                if (fac == (0.0_qp, 0.0_qp)) cycle
                a(i, k + 1:n) = a(i, k + 1:n) - fac * a(k, k + 1:n)
                z(i) = z(i) - fac * z(k)
            end do
        end do
        do k = n, 1, -1
            z(k) = (z(k) - sum(a(k, k + 1:n) * z(k + 1:n))) / a(k, k)
        end do
        if (.not. all(ieee_is_finite(real(z, dp)))) return
        dx = z * cs
        good = .true.
    end subroutine solve_scaled_qp

    pure subroutine velocity_qp(sys, x, c, dc, v, good)
        type(poly_system_t), intent(in) :: sys
        complex(qp), intent(in) :: x(:), c(:), dc(:)
        complex(qp), intent(out) :: v(:)
        logical, intent(out) :: good

        complex(qp) :: f(sys%n), jac(sys%n, sys%n), ft(sys%n)
        real(qp) :: rs(sys%n)

        call eval_qp(sys, x, c, dc, f, jac, ft, rs)
        call solve_scaled_qp(jac, rs, x, -ft, v, good)
    end subroutine velocity_qp

    ! Newton in quad precision with relative corrections; good when the
    ! correction falls below tol or stops contracting below accept.
    pure subroutine newton_qp(sys, x, c, maxit, tol, max_first, accept, good, &
            first, last)
        type(poly_system_t), intent(in) :: sys
        complex(qp), intent(inout) :: x(:)
        complex(qp), intent(in) :: c(:)
        integer, intent(in) :: maxit
        real(qp), intent(in) :: tol, max_first, accept
        logical, intent(out) :: good
        real(qp), intent(out) :: first, last

        complex(qp) :: f(sys%n), jac(sys%n, sys%n), ft(sys%n), dx(sys%n), zero(size(c))
        real(qp) :: rs(sys%n), nrm, prev
        logical :: sok
        integer :: it

        zero = (0.0_qp, 0.0_qp)
        good = .false.
        first = huge(1.0_qp)
        last = huge(1.0_qp)
        prev = huge(1.0_qp)
        do it = 1, maxit
            call eval_qp(sys, x, c, zero, f, jac, ft, rs)
            call solve_scaled_qp(jac, rs, x, -f, dx, sok)
            if (.not. sok) return
            nrm = maxval(abs(dx) / max(abs(x), tiny(1.0_qp)))
            x = x + dx
            if (it == 1) then
                first = nrm
                if (first > max_first) return
            end if
            last = nrm
            if (nrm < tol) then
                good = .true.
                return
            end if
            if (it > 1) then
                if (nrm > 0.5_qp * prev) then
                    good = (prev < accept)
                    return
                end if
            end if
            prev = nrm
        end do
    end subroutine newton_qp

    ! Track one solution x of F(x; ca) = 0 to F(x; cb) = 0 along the straight
    ! coefficient segment. acc is the final relative Newton correction.
    subroutine homotopy_track(sys, x, ca, cb, opts, ok, nsteps, acc)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(inout) :: x(:)
        complex(dp), intent(in) :: ca(:), cb(:)
        type(homotopy_options_t), intent(in) :: opts
        logical, intent(out) :: ok
        integer, intent(out) :: nsteps
        real(dp), intent(out), optional :: acc

        type(track_work_t) :: w
        real(dp) :: a

        call check_inputs(sys, size(x), size(ca), size(cb), "homotopy_track")
        call work_init(w, sys)
        call track_segment(sys, x, ca, cb, opts, w, ok, nsteps, a)
        if (present(acc)) acc = a
    end subroutine homotopy_track

    ! Random x0 on the unit torus and complex Gaussian coefficients c0; the
    ! coefficient of the first term of each equation is then solved so that
    ! F(x0; c0) = 0. That coefficient index must not occur in another
    ! equation.
    subroutine homotopy_plant(sys, x0, c0, seed)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(out) :: x0(:), c0(:)
        integer, intent(in) :: seed

        type(rng_t) :: g
        type(fortnum_status_t) :: st
        integer, allocatable :: owner(:)
        complex(dp) :: a, b, mono
        real(dp) :: u, re, im
        integer :: i, j, k, p

        call check_inputs(sys, size(x0), size(c0), size(c0), "homotopy_plant")
        allocate (owner(sys%ncoef))
        owner = 0
        do i = 1, sys%n
            do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                p = sys%coef(k)
                if (owner(p) == 0) then
                    owner(p) = i
                else if (owner(p) /= i) then
                    owner(p) = -1
                end if
            end do
        end do
        do i = 1, sys%n
            if (owner(sys%coef(sys%eq_ptr(i))) /= i) error stop &
                "homotopy_plant: first coefficient of an equation is shared"
        end do
        call rng_seed(g, int(seed, int64), st)
        do j = 1, sys%n
            call rng_uniform(g, u)
            x0(j) = exp(cmplx(0.0_dp, TWO_PI * u, dp))
        end do
        do k = 1, sys%ncoef
            call rng_normal(g, re)
            call rng_normal(g, im)
            c0(k) = cmplx(re, im, dp) / sqrt(2.0_dp)
        end do
        do i = 1, sys%n
            p = sys%coef(sys%eq_ptr(i))
            a = (0.0_dp, 0.0_dp)
            b = (0.0_dp, 0.0_dp)
            do k = sys%eq_ptr(i), sys%eq_ptr(i + 1) - 1
                mono = term_value(sys, k, x0)
                if (sys%coef(k) == p) then
                    a = a + mono
                else
                    b = b + c0(sys%coef(k)) * mono
                end if
            end do
            c0(p) = -b / a
        end do
        c0(sys%ncoef + 1:) = (0.0_dp, 0.0_dp)
    end subroutine homotopy_plant

    ! Track every column of sols from ca to cb in parallel. ok(i) reports
    ! whether path i reached cb, acc(i) its final relative Newton correction.
    subroutine homotopy_move(sys, sols, nsol, ca, cb, opts, out, ok, acc, stats)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: sols(:, :)
        integer, intent(in) :: nsol
        complex(dp), intent(in) :: ca(:), cb(:)
        type(homotopy_options_t), intent(in) :: opts
        complex(dp), intent(out) :: out(:, :)
        logical, intent(out) :: ok(:)
        real(dp), intent(out), optional :: acc(:)
        type(homotopy_stats_t), intent(out), optional :: stats

        real(dp), allocatable :: a(:), as(:)
        integer, allocatable :: steps(:), sub(:), ss(:)
        integer(int64), allocatable :: facts(:), fs(:)
        complex(dp), allocatable :: xs(:, :), os(:, :)
        logical, allocatable :: member(:), oks(:)
        type(homotopy_options_t) :: careful
        integer :: round, njump, i

        call check_inputs(sys, size(sols, 1), size(ca), size(cb), "homotopy_move")
        if (size(sols, 2) < nsol .or. size(out, 2) < nsol .or. size(ok) < nsol) &
            error stop "homotopy_move: arrays shorter than nsol"
        if (size(out, 1) /= sys%n) error stop "homotopy_move: bad size(out, 1)"
        allocate (a(nsol), steps(nsol), facts(nsol))
        !$omp parallel default(shared)
        call move_worker(sys, sols, nsol, ca, cb, opts, out, ok, a, steps, facts)
        !$omp end parallel
        ! Coinciding endpoints mean that a path jumped to another branch: the
        ! straight segment maps the start fibre one to one. Every member of a
        ! coinciding group is tracked again with a smaller step target.
        njump = 0
        do round = 1, JUMP_ROUNDS
            call coinciding(sys%n, out, ok, nsol, opts%dup_tol, member, njump)
            if (njump == 0) exit
            careful = opts
            careful%step_target = opts%step_target / 100.0_dp**round
            careful%max_first_correction = opts%max_first_correction / 10.0_dp**round
            careful%max_steps = opts%max_steps * 10
            sub = pack([(i, i=1, nsol)], member(1:nsol))
            allocate (xs(sys%n, size(sub)), os(sys%n, size(sub)), oks(size(sub)), &
                as(size(sub)), ss(size(sub)), fs(size(sub)))
            xs = sols(1:sys%n, sub)
            !$omp parallel default(shared)
            call move_worker(sys, xs, size(sub), ca, cb, careful, os, oks, as, ss, fs)
            !$omp end parallel
            out(1:sys%n, sub) = os
            ok(sub) = oks
            a(sub) = as
            steps(sub) = steps(sub) + ss
            facts(sub) = facts(sub) + fs
            deallocate (xs, os, oks, as, ss, fs)
        end do
        if (present(acc)) acc(1:nsol) = a
        if (present(stats)) then
            stats%loops = 0
            stats%tracks = nsol
            stats%failed = count(.not. ok(1:nsol))
            stats%steps = sum(int(steps, int64))
            stats%factorizations = sum(facts)
            stats%jumps = njump
        end if
    end subroutine homotopy_move

    ! member(i) is true when endpoint i coincides with another successful
    ! endpoint within the relative distance tol; njump counts such endpoints.
    subroutine coinciding(n, out, ok, nsol, tol, member, njump)
        integer, intent(in) :: n, nsol
        complex(dp), intent(in) :: out(:, :)
        logical, intent(in) :: ok(:)
        real(dp), intent(in) :: tol
        logical, allocatable, intent(out) :: member(:)
        integer, intent(out) :: njump

        real(dp), allocatable :: w1(:), keys(:)
        complex(dp), allocatable :: w2(:)
        integer, allocatable :: idx(:), perm(:)
        real(dp) :: window
        integer :: j, k, m

        allocate (member(nsol), w1(n), w2(n))
        member = .false.
        njump = 0
        ! Fixed weights in [1/2, 1] and fixed phases; any choice gives a valid
        ! key, and the window bounds its change for coinciding points.
        do j = 1, n
            w1(j) = 0.5_dp + 0.5_dp * modulo(0.6180339887498949_dp * j, 1.0_dp)
            w2(j) = exp(cmplx(0.0_dp, TWO_PI * modulo(0.4142135623730951_dp * j, 1.0_dp), dp))
        end do
        window = key_window(w1, tol)
        idx = pack([(j, j=1, nsol)], ok(1:nsol))
        m = size(idx)
        allocate (keys(m), perm(m))
        do j = 1, m
            keys(j) = point_key(out(1:n, idx(j)), w1, w2)
        end do
        call argsort(keys, perm)
        do j = 1, m
            do k = j + 1, m
                if (keys(perm(k)) - keys(perm(j)) > window) exit
                if (same_point(out(1:n, idx(perm(j))), out(1:n, idx(perm(k))), tol)) then
                    member(idx(perm(j))) = .true.
                    member(idx(perm(k))) = .true.
                end if
            end do
        end do
        njump = count(member)
    end subroutine coinciding

    subroutine move_worker(sys, sols, nsol, ca, cb, opts, out, ok, acc, steps, &
            facts)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: sols(:, :)
        integer, intent(in) :: nsol
        complex(dp), intent(in) :: ca(:), cb(:)
        type(homotopy_options_t), intent(in) :: opts
        complex(dp), intent(inout) :: out(:, :)
        logical, intent(inout) :: ok(:)
        real(dp), intent(inout) :: acc(:)
        integer, intent(inout) :: steps(:)
        integer(int64), intent(inout) :: facts(:)

        type(track_work_t) :: w
        complex(dp), allocatable :: x(:)
        integer(int64) :: f0
        integer :: i, nq
        real(dp) :: tr, sr

        call work_init(w, sys)
        allocate (x(sys%n))
        !$omp do schedule(dynamic)
        do i = 1, nsol
            x = sols(1:sys%n, i)
            f0 = w%nfact
            call track_segment(sys, x, ca, cb, opts, w, ok(i), steps(i), acc(i), &
                tr, sr)
            facts(i) = w%nfact - f0
            if (.not. ok(i) .and. opts%quad_fallback) then
                ! Resume from the last point accepted in double precision.
                call track_segment_qp(sys, x, ca, cb, opts, ok(i), nq, acc(i), tr, sr)
                steps(i) = steps(i) + nq
            end if
            out(1:sys%n, i) = x
        end do
        !$omp end do
    end subroutine move_worker

    ! Fill the fibre over c0 by monodromy loops c0 -> c1 -> c2 -> c0 with fresh
    ! random c1, c2, starting from one solution x0. Every loop tracks all
    ! current solutions in parallel; endpoints not yet known are added.
    subroutine homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status, stats)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: x0(:), c0(:)
        type(homotopy_options_t), intent(in) :: opts
        complex(dp), allocatable, intent(out) :: sols(:, :)
        integer, intent(out) :: nsol
        type(fortnum_status_t), intent(out) :: status
        type(homotopy_stats_t), intent(out), optional :: stats

        type(track_work_t) :: w
        type(rng_t) :: g, gl
        type(fortnum_status_t) :: rst
        type(homotopy_stats_t) :: st
        complex(dp), allocatable :: c1(:), c2(:), cand(:, :), wk2(:), x(:)
        logical, allocatable :: good(:), hit(:)
        real(dp), allocatable :: acc(:), wk1(:), skeys(:), ckeys(:), tkeys(:)
        integer, allocatable :: sidx(:), tidx(:), steps(:), nseg(:), perm(:)
        integer, allocatable :: cidx(:), newpos(:)
        integer(int64), allocatable :: facts(:)
        real(dp) :: u, re, im, window, first, last
        integer :: n, i, j, k, m, loop, stagnant, nnew, ngood
        logical :: cok

        call check_inputs(sys, size(x0), size(c0), size(c0), "homotopy_monodromy")
        n = sys%n
        nsol = 0
        call rng_seed(g, int(opts%seed, int64), rst)
        call rng_split(g, 1_int64, gl, rst)
        allocate (wk1(n), wk2(n), c1(sys%ncoef), c2(sys%ncoef))
        do j = 1, n
            call rng_uniform(gl, u)
            wk1(j) = 0.5_dp + 0.5_dp * u
            call rng_uniform(gl, u)
            wk2(j) = exp(cmplx(0.0_dp, TWO_PI * u, dp))
        end do
        window = key_window(wk1, opts%dup_tol)

        call work_init(w, sys)
        x = x0(1:n)
        call correct(sys, x, c0, opts%final_newton_iters, opts%final_tol, &
            huge(1.0_dp), w, cok, first, last)
        if (.not. cok) then
            allocate (sols(n, 0))
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "homotopy_monodromy: x0 is not a regular solution at c0")
            if (present(stats)) stats = st
            return
        end if
        allocate (sols(n, max(16, opts%target)), skeys(1), sidx(1))
        sols(:, 1) = x
        nsol = 1
        skeys(1) = point_key(x, wk1, wk2)
        sidx(1) = 1

        loop = 0
        stagnant = 0
        do
            if (opts%target > 0 .and. nsol >= opts%target) exit
            if (loop >= opts%max_loops) exit
            if (stagnant >= max(1, opts%stagnant_loops)) exit
            loop = loop + 1
            call rng_split(g, 2_int64 * loop, gl, rst)
            do k = 1, sys%ncoef
                call rng_normal(gl, re)
                call rng_normal(gl, im)
                c1(k) = cmplx(re, im, dp) / sqrt(2.0_dp)
            end do
            do k = 1, sys%ncoef
                call rng_normal(gl, re)
                call rng_normal(gl, im)
                c2(k) = cmplx(re, im, dp) / sqrt(2.0_dp)
            end do

            m = nsol
            allocate (cand(n, m), good(m), acc(m), steps(m), nseg(m), facts(m))
            !$omp parallel default(shared)
            call loop_worker(sys, sols, m, c0, c1, c2, opts, cand, good, acc, &
                steps, nseg, facts)
            !$omp end parallel
            st%tracks = st%tracks + sum(int(nseg, int64))
            st%failed = st%failed + count(.not. good)
            st%steps = st%steps + sum(int(steps, int64))
            st%factorizations = st%factorizations + sum(facts)

            ! Deduplicate in key order against the known set and among the
            ! new points; keys within the window are compared coordinatewise.
            ngood = count(good)
            allocate (ckeys(ngood), cidx(ngood), perm(ngood), newpos(ngood))
            j = 0
            do i = 1, m
                if (.not. good(i)) cycle
                j = j + 1
                cidx(j) = i
                ckeys(j) = point_key(cand(:, i), wk1, wk2)
            end do
            call argsort(ckeys, perm)
            allocate (hit(nsol))
            hit = .false.
            nnew = 0
            do j = 1, ngood
                i = cidx(perm(j))
                if (is_known(cand(:, i), ckeys(perm(j)))) cycle
                nnew = nnew + 1
                newpos(nnew) = perm(j)
            end do
            if (nsol + nnew > size(sols, 2)) call grow_sols(sols, nsol + nnew)
            allocate (tkeys(nsol + nnew), tidx(nsol + nnew))
            call merge_new()
            call move_alloc(tkeys, skeys)
            call move_alloc(tidx, sidx)
            if (opts%verbosity > 0) then
                write (output_unit, '(a,i0,a,i0,a,i0,a,i0,a,f0.1)') &
                    "homotopy_monodromy: loop ", loop, ", solutions ", nsol, &
                    " (new ", nnew, ", failed ", count(.not. good), &
                    "), mean steps/segment ", &
                    real(sum(int(steps, int64)), dp) / real(max(1, sum(nseg)), dp)
            end if
            if (nnew == 0) then
                stagnant = stagnant + 1
            else
                stagnant = 0
            end if
            deallocate (cand, good, acc, steps, nseg, facts, ckeys, cidx, perm, &
                newpos, hit)
        end do
        st%loops = loop
        call grow_sols(sols, nsol, exact=.true.)

        if (opts%target > 0) then
            if (nsol == opts%target) then
                call status_set(status, FORTNUM_OK, "")
            else if (nsol > opts%target) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "homotopy_monodromy: more solutions than the target count")
            else
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "homotopy_monodromy: target count not reached")
            end if
        else if (stagnant >= max(1, opts%stagnant_loops)) then
            call status_set(status, FORTNUM_OK, "")
        else
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "homotopy_monodromy: max_loops reached while still finding solutions")
        end if
        if (present(stats)) stats = st

    contains

        logical function is_known(z, key) result(known)
            complex(dp), intent(in) :: z(:)
            real(dp), intent(in) :: key

            integer :: p, q

            known = .true.
            p = lower_bound(skeys, nsol, key - window)
            do while (p <= nsol)
                if (skeys(p) > key + window) exit
                if (same_point(z, sols(:, sidx(p)), opts%dup_tol)) then
                    if (hit(sidx(p))) st%jumps = st%jumps + 1
                    hit(sidx(p)) = .true.
                    return
                end if
                p = p + 1
            end do
            do q = nnew, 1, -1
                if (ckeys(newpos(q)) < key - window) exit
                if (same_point(z, cand(:, cidx(newpos(q))), opts%dup_tol)) then
                    st%jumps = st%jumps + 1
                    return
                end if
            end do
            known = .false.
        end function is_known

        ! Append the new points and merge their sorted keys into skeys.
        subroutine merge_new()
            integer :: p, q, r

            p = 1
            q = 1
            do r = 1, nsol + nnew
                if (q > nnew) then
                    tkeys(r) = skeys(p)
                    tidx(r) = sidx(p)
                    p = p + 1
                    cycle
                end if
                if (p <= nsol) then
                    if (skeys(p) <= ckeys(newpos(q))) then
                        tkeys(r) = skeys(p)
                        tidx(r) = sidx(p)
                        p = p + 1
                        cycle
                    end if
                end if
                tkeys(r) = ckeys(newpos(q))
                tidx(r) = nsol + q
                sols(:, nsol + q) = cand(:, cidx(newpos(q)))
                q = q + 1
            end do
            nsol = nsol + nnew
        end subroutine merge_new

    end subroutine homotopy_monodromy

    subroutine loop_worker(sys, sols, m, c0, c1, c2, opts, cand, good, acc, &
            steps, nseg, facts)
        type(poly_system_t), intent(in) :: sys
        complex(dp), intent(in) :: sols(:, :)
        integer, intent(in) :: m
        complex(dp), intent(in) :: c0(:), c1(:), c2(:)
        type(homotopy_options_t), intent(in) :: opts
        complex(dp), intent(inout) :: cand(:, :)
        logical, intent(inout) :: good(:)
        real(dp), intent(inout) :: acc(:)
        integer, intent(inout) :: steps(:), nseg(:)
        integer(int64), intent(inout) :: facts(:)

        type(track_work_t) :: w
        complex(dp), allocatable :: x(:)
        integer(int64) :: f0
        integer :: i, ns
        logical :: ok

        call work_init(w, sys)
        allocate (x(sys%n))
        !$omp do schedule(dynamic)
        do i = 1, m
            x = sols(1:sys%n, i)
            f0 = w%nfact
            call track_segment(sys, x, c0, c1, opts, w, ok, steps(i), acc(i))
            nseg(i) = 1
            if (ok) then
                call track_segment(sys, x, c1, c2, opts, w, ok, ns, acc(i))
                steps(i) = steps(i) + ns
                nseg(i) = 2
            end if
            if (ok) then
                call track_segment(sys, x, c2, c0, opts, w, ok, ns, acc(i))
                steps(i) = steps(i) + ns
                nseg(i) = 3
            end if
            ! An endpoint must be accurate enough to be told apart from others.
            if (ok) ok = acc(i) <= opts%dup_tol
            good(i) = ok
            facts(i) = w%nfact - f0
            cand(1:sys%n, i) = x
        end do
        !$omp end do
    end subroutine loop_worker

    ! A real key that moves by at most key_window when every coordinate moves
    ! by a relative distance tol: log-moduli plus a projection of the phases,
    ! |log|a| - log|b|| <= tol / (1 - tol), |a/|a| - b/|b|| <= 2 tol.
    pure real(dp) function point_key(x, w1, w2) result(key)
        complex(dp), intent(in) :: x(:)
        real(dp), intent(in) :: w1(:)
        complex(dp), intent(in) :: w2(:)

        real(dp) :: a
        integer :: j

        key = 0.0_dp
        do j = 1, size(x)
            a = abs(x(j))
            if (a > TINY_SCALE) then
                key = key + w1(j) * log(a) + (real(w2(j), dp) * real(x(j), dp) &
                    - aimag(w2(j)) * aimag(x(j))) / a
            else
                key = key + w1(j) * log(TINY_SCALE)
            end if
        end do
    end function point_key

    pure real(dp) function key_window(w1, tol) result(window)
        real(dp), intent(in) :: w1(:), tol

        real(dp) :: lip, mag

        lip = sum(w1) / (1.0_dp - tol) + 2.0_dp * size(w1)
        ! Rounding in the key: each log term is at most 700 in magnitude.
        mag = sum(w1) * 700.0_dp + size(w1)
        window = 1.01_dp * lip * tol + 64.0_dp * epsilon(1.0_dp) * mag
    end function key_window

    pure logical function same_point(a, b, tol) result(same)
        complex(dp), intent(in) :: a(:), b(:)
        real(dp), intent(in) :: tol

        integer :: j

        same = .false.
        do j = 1, size(a)
            if (abs(a(j) - b(j)) > tol * max(abs(a(j)), abs(b(j)))) return
        end do
        same = .true.
    end function same_point

    ! First index p in 1..m with keys(p) >= v, or m + 1.
    pure integer function lower_bound(keys, m, v) result(p)
        real(dp), intent(in) :: keys(:)
        integer, intent(in) :: m
        real(dp), intent(in) :: v

        integer :: lo, hi, mid

        lo = 1
        hi = m + 1
        do while (lo < hi)
            mid = (lo + hi) / 2
            if (keys(mid) < v) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        p = lo
    end function lower_bound

    subroutine grow_sols(sols, need, exact)
        complex(dp), allocatable, intent(inout) :: sols(:, :)
        integer, intent(in) :: need
        logical, intent(in), optional :: exact

        complex(dp), allocatable :: tmp(:, :)
        integer :: cap, keep

        cap = max(need, 2 * size(sols, 2))
        if (present(exact)) then
            if (exact) cap = need
        end if
        keep = min(cap, size(sols, 2))
        allocate (tmp(size(sols, 1), cap))
        tmp(:, 1:keep) = sols(:, 1:keep)
        call move_alloc(tmp, sols)
    end subroutine grow_sols

    subroutine grow_int(a, cap)
        integer, allocatable, intent(inout) :: a(:)
        integer, intent(in) :: cap

        integer, allocatable :: tmp(:)

        allocate (tmp(cap))
        tmp(1:size(a)) = a
        call move_alloc(tmp, a)
    end subroutine grow_int

    subroutine check_inputs(sys, nx, nca, ncb, caller)
        type(poly_system_t), intent(in) :: sys
        integer, intent(in) :: nx, nca, ncb
        character(*), intent(in) :: caller

        if (.not. sys%finalized) error stop "homotopy: system not finalised"
        if (nx /= sys%n) then
            write (output_unit, '(a)') caller
            error stop "homotopy: point size differs from sys%n"
        end if
        if (nca < sys%ncoef .or. ncb < sys%ncoef) then
            write (output_unit, '(a)') caller
            error stop "homotopy: coefficient vector shorter than sys%ncoef"
        end if
    end subroutine check_inputs

end module fortnum_homotopy
