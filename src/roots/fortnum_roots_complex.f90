module fortnum_roots_complex
    ! Clean-room complex analytic-zero finder over an axis-aligned rectangle.
    !
    ! Replaces the ZEAL ICON=3 region search used by KIM wkb_dispersion
    ! (WKB_dispersion_solver='ZEAL').  Given an analytic f and a rectangle
    ! [ll, ur] in the complex plane, returns the distinct zeros inside, their
    ! function values, and their multiplicities.
    !
    ! METHOD (argument principle + formal orthogonal polynomials, NOT ZEAL
    ! source).  References:
    !   Kravanja & Van Barel, "Computing the Zeros of Analytic Functions",
    !     Lecture Notes in Mathematics 1727, Springer (2000).
    !   Kravanja, Van Barel, Ragos, Vrahatis, Zafiropoulos, "ZEAL: a
    !     mathematical software package ...", Comput. Phys. Commun. 124
    !     (2000) 212-232.
    !   Henrici, "Applied and Computational Complex Analysis", vol. 1 (1974),
    !     ch. 4 (argument principle, Newton sums of zeros).
    !
    ! The number of zeros counted with multiplicity in the box is the winding
    ! number N = (1/2 pi i) oint_C f'/f dz.  The ordinary moments
    !   s_p = (1/2 pi i) oint_C z^p f'/f dz = sum_k m_k z_k^p,  p = 0,1,2,...
    ! are the Newton power sums of the zeros z_k with multiplicities m_k.  The
    ! distinct zeros are the generalized eigenvalues of the Hankel pencil
    ! (H<, H) with H = [s_{i+j}], H< = [s_{i+j+1}] (Prony / formal orthogonal
    ! polynomials).  Multiplicities follow from the confluent Vandermonde
    ! system sum_k m_k z_k^p = s_p and are rounded to positive integers.  Each
    ! eigenvalue is polished by modified complex Newton.  When N exceeds m_max
    ! the box is bisected (ZEAL ICON=3 recursion).
    !
    ! The contour integrals run on the four rectangle edges with the package
    ! Gauss-Kronrod driver (fortnum_integrate_gk), real and imaginary parts
    ! integrated separately.  f'/f uses a complex central finite difference of
    ! f; the logarithmic derivative is never formed across a branch cut.
    !
    ! DERIVATIVE POLICY (ad.md): implicit_rule with differentiate_through=false.
    !   A zero satisfies f(z*, p) = 0, so the implicit function theorem gives
    !   dz*/dp = -f_p / f'(z*) without differentiating the contour integral or
    !   the eigensolve.  No consumer differentiates through this region search
    !   (differentiate_through=false), so no JVP/VJP is provided here; the
    !   scalar implicit rule in fortnum_roots covers the per-root sensitivity.
    !   Inactive: m_max, tolerances, max_refine, status, multiplicities.
    !
    ! The small generalized eigenproblem (size <= m_max) is solved with LAPACK
    ! ZGGEV.  No module-level save, no global procedure pointers; integrand
    ! state is carried by host association in internal procedures.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_integrate_gk, only: integrate_gk
    implicit none
    private

    public :: complex_root_fn_t, complex_region_roots

    ! 2 pi as a working constant for the 1/(2 pi i) contour normalisation.
    real(dp), parameter :: two_pi = 6.283185307179586476925286766559_dp

    ! Analytic function f whose zeros are sought.  kr in, fk = f(kr) out;
    ! ctx carries parameters (active for the implicit rule, never inspected
    ! by the finder).
    abstract interface
        subroutine complex_root_fn_t(kr, fk, ctx)
            import :: dp
            complex(dp), intent(in)  :: kr
            complex(dp), intent(out) :: fk
            class(*),    intent(in), optional :: ctx
        end subroutine complex_root_fn_t
    end interface

    ! Explicit interface for the LAPACK double-complex generalized eigensolver
    ! (the only LAPACK routine used here).  Declared so the call type-checks
    ! under -Werror=implicit-interface; LAPACK itself is linked, not vendored.
    interface
        subroutine zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, &
                vl, ldvl, vr, ldvr, work, lwork, rwork, info)
            import :: dp
            character,   intent(in)    :: jobvl, jobvr
            integer,     intent(in)    :: n, lda, ldb, ldvl, ldvr, lwork
            complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
            complex(dp), intent(out)   :: alpha(*), beta(*)
            complex(dp), intent(out)   :: vl(ldvl, *), vr(ldvr, *)
            complex(dp), intent(out)   :: work(*)
            real(dp),    intent(out)   :: rwork(*)
            integer,     intent(out)   :: info
        end subroutine zggev
    end interface

contains

    ! Find the distinct zeros of f in the rectangle [ll, ur].
    !
    ! ll, ur     lower-left / upper-right corners (complex(dp)).
    ! roots,fvals distinct zeros and f at them (allocatable out).
    ! mult       multiplicity of each distinct zero (allocatable out).
    ! nfound     number of distinct zeros.
    ! status     FORTNUM_OK, or an error code on failure.
    ! m_max      max zeros (counting multiplicity) handled per subregion
    !            before bisecting (ZEAL M); default 5.
    ! newtonz,newtonf  Newton stopping tolerances on |dz| and |f|.
    ! max_refine maximum Newton iterations per root; default 60.
    ! ctx        parameters forwarded to f.
    subroutine complex_region_roots(f, ll, ur, roots, fvals, mult, nfound, &
            status, m_max, newtonz, newtonf, max_refine, ctx)
        procedure(complex_root_fn_t)            :: f
        complex(dp),             intent(in)     :: ll, ur
        complex(dp), allocatable, intent(out)   :: roots(:), fvals(:)
        integer,     allocatable, intent(out)   :: mult(:)
        integer,                 intent(out)    :: nfound
        type(fortnum_status_t),  intent(out)    :: status
        integer,     intent(in), optional       :: m_max, max_refine
        real(dp),    intent(in), optional       :: newtonz, newtonf
        class(*),    intent(in), optional       :: ctx

        integer  :: mm, mref
        real(dp) :: nz, nf

        mm   = 5
        nz   = 5.0e-8_dp
        nf   = 1.0e-14_dp
        mref = 60
        if (present(m_max))     mm   = m_max
        if (present(newtonz))   nz   = newtonz
        if (present(newtonf))   nf   = newtonf
        if (present(max_refine)) mref = max_refine

        allocate(roots(0), fvals(0), mult(0))
        nfound = 0
        call status_set(status, FORTNUM_OK, "")

        if (real(ur) <= real(ll) .or. aimag(ur) <= aimag(ll)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "complex_region_roots: degenerate box (ur must exceed ll)")
            return
        end if
        if (mm < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "complex_region_roots: m_max must be >= 1")
            return
        end if

        call region_recurse(f, ll, ur, mm, nz, nf, mref, ctx, &
            roots, fvals, mult, nfound, status, 0)
    end subroutine complex_region_roots

    ! Recursive ICON=3 region search.  Counts zeros in the box; if the count
    ! exceeds m_max, bisects along the longer side and recurses; otherwise
    ! recovers the zeros from the moment pencil and appends them.
    recursive subroutine region_recurse(f, ll, ur, mm, nz, nf, mref, ctx, &
            roots, fvals, mult, nfound, status, depth)
        procedure(complex_root_fn_t)            :: f
        complex(dp),             intent(in)     :: ll, ur
        integer,                 intent(in)     :: mm, mref, depth
        real(dp),                intent(in)     :: nz, nf
        class(*),    intent(in), optional       :: ctx
        complex(dp), allocatable, intent(inout) :: roots(:), fvals(:)
        integer,     allocatable, intent(inout) :: mult(:)
        integer,                 intent(inout)  :: nfound
        type(fortnum_status_t),  intent(inout)  :: status

        integer, parameter :: max_depth = 40
        complex(dp) :: s(0:2*mm)
        integer     :: ntot

        call moments(f, ll, ur, 2*mm + 1, ctx, s, ntot, status)
        if (status%code /= FORTNUM_OK) return
        if (ntot == 0) return

        if (ntot > mm) then
            if (depth >= max_depth) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "complex_region_roots: subdivision depth limit; box has "// &
                    "more zeros than m_max in an irreducible cell")
                return
            end if
            call bisect_box(f, ll, ur, ntot, mm, nz, nf, mref, ctx, &
                roots, fvals, mult, nfound, status, depth)
            return
        end if

        call resolve_cell(f, ll, ur, s, ntot, mm, nz, nf, mref, ctx, &
            roots, fvals, mult, nfound, status)
    end subroutine region_recurse

    ! Split the box along its longer side and recurse into both halves.
    ! The cut must not pass through a zero (that puts a zero on the new shared
    ! edge and corrupts both winding numbers), so a small set of cut fractions
    ! near the centre is tried until both halves report integer winding numbers
    ! summing to ntot.  ZEAL perturbs ICON=3 cuts the same way.
    recursive subroutine bisect_box(f, ll, ur, ntot, mm, nz, nf, mref, ctx, &
            roots, fvals, mult, nfound, status, depth)
        procedure(complex_root_fn_t)            :: f
        complex(dp),             intent(in)     :: ll, ur
        integer,                 intent(in)     :: ntot, mm, mref, depth
        real(dp),                intent(in)     :: nz, nf
        class(*),    intent(in), optional       :: ctx
        complex(dp), allocatable, intent(inout) :: roots(:), fvals(:)
        integer,     allocatable, intent(inout) :: mult(:)
        integer,                 intent(inout)  :: nfound
        type(fortnum_status_t),  intent(inout)  :: status

        real(dp), parameter :: frac(5) = &
            [0.5_dp, 0.532_dp, 0.468_dp, 0.571_dp, 0.443_dp]
        logical     :: horiz
        complex(dp) :: lo1, hi1, lo2, hi2
        integer     :: n1, n2, t
        real(dp)    :: cut

        horiz = (real(ur) - real(ll) >= aimag(ur) - aimag(ll))
        do t = 1, size(frac)
            if (horiz) then
                cut = real(ll) + frac(t) * (real(ur) - real(ll))
                lo1 = ll;  hi1 = cmplx(cut, aimag(ur), dp)
                lo2 = cmplx(cut, aimag(ll), dp);  hi2 = ur
            else
                cut = aimag(ll) + frac(t) * (aimag(ur) - aimag(ll))
                lo1 = ll;  hi1 = cmplx(real(ur), cut, dp)
                lo2 = cmplx(real(ll), cut, dp);  hi2 = ur
            end if
            call count_only(f, lo1, hi1, ctx, n1, status)
            if (status%code /= FORTNUM_OK) cycle
            call count_only(f, lo2, hi2, ctx, n2, status)
            if (status%code /= FORTNUM_OK) cycle
            if (n1 + n2 /= ntot) cycle
            ! Both halves report integer winding numbers summing to ntot, so
            ! the cut is clean: recurse into each half and return.
            call region_recurse(f, lo1, hi1, mm, nz, nf, mref, ctx, &
                roots, fvals, mult, nfound, status, depth + 1)
            if (status%code /= FORTNUM_OK) return
            call region_recurse(f, lo2, hi2, mm, nz, nf, mref, ctx, &
                roots, fvals, mult, nfound, status, depth + 1)
            return
        end do
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "complex_region_roots: no clean subdivision cut found")
    end subroutine bisect_box

    ! Winding number alone (s_0) for a candidate sub-box, used to validate a
    ! subdivision cut before committing to the recursion.
    subroutine count_only(f, ll, ur, ctx, n, status)
        procedure(complex_root_fn_t)        :: f
        complex(dp),         intent(in)     :: ll, ur
        class(*),  intent(in), optional     :: ctx
        integer,             intent(out)    :: n
        type(fortnum_status_t), intent(out) :: status
        complex(dp) :: s(0:0)

        call moments(f, ll, ur, 1, ctx, s, n, status)
    end subroutine count_only

    ! Recover the <= mm zeros of a cell from its moments s(0:ntot..) and
    ! append the polished distinct zeros to the accumulators.
    subroutine resolve_cell(f, ll, ur, s, ntot, mm, nz, nf, mref, ctx, &
            roots, fvals, mult, nfound, status)
        procedure(complex_root_fn_t)            :: f
        complex(dp),             intent(in)     :: ll, ur
        complex(dp),             intent(in)     :: s(0:)
        integer,                 intent(in)     :: ntot, mm, mref
        real(dp),                intent(in)     :: nz, nf
        class(*),    intent(in), optional       :: ctx
        complex(dp), allocatable, intent(inout) :: roots(:), fvals(:)
        integer,     allocatable, intent(inout) :: mult(:)
        integer,                 intent(inout)  :: nfound
        type(fortnum_status_t),  intent(inout)  :: status

        complex(dp) :: zk(mm), fk
        integer     :: mk(mm), ndist, k

        call distinct_zeros(s, ntot, mm, zk, ndist, status)
        if (status%code /= FORTNUM_OK) return

        call multiplicities(s, zk, ndist, ntot, mk, status)
        if (status%code /= FORTNUM_OK) return

        do k = 1, ndist
            call newton_polish(f, zk(k), ll, ur, nz, nf, mref, mk(k), ctx, fk, status)
            if (status%code /= FORTNUM_OK) return
            call append_root(roots, fvals, mult, nfound, zk(k), fk, mk(k))
        end do
    end subroutine resolve_cell

    ! Ordinary moments s_p = (1/2 pi i) oint_C z^p f'/f dz for p = 0..np-1.
    ! ntot returns s_0 rounded to the nearest integer (total zero count).
    subroutine moments(f, ll, ur, np, ctx, s, ntot, status)
        procedure(complex_root_fn_t)        :: f
        complex(dp),         intent(in)     :: ll, ur
        integer,             intent(in)     :: np
        class(*),  intent(in), optional     :: ctx
        complex(dp),         intent(out)    :: s(0:)
        integer,             intent(out)    :: ntot
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: corner(4), edge
        integer     :: p, e
        real(dp)    :: re_n

        call status_set(status, FORTNUM_OK, "")
        corner(1) = ll
        corner(2) = cmplx(real(ur), aimag(ll), dp)
        corner(3) = ur
        corner(4) = cmplx(real(ll), aimag(ur), dp)

        do p = 0, np - 1
            s(p) = (0.0_dp, 0.0_dp)
            do e = 1, 4
                call edge_integral(f, corner(e), corner(mod(e, 4) + 1), p, &
                    ctx, edge, status)
                if (status%code /= FORTNUM_OK) return
                s(p) = s(p) + edge
            end do
            s(p) = s(p) / cmplx(0.0_dp, two_pi, dp)
        end do

        re_n = real(s(0))
        ntot = nint(re_n)
        if (ntot < 0 .or. abs(re_n - real(ntot, dp)) > 0.25_dp) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "complex_region_roots: winding number not near an integer; "// &
                "f may have a zero or pole on the contour")
            ntot = max(ntot, 0)
        end if
    end subroutine moments

    ! Integral of z^p f'(z)/f(z) along the straight edge za -> zb.
    ! Parametrise z(t) = za + t (zb - za), t in [0,1]; the integrand carries
    ! the dz = (zb - za) dt factor.  Real and imaginary parts integrate
    ! separately with the package Gauss-Kronrod driver.
    subroutine edge_integral(f, za, zb, p, ctx, val, status)
        procedure(complex_root_fn_t)        :: f
        complex(dp),         intent(in)     :: za, zb
        integer,             intent(in)     :: p
        class(*),  intent(in), optional     :: ctx
        complex(dp),         intent(out)    :: val
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: dzdt
        real(dp)    :: rr, ri, aerr
        integer     :: ierr

        call status_set(status, FORTNUM_OK, "")
        dzdt = zb - za

        call integrate_gk(g_re, 0.0_dp, 1.0_dp, 1.0e-12_dp, 1.0e-11_dp, &
            rr, aerr, ierr, key=31, limit=400)
        if (ierr /= 0 .and. ierr /= 2) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "complex_region_roots: edge integral (real part) failed")
            return
        end if
        call integrate_gk(g_im, 0.0_dp, 1.0_dp, 1.0e-12_dp, 1.0e-11_dp, &
            ri, aerr, ierr, key=31, limit=400)
        if (ierr /= 0 .and. ierr /= 2) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "complex_region_roots: edge integral (imag part) failed")
            return
        end if
        val = cmplx(rr, ri, dp)

    contains

        ! Full complex integrand value at parameter t.
        function integrand(t) result(gz)
            real(dp), intent(in) :: t
            complex(dp) :: gz, z, lfp, zp
            z   = za + cmplx(t, 0.0_dp, dp) * dzdt
            lfp = log_deriv(f, z, ctx)
            if (p == 0) then
                zp = (1.0_dp, 0.0_dp)
            else
                zp = z ** p
            end if
            gz  = zp * lfp * dzdt
        end function integrand

        function g_re(t) result(y)
            real(dp), intent(in) :: t
            real(dp) :: y
            y = real(integrand(t))
        end function g_re

        function g_im(t) result(y)
            real(dp), intent(in) :: t
            real(dp) :: y
            y = aimag(integrand(t))
        end function g_im
    end subroutine edge_integral

    ! Logarithmic derivative f'(z)/f(z) by a complex central finite
    ! difference of f.  The quotient f'/f is formed (never log f), so branch
    ! cuts of f play no role.  Step scales with |z| for relative accuracy.
    function log_deriv(f, z, ctx) result(lfp)
        procedure(complex_root_fn_t)    :: f
        complex(dp), intent(in)         :: z
        class(*),  intent(in), optional :: ctx
        complex(dp) :: lfp, fp, fm, fz, fpr
        real(dp)    :: h

        h = 1.0e-6_dp * max(abs(z), 1.0_dp)
        call f(z + cmplx(h, 0.0_dp, dp), fp, ctx)
        call f(z - cmplx(h, 0.0_dp, dp), fm, ctx)
        call f(z, fz, ctx)
        fpr = (fp - fm) / cmplx(2.0_dp*h, 0.0_dp, dp)
        if (fz == (0.0_dp, 0.0_dp)) then
            lfp = (0.0_dp, 0.0_dp)
        else
            lfp = fpr / fz
        end if
    end function log_deriv

    ! Distinct zeros from the Hankel moment pencil (H<, H).
    !
    ! ndist = rank of the Hankel matrix [s_{i+j}] (i,j = 0..ntot-1), found by
    ! the largest leading principal block that stays well conditioned.  The
    ! ndist distinct zeros are the generalized eigenvalues of the ndist-by-ndist
    ! pencil H< v = z H v, solved with LAPACK ZGGEV.
    subroutine distinct_zeros(s, ntot, mm, zk, ndist, status)
        complex(dp),         intent(in)     :: s(0:)
        integer,             intent(in)     :: ntot, mm
        complex(dp),         intent(out)    :: zk(:)
        integer,             intent(out)    :: ndist
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: hmat(mm, mm), hlt(mm, mm)
        integer     :: n, i, j

        call status_set(status, FORTNUM_OK, "")
        n = rank_hankel(s, ntot, mm)
        ndist = n
        if (n == 0) return

        do j = 1, n
            do i = 1, n
                hmat(i, j) = s(i + j - 2)
                hlt(i, j)  = s(i + j - 1)
            end do
        end do
        call solve_pencil(hlt(1:n, 1:n), hmat(1:n, 1:n), n, zk, status)
    end subroutine distinct_zeros

    ! Numerical rank of the leading Hankel blocks [s_{i+j}], capped at mm and
    ! at ntot.  The block of size n+1 is declared singular when its scaled
    ! determinant magnitude collapses relative to the size-n block.
    function rank_hankel(s, ntot, mm) result(n)
        complex(dp), intent(in) :: s(0:)
        integer,     intent(in) :: ntot, mm
        integer :: n, k, nmax
        complex(dp) :: blk(mm, mm)
        real(dp)    :: dcur, scal0
        integer     :: i, j

        nmax = min(ntot, mm)
        scal0 = abs(s(0))
        if (scal0 == 0.0_dp) scal0 = 1.0_dp
        n = 0
        do k = 1, nmax
            do j = 1, k
                do i = 1, k
                    blk(i, j) = s(i + j - 2)
                end do
            end do
            dcur = abs(det_lu(blk(1:k, 1:k), k))
            ! A genuine rank-k Hankel block keeps a determinant well above the
            ! roundoff floor set by the moment scale; once it collapses, the
            ! distinct-zero count is the previous k.
            if (dcur <= 1.0e-9_dp * scal0**k) exit
            n = k
        end do
        if (n == 0 .and. nmax >= 1) n = 1
    end function rank_hankel

    ! Confluent-Vandermonde multiplicities: solve sum_k m_k z_k^p = s_p for
    ! p = 0..ndist-1 by Gaussian elimination, round to positive integers,
    ! then rescale so they sum to the total winding number ntot.
    subroutine multiplicities(s, zk, ndist, ntot, mk, status)
        complex(dp),         intent(in)     :: s(0:)
        complex(dp),         intent(in)     :: zk(:)
        integer,             intent(in)     :: ndist, ntot
        integer,             intent(out)    :: mk(:)
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: vmat(ndist, ndist), rhs(ndist), msol(ndist)
        integer     :: i, k, ssum

        call status_set(status, FORTNUM_OK, "")
        if (ndist == 0) return
        do i = 1, ndist
            do k = 1, ndist
                vmat(i, k) = zk(k) ** (i - 1)
            end do
            rhs(i) = s(i - 1)
        end do
        call solve_lin(vmat, rhs, ndist, msol, status)
        if (status%code /= FORTNUM_OK) return

        ssum = 0
        do k = 1, ndist
            mk(k) = max(1, nint(real(msol(k))))
            ssum  = ssum + mk(k)
        end do
        ! Reconcile rounded multiplicities with the integer winding number:
        ! adjust the largest entry to absorb a unit residual when consistent.
        if (ssum /= ntot .and. ndist == 1) mk(1) = ntot
        if (ssum /= ntot .and. ndist > 1 .and. abs(ssum - ntot) <= ndist) then
            k = maxloc(mk(1:ndist), 1)
            mk(k) = mk(k) + (ntot - ssum)
            if (mk(k) < 1) mk(k) = 1
        end if
    end subroutine multiplicities

    ! Modified complex Newton on f, keeping the iterate inside [ll, ur].
    ! Uses the central-difference derivative; the multiplicity factor m
    ! restores quadratic convergence at a multiplicity-m zero, where the
    ! plain Newton step is only linearly convergent. The damping shrinks the
    ! step when a full step would leave the box.
    subroutine newton_polish(f, z, ll, ur, nz, nf, mref, m, ctx, fz, status)
        procedure(complex_root_fn_t)        :: f
        complex(dp),         intent(inout)  :: z
        complex(dp),         intent(in)     :: ll, ur
        real(dp),            intent(in)     :: nz, nf
        integer,             intent(in)     :: mref
        integer,             intent(in)     :: m
        class(*),  intent(in), optional     :: ctx
        complex(dp),         intent(out)    :: fz
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: zc, fpr, fm, fp, step, znew
        real(dp)    :: h, lam, refl, imfl, reup, imup, pad
        integer     :: it

        call status_set(status, FORTNUM_OK, "")
        zc = z
        refl = real(ll); reup = real(ur)
        imfl = aimag(ll); imup = aimag(ur)
        pad  = 0.05_dp * max(reup - refl, imup - imfl)

        do it = 1, mref
            call f(zc, fz, ctx)
            if (abs(fz) <= nf) then
                z = zc
                return
            end if
            h = 1.0e-7_dp * max(abs(zc), 1.0_dp)
            call f(zc + cmplx(h, 0.0_dp, dp), fp, ctx)
            call f(zc - cmplx(h, 0.0_dp, dp), fm, ctx)
            fpr = (fp - fm) / cmplx(2.0_dp*h, 0.0_dp, dp)
            if (fpr == (0.0_dp, 0.0_dp)) then
                z = zc
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "complex_region_roots: zero derivative during Newton polish")
                return
            end if
            step = real(max(m, 1), dp) * fz / fpr
            lam  = 1.0_dp
            do
                znew = zc - cmplx(lam, 0.0_dp, dp) * step
                if (real(znew) >= refl - pad .and. real(znew) <= reup + pad &
                    .and. aimag(znew) >= imfl - pad &
                    .and. aimag(znew) <= imup + pad) exit
                lam = 0.5_dp * lam
                if (lam < 1.0e-3_dp) exit
            end do
            znew = zc - cmplx(lam, 0.0_dp, dp) * step
            if (abs(znew - zc) <= nz) then
                zc = znew
                call f(zc, fz, ctx)
                z = zc
                return
            end if
            zc = znew
        end do
        z = zc
        call f(zc, fz, ctx)
        if (abs(fz) > max(nf, 1.0e-8_dp)) &
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "complex_region_roots: Newton polish did not converge")
    end subroutine newton_polish

    ! Append one distinct zero to the growing result arrays.
    subroutine append_root(roots, fvals, mult, nfound, z, fz, m)
        complex(dp), allocatable, intent(inout) :: roots(:), fvals(:)
        integer,     allocatable, intent(inout) :: mult(:)
        integer,                 intent(inout)  :: nfound
        complex(dp),             intent(in)     :: z, fz
        integer,                 intent(in)     :: m

        complex(dp), allocatable :: tr(:), tf(:)
        integer,     allocatable :: tm(:)

        allocate(tr(nfound + 1), tf(nfound + 1), tm(nfound + 1))
        if (nfound > 0) then
            tr(1:nfound) = roots
            tf(1:nfound) = fvals
            tm(1:nfound) = mult
        end if
        tr(nfound + 1) = z
        tf(nfound + 1) = fz
        tm(nfound + 1) = m
        call move_alloc(tr, roots)
        call move_alloc(tf, fvals)
        call move_alloc(tm, mult)
        nfound = nfound + 1
    end subroutine append_root

    ! Generalized eigenvalues of A v = z B v via LAPACK ZGGEV (alpha/beta).
    subroutine solve_pencil(a, b, n, ev, status)
        complex(dp),         intent(in)     :: a(n, n), b(n, n)
        integer,             intent(in)     :: n
        complex(dp),         intent(out)    :: ev(:)
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: ac(n, n), bc(n, n)
        complex(dp) :: alpha(n), beta(n), vl(1, 1), vr(1, 1)
        complex(dp) :: work(8*n)
        real(dp)    :: rwork(8*n)
        integer     :: info, k

        call status_set(status, FORTNUM_OK, "")
        ac = a
        bc = b
        call zggev('N', 'N', n, ac, n, bc, n, alpha, beta, vl, 1, vr, 1, &
            work, 8*n, rwork, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "complex_region_roots: ZGGEV failed on the moment pencil")
            return
        end if
        do k = 1, n
            if (beta(k) == (0.0_dp, 0.0_dp)) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "complex_region_roots: infinite generalized eigenvalue")
                return
            end if
            ev(k) = alpha(k) / beta(k)
        end do
    end subroutine solve_pencil

    ! Solve the dense complex system A x = b by Gaussian elimination with
    ! partial pivoting (small n; the confluent Vandermonde of distinct zeros).
    subroutine solve_lin(a, b, n, x, status)
        complex(dp),         intent(in)     :: a(n, n), b(n)
        integer,             intent(in)     :: n
        complex(dp),         intent(out)    :: x(n)
        type(fortnum_status_t), intent(out) :: status

        complex(dp) :: m(n, n), rhs(n), factor, tmp(n), tval
        integer     :: i, j, k, piv
        real(dp)    :: amax

        call status_set(status, FORTNUM_OK, "")
        m = a
        rhs = b
        do k = 1, n
            piv = k
            amax = abs(m(k, k))
            do i = k + 1, n
                if (abs(m(i, k)) > amax) then
                    amax = abs(m(i, k))
                    piv = i
                end if
            end do
            if (amax == 0.0_dp) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "complex_region_roots: singular Vandermonde in multiplicities")
                return
            end if
            if (piv /= k) then
                tmp = m(k, :); m(k, :) = m(piv, :); m(piv, :) = tmp
                tval = rhs(k); rhs(k) = rhs(piv); rhs(piv) = tval
            end if
            do i = k + 1, n
                factor = m(i, k) / m(k, k)
                do j = k, n
                    m(i, j) = m(i, j) - factor * m(k, j)
                end do
                rhs(i) = rhs(i) - factor * rhs(k)
            end do
        end do
        do i = n, 1, -1
            x(i) = rhs(i)
            do j = i + 1, n
                x(i) = x(i) - m(i, j) * x(j)
            end do
            x(i) = x(i) / m(i, i)
        end do
    end subroutine solve_lin

    ! Determinant of a small complex matrix by LU with partial pivoting.
    function det_lu(a, n) result(d)
        complex(dp), intent(in) :: a(n, n)
        integer,     intent(in) :: n
        complex(dp) :: d, m(n, n), factor, tmp(n)
        integer     :: i, j, k, piv
        real(dp)    :: amax

        m = a
        d = (1.0_dp, 0.0_dp)
        do k = 1, n
            piv = k
            amax = abs(m(k, k))
            do i = k + 1, n
                if (abs(m(i, k)) > amax) then
                    amax = abs(m(i, k))
                    piv = i
                end if
            end do
            if (amax == 0.0_dp) then
                d = (0.0_dp, 0.0_dp)
                return
            end if
            if (piv /= k) then
                tmp = m(k, :); m(k, :) = m(piv, :); m(piv, :) = tmp
                d = -d
            end if
            do i = k + 1, n
                factor = m(i, k) / m(k, k)
                do j = k, n
                    m(i, j) = m(i, j) - factor * m(k, j)
                end do
            end do
            d = d * m(k, k)
        end do
    end function det_lu

end module fortnum_roots_complex
