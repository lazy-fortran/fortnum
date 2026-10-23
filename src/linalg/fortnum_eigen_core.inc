    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: krylov_schur, lanczos, real_krylov_schur, real_lanczos
    public :: eigen_bounds, inertia_count, real_operator, complex_operator
    public :: real_solve, complex_solve
    abstract interface
        subroutine complex_operator(x, y)
            import rk
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
        end subroutine
        subroutine real_operator(x, y)
            import rk
            real(rk), intent(in) :: x(:)
            real(rk), intent(out) :: y(:)
        end subroutine
        subroutine complex_solve(x, y, info)
            import rk
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
            integer, intent(out) :: info
        end subroutine
        subroutine real_solve(x, y, info)
            import rk
            real(rk), intent(in) :: x(:)
            real(rk), intent(out) :: y(:)
            integer, intent(out) :: info
        end subroutine
        function spectral_map(theta) result(value)
            import rk
            complex(rk), intent(in) :: theta
            complex(rk) :: value
        end function
        subroutine inertia_provider(shift, negative, zero, positive, info)
            import rk
            real(rk), intent(in) :: shift
            integer, intent(out) :: negative, zero, positive, info
        end subroutine
    end interface
    interface
        subroutine GEMM(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            import rk
            character, intent(in) :: ta, tb
            integer, intent(in) :: m, n, k, lda, ldb, ldc
            complex(rk), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
            complex(rk), intent(inout) :: c(ldc, *)
        end subroutine
        subroutine GEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            import rk
            character, intent(in) :: trans
            integer, intent(in) :: m, n, lda, incx, incy
            complex(rk), intent(in) :: alpha, beta, a(lda, *), x(*)
            complex(rk), intent(inout) :: y(*)
        end subroutine
        subroutine GEES(job, sort, select, n, a, lda, sdim, w, vs, ldvs, &
                work, lwork, rwork, bwork, info)
            import rk
            character, intent(in) :: job, sort
            logical, external :: select
            integer, intent(in) :: n, lda, ldvs, lwork
            integer, intent(out) :: sdim, info
            complex(rk), intent(inout) :: a(lda, *)
            complex(rk), intent(out) :: w(*), vs(ldvs, *), work(*)
            real(rk), intent(out) :: rwork(*)
            logical, intent(out) :: bwork(*)
        end subroutine
        subroutine TREXC(compq, n, t, ldt, q, ldq, ifst, ilst, info)
            import rk
            character, intent(in) :: compq
            integer, intent(in) :: n, ldt, ldq, ifst, ilst
            integer, intent(out) :: info
            complex(rk), intent(inout) :: t(ldt, *), q(ldq, *)
        end subroutine
        subroutine TREVC(side, how, select, n, t, ldt, vl, ldvl, vr, ldvr, &
                mm, m, work, rwork, info)
            import rk
            character, intent(in) :: side, how
            logical, intent(inout) :: select(*)
            integer, intent(in) :: n, ldt, ldvl, ldvr, mm
            integer, intent(out) :: m, info
            complex(rk), intent(in) :: t(ldt, *)
            complex(rk), intent(inout) :: vl(ldvl, *), vr(ldvr, *)
            complex(rk), intent(out) :: work(*)
            real(rk), intent(out) :: rwork(*)
        end subroutine
        subroutine HEEV(job, uplo, n, a, lda, w, work, lwork, rwork, info)
            import rk
            character, intent(in) :: job, uplo
            integer, intent(in) :: n, lda, lwork
            integer, intent(out) :: info
            complex(rk), intent(inout) :: a(lda, *)
            real(rk), intent(out) :: w(*), rwork(*)
            complex(rk), intent(out) :: work(*)
        end subroutine
    end interface
    type, abstract, public :: real_eigen_operator_t
    contains
        procedure(real_eigen_operator_t_binding), deferred :: apply
    end type
    type, extends(real_eigen_operator_t) :: real_operator_adapter_t
        procedure(real_operator), pointer, nopass :: callback => null()
    contains
        procedure :: apply => real_operator_adapter_t_call
    end type
    abstract interface
        subroutine real_eigen_operator_t_binding(self, x, y)
            import :: rk, real_eigen_operator_t
            class(real_eigen_operator_t), intent(inout) :: self
            real(rk), intent(in) :: x(:)
            real(rk), intent(out) :: y(:)
        end subroutine
    end interface
    type, abstract, public :: real_eigen_solve_t
    contains
        procedure(real_eigen_solve_t_binding), deferred :: solve
    end type
    type, extends(real_eigen_solve_t) :: real_solve_adapter_t
        procedure(real_solve), pointer, nopass :: callback => null()
    contains
        procedure :: solve => real_solve_adapter_t_call
    end type
    abstract interface
        subroutine real_eigen_solve_t_binding(self, x, y, info)
            import :: rk, real_eigen_solve_t
            class(real_eigen_solve_t), intent(inout) :: self
            real(rk), intent(in) :: x(:)
            real(rk), intent(out) :: y(:)
            integer, intent(out) :: info
        end subroutine
    end interface
    type, abstract, public :: complex_eigen_operator_t
    contains
        procedure(complex_eigen_operator_t_binding), deferred :: apply
    end type
    type, extends(complex_eigen_operator_t) :: complex_operator_adapter_t
        procedure(complex_operator), pointer, nopass :: callback => null()
    contains
        procedure :: apply => complex_operator_adapter_t_call
    end type
    abstract interface
        subroutine complex_eigen_operator_t_binding(self, x, y)
            import :: rk, complex_eigen_operator_t
            class(complex_eigen_operator_t), intent(inout) :: self
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
        end subroutine
    end interface
    type, abstract, public :: complex_eigen_solve_t
    contains
        procedure(complex_eigen_solve_t_binding), deferred :: solve
    end type
    type, extends(complex_eigen_solve_t) :: complex_solve_adapter_t
        procedure(complex_solve), pointer, nopass :: callback => null()
    contains
        procedure :: solve => complex_solve_adapter_t_call
    end type
    abstract interface
        subroutine complex_eigen_solve_t_binding(self, x, y, info)
            import :: rk, complex_eigen_solve_t
            class(complex_eigen_solve_t), intent(inout) :: self
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
            integer, intent(out) :: info
        end subroutine
    end interface
    public :: krylov_schur_object
    public :: lanczos_object
    public :: real_krylov_schur_object
    public :: real_lanczos_object
    type, extends(complex_eigen_operator_t) :: real_complex_operator_t
        class(real_eigen_operator_t), pointer :: operator => null()
        real(rk), allocatable :: xr(:), xi(:), yr(:), yi(:)
    contains
        procedure :: apply => real_complex_operator_t_call
    end type
    type, extends(complex_eigen_solve_t) :: real_complex_solve_t
        class(real_eigen_solve_t), pointer :: operator => null()
        real(rk), allocatable :: xr(:), xi(:), yr(:), yi(:)
    contains
        procedure :: solve => real_complex_solve_t_call
    end type

contains
    subroutine krylov_schur_object(apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma, transform, map)
        class(complex_eigen_operator_t), target, intent(inout) :: apply
        complex(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(complex_eigen_operator_t), target, intent(inout), optional :: mass
        class(complex_eigen_solve_t), target, intent(inout), optional :: solve_mass, solve_shift, transform
        complex(rk), intent(in), optional :: sigma
        procedure(spectral_map), optional :: map
        call iterate(.false., apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma, transform, map)
    end subroutine

    subroutine lanczos_object(apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma)
        class(complex_eigen_operator_t), target, intent(inout) :: apply
        complex(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: vectors(:, :)
        real(rk), intent(out) :: values(:), residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(complex_eigen_operator_t), target, intent(inout), optional :: mass
        class(complex_eigen_solve_t), target, intent(inout), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        complex(rk) :: cv(size(values))
        if (present(sigma)) then
            call iterate(.true., apply, start, cv, vectors, residuals, tolerance, &
          options, info, stats, mass, solve_mass, solve_shift, cmplx(sigma, 0.0_rk, rk))
        else
            call iterate(.true., apply, start, cv, vectors, residuals, tolerance, &
                options, info, stats, mass, solve_mass, solve_shift)
        end if
        values = real(cv, rk)
    end subroutine

    subroutine iterate(hermitian, apply, start, values, vectors, residuals, &
  tolerance, options, info, stats, mass, solve_mass, solve_shift, sigma, transform, map)
        logical, intent(in) :: hermitian
        class(complex_eigen_operator_t), target, intent(inout) :: apply
        complex(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(complex_eigen_operator_t), target, intent(inout), optional :: mass
        class(complex_eigen_solve_t), target, intent(inout), optional :: solve_mass, solve_shift, transform
        complex(rk), intent(in), optional :: sigma
        procedure(spectral_map), optional :: map
        class(complex_eigen_operator_t), pointer :: mass_call
        class(complex_eigen_solve_t), pointer :: mass_solve_call, shift_call, transform_call
        procedure(spectral_map), pointer :: map_call
        complex(rk), allocatable :: q(:, :), mq(:, :), bq(:, :), mbq(:, :)
        complex(rk), allocatable :: nq(:, :), nmq(:, :), nbq(:, :), nmbq(:, :)
        complex(rk), allocatable :: h(:, :), u(:, :), ev(:, :), work(:), theta(:)
      complex(rk), allocatable :: w(:), mw(:), tmp(:), kx(:), mx(:), saved(:), msaved(:)
        complex(rk), parameter :: one=cmplx(1.0_rk,0.0_rk,rk), zero=cmplx(0.0_rk, &
            0.0_rk,rk)
        complex(rk) :: coeff, lambda, dummy(1, 1), shift
        real(rk), allocatable :: rw(:), hw(:)
        logical, allocatable :: bw(:)
        real(rk) :: beta, before, score, best, denom, defect
        real(rk) :: kscale, mscale, xnorm
        integer :: n, m, nev, k, keep, cycle, i, j, p, ierr, count, seed
        logical :: found
        info = EIGEN_INVALID
        stats = eigen_stats_t()
        values = 0
        vectors = 0
        residuals = huge(1.0_rk)
        n = size(start)
        nev = size(values)
        m = min(n, options%max_basis)
        if (n < 1 .or. nev < 1 .or. nev > m) return
        if (m < n .and. m <= nev) return
        if (size(vectors, 1) /= n .or. size(vectors, 2) /= nev) return
        if (size(residuals) /= nev) return
        if (options%max_restarts < 0 .or. options%retain < 0) return
        if (options%which < 1 .or. options%which > 4) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0 .or. .not. finite(start)) return
        if (.not. ieee_is_finite(options%operator_scale)) return
        if (.not. ieee_is_finite(options%mass_scale)) return
        if (options%operator_scale < 0 .or. options%mass_scale < 0) return
        if (options%operator_scale > real(huge(1.0_rk), kind(options%operator_scale))) &
            return
        if (options%mass_scale > real(huge(1.0_rk), kind(options%mass_scale))) return
        kscale = real(options%operator_scale, rk)
        mscale = real(options%mass_scale, rk)
        if (.not. present(mass) .and. options%mass_scale == 0) mscale = 1
        if (present(solve_shift) .neqv. present(sigma)) return
        if (present(transform) .neqv. present(map)) return
        if (present(transform) .and. present(solve_shift)) return
        if (present(mass)) then
            if (.not. present(solve_shift) .and. .not. present(transform)) then
                if (.not. present(solve_mass)) return
            end if
        else
            if (present(solve_mass)) return
        end if
        shift = 0
        if (present(sigma)) shift = sigma
        if (.not. finite([shift])) return
        keep = max(nev, min(m - 1, max(nev + 2, m/2)))
        if (options%retain > 0) keep = options%retain
        if (m < n) then
            if (keep < nev .or. keep >= m) return
        end if
        allocate (q(n, m), mq(n, m), bq(n, m), mbq(n, m), &
            nq(n, m), nmq(n, m), nbq(n, m), nmbq(n, m), &
            h(m, m), u(m, m), ev(m, m), work(max(4*m, 1)), theta(m), &
            w(n), mw(n), tmp(n), kx(n), mx(n), saved(n), msaved(n), &
            rw(max(3*m, 1)), hw(m), bw(m))
        nullify (mass_call, mass_solve_call, shift_call, transform_call, map_call)
        if (present(mass)) mass_call => mass
        if (present(solve_mass)) mass_solve_call => solve_mass
        if (present(solve_shift)) shift_call => solve_shift
        if (present(transform)) transform_call => transform
        if (present(map)) map_call => map
        k = 0
        seed = 0
        w = start
        call metric(w, mw)
        if (info < 0) return
        call orthogonalize(w, mw, k, beta)
        if (info < 0) return
        if (beta <= tiny(1.0_rk)) then
            info = EIGEN_INVALID
            return
        end if
        do cycle = 0, options%max_restarts
            stats%restarts = cycle
            do while (k < m)
                if (beta == 0) then
                    ! Happy breakdown: deterministic injections expose missing
                    ! invariant subspaces and repeated eigenvalues.
                    found = .false.
                    do p = 1, n
                        seed = mod(seed, n) + 1
                        w = 0
                        w(seed) = 1
                        call metric(w, mw)
                        if (info < 0) return
                        call orthogonalize(w, mw, k, beta)
                        if (info < 0) return
                        if (beta > sqrt(epsilon(1.0_rk))) then
                            found = .true.
                            exit
                        end if
                    end do
                    if (.not. found) then
                        info = EIGEN_METRIC_FAILED
                        return
                    end if
                end if
                k = k + 1
                q(:, k) = w/beta
                mq(:, k) = mw/beta
                if (options%mass_scale == 0) &
                    mscale = max(mscale, cnorm(mq(:, k))/cnorm(q(:, k)))
                call operator(q(:, k), bq(:, k))
                if (info < 0) return
                call metric(bq(:, k), mbq(:, k))
                if (info < 0) return
                w = bq(:, k)
                mw = mbq(:, k)
                before = cnorm(w)
                call orthogonalize(w, mw, k, beta)
                if (info < 0) return
                if (beta <= 64*epsilon(1.0_rk)*max(before, tiny(1.0_rk))) beta = 0
            end do
            call GEMM('C', 'N', m, m, n, one, q, n, mbq, n, zero, h, m)
            if (hermitian) then
                defect = maxval(abs(h - conjg(transpose(h))))
                if (defect > 100*sqrt(epsilon(1.0_rk))* &
                        max(maxval(abs(h)), tiny(1.0_rk))) then
                    info = EIGEN_METRIC_FAILED
                    return
                end if
                u = (h + conjg(transpose(h)))/2
                call HEEV('V', 'U', m, u, m, hw, work, size(work), rw, ierr)
                h = 0
                do i = 1, m
                    h(i, i) = hw(i)
                end do
            else
                call GEES('V', 'N', unused_select, m, h, m, count, theta, u, m, &
                    work, size(work), rw, bw, ierr)
            end if
            if (ierr /= 0) then
                info = EIGEN_DENSE_FAILED
                return
            end if
            do i = 1, m
                p = i
                best = huge(1.0_rk)
                do j = i, m
                    lambda = physical(h(j, j))
                    score = ordering(lambda)
                    if (score < best) then
                        best = score
                        p = j
                    end if
                end do
                if (p == i) cycle
                if (hermitian) then
                    coeff = h(i, i)
                    h(i, i) = h(p, p)
                    h(p, p) = coeff
                    work(1:m) = u(:, i)
                    u(:, i) = u(:, p)
                    u(:, p) = work(1:m)
                else
                    call TREXC('V', m, h, m, u, m, p, i, ierr)
                    if (ierr /= 0) then
                        info = EIGEN_DENSE_FAILED
                        return
                    end if
                end if
            end do
            ev = u
            if (.not. hermitian) then
            call TREVC('R', 'B', bw, m, h, m, dummy, 1, ev, m, m, count, work, rw, ierr)
                if (ierr /= 0) then
                    info = EIGEN_DENSE_FAILED
                    return
                end if
            end if
            stats%converged = 0
            do j = 1, nev
                call GEMV('N', n, m, one, q, n, ev(:, j), 1, zero, vectors(:, j), 1)
                call metric(vectors(:, j), mx)
                if (info < 0) return
                beta = sqrt(max(0.0_rk, real(dot_product(vectors(:, j), mx), rk)))
                if (beta <= tiny(1.0_rk)) then
                    info = EIGEN_METRIC_FAILED
                    return
                end if
                vectors(:, j) = vectors(:, j)/beta
                mx = mx/beta
                lambda = physical(h(j, j))
                values(j) = lambda
                call apply%apply(vectors(:, j), kx)
                stats%applications = stats%applications + 1
                if (.not. finite(kx)) then
                    info = EIGEN_OPERATOR_FAILED
                    return
                end if
                if (.not. finite([lambda])) then
                    info = EIGEN_OPERATOR_FAILED
                    return
                end if
                w = kx - lambda*mx
                if (.not. finite(w)) then
                    info = EIGEN_OPERATOR_FAILED
                    return
                end if
                residuals(j) = cnorm(w)
                xnorm = cnorm(vectors(:, j))
                if (options%operator_scale == 0) &
                    kscale = max(kscale, cnorm(kx)/xnorm)
                if (options%mass_scale == 0) mscale = max(mscale, cnorm(mx)/xnorm)
                denom = (kscale + abs(lambda)*mscale)*xnorm
                if (.not. ieee_is_finite(denom)) then
                    info = EIGEN_OPERATOR_FAILED
                    return
                end if
                if (residuals(j) <= tolerance*denom) &
                    stats%converged = stats%converged + 1
            end do
            if (stats%converged == nev) then
                stats%locked = nev
                info = EIGEN_OK
                return
            end if
            if (cycle == options%max_restarts .or. m == n) exit
            ! Soft locking retains Schur directions and their operator images.
            call GEMM('N', 'N', n, keep, m, one, q, n, u, m, zero, nq, n)
            call GEMM('N', 'N', n, keep, m, one, mq, n, u, m, zero, nmq, n)
            call GEMM('N', 'N', n, keep, m, one, bq, n, u, m, zero, nbq, n)
            call GEMM('N', 'N', n, keep, m, one, mbq, n, u, m, zero, nmbq, n)
            stats%locked = 0
            do j = 1, min(nev, keep)
                w = nbq(:, j)
                call GEMV('N', n, keep, -one, nq, n, h(:, j), 1, one, w, 1)
                if (cnorm(w) > tolerance*max(1.0_rk, cnorm(nbq(:, j)))) exit
                stats%locked = stats%locked + 1
            end do
            q(:, 1:keep) = nq(:, 1:keep)
            mq(:, 1:keep) = nmq(:, 1:keep)
            bq(:, 1:keep) = nbq(:, 1:keep)
            mbq(:, 1:keep) = nmbq(:, 1:keep)
            k = keep
            ! Residual columns are collinear in exact Krylov-Schur arithmetic.
            ! Choose the largest to avoid expanding from a locked column.
            best = -1
            before = 0
            do j = 1, keep
                w = bq(:, j)
                mw = mbq(:, j)
                before = max(before, cnorm(w))
                call orthogonalize(w, mw, k, beta)
                if (info < 0) return
                if (beta > best) then
                    best = beta
                    saved = w
                    msaved = mw
                end if
            end do
            w = saved
            mw = msaved
            beta = best
            if (beta <= 64*epsilon(1.0_rk)*max(before, tiny(1.0_rk))) beta = 0
        end do
        info = EIGEN_NOT_CONVERGED
    contains
        subroutine metric(x, y)
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
            info = EIGEN_NOT_CONVERGED
            if (associated(mass_call)) then
                call mass_call%apply(x, y)
            else
                y = x
            end if
            if (.not. finite(y)) info = EIGEN_METRIC_FAILED
        end subroutine
        subroutine operator(x, y)
            complex(rk), intent(in) :: x(:)
            complex(rk), intent(out) :: y(:)
            integer :: code
            code = 0
            ! Shifted/transformed images do not measure the original K norm.
            ! Sample K on the available basis vector only in automatic mode.
            if (options%operator_scale == 0) then
                if (associated(transform_call) .or. associated(shift_call)) then
                    call apply%apply(x, kx)
                    stats%applications = stats%applications + 1
                    if (.not. finite(kx)) then
                        info = EIGEN_OPERATOR_FAILED
                        return
                    end if
                    kscale = max(kscale, cnorm(kx)/cnorm(x))
                end if
            end if
            if (associated(transform_call)) then
                call transform_call%solve(x, y, code)
            else if (associated(shift_call)) then
                call metric(x, tmp)
                if (info < 0) return
                call shift_call%solve(tmp, y, code)
            else
                call apply%apply(x, tmp)
                if (.not. finite(tmp)) then
                    info = EIGEN_OPERATOR_FAILED
                    return
                end if
                if (options%operator_scale == 0) &
                    kscale = max(kscale, cnorm(tmp)/cnorm(x))
                if (associated(mass_call)) then
                    call mass_solve_call%solve(tmp, y, code)
                else
                    y = tmp
                end if
            end if
            stats%applications = stats%applications + 1
            if (code /= 0) then
                info = EIGEN_OPERATOR_FAILED
                return
            end if
            if (.not. finite(y)) info = EIGEN_OPERATOR_FAILED
        end subroutine
        subroutine orthogonalize(x, mx, k, beta)
            complex(rk), intent(inout) :: x(:), mx(:)
            integer, intent(in) :: k
            real(rk), intent(out) :: beta
            complex(rk) :: c, energy
            integer :: pass, j, l
            do pass = 1, 2
                call metric(x, mx)
                if (info < 0) return
                ! Keep one team for the dependent projection/update sequence.
                ! Every thread advances j in order; worksharing barriers make
                ! the reduced coefficient and updated vectors visible to all.
                !$omp parallel if(n >= 8192) private(j, l)
                do j = 1, k
                    !$omp single
                    c = 0
                    !$omp end single
                    !$omp do simd reduction(+:c)
                    do l = 1, n
                        c = c + conjg(q(l, j))*mx(l)
                    end do
                    !$omp end do simd
                    !$omp do simd
                    do l = 1, n
                        x(l) = x(l) - c*q(l, j)
                        mx(l) = mx(l) - c*mq(l, j)
                    end do
                    !$omp end do simd
                end do
                !$omp end parallel
            end do
            ! Recompute the metric after cancellation; recursively subtracting
            ! Mq images can lose positivity near an invariant subspace.
            call metric(x, mx)
            if (info < 0) return
            energy = dot_product(x, mx)
            beta = 0
            if (.not. finite([energy])) then
                info = EIGEN_METRIC_FAILED
                return
            end if
            if (real(energy, rk) < 0 .or. &
                    abs(aimag(energy)) > 100*sqrt(epsilon(1.0_rk))* &
                    max(abs(energy), tiny(1.0_rk))) then
                info = EIGEN_METRIC_FAILED
                return
            end if
            beta = sqrt(real(energy, rk))
        end subroutine
        function physical(t) result(v)
            complex(rk), intent(in) :: t
            complex(rk) :: v
            v = t
            if (associated(map_call)) then
                v = map_call(t)
            else if (associated(shift_call)) then
                if (abs(t) <= tiny(1.0_rk)) then
                    v = huge(1.0_rk)
                else
                    v = shift + 1/t
                end if
            end if
        end function
        function ordering(v) result(s)
            complex(rk), intent(in) :: v
            real(rk) :: s
            s = huge(1.0_rk)
            if (.not. finite([v])) return
            select case (options%which)
            case (EIGEN_SMALLEST); s = real(v, rk)
            case (EIGEN_LARGEST); s = -real(v, rk)
            case (EIGEN_NEAREST); s = abs(v - shift)
            case (EIGEN_LARGEST_MAGNITUDE); s = -abs(v)
            end select
        end function
    end subroutine
    pure logical function finite(x)
        complex(rk), intent(in) :: x(:)
        finite = all(ieee_is_finite(real(x, rk))) .and. all(ieee_is_finite(aimag(x)))
    end function
    real(rk) function cnorm(x) result(v)
        complex(rk), intent(in) :: x(:)
        real(rk) :: scale
        integer :: i
        scale = 0
        !$omp parallel do simd if(size(x) >= 8192) reduction(max:scale)
        do i = 1, size(x)
            scale = max(scale, abs(x(i)))
        end do
        v = 0
        if (scale > 0) then
            !$omp parallel do simd if(size(x) >= 8192) reduction(+:v)
            do i = 1, size(x)
                v = v + abs(x(i)/scale)**2
            end do
            v = scale*sqrt(v)
        end if
    end function
    logical function unused_select(z)
        complex(rk), intent(in) :: z
        unused_select = .false.
    end function
    subroutine real_krylov_schur_object(apply, start, values, vectors, residuals, &
            tolerance, options, info, stats, mass, solve_mass, solve_shift, sigma)
        class(real_eigen_operator_t), target, intent(inout) :: apply
        real(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(real_eigen_operator_t), target, intent(inout), optional :: mass
        class(real_eigen_solve_t), target, intent(inout), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        call real_bridge(.false., apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma)
    end subroutine
    subroutine real_lanczos_object(apply, start, values, vectors, residuals, &
            tolerance, options, info, stats, mass, solve_mass, solve_shift, sigma)
        class(real_eigen_operator_t), target, intent(inout) :: apply
        real(rk), intent(in) :: start(:)
        real(rk), intent(out) :: values(:), vectors(:, :), residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(real_eigen_operator_t), target, intent(inout), optional :: mass
        class(real_eigen_solve_t), target, intent(inout), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        complex(rk) :: cv(size(values)), cx(size(vectors, 1), size(vectors, 2))
        call real_bridge(.true., apply, start, cv, cx, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma)
        values = real(cv, rk)
        vectors = real(cx, rk)
    end subroutine
    subroutine real_bridge(hermitian, apply, start, values, vectors, residuals, &
   tolerance, options, info, stats, real_mass, real_solve_mass, real_solve_shift, sigma)
        logical, intent(in) :: hermitian
        class(real_eigen_operator_t), target, intent(inout) :: apply
        real(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        class(real_eigen_operator_t), target, intent(inout), optional :: real_mass
        class(real_eigen_solve_t), target, intent(inout), optional :: real_solve_mass, real_solve_shift
        real(rk), intent(in), optional :: sigma
        type(real_complex_operator_t), target :: op, mop
        type(real_complex_solve_t), target :: msolve, shifted
        class(complex_eigen_operator_t), pointer :: pm
        class(complex_eigen_solve_t), pointer :: ps, pi
        nullify(pm, ps, pi)
        op%operator => apply
        allocate(op%xr(size(start)), op%xi(size(start)), &
                 op%yr(size(start)), op%yi(size(start)))
        if (present(real_mass)) then
            mop%operator => real_mass
            allocate(mop%xr(size(start)), mop%xi(size(start)), &
                     mop%yr(size(start)), mop%yi(size(start)))
            pm => mop
        end if
        if (present(real_solve_mass)) then
            msolve%operator => real_solve_mass
            allocate(msolve%xr(size(start)), msolve%xi(size(start)), &
                     msolve%yr(size(start)), msolve%yi(size(start)))
            ps => msolve
        end if
        if (present(real_solve_shift)) then
            shifted%operator => real_solve_shift
            allocate(shifted%xr(size(start)), shifted%xi(size(start)), &
                     shifted%yr(size(start)), shifted%yi(size(start)))
            pi => shifted
        end if
        if (present(sigma)) then
            call iterate(hermitian, op, cmplx(start, 0.0_rk, rk), values, vectors, &
                residuals, tolerance, options, info, stats, pm, ps, pi, &
                cmplx(sigma, 0.0_rk, rk))
        else
            call iterate(hermitian, op, cmplx(start, 0.0_rk, rk), values, vectors, &
                residuals, tolerance, options, info, stats, pm, ps, pi)
        end if
    end subroutine

    subroutine eigen_bounds(residual, mass_lower, metric_norm, radius, info, &
            eigenvector_condition, rayleigh, gap_lower, gap_upper, lower, upper)
        !! Floating-point estimates, not outward-rounded certificates.
        real(rk), intent(in) :: residual, mass_lower, metric_norm
        real(rk), intent(out) :: radius
        integer, intent(out) :: info
        real(rk), intent(in), optional :: eigenvector_condition, rayleigh
        real(rk), intent(in), optional :: gap_lower, gap_upper
        real(rk), intent(out), optional :: lower, upper
        real(rk) :: rho, lo, hi
        info = EIGEN_INVALID
        radius = huge(1.0_rk)
        if (present(lower)) lower = -huge(1.0_rk)
        if (present(upper)) upper = huge(1.0_rk)
        if (.not. all(ieee_is_finite([residual, mass_lower, metric_norm]))) return
        if (residual < 0 .or. mass_lower <= 0 .or. metric_norm <= 0) return
        if (.not. ieee_is_finite(sqrt(mass_lower)*metric_norm)) return
        if (sqrt(mass_lower)*metric_norm <= tiny(1.0_rk)) return
        call bound_radius(residual, mass_lower, metric_norm, rho)
        radius = rho
        if (present(eigenvector_condition)) then
            if (.not. ieee_is_finite(eigenvector_condition)) return
            if (eigenvector_condition < 1) return
            radius = rho*eigenvector_condition
        end if
        if (present(rayleigh) .or. present(gap_lower) .or. present(gap_upper)) then
            if (.not. present(rayleigh)) return
            if (.not. present(gap_lower)) return
            if (.not. present(gap_upper)) return
            if (present(eigenvector_condition)) return
            if (.not. all(ieee_is_finite([rayleigh, gap_lower, gap_upper]))) return
            if (rayleigh <= gap_lower .or. rayleigh >= gap_upper) return
            if (rho >= sqrt(rayleigh - gap_lower)*sqrt(gap_upper - rayleigh)) return
            call bound_temple(rho, rayleigh, gap_lower, gap_upper, lo, hi)
            if (.not. all(ieee_is_finite([lo, hi]))) return
            if (present(lower)) lower = lo
            if (present(upper)) upper = hi
        end if
        if (.not. ieee_is_finite(radius)) return
        info = EIGEN_OK
    end subroutine

    subroutine inertia_count(factorization, n, left, right, count, info)
        !! Negative inertia of K-s M counts eigenvalues strictly below s.
        !! The difference counts [left,right), including endpoint multiplicity.
        procedure(inertia_provider) :: factorization
        integer, intent(in) :: n
        real(rk), intent(in) :: left, right
        integer, intent(out) :: count, info
        integer :: nl, nr, z, p, code
        info = EIGEN_INVALID
        count = 0
        if (n < 1 .or. left >= right) return
        if (.not. all(ieee_is_finite([left, right]))) return
        call factorization(left, nl, z, p, code)
        if (code /= 0) then
            info = EIGEN_OPERATOR_FAILED
            return
        end if
        if (min(nl, z, p) < 0 .or. nl + z + p /= n) return
        call factorization(right, nr, z, p, code)
        if (code /= 0) then
            info = EIGEN_OPERATOR_FAILED
            return
        end if
        if (min(nr, z, p) < 0 .or. nr + z + p /= n .or. nr < nl) return
        count = nr - nl
        info = EIGEN_OK
    end subroutine

    subroutine real_operator_adapter_t_call(self, x, y)
        class(real_operator_adapter_t), intent(inout) :: self
        real(rk), intent(in) :: x(:)
        real(rk), intent(out) :: y(:)
        call self%callback(x, y)
    end subroutine

    subroutine real_solve_adapter_t_call(self, x, y, info)
        class(real_solve_adapter_t), intent(inout) :: self
        real(rk), intent(in) :: x(:)
        real(rk), intent(out) :: y(:)
        integer, intent(out) :: info
        call self%callback(x, y, info)
    end subroutine

    subroutine complex_operator_adapter_t_call(self, x, y)
        class(complex_operator_adapter_t), intent(inout) :: self
        complex(rk), intent(in) :: x(:)
        complex(rk), intent(out) :: y(:)
        call self%callback(x, y)
    end subroutine

    subroutine complex_solve_adapter_t_call(self, x, y, info)
        class(complex_solve_adapter_t), intent(inout) :: self
        complex(rk), intent(in) :: x(:)
        complex(rk), intent(out) :: y(:)
        integer, intent(out) :: info
        call self%callback(x, y, info)
    end subroutine

    subroutine krylov_schur(apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma, transform, map)
        procedure(complex_operator) :: apply
        complex(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        procedure(complex_operator), optional :: mass
        procedure(complex_solve), optional :: solve_mass, solve_shift, transform
        complex(rk), intent(in), optional :: sigma
        procedure(spectral_map), optional :: map
        type(complex_operator_adapter_t), target :: apply_adapter
        type(complex_operator_adapter_t), target :: mass_adapter
        class(complex_eigen_operator_t), pointer :: mass_object
        type(complex_solve_adapter_t), target :: solve_mass_adapter
        class(complex_eigen_solve_t), pointer :: solve_mass_object
        type(complex_solve_adapter_t), target :: solve_shift_adapter
        class(complex_eigen_solve_t), pointer :: solve_shift_object
        type(complex_solve_adapter_t), target :: transform_adapter
        class(complex_eigen_solve_t), pointer :: transform_object
        nullify(mass_object, solve_mass_object, solve_shift_object, transform_object)

        apply_adapter%callback => apply
        if (present(mass)) then
            mass_adapter%callback => mass
            mass_object => mass_adapter
        end if
        if (present(solve_mass)) then
            solve_mass_adapter%callback => solve_mass
            solve_mass_object => solve_mass_adapter
        end if
        if (present(solve_shift)) then
            solve_shift_adapter%callback => solve_shift
            solve_shift_object => solve_shift_adapter
        end if
        if (present(transform)) then
            transform_adapter%callback => transform
            transform_object => transform_adapter
        end if
        call krylov_schur_object(apply_adapter, &
            start, &
            values, &
            vectors, &
            residuals, &
            tolerance, &
            options, &
            info, &
            stats, &
            mass_object, &
            solve_mass_object, &
            solve_shift_object, &
            sigma, &
            transform_object, &
            map)
    end subroutine

    subroutine lanczos(apply, start, values, vectors, residuals, tolerance, &
            options, info, stats, mass, solve_mass, solve_shift, sigma)
        procedure(complex_operator) :: apply
        complex(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: vectors(:, :)
        real(rk), intent(out) :: values(:), residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        procedure(complex_operator), optional :: mass
        procedure(complex_solve), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        type(complex_operator_adapter_t), target :: apply_adapter
        type(complex_operator_adapter_t), target :: mass_adapter
        class(complex_eigen_operator_t), pointer :: mass_object
        type(complex_solve_adapter_t), target :: solve_mass_adapter
        class(complex_eigen_solve_t), pointer :: solve_mass_object
        type(complex_solve_adapter_t), target :: solve_shift_adapter
        class(complex_eigen_solve_t), pointer :: solve_shift_object
        nullify(mass_object, solve_mass_object, solve_shift_object)

        apply_adapter%callback => apply
        if (present(mass)) then
            mass_adapter%callback => mass
            mass_object => mass_adapter
        end if
        if (present(solve_mass)) then
            solve_mass_adapter%callback => solve_mass
            solve_mass_object => solve_mass_adapter
        end if
        if (present(solve_shift)) then
            solve_shift_adapter%callback => solve_shift
            solve_shift_object => solve_shift_adapter
        end if
        call lanczos_object(apply_adapter, &
            start, &
            values, &
            vectors, &
            residuals, &
            tolerance, &
            options, &
            info, &
            stats, &
            mass_object, &
            solve_mass_object, &
            solve_shift_object, &
            sigma)
    end subroutine

    subroutine real_krylov_schur(apply, start, values, vectors, residuals, &
            tolerance, options, info, stats, mass, solve_mass, solve_shift, sigma)
        procedure(real_operator) :: apply
        real(rk), intent(in) :: start(:)
        complex(rk), intent(out) :: values(:), vectors(:, :)
        real(rk), intent(out) :: residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        procedure(real_operator), optional :: mass
        procedure(real_solve), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        type(real_operator_adapter_t), target :: apply_adapter
        type(real_operator_adapter_t), target :: mass_adapter
        class(real_eigen_operator_t), pointer :: mass_object
        type(real_solve_adapter_t), target :: solve_mass_adapter
        class(real_eigen_solve_t), pointer :: solve_mass_object
        type(real_solve_adapter_t), target :: solve_shift_adapter
        class(real_eigen_solve_t), pointer :: solve_shift_object
        nullify(mass_object, solve_mass_object, solve_shift_object)

        apply_adapter%callback => apply
        if (present(mass)) then
            mass_adapter%callback => mass
            mass_object => mass_adapter
        end if
        if (present(solve_mass)) then
            solve_mass_adapter%callback => solve_mass
            solve_mass_object => solve_mass_adapter
        end if
        if (present(solve_shift)) then
            solve_shift_adapter%callback => solve_shift
            solve_shift_object => solve_shift_adapter
        end if
        call real_krylov_schur_object(apply_adapter, &
            start, &
            values, &
            vectors, &
            residuals, &
            tolerance, &
            options, &
            info, &
            stats, &
            mass_object, &
            solve_mass_object, &
            solve_shift_object, &
            sigma)
    end subroutine

    subroutine real_lanczos(apply, start, values, vectors, residuals, &
            tolerance, options, info, stats, mass, solve_mass, solve_shift, sigma)
        procedure(real_operator) :: apply
        real(rk), intent(in) :: start(:)
        real(rk), intent(out) :: values(:), vectors(:, :), residuals(:)
        real(rk), intent(in) :: tolerance
        type(eigen_options_t), intent(in) :: options
        integer, intent(out) :: info
        type(eigen_stats_t), intent(out) :: stats
        procedure(real_operator), optional :: mass
        procedure(real_solve), optional :: solve_mass, solve_shift
        real(rk), intent(in), optional :: sigma
        type(real_operator_adapter_t), target :: apply_adapter
        type(real_operator_adapter_t), target :: mass_adapter
        class(real_eigen_operator_t), pointer :: mass_object
        type(real_solve_adapter_t), target :: solve_mass_adapter
        class(real_eigen_solve_t), pointer :: solve_mass_object
        type(real_solve_adapter_t), target :: solve_shift_adapter
        class(real_eigen_solve_t), pointer :: solve_shift_object
        nullify(mass_object, solve_mass_object, solve_shift_object)

        apply_adapter%callback => apply
        if (present(mass)) then
            mass_adapter%callback => mass
            mass_object => mass_adapter
        end if
        if (present(solve_mass)) then
            solve_mass_adapter%callback => solve_mass
            solve_mass_object => solve_mass_adapter
        end if
        if (present(solve_shift)) then
            solve_shift_adapter%callback => solve_shift
            solve_shift_object => solve_shift_adapter
        end if
        call real_lanczos_object(apply_adapter, &
            start, &
            values, &
            vectors, &
            residuals, &
            tolerance, &
            options, &
            info, &
            stats, &
            mass_object, &
            solve_mass_object, &
            solve_shift_object, &
            sigma)
    end subroutine

    subroutine real_complex_operator_t_call(self, x, y)
        class(real_complex_operator_t), intent(inout) :: self
        complex(rk), intent(in) :: x(:)
        complex(rk), intent(out) :: y(:)
        self%xr = real(x, rk)
        self%xi = aimag(x)
        call self%operator%apply(self%xr, self%yr)
        call self%operator%apply(self%xi, self%yi)
        y = cmplx(self%yr, self%yi, rk)
    end subroutine

    subroutine real_complex_solve_t_call(self, x, y, info)
        class(real_complex_solve_t), intent(inout) :: self
        complex(rk), intent(in) :: x(:)
        complex(rk), intent(out) :: y(:)
        integer, intent(out) :: info
        integer :: second
        self%xr = real(x, rk)
        self%xi = aimag(x)
        call self%operator%solve(self%xr, self%yr, info)
        if (info /= 0) return
        call self%operator%solve(self%xi, self%yi, second)
        info = second
        if (info /= 0) return
        y = cmplx(self%yr, self%yi, rk)
    end subroutine
