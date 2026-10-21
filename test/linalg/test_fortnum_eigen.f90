program test_fortnum_eigen
    use fortnum_kinds, only: dp, sp
    use fortnum_eigen
    use fortnum_eigen_temple_products_dp_kernel, only: fortnum_eigen_temple_products_dp
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: n = 24, nev = 3
    real(dp) :: start(n), val(nev), vec(n, nev), res(nev), s(n), exact(nev), z(n)
    real(dp) :: radius, lower, upper, theta, rho, shift, gram(nev, nev), scale
    complex(dp) :: cstart(n), cv(nev), cx(n, nev), cz(n)
    real(dp) :: tl,tu,ar,at,al,ab
    real(sp) :: ss(8), sv(2), sx(8, 2), sr(2), srad
    complex(sp) :: cs(8), csv(2), csx(8, 2)
    type(eigen_options_t) :: opt
    type(eigen_stats_t) :: stats
    integer :: i, j, info, count
    logical :: fail_solve
    start = [(sin(real(i, dp)*1.234_dp), i=1, n)]
    s = [(1 + real(i, dp)/n, i=1, n)]
    exact = [(2 - 2*cos(real(i, dp)*acos(-1.0_dp)/(n + 1)), i=1, nev)]
    opt = eigen_options_t(max_basis=12, max_restarts=200, retain=6)
    call lanczos(toeplitz, start, val, vec, res, 1.e-9_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'Toeplitz restarted convergence')
    call check(stats%restarts > 0, 'restart exercised')
    call check(maxval(abs(val - exact)) < 1.e-7_dp, 'Toeplitz sine spectrum')
    do j = 1, nev
        call toeplitz(vec(:, j), z)
        call check(abs(norm2(z - val(j)*vec(:, j)) - res(j)) < 1.e-13_dp, &
                   'independent original residual')
    end do
    call lanczos(pencil, start, val, vec, res, 1.e-9_dp, opt, info, stats, &
                 mass=mass, solve_mass=mass_solve)
    call check(info == EIGEN_OK, 'generalized SPD convergence')
    call check(maxval(abs(val - exact)) < 1.e-7_dp, 'congruent pencil spectrum')
    do j = 1, nev
        do i = 1, nev
            gram(i, j) = sum(vec(:, i)*s*s*vec(:, j))
        end do
        gram(j, j) = gram(j, j) - 1
    end do
    call check(maxval(abs(gram)) < 1.e-10_dp, 'mass orthogonality')
    ! Manufactured diagonal generalized pencil K_ii=(i+1)i^2, M_ii=i+1.
    shift = 30
    fail_solve = .false.
    opt%which = EIGEN_NEAREST
    opt%max_basis = 10
    call lanczos(diag_k, start, val, vec, res, 1.e-10_dp, opt, info, stats, &
                 mass=diag_mass, solve_shift=shift_solve, sigma=shift)
    call check(info == EIGEN_OK, 'shift-invert original-pencil convergence')
    call check(maxval(abs(val - [25.0_dp, 36.0_dp, 16.0_dp])) < 1.e-7_dp, &
               'shift-invert physical ordering')
    fail_solve = .true.
    call lanczos(diag_k, start, val, vec, res, 1.e-10_dp, opt, info, stats, &
                 mass=diag_mass, solve_shift=shift_solve, sigma=shift)
    call check(info == EIGEN_OPERATOR_FAILED, 'supplied solve failure')
    fail_solve = .false.
    ! Upper triangular, nonnormal resistive/rotating model: exact diagonal roots.
    opt%which = EIGEN_LARGEST
    opt%max_basis = 12
    cstart = cmplx(start, 0.3_dp, dp)
    call krylov_schur(nonnormal, cstart, cv, cx, res, 1.e-9_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'complex nonnormal restarted convergence')
    do j = 1, nev
        call check(abs(cv(j) - cmplx(n - j + 1, 0.1_dp*(n - j + 1), dp)) < 1.e-6_dp, &
                   'triangular complex exact spectrum')
        call nonnormal(cx(:, j), cz)
        call check(norm2(abs(cz - cv(j)*cx(:, j))) < 1.e-6_dp, &
                   'nonnormal true residual')
    end do
    call krylov_schur(real_nonnormal, start, cv, cx, res, 1.e-9_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'real Arnoldi callback')
    call check(abs(cv(1) - n) < 1.e-6_dp, 'real Arnoldi exact spectrum')
    ! Custom spectral transformation B=2K+I, with inverse spectral map.
    call krylov_schur(nonnormal, cstart, cv, cx, res, 1.e-9_dp, opt, info, stats, &
                      transform=transform, map=map)
    call check(info == EIGEN_OK, 'custom transform convergence')
    call check(abs(cv(1) - cmplx(n, 0.1_dp*n, dp)) < 1.e-6_dp, 'custom map')
    ! Happy breakdown and deflation: start sees one eigenspace only.
    opt%max_basis = n
    opt%which = EIGEN_SMALLEST
    start = 0
    start(n) = 1
    call lanczos(repeated, start, val, vec, res, 1.e-11_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'happy breakdown injection')
    call check(maxval(abs(val - [1.0_dp, 1.0_dp, 2.0_dp])) < 1.e-10_dp, &
               'repeated eigenvalues retained')
    ! A full-basis complex Hermitian path.
    call lanczos(complex_hermitian, cstart, val, cx, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'complex Hermitian convergence')
    call check(maxval(abs(val - exact)) < 1.e-9_dp, 'unitary Toeplitz spectrum')
    start = 0
    call lanczos(toeplitz, start, val, vec, res, 1.e-9_dp, opt, info, stats)
    call check(info == EIGEN_INVALID, 'zero start rejected')
    start = 1
    call lanczos(pencil, start, val, vec, res, 1.e-9_dp, opt, info, stats, mass=mass)
    call check(info == EIGEN_INVALID, 'missing mass solve rejected')
    opt%max_basis = 6
    opt%retain = 3
    opt%max_restarts = 0
    call lanczos(toeplitz, start, val, vec, res, 1.e-15_dp, opt, info, stats)
    call check(info == EIGEN_NOT_CONVERGED, 'nonconvergence explicit')
    ! A tiny operator must not make an unresolved Ritz pair pass an absolute test.
    scale = 1.e-12_dp
    call lanczos(scaled_toeplitz, start, val, vec, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_NOT_CONVERGED, 'small operator remains unconverged')
    call check(any(res > 1.e-10_dp*(4*scale + abs(val))), &
               'independent Toeplitz norm bound rejects unresolved Ritz pairs')
    opt%operator_scale = 4*scale
    opt%mass_scale = 1
    call lanczos(scaled_toeplitz, start, val, vec, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_NOT_CONVERGED, 'supplied norm rejects unresolved pairs')
    opt%max_basis = n
    call lanczos(scaled_toeplitz, start, val, vec, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_OK, 'small operator full basis converges')
    call check(maxval(abs(val/scale - exact)) < 1.e-10_dp, 'scaled sine spectrum')
    do j = 1, nev
        call scaled_toeplitz(vec(:, j), z)
        call check(norm2(z - val(j)*vec(:, j)) <= &
                   1.e-10_dp*(4*scale + abs(val(j)))*norm2(vec(:, j)), &
                   'independent supplied-scale backward error')
    end do
    opt%operator_scale = -1
    call lanczos(toeplitz, start, val, vec, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_INVALID, 'negative operator scale rejected')
    opt%operator_scale = 0
    opt%mass_scale = ieee_value(1.0_dp, ieee_quiet_nan)
    call lanczos(toeplitz, start, val, vec, res, 1.e-10_dp, opt, info, stats)
    call check(info == EIGEN_INVALID, 'nonfinite mass scale rejected')
    ! Exact 2x2 Hermitian Rayleigh quotient: eigenvalues 1 and 3.
    theta = 1.02_dp
    rho = sqrt(0.0396_dp)
    call eigen_bounds(rho, 1.0_dp, 1.0_dp, radius, info, &
           rayleigh=theta, gap_lower=0.0_dp, gap_upper=2.5_dp, lower=lower, upper=upper)
    call check(info == EIGEN_OK, 'Temple applicability')
    call check(lower <= 1 .and. upper >= 1, 'Temple contains exact eigenvalue')
    call check(abs(theta - 1) <= radius, 'Hermitian residual inclusion')
    call eigen_bounds(0.125_dp, 4.0_dp, 2.0_dp, radius, info, &
                      eigenvector_condition=8.0_dp)
    call check(info == EIGEN_OK .and. abs(radius - 0.25_dp) < 1.e-15_dp, &
               'generalized Bauer-Fike scaling')
    call eigen_bounds(rho, 1.0_dp, 1.0_dp, radius, info, &
                      rayleigh=theta, gap_lower=theta, gap_upper=3.0_dp)
    call check(info == EIGEN_INVALID, 'invalid Temple gap')
    call inertia_count(inertia, n, 4.0_dp, 25.0_dp, count, info)
    call check(info == EIGEN_OK .and. count == 3, 'inertia endpoint count [4,25)')
    ! Generated products checked by the independent analytic r derivative.
    call fortnum_eigen_temple_products_dp(0.25_dp,1.5_dp,0.0_dp,3.0_dp, &
        1.0_dp,0.0_dp,0.0_dp,0.0_dp,2.0_dp,3.0_dp, &
        lower,upper,tl,tu,ar,at,al,ab)
    call check(abs(tl+1.0_dp/3) < 1.e-14_dp,'Temple generated JVP lower')
    call check(abs(tu-1.0_dp/3) < 1.e-14_dp,'Temple generated JVP upper')
    call check(abs(2*tl+3*tu-ar) < 1.e-14_dp,'Temple generated adjoint identity')
    ! Both precision families and both scalar callback families compile and run.
    opt = eigen_options_t(max_basis=8)
    ss = 1
    call lanczos(single_diag, ss, sv, sx, sr, 1.e-5_sp, opt, info, stats)
    call check(info == EIGEN_OK, 'FP32 Lanczos')
    call check(maxval(abs(sv - [1.0_sp, 2.0_sp])) < 1.e-4_sp, 'FP32 exact values')
    cs = cmplx(ss, 0.0_sp, sp)
    call krylov_schur(single_complex, cs, csv, csx, sr, 1.e-5_sp, opt, info, stats)
    call check(info == EIGEN_OK, 'FP32 complex Arnoldi')
    call eigen_bounds(0.125_sp, 1.0_sp, 1.0_sp, srad, info)
    call check(info == EIGEN_OK .and. srad == 0.125_sp, 'FP32 bounds')
    call large_complex_projection()
    print *, 'matrix-free eigen: analytic spectra, pencils, restarts and failures PASS'
contains
    subroutine large_complex_projection()
        ! Three distinct complex diagonal roots span the full Krylov space.
        ! This size exercises the threaded conjugate projections in both kinds.
        integer, parameter :: large_n = 8192
        complex(dp), allocatable :: v64(:), x64(:, :), r64(:)
        complex(sp), allocatable :: v32(:), x32(:, :), r32(:)
        complex(dp) :: lambda64(2)
        complex(sp) :: lambda32(2)
        real(dp) :: residual64(2)
        real(sp) :: residual32(2)
        integer :: k, mode

        allocate (v64(large_n), x64(large_n, 2), r64(large_n), &
                  v32(large_n), x32(large_n, 2), r32(large_n))
        v64 = [(cmplx(1.0_dp, sin(real(k, dp)), dp), k=1, large_n)]
        v32 = cmplx(v64, kind=sp)
        opt = eigen_options_t(max_basis=3, retain=2, max_restarts=0)
        call krylov_schur(three_roots_dp, v64, lambda64, x64, residual64, &
                         1.e-10_dp, opt, info, stats)
        call check(info == EIGEN_OK, 'large complex FP64 projection')
        call krylov_schur(three_roots_sp, v32, lambda32, x32, residual32, &
                         1.e-4_sp, opt, info, stats)
        call check(info == EIGEN_OK, 'large complex FP32 projection')
        do mode = 1, 2
            call check(abs(lambda64(mode) - cmplx(mode, 0.2_dp*mode, dp)) &
                       < 1.e-9_dp, 'large FP64 exact diagonal spectrum')
            call check(abs(lambda32(mode) - cmplx(mode, 0.2_sp*mode, sp)) &
                       < 1.e-3_sp, 'large FP32 exact diagonal spectrum')
            call three_roots_dp(x64(:, mode), r64)
            call three_roots_sp(x32(:, mode), r32)
            call check(norm2(abs(r64 - lambda64(mode)*x64(:, mode))) &
                       < 1.e-9_dp, 'large FP64 independent residual')
            call check(norm2(abs(r32 - lambda32(mode)*x32(:, mode))) &
                       < 1.e-3_sp, 'large FP32 independent residual')
        end do
    end subroutine
    subroutine three_roots_dp(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer :: k, root
        do k = 1, size(x)
            root = 1 + mod(k - 1, 3)
            y(k) = cmplx(root, 0.2_dp*root, dp)*x(k)
        end do
    end subroutine
    subroutine three_roots_sp(x, y)
        complex(sp), intent(in) :: x(:)
        complex(sp), intent(out) :: y(:)
        integer :: k, root
        do k = 1, size(x)
            root = 1 + mod(k - 1, 3)
            y(k) = cmplx(root, 0.2_sp*root, sp)*x(k)
        end do
    end subroutine
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (ok) return
        print *, 'FAIL ', label, ' info=', info, ' restarts=', stats%restarts
        print *, 'values=', val, ' residuals=', res
        error stop 1
    end subroutine
    subroutine toeplitz(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = 2*x
        do k = 1, size(x) - 1
            y(k) = y(k) - x(k + 1)
            y(k + 1) = y(k + 1) - x(k)
        end do
    end subroutine
    subroutine scaled_toeplitz(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        call toeplitz(x, y)
        y = scale*y
    end subroutine
    subroutine pencil(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        call toeplitz(s*x, y)
        y = s*y
    end subroutine
    subroutine mass(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = s*s*x
    end subroutine
    subroutine mass_solve(x, y, code)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: code
        y = x/(s*s)
        code = 0
    end subroutine
    subroutine diag_k(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = [(real((k + 1)*k*k, dp)*x(k), k=1, size(x))]
    end subroutine
    subroutine diag_mass(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = [(real(k + 1, dp)*x(k), k=1, size(x))]
    end subroutine
    subroutine shift_solve(x, y, code)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: code
        integer :: k
        code = 1
        if (fail_solve) return
        y = [(x(k)/(real(k + 1, dp)*(real(k*k, dp) - shift)), k=1, size(x))]
        code = 0
    end subroutine
    subroutine nonnormal(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer :: k
        y = [(cmplx(k, 0.1_dp*k, dp)*x(k), k=1, size(x))]
        y(:size(x) - 1) = y(:size(x) - 1) + cmplx(0.2_dp, 0.15_dp, dp)*x(2:)
    end subroutine
    subroutine real_nonnormal(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = [(real(k, dp)*x(k), k=1, size(x))]
        y(:size(x) - 1) = y(:size(x) - 1) + 0.2_dp*x(2:)
    end subroutine
    subroutine transform(x, y, code)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer, intent(out) :: code
        call nonnormal(x, y)
        y = 2*y + x
        code = 0
    end subroutine
    function map(t) result(v)
        complex(dp), intent(in) :: t
        complex(dp) :: v
        v = (t - 1)/2
    end function
    subroutine repeated(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = [(real(max(1, k - 1), dp)*x(k), k=1, size(x))]
    end subroutine
    subroutine complex_hermitian(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        integer :: k
        y = 2*x
        do k = 1, size(x) - 1
            y(k) = y(k) + cmplx(0.0_dp, 1.0_dp, dp)*x(k + 1)
            y(k + 1) = y(k + 1) - cmplx(0.0_dp, 1.0_dp, dp)*x(k)
        end do
    end subroutine
    subroutine inertia(t, negative, zero, positive, code)
        real(dp), intent(in) :: t
        integer, intent(out) :: negative, zero, positive, code
        integer :: k
        negative = 0
        zero = 0
        do k = 1, n
            if (real(k*k, dp) < t) negative = negative + 1
            if (real(k*k, dp) == t) zero = zero + 1
        end do
        positive = n - negative - zero
        code = 0
    end subroutine
    subroutine single_diag(x, y)
        real(sp), intent(in) :: x(:)
        real(sp), intent(out) :: y(:)
        integer :: k
        y = [(real(k, sp)*x(k), k=1, size(x))]
    end subroutine
    subroutine single_complex(x, y)
        complex(sp), intent(in) :: x(:)
        complex(sp), intent(out) :: y(:)
        integer :: k
        y = [(cmplx(k, k, sp)*x(k), k=1, size(x))]
    end subroutine
end program
