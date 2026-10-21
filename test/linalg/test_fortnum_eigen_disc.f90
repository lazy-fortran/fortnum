program test_fortnum_eigen_disc
    use fortnum_kinds, only: dp
    use fortnum_eigen
    implicit none
    real(dp), allocatable :: mass_diagonal(:), diagonal(:), factor(:), work(:)
    real(dp), allocatable :: start(:), vectors(:, :)
    real(dp) :: values(3), residuals(3), roots(3), lo, hi, mid, error, h, tol
    integer :: n, i, j, info
    character(32) :: argument
    type(eigen_options_t) :: options
    type(eigen_stats_t) :: stats
    n = 64
    if (command_argument_count() > 0) then
        call get_command_argument(1, argument)
        read (argument, *) n
    end if
    if (n < 32) error stop 'disc needs n>=32'
   allocate (mass_diagonal(n), diagonal(n), factor(n), work(n), start(n), vectors(n, 3))
    h = 1.0_dp/n
    ! Radial finite-volume weak form integral r u' v' dr = lambda integral r u v dr.
    ! Cell centers (i-1/2)h; zero flux at r=0, u(1)=0 by a half-cell face.
    mass_diagonal = [(real(i, dp) - 0.5_dp, i=1, n)]*h*h
    diagonal = [(2.0_dp*i - 1, i=1, n)]
    diagonal(n) = 3.0_dp*n - 1
    factor(1) = diagonal(1)
    do i = 2, n
        factor(i) = diagonal(i) - real(i - 1, dp)**2/factor(i - 1)
    end do
    start = 1
    tol = 1.e-10_dp
    options = eigen_options_t(max_basis=20, retain=10, max_restarts=100, &
                              which=EIGEN_NEAREST)
    call lanczos(stiffness, start, values, vectors, residuals, tol, &
                 options, info, stats, mass=mass, solve_shift=solve, sigma=0.0_dp)
    if (info /= EIGEN_OK) then
        print *, 'disc failure', info, stats
        error stop 1
    end if
    ! Independent symmetric row-sum bounds: ||K||_2 <= 4n, ||M||_2 <= 1/n.
    do j = 1, 3
        call stiffness(vectors(:, j), work)
        error = norm2(work - values(j)*mass_diagonal*vectors(:, j))
        if (error > tol*(4*n + abs(values(j))/n)*norm2(vectors(:, j))) &
            error stop 'disc backward error exceeds analytic norm bound'
    end do
    ! Independent roots of the exact entire J0 series, bracketed per radial mode.
    do j = 1, 3
        lo = real(3*j - 1, dp)
        hi = real(3*j, dp)
        do i = 1, 55
            mid = (lo + hi)/2
            if (series(lo)*series(mid) > 0) then
                lo = mid
            else
                hi = mid
            end if
        end do
        roots(j) = (lo + hi)/2
    end do
    error = maxval(abs(values - roots*roots)/(roots*roots))
    print *, 'disc n,cycles,applications,max-relative-error:', &
        n, stats%restarts, stats%applications, error
    print *, 'values:', values, ' residuals:', residuals
    if (error > 0.01_dp) error stop 'disc analytic Bessel limit'
contains
    function series(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y, term
        integer :: k
        y = 1
        term = 1
        do k = 1, 100
            term = -term*(x/(2*k))**2
            y = y + term
            if (abs(term) < epsilon(y)*max(1.0_dp, abs(y))) exit
        end do
    end function
    subroutine stiffness(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer :: k
        y = diagonal*x
        do k = 1, n - 1
            y(k) = y(k) - k*x(k + 1)
            y(k + 1) = y(k + 1) - k*x(k)
        end do
    end subroutine
    subroutine mass(x, y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        y = mass_diagonal*x
    end subroutine
    subroutine solve(x, y, code)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: y(:)
        integer, intent(out) :: code
        integer :: k
        work(1) = x(1)
        do k = 2, n
            work(k) = x(k) + real(k - 1, dp)*work(k - 1)/factor(k - 1)
        end do
        y(n) = work(n)/factor(n)
        do k = n - 1, 1, -1
            y(k) = (work(k) + real(k, dp)*y(k + 1))/factor(k)
        end do
        code = 0
    end subroutine
end program
