program test_fortnum_nuts
    ! Behavioral tests for fortnum_nuts on targets with known moments:
    ! a correlated 5-D Gaussian with an identity metric (step adaptation must
    ! find a usable step), an ill-conditioned Gaussian with the exact metric
    ! (short trajectories), and a Gaussian truncated by a callback that fails
    ! outside the support (failures become divergences, never draws).
    ! Tolerances are statistical: several Monte Carlo standard errors at the
    ! fixed seeds, so a correct sampler passes with a wide margin.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, status_set
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_nuts, only: nuts_options_t, nuts_result_t, nuts_sample
    implicit none

    type :: gauss_t
        real(dp), allocatable :: mean(:), prec(:, :)
        real(dp) :: upper = huge(1.0_dp)     !! fail when x(1) > upper
        integer :: calls = 0
    end type gauss_t

    integer :: nfail

    nfail = 0
    call test_correlated_identity_metric(nfail)
    call test_exact_metric(nfail)
    call test_failing_region(nfail)
    call test_bad_input(nfail)
    if (nfail > 0) then
        write(error_unit, '(I0,A)') nfail, ' NUTS check(s) failed'
        error stop 1
    end if
    print '(A)', 'fortnum_nuts: all checks passed'

contains

    subroutine logp(n, x, f, g, context, status)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n)
        real(dp), intent(out) :: f, g(n)
        class(*), intent(inout) :: context
        type(fortnum_status_t), intent(out) :: status
        select type (c => context)
        type is (gauss_t)
            c%calls = c%calls + 1
            if (x(1) > c%upper) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'outside support')
                f = 0.0_dp
                g = 0.0_dp
                return
            end if
            g = -matmul(c%prec, x - c%mean)
            f = 0.5_dp * dot_product(x - c%mean, g)
        end select
    end subroutine logp

    subroutine covariance_5(cov)
        real(dp), intent(out) :: cov(5, 5)
        real(dp) :: sd(5)
        integer :: i, j
        sd = [1.0_dp, 2.0_dp, 0.5_dp, 3.0_dp, 1.0_dp]
        do j = 1, 5
            do i = 1, 5
                cov(i, j) = sd(i) * sd(j) * 0.6_dp**abs(i - j)
            end do
        end do
    end subroutine covariance_5

    subroutine invert_spd(a, inv)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: inv(:, :)
        real(dp) :: l(size(a, 1), size(a, 1)), e(size(a, 1)), y(size(a, 1))
        integer :: n, j, i
        n = size(a, 1)
        call cholesky(a, l)
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            do i = 1, n
                y(i) = (e(i) - dot_product(l(i, 1:i - 1), y(1:i - 1))) / l(i, i)
            end do
            do i = n, 1, -1
                inv(i, j) = (y(i) - dot_product(l(i + 1:n, i), inv(i + 1:n, j))) / l(i, i)
            end do
        end do
    end subroutine invert_spd

    subroutine cholesky(a, l)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: l(:, :)
        integer :: i, j, n
        n = size(a, 1)
        l = 0.0_dp
        do j = 1, n
            l(j, j) = sqrt(a(j, j) - dot_product(l(j, 1:j - 1), l(j, 1:j - 1)))
            do i = j + 1, n
                l(i, j) = (a(i, j) - dot_product(l(i, 1:j - 1), l(j, 1:j - 1))) / l(j, j)
            end do
        end do
    end subroutine cholesky

    subroutine moments(x, mean, cov)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: mean(size(x, 1)), cov(size(x, 1), size(x, 1))
        integer :: k, m
        m = size(x, 2)
        mean = sum(x, dim=2) / m
        cov = 0.0_dp
        do k = 1, m
            cov = cov + spread(x(:, k) - mean, 2, size(x, 1)) * spread(x(:, k) - mean, 1, size(x, 1))
        end do
        cov = cov / (m - 1)
    end subroutine moments

    subroutine check(ok, what, nfail)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        integer, intent(inout) :: nfail
        if (.not. ok) then
            nfail = nfail + 1
            write(error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine check

    subroutine test_correlated_identity_metric(nfail)
        integer, intent(inout) :: nfail
        type(gauss_t) :: c
        type(nuts_options_t) :: opts
        type(nuts_result_t) :: res
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp) :: cov(5, 5), eye(5, 5), mean(5), est(5, 5), x0(5)
        integer :: i

        call covariance_5(cov)
        c%mean = [1.0_dp, -2.0_dp, 0.0_dp, 5.0_dp, 0.5_dp]
        allocate(c%prec(5, 5))
        call invert_spd(cov, c%prec)
        eye = 0.0_dp
        do i = 1, 5
            eye(i, i) = 1.0_dp
        end do
        x0 = 0.0_dp
        opts%n_warmup = 1000
        opts%n_samples = 4000
        call rng_seed(rng, 20261026_int64, st)
        call nuts_sample(logp, 5, x0, eye, opts, rng, c, res, st)
        call check(st%code == FORTNUM_OK, 'identity metric: status', nfail)
        call moments(res%x, mean, est)
        do i = 1, 5
            ! effective sample size is at least a few hundred here; 0.2 sd is > 4 MCSE
            call check(abs(mean(i) - c%mean(i)) < 0.2_dp * sqrt(cov(i, i)), &
                'identity metric: mean', nfail)
            call check(abs(est(i, i) / cov(i, i) - 1.0_dp) < 0.25_dp, &
                'identity metric: variance', nfail)
        end do
        call check(abs(est(1, 2) / sqrt(est(1, 1) * est(2, 2)) - 0.6_dp) < 0.1_dp, &
            'identity metric: correlation', nfail)
        call check(abs(sum(res%accept_stat) / size(res%accept_stat) - opts%target_accept) < 0.1_dp, &
            'identity metric: adapted acceptance', nfail)
        call check(res%n_divergent == 0, 'identity metric: no divergences', nfail)
    end subroutine test_correlated_identity_metric

    subroutine test_exact_metric(nfail)
        ! condition number 1e6: hopeless with the identity metric at this
        ! budget, trivial with the exact Cholesky factor of the covariance
        integer, intent(inout) :: nfail
        type(gauss_t) :: c
        type(nuts_options_t) :: opts
        type(nuts_result_t) :: res
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp) :: cov(3, 3), l(3, 3), mean(3), est(3, 3)
        integer :: i

        cov = reshape([1.0e-6_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.99_dp * 1.0e1_dp, &
            0.0_dp, 0.99_dp * 1.0e1_dp, 1.0e2_dp], [3, 3])
        c%mean = [0.0_dp, 1.0_dp, 2.0_dp]
        allocate(c%prec(3, 3))
        call invert_spd(cov, c%prec)
        call cholesky(cov, l)
        opts%n_warmup = 300
        opts%n_samples = 2000
        call rng_seed(rng, 7_int64, st)
        call nuts_sample(logp, 3, c%mean + 1.0_dp * [1.0e-3_dp, 1.0_dp, 10.0_dp], l, opts, rng, c, res, st)
        call check(st%code == FORTNUM_OK, 'exact metric: status', nfail)
        call moments(res%x, mean, est)
        do i = 1, 3
            call check(abs(mean(i) - c%mean(i)) < 0.15_dp * sqrt(cov(i, i)), 'exact metric: mean', nfail)
            call check(abs(est(i, i) / cov(i, i) - 1.0_dp) < 0.15_dp, 'exact metric: variance', nfail)
        end do
        call check(maxval(res%depth) <= 4, 'exact metric: short trajectories', nfail)
        call check(res%step_size > 0.3_dp, 'exact metric: step size of order one', nfail)
    end subroutine test_exact_metric

    subroutine test_failing_region(nfail)
        ! standard normal in 2-D, model fails for x(1) > 1: draws stay in the
        ! support and follow the truncated normal (mean of x1 = -phi(1)/Phi(1))
        integer, intent(inout) :: nfail
        type(gauss_t) :: c
        type(nuts_options_t) :: opts
        type(nuts_result_t) :: res
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp) :: eye(2, 2), truth

        c%mean = [0.0_dp, 0.0_dp]
        c%prec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        c%upper = 1.0_dp
        eye = c%prec
        opts%n_warmup = 500
        opts%n_samples = 4000
        call rng_seed(rng, 11_int64, st)
        call nuts_sample(logp, 2, [-0.5_dp, 0.0_dp], eye, opts, rng, c, res, st)
        call check(st%code == FORTNUM_OK, 'failing region: status', nfail)
        call check(all(res%x(1, :) <= 1.0_dp), 'failing region: draws inside support', nfail)
        call check(res%n_divergent > 0, 'failing region: failures reported as divergences', nfail)
        truth = -exp(-0.5_dp) / sqrt(2.0_dp * acos(-1.0_dp)) / (0.5_dp * erfc(-1.0_dp / sqrt(2.0_dp)))
        call check(abs(sum(res%x(1, :)) / size(res%x, 2) - truth) < 0.1_dp, &
            'failing region: truncated mean', nfail)
    end subroutine test_failing_region

    subroutine test_bad_input(nfail)
        integer, intent(inout) :: nfail
        type(gauss_t) :: c
        type(nuts_options_t) :: opts
        type(nuts_result_t) :: res
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp) :: l(2, 2)

        c%mean = [0.0_dp, 0.0_dp]
        c%prec = reshape([1.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [2, 2])
        l = 0.0_dp
        call rng_seed(rng, 1_int64, st)
        call nuts_sample(logp, 2, [0.0_dp, 0.0_dp], l, opts, rng, c, res, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, 'bad input: singular metric rejected', nfail)
        c%upper = -1.0_dp
        l = c%prec
        call nuts_sample(logp, 2, [0.0_dp, 0.0_dp], l, opts, rng, c, res, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, 'bad input: failing start rejected', nfail)
    end subroutine test_bad_input

end program test_fortnum_nuts
