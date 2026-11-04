program test_multifidelity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_multifidelity, only: autoregressive_gp_t
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(autoregressive_gp_t) :: gp
    type(fortnum_status_t) :: status
    integer, parameter :: n = 17
    real(dp) :: x(1, n), lo(1, n), hi(1, n), y(1), j(1, 1), v, &
        yp(1), ym(1), dummy(1, 1), dot(1), bar(1), error, train_error
    integer :: i
    call noisy_ar1()
    call dense_ar1()
    call single_site()
    call parallel_control()
    do i = 1, n
        x(1, i) = -2 + 4*real(i - 1, dp)/(n - 1)
    end do
    lo = sin(x); hi = 1.2_dp*lo + 0.1_dp*cos(2*x)
    call gp%fit(x, lo, hi, [.7_dp], .2_dp, 1.2_dp, 1.e-7_dp, status)
    call require(status%code == 0, 'fit')
    train_error = 0
    do i = 1, n
        call gp%predict(x(:, i), lo(:, i), reshape(cos(x(:, i)), [1, 1]), &
            y, j, v, status)
        call require(status%code == 0, 'predict')
        train_error = max(train_error, abs(y(1) - hi(1, i)))
    end do
    call require(train_error < 2.e-9_dp, 'training interpolation')
    call gp%predict([.37_dp], [sin(.37_dp)], reshape([cos(.37_dp)], [1, 1]), &
        y, j, v, status)
    error = abs(y(1) - (1.2_dp*sin(.37_dp) + .1_dp*cos(.74_dp)))
    print '(A,ES14.6)', 'held-out analytic error=', error
    print '(A,ES14.6)', 'training analytic error=', train_error
    call require(error < 2.e-5_dp, 'held-out reference')
    call require(v >= 0 .and. v < 1.e-5_dp, 'posterior variance')
    call gp%predict([.3701_dp], [sin(.3701_dp)], &
        reshape([cos(.3701_dp)], [1, 1]), yp, dummy, v, status)
    call gp%predict([.3699_dp], [sin(.3699_dp)], &
        reshape([cos(.3699_dp)], [1, 1]), ym, dummy, v, status)
    error = abs((yp(1) - ym(1))/2.e-4_dp - j(1, 1))
    print '(A,ES14.6)', 'finite-difference derivative error=', error
    call require(error < 1.e-7_dp, 'derivative finite differences')
    call gp%jvp([.37_dp], [sin(.37_dp)], reshape([cos(.37_dp)], [1, 1]), &
        [.3_dp], dot, status)
    call gp%vjp([.37_dp], [sin(.37_dp)], reshape([cos(.37_dp)], [1, 1]), &
        [.8_dp], bar, status)
    call require(abs(.8_dp*dot(1) - .3_dp*bar(1)) < 1.e-14_dp, 'transpose')
    call gp%fit(x, lo, 1.2_dp*lo, [.7_dp], .2_dp, 1.2_dp, 1.e-7_dp, status)
    call gp%predict([.37_dp], [sin(.37_dp)], reshape([cos(.37_dp)], [1, 1]), &
        y, j, v, status)
    call require(abs(y(1) - 1.2_dp*sin(.37_dp)) < 1.e-15_dp, 'zero discrepancy')
    call gp%fit(x, lo, hi, [-1._dp], .2_dp, 1.2_dp, 0._dp, status)
    call require(status%code /= 0, 'invalid length rejected')
    call gp%predict([.37_dp], [sin(.37_dp)], reshape([cos(.37_dp)], [1, 1]), &
        y, j, v, status)
    call require(status%code /= 0, 'failed fit invalidates state')
contains
    subroutine noisy_ar1()
        use fortnum_ar1, only: ar1_gp_t, ar1_covariance, rbf_covariance, covariance_factor
        use fortnum_gaussian_constraints, only: constraint_nullspace
        type(ar1_gp_t) :: model
        real(dp) :: low(1, 1), cross(1, 1), high(1, 1), delta(1, 1), nl(1, 1), nh(1, 1)
        real(dp) :: prediction(1), variance(1, 1), prior(1, 1), cl(1, 1), ch(1, 1)
        real(dp) :: density, grad(1), score(2), expected(2), det, fp, fm, h
        real(dp) :: derivative(2, 2, 1), points(3), k(3, 3), dk(3, 3, 2)
        real(dp) :: l(3, 3), dl(3, 3, 2), kp(3, 3), km(3, 3), scratch(3, 3, 2)
        real(dp) :: lp(3, 3), lm(3, 3), dummy_deriv(3, 3, 2)
        real(dp), allocatable :: covariance(:, :), dr(:, :), q(:, :)
        integer :: j, code
        low = 1.3_dp; cross = .7_dp; high = 1.1_dp; delta = .2_dp
        nl = .4_dp; nh = .1_dp; h = 1.e-5_dp
        call model%condition(low, cross, high, delta, nl, nh, 1.2_dp, [.4_dp], [1._dp], status)
        call require(status%code == 0, 'noisy nonnested AR1 condition')
        call ar1_covariance(low, cross, high, delta, nl, nh, 1.2_dp, covariance, dr, status)
        det = covariance(1, 1)*covariance(2, 2) - covariance(1, 2)**2
        expected = [(covariance(2, 2)*.4_dp - covariance(1, 2))/det, &
            (covariance(1, 1) - covariance(1, 2)*.4_dp)/det]
        call model%log_likelihood(density, status)
        call require(abs(density + (dot_product([.4_dp, 1._dp], expected) &
            + log(det) + 2*log(2*acos(-1._dp)))/2) < 1.e-14_dp, 'AR1 density determinant')
        derivative(:, :, 1) = dr
        call model%gradient(derivative, grad, status, score)
        call require(status%code == 0, 'AR1 gradient status')
        call require(maxval(abs(score - expected)) < 1.e-14_dp, 'AR1 mean score')
        prior = 2; cl = .5_dp; ch = .8_dp
        call model%predict(prior, cl, ch, prediction, variance, status)
        call require(status%code == 0, 'AR1 prediction status')
        call require(abs(prediction(1) - dot_product([.5_dp, .8_dp], expected)) &
            < 1.e-14_dp, 'AR1 conditional mean')
        call require(abs(variance(1, 1) - (2 - (.25_dp*covariance(2, 2) &
            - .8_dp*covariance(1, 2) + .64_dp*covariance(1, 1))/det)) &
            < 1.e-14_dp, 'AR1 conditional variance retains low noise')
        call model%condition(low, cross, high, delta, nl, nh, 1.2_dp + h, [.4_dp], [1._dp], status)
        call model%log_likelihood(fp, status)
        call model%condition(low, cross, high, delta, nl, nh, 1.2_dp - h, [.4_dp], [1._dp], status)
        call model%log_likelihood(fm, status)
        call require(abs((fp - fm)/(2*h) - grad(1)) < 1.e-9_dp, 'AR1 rho gradient FD')
        points = [0._dp, .3_dp, 1._dp]
        call rbf_covariance(points, points, 1.1_dp, .4_dp, k, dk, status)
        call require(status%code == 0, 'RBF covariance')
        call covariance_factor(k, dk, l, dl, status)
        call require(status%code == 0, 'covariance factor')
        call require(maxval(abs(matmul(l, transpose(l)) - k)) < 1.e-14_dp, 'factor residual')
        do j = 1, 2
            call rbf_covariance(points, points, 1.1_dp*exp(h*real(2 - j, dp)), &
                .4_dp*exp(h*real(j - 1, dp)), kp, scratch, status)
            call rbf_covariance(points, points, 1.1_dp*exp(-h*real(2 - j, dp)), &
                .4_dp*exp(-h*real(j - 1, dp)), km, scratch, status)
            call covariance_factor(kp, dk, lp, dummy_deriv, status)
            call covariance_factor(km, dk, lm, dummy_deriv, status)
            call require(maxval(abs((kp - km)/(2*h) - dk(:, :, j))) < 1.e-9_dp, 'kernel FD')
            call require(maxval(abs((lp - lm)/(2*h) - dl(:, :, j))) < 1.e-9_dp, 'Cholesky FD')
        end do
        call constraint_nullspace(reshape([1._dp, 2._dp, 3._dp, &
            2._dp, 4._dp, 6._dp], [3, 2]), q, code)
        call require(code == 0, 'constraint nullspace')
        call require(size(q, 2) == 2, 'dependent protected columns')
        call require(maxval(abs(matmul([1._dp, 2._dp, 3._dp], q))) &
            < 2.e-15_dp, 'constraint orthogonality')
        call model%condition(low, cross, high, delta, nl, nh, 0._dp, [.4_dp], [1._dp], status)
        call model%predict(prior, 0*cl, ch, prediction, variance, status)
        call require(abs(prediction(1) - .8_dp/.3_dp) < 1.e-14_dp, 'rho zero independent levels')
        nl = -10
        call model%condition(low, cross, high, delta, nl, nh, 1._dp, [.4_dp], [1._dp], status)
        call require(status%code /= 0, 'indefinite training covariance')
        call model%predict(prior, cl, ch, prediction, variance, status)
        call require(status%code /= 0, 'failed condition invalidates posterior')
        call rbf_covariance(points, points, 1._dp, 0._dp, k, dk, status)
        call require(status%code /= 0, 'RBF rejects zero length')
        call rbf_covariance(points, points, huge(1._dp), .4_dp, k, dk, status)
        call require(status%code /= 0, 'RBF rejects overflowing covariance')
    end subroutine

    subroutine dense_ar1()
        use fortnum_ar1, only: ar1_gp_t, ar1_covariance, rbf_covariance
        type(ar1_gp_t) :: model
        real(dp) :: xl(2), xh(3), xt(2), rho, amp, ell, h, fp, fm, density, pert_amp, pert_ell
        real(dp) :: low(2, 2), cross(2, 3), high(3, 3), delta(3, 3)
        real(dp) :: nl(2, 2), nh(3, 3), dl(2, 2, 2), dc(2, 3, 2), dh(3, 3, 2)
        real(dp) :: scratch(3, 3, 2), k(5, 5), chol(5, 5), rhs(5, 3)
        real(dp) :: target(2, 2), cl(2, 2), ch(3, 2), prediction(2), cov(2, 2)
        real(dp) :: full_cross(5, 2), expected(2, 2), y(5), deriv(5, 5, 3), grad(3)
        real(dp), allocatable :: assembled(:, :), dr(:, :)
        integer :: i, j, p, info
        external :: dposv
        xl = [.1_dp, .8_dp]; xh = [0._dp, .4_dp, 1._dp]; xt = [.2_dp, .6_dp]
        y = [.3_dp, -.5_dp, .7_dp, -.1_dp, .2_dp]
        amp = 1.3_dp; ell = .35_dp; rho = .9_dp; h = 1.e-5_dp
        nl = reshape([.04_dp, .01_dp, .01_dp, .09_dp], [2, 2])
        nh = 0
        do i = 1, 3
            nh(i, i) = .02_dp*i
        end do
        call rbf_covariance(xl, xl, amp, ell, low, dl, status)
        call rbf_covariance(xl, xh, amp, ell, cross, dc, status)
        call rbf_covariance(xh, xh, amp, ell, high, dh, status)
        call rbf_covariance(xh, xh, .4_dp, .2_dp, delta, scratch, status)
        call model%condition(low, cross, high, delta, nl, nh, rho, y(:2), y(3:), status)
        call require(status%code == 0, 'unequal nonnested AR1 design')
        ! Independent dense Gaussian conditioning with directly evaluated kernels.
        do j = 1, 5
            do i = 1, 5
                if (i <= 2 .and. j <= 2) then
                    k(i, j) = amp**2*exp(-.5_dp*((xl(i) - xl(j))/ell)**2) + nl(i, j)
                else if (i <= 2) then
                    k(i, j) = rho*amp**2*exp(-.5_dp*((xl(i) - xh(j - 2))/ell)**2)
                else if (j <= 2) then
                    k(i, j) = rho*amp**2*exp(-.5_dp*((xh(i - 2) - xl(j))/ell)**2)
                else
                    k(i, j) = rho**2*amp**2*exp(-.5_dp*((xh(i - 2) - xh(j - 2))/ell)**2) &
                        + .4_dp**2*exp(-.5_dp*((xh(i - 2) - xh(j - 2))/.2_dp)**2) + nh(i - 2, j - 2)
                end if
            end do
        end do
        do j = 1, 2
            cl(:, j) = rho*amp**2*exp(-.5_dp*((xl - xt(j))/ell)**2)
            ch(:, j) = rho**2*amp**2*exp(-.5_dp*((xh - xt(j))/ell)**2) &
                + .4_dp**2*exp(-.5_dp*((xh - xt(j))/.2_dp)**2)
            target(:, j) = rho**2*amp**2*exp(-.5_dp*((xt - xt(j))/ell)**2) &
                + .4_dp**2*exp(-.5_dp*((xt - xt(j))/.2_dp)**2)
        end do
        full_cross(:2, :) = cl; full_cross(3:, :) = ch
        chol = k; rhs(:, 1) = y; rhs(:, 2:) = full_cross
        call dposv('U', 5, 3, chol, 5, rhs, 5, info)
        call require(info == 0, 'dense conditioning oracle')
        call model%predict(target, cl, ch, prediction, cov, status)
        call require(status%code == 0, 'multivariate AR1 prediction')
        expected = target - matmul(transpose(full_cross), rhs(:, 2:))
        call require(maxval(abs(cov - expected)) < 2.e-14_dp, 'dense full covariance oracle')
        call require(maxval(abs(prediction - matmul(transpose(full_cross), rhs(:, 1)))) &
            < 2.e-14_dp, 'dense mean oracle')
        call model%log_likelihood(density, status)
        call require(abs(density + .5_dp*(dot_product(y, rhs(:, 1)) &
            + 2*sum([(log(chol(i, i)), i=1, 5)]) + 5*log(2*acos(-1._dp)))) &
            < 2.e-14_dp, 'dense normalized density oracle')
        call ar1_covariance(low, cross, high, delta, nl, nh, rho, assembled, dr, status)
        deriv(:, :, 1) = dr
        do p = 1, 2
            deriv(:2, :2, p + 1) = dl(:, :, p)
            deriv(:2, 3:, p + 1) = rho*dc(:, :, p)
            deriv(3:, :2, p + 1) = transpose(rho*dc(:, :, p))
            deriv(3:, 3:, p + 1) = rho**2*dh(:, :, p)
        end do
        call model%gradient(deriv, grad, status)
        call require(status%code == 0, 'multiple covariance gradients')
        do p = 1, 3
            pert_amp = amp; pert_ell = ell
            if (p == 2) pert_amp = amp*exp(h)
            if (p == 3) pert_ell = ell*exp(h)
            call rbf_covariance(xl, xl, pert_amp, pert_ell, low, dl, status)
            call rbf_covariance(xl, xh, pert_amp, pert_ell, cross, dc, status)
            call rbf_covariance(xh, xh, pert_amp, pert_ell, high, dh, status)
            call rbf_covariance(xh, xh, .4_dp, .2_dp, delta, scratch, status)
            call model%condition(low, cross, high, delta, nl, nh, &
                rho + h*real(merge(1, 0, p == 1), dp), y(:2), y(3:), status)
            call model%log_likelihood(fp, status)
            pert_amp = amp; pert_ell = ell
            if (p == 2) pert_amp = amp*exp(-h)
            if (p == 3) pert_ell = ell*exp(-h)
            call rbf_covariance(xl, xl, pert_amp, pert_ell, low, dl, status)
            call rbf_covariance(xl, xh, pert_amp, pert_ell, cross, dc, status)
            call rbf_covariance(xh, xh, pert_amp, pert_ell, high, dh, status)
            call rbf_covariance(xh, xh, .4_dp, .2_dp, delta, scratch, status)
            call model%condition(low, cross, high, delta, nl, nh, &
                rho - h*real(merge(1, 0, p == 1), dp), y(:2), y(3:), status)
            call model%log_likelihood(fm, status)
            call require(abs((fp - fm)/(2*h) - grad(p)) < 2.e-8_dp, &
                'nonnested rho and log kernel marginal-likelihood gradients')
        end do
        ! With delta off both instruments observe one scaled GP.
        call rbf_covariance(xl, xl, amp, ell, low, dl, status)
        call rbf_covariance(xl, xh, amp, ell, cross, dc, status)
        call rbf_covariance(xh, xh, amp, ell, high, dh, status)
        delta = 0
        call model%condition(low, cross, high, delta, nl, nh, rho, y(:2), y(3:), status)
        call require(status%code == 0, 'single GP condition')
        do j = 1,3
            k(3:,j+2) = k(3:,j+2)-.4_dp**2*exp(-.5_dp*((xh-xh(j))/.2_dp)**2)
        end do
        do j = 1,2
            ch(:,j) = rho**2*amp**2*exp(-.5_dp*((xh-xt(j))/ell)**2)
            target(:,j) = rho**2*amp**2*exp(-.5_dp*((xt-xt(j))/ell)**2)
        end do
        full_cross(3:,:) = ch
        chol = k; rhs(:,1) = y; rhs(:,2:) = full_cross
        call dposv('U',5,3,chol,5,rhs,5,info)
        call require(info == 0, 'single GP dense condition')
        call model%predict(target,cl,ch,prediction,cov,status)
        call require(status%code == 0, 'single GP prediction')
        expected = target-matmul(transpose(full_cross),rhs(:,2:))
        call require(maxval(abs(cov-expected)) < 2.e-14_dp, 'single GP posterior covariance')
        call require(maxval(abs(prediction-matmul(transpose(full_cross),rhs(:,1)))) &
            < 2.e-14_dp, 'single GP posterior mean')
    end subroutine

    subroutine single_site()
        type(autoregressive_gp_t) :: control
        real(dp) :: value(1), row(1, 1), variance, expected
        ! Conditional Gaussian moments: Kennedy and O'Hagan (2000),
        ! Biometrika 87, 1-13, https://doi.org/10.1093/biomet/87.1.1.
        call control%fit(reshape([0._dp], [1, 1]), reshape([0._dp], [1, 1]), &
            reshape([1._dp], [1, 1]), [.7_dp], 2._dp, 0._dp, 0._dp, status)
        call require(status%code == 0, 'single-site fit')
        call control%predict([.5_dp], [0._dp], reshape([0._dp], [1, 1]), &
            value, row, variance, status)
        expected = exp(-.5_dp*(.5_dp/.7_dp)**2)
        call require(abs(value(1) - expected) < 1.e-15_dp, 'single-site mean')
        call require(abs(row(1, 1) + .5_dp/.7_dp**2*expected) < 1.e-15_dp, &
            'single-site derivative')
        call require(abs(variance - 4*(1 - expected**2)) < 1.e-15_dp, &
            'single-site variance')
    end subroutine

    subroutine parallel_control()
        use, intrinsic :: iso_fortran_env, only: int64
        type(autoregressive_gp_t) :: control
        real(dp) :: sites(2, 256), low(1, 256), value(1), row(1, 2), variance
        integer(int64) :: first, last, rate
        integer :: k
        do k = 1, 256
            sites(:, k) = [real(mod(k - 1, 16), dp), real((k - 1)/16, dp)]/15
        end do
        low(1, :) = sin(sites(1, :)) + cos(sites(2, :))
        call control%fit(sites, low, .7_dp*low, [.2_dp, .2_dp], 1._dp, .7_dp, &
            1.e-4_dp, status)
        call require(status%code == 0, 'parallel training')
        call system_clock(first, rate)
        do k = 1, 1000
            call control%predict([.3_dp, .4_dp], [sin(.3_dp) + cos(.4_dp)], &
                reshape([cos(.3_dp), -sin(.4_dp)], [1, 2]), value, row, &
                variance, status)
        end do
        call system_clock(last)
        call require(status%code == 0, 'parallel prediction')
        call require(abs(value(1) - .7_dp*(sin(.3_dp) + cos(.4_dp))) < 1.e-14_dp, &
            'parallel zero correction')
        print '(A,ES14.6)', '256-site GP wall seconds per call=', &
            real(last - first, dp)/real(rate, dp)/1000
    end subroutine
    subroutine require(ok, msg)
        logical, intent(in)::ok
        character(*), intent(in)::msg
        if (.not. ok) then
            print *, msg
            error stop 1
        end if
    end subroutine
end program
