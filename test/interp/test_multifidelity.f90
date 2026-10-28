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
