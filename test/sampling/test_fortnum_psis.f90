program test_fortnum_psis
    ! Behavioral tests for fortnum_psis with independent oracles:
    ! - gpd_fit recovers the shape and scale of generalized Pareto draws
    !   produced by the inverse CDF;
    ! - importance sampling N(0, 1) with a wider normal proposal has bounded
    !   weights (true k <= 0) and k-hat < 0.5, and the smoothed weights
    !   reproduce E[x^2] = 1;
    ! - a narrower proposal N(0, 0.5^2) has weights with Pareto tail index
    !   1 / k = 1 / (1 - 0.5^2), i.e. k = 0.75 (density ratio
    !   exp(3 x^2 / 2) under x ~ N(0, 1/4)). At S = 4000 the estimator is
    !   pre-asymptotic: an independent numpy port of R loo's psis gives a mean
    !   k-hat of 0.68 (5-95 %: 0.51-0.87) over 200 replicates, so the test
    !   averages 40 replicates and checks that the mean lies in [0.6, 0.8];
    ! - cross-implementation values: gpd_fit and the tail k-hat on fixed
    !   deterministic inputs agree with that numpy port to 1e-12.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform, rng_normal
    use fortnum_psis, only: psis_smooth, gpd_fit, psis_khat_threshold
    implicit none

    integer :: nfail

    nfail = 0
    call test_gpd_recovery(nfail)
    call test_good_proposal(nfail)
    call test_bad_proposal(nfail)
    call test_input_checks(nfail)
    call test_reference_values(nfail)
    if (nfail > 0) then
        write(error_unit, '(I0,A)') nfail, ' PSIS check(s) failed'
        error stop 1
    end if
    print '(A)', 'fortnum_psis: all checks passed'

contains

    subroutine check(ok, what, nfail)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        integer, intent(inout) :: nfail
        if (.not. ok) then
            nfail = nfail + 1
            write(error_unit, '(A)') 'FAIL: '//what
        end if
    end subroutine check

    subroutine test_gpd_recovery(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        real(dp) :: x(4000), u, k, sigma, kk
        real(dp), parameter :: shapes(3) = [0.3_dp, 0.7_dp, -0.2_dp]
        integer :: i, t

        do t = 1, 3
            kk = shapes(t)
            call rng_seed(rng, int(100 + t, int64), st)
            do i = 1, size(x)
                call rng_uniform(rng, u)
                x(i) = 2.0_dp * ((1.0_dp - u)**(-kk) - 1.0_dp) / kk
            end do
            call sort(x)
            call gpd_fit(x, k, sigma, st)
            call check(st%code == FORTNUM_OK, 'gpd: status', nfail)
            ! standard error of k is about (1 + k) / sqrt(n) = 0.03
            call check(abs(k - kk) < 0.1_dp, 'gpd: shape recovered', nfail)
            call check(abs(sigma / 2.0_dp - 1.0_dp) < 0.15_dp, 'gpd: scale recovered', nfail)
        end do
    end subroutine test_gpd_recovery

    subroutine draws(seed, scale, x, lw)
        integer, intent(in) :: seed
        real(dp), intent(in) :: scale
        real(dp), intent(out) :: x(:), lw(:)
        type(rng_t) :: rng
        type(fortnum_status_t) :: st
        integer :: i
        call rng_seed(rng, int(seed, int64), st)
        do i = 1, size(x)
            call rng_normal(rng, x(i))
            x(i) = scale * x(i)
        end do
        ! log N(x; 0, 1) - log N(x; 0, scale^2) up to a constant
        lw = -0.5_dp * x**2 + 0.5_dp * (x / scale)**2
    end subroutine draws

    subroutine test_good_proposal(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: st
        real(dp) :: x(4000), lw(4000), sm(4000), w(4000), khat

        call draws(5, 1.5_dp, x, lw)
        call psis_smooth(lw, sm, khat, st)
        call check(st%code == FORTNUM_OK, 'good proposal: status', nfail)
        call check(khat < 0.5_dp, 'good proposal: k-hat < 0.5', nfail)
        call check(abs(maxval(sm) - maxval(lw)) < 1.0e-12_dp, 'good proposal: max preserved', nfail)
        w = exp(sm - maxval(sm))
        w = w / sum(w)
        call check(abs(sum(w * x**2) - 1.0_dp) < 0.1_dp, 'good proposal: E[x^2]', nfail)
    end subroutine test_good_proposal

    subroutine test_bad_proposal(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: st
        real(dp) :: x(4000), lw(4000), sm(4000), khat, mean_k
        integer :: r

        mean_k = 0.0_dp
        do r = 1, 40
            call draws(1000 + r, 0.5_dp, x, lw)
            call psis_smooth(lw, sm, khat, st)
            call check(st%code == FORTNUM_OK, 'bad proposal: status', nfail)
            call check(all(sm <= maxval(lw) + 1.0e-12_dp), 'bad proposal: truncated at max', nfail)
            mean_k = mean_k + khat / 40
        end do
        call check(mean_k > 0.6_dp .and. mean_k < 0.8_dp, 'bad proposal: mean k-hat', nfail)
        ! min(1 - 1/log10(S), 0.7): 0.7 at S = 4000, 0.5 at S = 100
        call check(abs(psis_khat_threshold(4000) - 0.7_dp) < 1.0e-15_dp, 'threshold S = 4000', nfail)
        call check(abs(psis_khat_threshold(100) - 0.5_dp) < 1.0e-15_dp, 'threshold S = 100', nfail)
    end subroutine test_bad_proposal

    subroutine test_reference_values(nfail)
        ! numpy port of R loo (scratch psis_ref.py, see PR); inputs are exact
        ! GPD quantiles (k = 0.5, sigma = 2) and a smooth log-weight ramp
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: st
        real(dp) :: x(200), lw(100), sm(100), k, sigma, u
        integer :: i

        do i = 1, 200
            u = (i - 0.5_dp) / 200
            x(i) = 2.0_dp * ((1.0_dp - u)**(-0.5_dp) - 1.0_dp) / 0.5_dp
        end do
        call gpd_fit(x, k, sigma, st)
        call check(abs(k - 0.5021985369632893_dp) < 1.0e-12_dp, 'reference: gpd k', nfail)
        call check(abs(sigma - 1.9902342358721519_dp) < 1.0e-12_dp, 'reference: gpd sigma', nfail)
        lw = [(log(1.0_dp + i / 7.0_dp)**2, i = 1, 100)]
        call psis_smooth(lw, sm, k, st)
        call check(abs(k - (-0.19733350257757673_dp)) < 1.0e-12_dp, 'reference: psis k-hat', nfail)
    end subroutine test_reference_values

    subroutine test_input_checks(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: st
        real(dp) :: lw(5), sm(5), khat
        lw = 0.0_dp
        call psis_smooth(lw, sm, khat, st)
        call check(st%code == FORTNUM_DOMAIN_ERROR, 'input: too few draws rejected', nfail)
    end subroutine test_input_checks

    subroutine sort(a)
        real(dp), intent(inout) :: a(:)
        real(dp) :: t
        integer :: i, j
        do i = 2, size(a)
            t = a(i)
            j = i - 1
            do while (j >= 1)
                if (a(j) <= t) exit
                a(j + 1) = a(j)
                j = j - 1
            end do
            a(j + 1) = t
        end do
    end subroutine sort

end program test_fortnum_psis
