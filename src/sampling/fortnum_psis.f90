module fortnum_psis
    !! Pareto-smoothed importance sampling (Vehtari, Simpson, Gelman, Yao and
    !! Gabry 2024, JMLR 25(72)) and the generalized Pareto fit it uses (Zhang
    !! and Stephens 2009, Technometrics 51, 316, with the weakly informative
    !! shrinkage of k towards 0.5 used by Vehtari et al.).
    !!
    !! psis_smooth replaces the M = ceiling(min(S/5, 3 sqrt(S))) largest of S
    !! importance log weights by the expected order statistics of a generalized
    !! Pareto distribution fitted to their exceedances over the next-largest
    !! weight, truncates at the largest raw weight and returns the fitted shape
    !! k-hat. k-hat < 0.5 indicates finite-variance weights; with
    !! k-hat < min(1 - 1/log10(S), 0.7) the smoothed estimate is reliable
    !! (Vehtari et al. 2024, Sec. 3); larger values reject the proposal.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    private

    public :: psis_smooth, gpd_fit, psis_khat_threshold

contains

    subroutine psis_smooth(log_w, smoothed, khat, status)
        !! log_w: raw log importance ratios (any additive constant). smoothed:
        !! the Pareto-smoothed log weights, normalized so that their maximum
        !! is that of log_w. khat: fitted tail shape (huge if the tail is
        !! degenerate and no fit was possible, which callers should reject).
        real(dp), intent(in) :: log_w(:)
        real(dp), intent(out) :: smoothed(size(log_w))
        real(dp), intent(out) :: khat
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: sorted(:), tail(:)
        integer, allocatable :: order(:)
        real(dp) :: lw_max, cutoff, sigma, p
        integer :: s, m, i

        call status_set(status, FORTNUM_OK, '')
        s = size(log_w)
        khat = huge(1.0_dp)
        smoothed = log_w
        if (s < 10) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'psis_smooth: need at least 10 draws')
            return
        end if
        if (.not. all(ieee_is_finite(log_w))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'psis_smooth: non-finite log weight')
            return
        end if
        lw_max = maxval(log_w)
        smoothed = log_w - lw_max
        m = ceiling(min(0.2_dp * s, 3.0_dp * sqrt(real(s, dp))))
        order = sort_index(smoothed)
        sorted = smoothed(order)
        cutoff = sorted(s - m)
        if (sorted(s) - sorted(s - m + 1) <= epsilon(1.0_dp) / 100.0_dp) then
            smoothed = smoothed + lw_max
            return
        end if
        tail = exp(sorted(s - m + 1:s)) - exp(cutoff)
        call gpd_fit(tail, khat, sigma, status)
        if (status%code /= FORTNUM_OK) return
        if (ieee_is_finite(khat)) then
            do i = 1, m
                p = (i - 0.5_dp) / m
                smoothed(order(s - m + i)) = log(gpd_quantile(p, khat, sigma) + exp(cutoff))
            end do
        end if
        smoothed = min(smoothed, 0.0_dp) + lw_max
    end subroutine psis_smooth

    subroutine gpd_fit(x, k, sigma, status)
        !! Generalized Pareto shape k and scale sigma for exceedances x > 0
        !! (sorted ascending), by the profile-likelihood posterior mean of
        !! Zhang and Stephens (2009) with their grid of 30 + floor(sqrt(n))
        !! points and prior parameter 3, then k <- (n k + 5) / (n + 10).
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: k, sigma
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: theta(:), l_theta(:), w(:)
        real(dp) :: xstar, theta_hat, lmax
        integer :: n, m, j

        call status_set(status, FORTNUM_OK, '')
        n = size(x)
        if (n < 2) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'gpd_fit: need at least 2 exceedances')
            k = huge(1.0_dp)
            sigma = 0.0_dp
            return
        end if
        m = 30 + int(sqrt(real(n, dp)))
        xstar = x(max(1, floor(n / 4.0_dp + 0.5_dp)))
        allocate(theta(m), l_theta(m), w(m))
        do j = 1, m
            theta(j) = 1.0_dp / x(n) + (1.0_dp - sqrt(m / (j - 0.5_dp))) / 3.0_dp / xstar
            k = sum(log(1.0_dp - theta(j) * x)) / n
            l_theta(j) = n * (log(-theta(j) / k) - k - 1.0_dp)
        end do
        lmax = maxval(l_theta)
        w = exp(l_theta - lmax)
        w = w / sum(w)
        theta_hat = sum(theta * w)
        k = sum(log(1.0_dp - theta_hat * x)) / n
        sigma = -k / theta_hat
        k = (n * k + 0.5_dp * 10.0_dp) / (n + 10.0_dp)
    end subroutine gpd_fit

    pure real(dp) function psis_khat_threshold(s)
        !! Sample-size dependent k-hat threshold min(1 - 1/log10(S), 0.7).
        integer, intent(in) :: s
        psis_khat_threshold = min(1.0_dp - 1.0_dp / log10(real(max(s, 11), dp)), 0.7_dp)
    end function psis_khat_threshold

    pure real(dp) function gpd_quantile(p, k, sigma)
        real(dp), intent(in) :: p, k, sigma
        if (abs(k) < 1.0e-12_dp) then
            gpd_quantile = -sigma * log(1.0_dp - p)
        else
            gpd_quantile = sigma * (exp(-k * log(1.0_dp - p)) - 1.0_dp) / k
        end if
    end function gpd_quantile

    function sort_index(a) result(idx)
        !! Indices that sort a ascending (heap sort; stable order not needed).
        real(dp), intent(in) :: a(:)
        integer :: idx(size(a))
        integer :: n, i, root, child, end, t

        n = size(a)
        idx = [(i, i = 1, n)]
        do i = n / 2, 1, -1
            call sift(i, n)
        end do
        do end = n, 2, -1
            t = idx(1)
            idx(1) = idx(end)
            idx(end) = t
            call sift(1, end - 1)
        end do
    contains
        subroutine sift(start, last)
            integer, intent(in) :: start, last
            root = start
            do while (2 * root <= last)
                child = 2 * root
                if (child < last) then
                    if (a(idx(child + 1)) > a(idx(child))) child = child + 1
                end if
                if (a(idx(root)) >= a(idx(child))) return
                t = idx(root)
                idx(root) = idx(child)
                idx(child) = t
                root = child
            end do
        end subroutine sift
    end function sort_index

end module fortnum_psis
