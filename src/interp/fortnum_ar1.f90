module fortnum_ar1
    !! Noisy two-level AR(1) conditioning by a Schur complement. Covariances
    !! and their derivatives are supplied by the caller, including constraints.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use fortnum_cholesky, only: cholesky_factorization_t
    use fortnum_generated_mf_distance, only: mf_distance
    use fortnum_generated_mf_covariance, only: mf_covariance
    implicit none
    private
    public :: ar1_covariance, rbf_covariance, covariance_factor

    type, public :: ar1_gp_t
        private
        type(cholesky_factorization_t) :: low_factor, conditional_factor
        real(dp), allocatable :: bridge(:, :), alpha(:), low_alpha(:)
        real(dp) :: log_density = 0
    contains
        procedure :: condition => ar1_condition
        procedure :: predict => ar1_predict
        procedure :: log_likelihood => ar1_log_likelihood
        procedure :: gradient => ar1_gradient
    end type
contains
    subroutine rbf_covariance(x, y, amplitude, length, k, derivatives, status)
        !! derivatives(:,:,1:2) are with respect to log amplitude/log length.
        real(dp), intent(in) :: x(:), y(:), amplitude, length
        real(dp), intent(out) :: k(:, :), derivatives(:, :, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: distance, dx, dc, dl, dd, da
        integer :: i, j
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'RBF: invalid arguments')
        if (any(shape(k) /= [size(x), size(y)])) return
        if (any(shape(derivatives) /= [size(x), size(y), 2])) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) return
        if (.not. all(ieee_is_finite([amplitude, length]))) return
        if (amplitude < 0 .or. length <= 0) return
        !$omp parallel do default(shared) private(i,j,distance,dx,dc,dl,dd,da) &
        !$omp if(size(x)*size(y)>16384)
        do j = 1, size(y)
            do i = 1, size(x)
                call mf_distance(x(i), y(j), length, distance, dx, dc, dl)
                call mf_covariance(distance, amplitude, k(i, j), dd, da)
                derivatives(i, j, :) = [da*amplitude, dd*dl*length]
            end do
        end do
        !$omp end parallel do
        if (.not. all(ieee_is_finite(k))) return
        if (.not. all(ieee_is_finite(derivatives))) return
        call status_set(status, 0, '')
    end subroutine

    subroutine covariance_factor(k, dk, l, dl, status)
        !! Differentiate the Cholesky recurrence; no numerical jitter is added.
        real(dp), intent(in) :: k(:, :), dk(:, :, :)
        real(dp), intent(out) :: l(:, :), dl(:, :, :)
        type(fortnum_status_t), intent(out) :: status
        type(cholesky_factorization_t) :: factor
        real(dp) :: v
        integer :: n, i, j, p
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'factor: invalid shape or covariance')
        n = size(k, 1)
        if (any(shape(k) /= [n, n]) .or. any(shape(l) /= [n, n])) return
        if (any(shape(dk) /= [n, n, size(dk, 3)])) return
        if (any(shape(dl) /= shape(dk))) return
        if (.not. symmetric_finite(k) .or. .not. all(ieee_is_finite(dk))) return
        call factor%factorize(k, status)
        if (status%code /= 0) return
        l = factor%lower
        dl = 0
        do p = 1, size(dk, 3)
            do j = 1, n
                v = dk(j, j, p) - 2*dot_product(dl(j, :j - 1, p), l(j, :j - 1))
                dl(j, j, p) = v/(2*l(j, j))
                do i = j + 1, n
                    v = dk(i, j, p) - dot_product(dl(i, :j - 1, p), l(j, :j - 1)) &
                        - dot_product(l(i, :j - 1), dl(j, :j - 1, p))
                    dl(i, j, p) = (v - l(i, j)*dl(j, j, p))/l(j, j)
                end do
            end do
        end do
    end subroutine

    subroutine ar1_covariance(low, cross, high, delta, noise_low, noise_high, &
            rho, k, rho_derivative, status)
        real(dp), intent(in) :: low(:, :), cross(:, :), high(:, :), delta(:, :)
        real(dp), intent(in) :: noise_low(:, :), noise_high(:, :), rho
        real(dp), allocatable, intent(out) :: k(:, :), rho_derivative(:, :)
        type(fortnum_status_t), intent(out) :: status
        integer :: nl, nh
        nl = size(low, 1); nh = size(high, 1)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'AR1: invalid covariance blocks')
        if (nl < 1 .or. nh < 1) return
        if (any(shape(low) /= [nl, nl]) .or. any(shape(high) /= [nh, nh])) return
        if (any(shape(cross) /= [nl, nh]) .or. any(shape(delta) /= [nh, nh])) return
        if (any(shape(noise_low) /= [nl, nl])) return
        if (any(shape(noise_high) /= [nh, nh])) return
        if (.not. symmetric_finite(low) .or. .not. symmetric_finite(high)) return
        if (.not. symmetric_finite(delta) .or. .not. symmetric_finite(noise_low)) return
        if (.not. symmetric_finite(noise_high)) return
        if (.not. all(ieee_is_finite(cross)) .or. .not. ieee_is_finite(rho)) return
        allocate (k(nl + nh, nl + nh), rho_derivative(nl + nh, nl + nh))
        k(:nl, :nl) = low + noise_low
        k(:nl, nl + 1:) = rho*cross
        k(nl + 1:, :nl) = transpose(k(:nl, nl + 1:))
        k(nl + 1:, nl + 1:) = rho**2*high + delta + noise_high
        rho_derivative = 0
        rho_derivative(:nl, nl + 1:) = cross
        rho_derivative(nl + 1:, :nl) = transpose(cross)
        rho_derivative(nl + 1:, nl + 1:) = 2*rho*high
        if (.not. all(ieee_is_finite(k))) return
        call status_set(status, 0, '')
    end subroutine

    subroutine ar1_condition(self, low, cross, high, delta, noise_low, noise_high, &
            rho, y_low, y_high, status)
        class(ar1_gp_t), intent(inout) :: self
        real(dp), intent(in) :: low(:, :), cross(:, :), high(:, :), delta(:, :)
        real(dp), intent(in) :: noise_low(:, :), noise_high(:, :), rho, y_low(:), y_high(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: k(:, :), dr(:, :), b(:, :), ah(:), al(:)
        real(dp) :: ld_low, ld_high
        integer :: nl, nh
        if (allocated(self%alpha)) deallocate (self%alpha, self%low_alpha, self%bridge)
        call ar1_covariance(low, cross, high, delta, noise_low, noise_high, rho, k, dr, status)
        if (status%code /= 0) return
        nl = size(low, 1); nh = size(high, 1)
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'AR1: invalid observations')
        if (size(y_low) /= nl .or. size(y_high) /= nh) return
        if (.not. all(ieee_is_finite(y_low)) .or. .not. all(ieee_is_finite(y_high))) return
        call self%low_factor%factorize(k(:nl, :nl), status)
        if (status%code /= 0) return
        b = k(:nl, nl + 1:)
        call self%low_factor%solve(b, status)
        if (status%code /= 0) return
        call self%conditional_factor%factorize( &
            k(nl + 1:, nl + 1:) - matmul(k(nl + 1:, :nl), b), status)
        if (status%code /= 0) return
        al = y_low
        call self%low_factor%solve(al, status)
        if (status%code /= 0) return
        ah = y_high - matmul(k(nl + 1:, :nl), al)
        call self%conditional_factor%solve(ah, status)
        if (status%code /= 0) return
        self%alpha = [al - matmul(b, ah), ah]
        self%low_alpha = al; self%bridge = b
        call self%low_factor%log_determinant(ld_low, status)
        if (status%code /= 0) return
        call self%conditional_factor%log_determinant(ld_high, status)
        if (status%code /= 0) return
        self%log_density = -0.5_dp*(dot_product([y_low, y_high], self%alpha) &
            + ld_low + ld_high &
            + (nl + nh)*log(2*acos(-1._dp)))
    end subroutine

    subroutine ar1_predict(self, prior, cross_low, cross_high, mean, covariance, status)
        class(ar1_gp_t), intent(in) :: self
        real(dp), intent(in) :: prior(:, :), cross_low(:, :), cross_high(:, :)
        real(dp), intent(out) :: mean(:), covariance(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: bl(:, :), bh(:, :), conditional_cross(:, :)
        integer :: nl, nh, np
        mean = 0; covariance = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'AR1: invalid prediction')
        if (.not. allocated(self%alpha)) return
        nl = size(self%bridge, 1); nh = size(self%bridge, 2); np = size(mean)
        if (any(shape(prior) /= [np, np])) return
        if (any(shape(covariance) /= [np, np])) return
        if (any(shape(cross_low) /= [nl, np])) return
        if (any(shape(cross_high) /= [nh, np])) return
        if (.not. symmetric_finite(prior)) return
        if (.not. all(ieee_is_finite(cross_low))) return
        if (.not. all(ieee_is_finite(cross_high))) return
        bl = cross_low
        call self%low_factor%solve(bl, status)
        if (status%code /= 0) return
        conditional_cross = cross_high - matmul(transpose(self%bridge), cross_low)
        bh = conditional_cross
        call self%conditional_factor%solve(bh, status)
        if (status%code /= 0) return
        mean = matmul(transpose(cross_low), self%low_alpha) &
            + matmul(transpose(conditional_cross), self%alpha(nl + 1:))
        covariance = prior - matmul(transpose(cross_low), bl) &
            - matmul(transpose(conditional_cross), bh)
        covariance = (covariance + transpose(covariance))/2
    end subroutine

    subroutine ar1_log_likelihood(self, value, status)
        class(ar1_gp_t), intent(in) :: self
        real(dp), intent(out) :: value
        type(fortnum_status_t), intent(out) :: status
        value = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'AR1: not conditioned')
        if (.not. allocated(self%alpha)) return
        value = self%log_density
        call status_set(status, 0, '')
    end subroutine

    subroutine ar1_gradient(self, derivatives, gradient, status, score)
        !! Exact log marginal likelihood gradient; derivatives include both
        !! noise and kernel terms. score is the derivative w.r.t. the means.
        class(ar1_gp_t), intent(in) :: self
        real(dp), intent(in) :: derivatives(:, :, :)
        real(dp), intent(out) :: gradient(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(out), optional :: score(:)
        real(dp), allocatable :: inverse(:, :), a(:, :), s(:, :)
        integer :: nl, nh, i, p, n
        gradient = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'AR1: invalid gradient')
        if (.not. allocated(self%alpha)) return
        nl = size(self%bridge, 1); nh = size(self%bridge, 2); n = nl + nh
        if (any(shape(derivatives) /= [n, n, size(gradient)])) return
        if (.not. all(ieee_is_finite(derivatives))) return
        if (present(score)) then
            if (size(score) /= n) return
        end if
        allocate (a(nl, nl), s(nh, nh), inverse(n, n))
        a = 0; s = 0
        do i = 1, nl
            a(i, i) = 1
        end do
        do i = 1, nh
            s(i, i) = 1
        end do
        call self%low_factor%solve(a, status)
        if (status%code /= 0) return
        call self%conditional_factor%solve(s, status)
        if (status%code /= 0) return
        inverse(nl + 1:, nl + 1:) = s
        inverse(:nl, nl + 1:) = -matmul(self%bridge, s)
        inverse(nl + 1:, :nl) = transpose(inverse(:nl, nl + 1:))
        inverse(:nl, :nl) = a - matmul(inverse(:nl, nl + 1:), transpose(self%bridge))
        do p = 1, size(gradient)
            gradient(p) = (dot_product(self%alpha, &
                matmul(derivatives(:, :, p), self%alpha)) &
                - sum(inverse*derivatives(:, :, p)))/2
        end do
        if (present(score)) score = self%alpha
    end subroutine

    logical function symmetric_finite(a) result(ok)
        real(dp), intent(in) :: a(:, :)
        ok = .false.
        if (size(a, 1) /= size(a, 2)) return
        if (.not. all(ieee_is_finite(a))) return
        ok = maxval(abs(a - transpose(a))) <= 64*epsilon(1._dp)*max(1._dp, maxval(abs(a)))
    end function
end module fortnum_ar1
