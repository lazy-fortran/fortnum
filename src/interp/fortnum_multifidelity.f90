module fortnum_multifidelity
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_DOMAIN_ERROR
    use fortnum_generated_mf_distance, only: mf_distance
    use fortnum_generated_mf_covariance, only: mf_covariance
    use fortnum_generated_mf_autoregressive, only: mf_autoregressive
    implicit none
    private
    type, public :: autoregressive_gp_t
        private
        real(dp), allocatable :: x(:, :), length(:), alpha(:, :), factor(:, :)
        real(dp) :: amplitude = 0, rho = 1
    contains
        procedure :: fit => gp_fit
        procedure :: predict => gp_predict
        procedure :: jvp => gp_jvp
        procedure :: vjp => gp_vjp
    end type
    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine
        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import dp
            character, intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine gp_fit(self, x, low, high, length, amplitude, rho, nugget, status)
        class(autoregressive_gp_t), intent(inout) :: self
        real(dp), intent(in) :: x(:, :), low(:, :), high(:, :), length(:)
        real(dp), intent(in) :: amplitude, rho, nugget
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: k(:, :), a(:, :)
        real(dp) :: unused(size(x, 1))
        integer :: n, i, j, info
        if (allocated(self%x)) deallocate (self%x, self%length, self%alpha, self%factor)
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        'multifidelity: invalid training data')
        n = size(x, 2)
        if (n < 1 .or. size(x, 1) < 1 .or. size(low, 1) < 1) return
        if (any(shape(low) /= shape(high)) .or. size(low, 2) /= n) return
        if (size(length) /= size(x, 1)) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(low))) return
        if (.not. all(ieee_is_finite(high)) .or. &
            .not. all(ieee_is_finite(length))) return
        if (.not. all(ieee_is_finite([amplitude, rho, nugget]))) return
        if (any(length <= 0) .or. amplitude <= 0 .or. nugget < 0) return
        allocate (k(n, n))
        !$omp parallel do default(shared) private(i,j,unused) if(n>128)
        do j = 1, n
            do i = 1, n
                call covariance(x(:, i), x(:, j), length, amplitude, k(i, j), unused)
            end do
            k(j, j) = k(j, j) + nugget**2
        end do
        !$omp end parallel do
        call dpotrf('L', n, k, n, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            'multifidelity: singular covariance')
            return
        end if
        a = transpose(high - rho*low)
        call dpotrs('L', n, size(low, 1), k, n, a, n, info)
        if (info /= 0) return
        self%x = x; self%length = length; self%amplitude = amplitude; self%rho = rho
        self%factor = k; self%alpha = a
        call status_set(status, 0, '')
    end subroutine

    subroutine gp_predict(self, x, low, low_jac, value, jac, variance, status)
        class(autoregressive_gp_t), intent(in) :: self
        real(dp), intent(in) :: x(:), low(:), low_jac(:, :)
        real(dp), intent(out) :: value(:), jac(:, :), variance
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: k(:), dk(:, :), solved(:, :)
        real(dp) :: correction, dr, dl, dc
        integer :: i, j, n, info
        value = 0; jac = 0; variance = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
                        'multifidelity: invalid prediction input')
        if (.not. allocated(self%x)) return
        if (size(x) /= size(self%x, 1) .or. size(low) /= size(self%alpha, 2)) return
        if (any(shape(low_jac) /= [size(low), size(x)])) return
        if (size(value) /= size(low) .or. any(shape(jac) /= shape(low_jac))) return
        if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(low))) return
        if (.not. all(ieee_is_finite(low_jac))) return
        n = size(self%x, 2)
        allocate (k(n), dk(size(x), n), solved(n, 1))
        !$omp parallel do default(shared) private(i) if(n>128)
        do i = 1, n
            call covariance(x, self%x(:, i), self%length, self%amplitude, &
                            k(i), dk(:, i))
        end do
        !$omp end parallel do
        do j = 1, size(low)
            correction = dot_product(k, self%alpha(:, j))
            call mf_autoregressive(self%rho, low(j), correction, value(j), dr, dl, dc)
            jac(j, :) = dl*low_jac(j, :) + dc*matmul(dk, self%alpha(:, j))
        end do
        solved(:, 1) = k
        call dpotrs('L', n, 1, self%factor, n, solved, n, info)
        if (info /= 0) return
        variance = max(0.0_dp, self%amplitude**2 - dot_product(k, solved(:, 1)))
        call status_set(status, 0, '')
    end subroutine

    subroutine covariance(x, center, length, amplitude, value, row)
        real(dp), intent(in) :: x(:), center(:), length(:), amplitude
        real(dp), intent(out) :: value, row(:)
        real(dp) :: distance, term, dx, dc, dl, dd, da
        integer :: i
        distance = 0
        do i = 1, size(x)
            call mf_distance(x(i), center(i), length(i), term, dx, dc, dl)
            distance = distance + term
            row(i) = dx
        end do
        call mf_covariance(distance, amplitude, value, dd, da)
        row = dd*row
    end subroutine

    subroutine gp_jvp(self, x, low, low_jac, direction, value_dot, status)
        class(autoregressive_gp_t), intent(in) :: self
        real(dp), intent(in) :: x(:), low(:), low_jac(:, :), direction(:)
        real(dp), intent(out) :: value_dot(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: value(size(low)), jac(size(low), size(x)), variance
        call self%predict(x, low, low_jac, value, jac, variance, status)
        value_dot = 0
        if (status%code /= 0) return
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'multifidelity: direction size')
        if (size(direction) /= size(x) .or. size(value_dot) /= size(low)) return
        value_dot = matmul(jac, direction)
        call status_set(status, 0, '')
    end subroutine

    subroutine gp_vjp(self, x, low, low_jac, weight, x_bar, status)
        class(autoregressive_gp_t), intent(in) :: self
        real(dp), intent(in) :: x(:), low(:), low_jac(:, :), weight(:)
        real(dp), intent(out) :: x_bar(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp) :: value(size(low)), jac(size(low), size(x)), variance
        call self%predict(x, low, low_jac, value, jac, variance, status)
        x_bar = 0
        if (status%code /= 0) return
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'multifidelity: cotangent size')
        if (size(weight) /= size(low) .or. size(x_bar) /= size(x)) return
        x_bar = matmul(transpose(jac), weight)
        call status_set(status, 0, '')
    end subroutine
end module
