module fortnum_psd_factor
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: psd_factor
    interface
        subroutine dpstrf(uplo, n, a, lda, piv, rank, tol, work, info)
            import :: dp
            character(1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: piv(*), rank, info
            real(dp), intent(in) :: tol
            real(dp), intent(out) :: work(*)
        end subroutine
    end interface
contains
    subroutine psd_factor(a, b, status, rank, backward_error)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: b(:, :)
        type(fortnum_status_t), intent(out) :: status
        integer, intent(out), optional :: rank
        real(dp), intent(out), optional :: backward_error
        real(dp), allocatable :: lower(:, :), work(:), scaled(:, :)
        integer, allocatable :: piv(:)
        integer :: n, numerical_rank, info, i, j
        real(dp) :: scale, root_scale, tolerance, residual

        call status_set(status, FORTNUM_DOMAIN_ERROR, 'PSD factor: invalid shape or data')
        b = 0
        if (present(rank)) rank = 0
        if (present(backward_error)) backward_error = huge(1.0_dp)
        n = size(a, 1)
        if (n < 1) return
        if (size(a, 2) /= n) return
        if (size(b, 1) /= n .or. size(b, 2) /= n) return
        if (.not. all(ieee_is_finite(a))) return
        scale = maxval(abs(a))
        if (scale == 0) then
            if (present(backward_error)) backward_error = 0
            call status_set(status, FORTNUM_OK, '')
            return
        end if
        tolerance = 128.0_dp*epsilon(1.0_dp)*real(n, dp)
        if (maxval(abs(a/scale-transpose(a)/scale)) > tolerance) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'PSD factor: nonsymmetric matrix')
            return
        end if
        allocate(lower(n, n), work(2*n), piv(n), scaled(n, n))
        lower = a
        ! Zero tolerance retains every positive pivot, regardless of global scale.
        call dpstrf('L', n, lower, n, piv, numerical_rank, 0.0_dp, work, info)
        if (info < 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'PSD factor: LAPACK argument error')
            return
        end if
        do j = 1, numerical_rank
            do i = j, n
                b(piv(i), j) = lower(i, j)
            end do
        end do
        if (.not. all(ieee_is_finite(b))) then
            b = 0
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'PSD factor: nonfinite factor')
            return
        end if
        ! Check the original, unprojected matrix without squaring its scale.
        root_scale = sqrt(scale)
        scaled = b/root_scale
        lower = matmul(scaled, transpose(scaled))
        residual = maxval(abs(lower-a/scale))
        if (present(backward_error)) backward_error = residual
        if (.not. ieee_is_finite(residual) .or. residual > tolerance) then
            b = 0
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                'PSD factor: matrix fails positive semidefinite reconstruction')
            return
        end if
        if (present(rank)) rank = numerical_rank
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module
