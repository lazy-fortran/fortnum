module fortnum_gram_factor
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    type, public :: gram_factor_workspace_t
        private
        real(dp), allocatable :: a(:, :), tau(:), work(:)
    contains
        procedure :: factor => gram_factor
    end type
    interface
        subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
            import :: dp
            integer, intent(in) :: m, n, lda, lwork
            real(dp), intent(inout) :: a(lda, *)
            real(dp), intent(out) :: tau(*), work(*)
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine gram_factor(self, w, b, status)
        class(gram_factor_workspace_t), intent(inout) :: self
        real(dp), intent(in) :: w(:, :)
        real(dp), intent(out) :: b(:, :)
        type(fortnum_status_t), intent(out) :: status
        integer :: n, m, i, j, info, lwork
        real(dp) :: query(1)
        logical :: resize

        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Gram factor: invalid shape or data')
        b = 0
        n = size(w, 1)
        m = size(w, 2)
        if (n < 1) return
        if (size(b, 1) /= n .or. size(b, 2) /= n) return
        if (.not. all(ieee_is_finite(w))) return
        if (m == 0) then
            call status_set(status, FORTNUM_OK, '')
            return
        end if
        resize = .true.
        if (allocated(self%a)) then
            resize = size(self%a, 1) /= m .or. size(self%a, 2) /= n
        end if
        if (resize) then
            if (allocated(self%a)) deallocate(self%a, self%tau, self%work)
            allocate(self%a(m, n), self%tau(min(m, n)))
            call dgeqrf(m, n, self%a, m, self%tau, query, -1, info)
            if (info /= 0) then
                deallocate(self%a, self%tau)
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    'Gram factor: LAPACK workspace query failed')
                return
            end if
            lwork = max(1, n, ceiling(query(1)))
            allocate(self%work(lwork))
        end if
        self%a = transpose(w)
        call dgeqrf(m, n, self%a, m, self%tau, self%work, size(self%work), info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'Gram factor: LAPACK QR failed')
            return
        end if
        do j = 1, min(m, n)
            do i = j, n
                b(i, j) = self%a(j, i)
            end do
        end do
        if (.not. all(ieee_is_finite(b))) then
            b = 0
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                'Gram factor: result outside finite range')
            return
        end if
        call status_set(status, FORTNUM_OK, '')
    end subroutine
end module
