module fortnum_verified_matvec
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_support_denormal
    use fortnum_interval, only: interval_t, ipoint
    use fortnum_rounding, only: add_up, sub_up, sub_down, add_down, mul_up, gamma_up
    implicit none
    private
    public :: matvec_bound_t, prepare_matvec_bound, enclose_matvec, matvec_row_bound

    ! A plan is valid only for the immutable binary64 matrix used at preparation.
    ! It owns row magnitude bounds, never a second matrix or per-call allocation.
    type :: matvec_bound_t
        private
        real(dp), allocatable :: rows(:)
        real(dp) :: gamma = 0, underflow = 0
        integer :: columns = 0
        logical :: ready = .false.
    end type
contains
    real(dp) function matvec_row_bound(plan) result(bound)
        type(matvec_bound_t), intent(in) :: plan
        bound = huge(0._dp)
        if (plan%ready) bound = maxval(plan%rows)
    end function
    subroutine prepare_matvec_bound(a, plan, ok)
        real(dp), intent(in) :: a(:,:)
        type(matvec_bound_t), intent(out) :: plan
        logical, intent(out) :: ok
        integer :: i,j,n
        ok = .false.
        n = size(a,2)
        if (n < 1 .or. size(a,1) < 1) return
        if (n > huge(n)/2) return
        if (.not. all(ieee_is_finite(a))) return
        if (.not. ieee_support_denormal(0._dp)) return
        allocate(plan%rows(size(a,1)))
        plan%rows = 0
        do j = 1,n
            do i = 1,size(a,1)
                plan%rows(i) = add_up(plan%rows(i),abs(a(i,j)))
            end do
        end do
        plan%gamma = gamma_up(2*n)
        ! One absolute tiny allowance per product/add, amplified by gamma.
        ! tiny (smallest NORMAL) deliberately dominates half a subnormal ulp.
        plan%underflow = mul_up(mul_up(real(2*n,dp),tiny(0._dp)), &
            add_up(1._dp,plan%gamma))
        if (.not. all(ieee_is_finite(plan%rows))) return
        if (.not. ieee_is_finite(plan%gamma)) return
        plan%columns = n; plan%ready = .true.; ok = .true.
    end subroutine

    subroutine enclose_matvec(a, plan, x, ax, ok, accumulate, sign)
        real(dp), intent(in) :: a(:,:)
        type(matvec_bound_t), intent(in) :: plan
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(inout) :: ax(:)
        logical, intent(in), optional :: accumulate
        integer, intent(in), optional :: sign
        logical, intent(out) :: ok
        real(dp) :: width, xmax, error, magnitude, value, lower, upper
        logical :: adding
        integer :: factor
        integer :: i,j
        ok = .false.
        adding = .false.; factor = 1
        if (present(accumulate)) adding = accumulate
        if (present(sign)) factor = sign
        if (abs(factor) /= 1) return
        if (.not. adding) ax = ipoint(0._dp)
        if (.not. plan%ready) return
        if (size(a,2) /= plan%columns) return
        if (size(a,1) /= size(plan%rows)) return
        if (size(x) /= plan%columns) return
        if (size(ax) /= size(a,1)) return
        if (.not. all(ieee_is_finite(x%lo))) return
        if (.not. all(ieee_is_finite(x%hi))) return
        if (any(x%lo > x%hi)) return
        ! Lower endpoints are exact point centers; no midpoint-rounding premise.
        xmax = maxval(abs(x%lo)); width = 0
        do j = 1,size(x)
            width = max(width,sub_up(x(j)%hi,x(j)%lo))
        end do
        magnitude = add_up(mul_up(plan%gamma,xmax),width)
        do i = 1,size(a,1)
            value = 0
            do j = 1,size(x)
                value = value + a(i,j)*x(j)%lo
            end do
            value = factor*value
            error = add_up(mul_up(plan%rows(i),magnitude),plan%underflow)
            lower = sub_down(value,error); upper = add_up(value,error)
            if (adding) then
                lower = add_down(ax(i)%lo,lower)
                upper = add_up(ax(i)%hi,upper)
            end if
            ax(i)%lo = lower; ax(i)%hi = upper
        end do
        ok = all(ieee_is_finite(ax%lo)) .and. all(ieee_is_finite(ax%hi))
    end subroutine
end module
