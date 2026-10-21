module fortnum_precision_vector
    use, intrinsic :: iso_fortran_env, only: int16, int64, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortgen_storage, only: pack_fp16, unpack_fp16, pack_bf16, unpack_bf16
    use fortnum_precision, only: FP64, FP32, FP16, BF16, valid_storage
    implicit none
    private
    type, public :: precision_vector_t
        private
        integer :: format=FP64, count=0
        real(real64) :: scaling=1
        real(real64), allocatable :: wide(:)
        real(real32), allocatable :: single(:)
        integer(int16), allocatable :: half(:)
    contains
        procedure :: store => vector_store
        procedure :: at => vector_at
        procedure :: load64 => vector_load64
        procedure :: load32 => vector_load32
        generic :: load => load64, load32
        procedure :: bytes => vector_bytes
        procedure :: size => vector_size
        procedure :: at_scaled => vector_at_scaled
        procedure :: scale_factor => vector_scale
        procedure :: is_fp64 => vector_is_fp64
    end type
contains
    subroutine vector_store(self, values, storage, ok)
        class(precision_vector_t), intent(inout) :: self
        real(real64), intent(in) :: values(:)
        integer, intent(in) :: storage
        logical, intent(out) :: ok
        type(precision_vector_t) :: next
        real(real64) :: maximum
        integer :: i
        ok = .false.
        if (.not. valid_storage(storage)) return
        if (.not. all(ieee_is_finite(values))) return
        maximum = maxval(abs(values), mask=.true.)
        if (size(values) == 0) maximum = 0
        next%count = size(values)
        next%format = storage
        ! One exact power-of-two scale protects range in all narrowed formats.
        if (maximum > 0 .and. storage /= FP64) &
            next%scaling = scale(1.0_real64, max(-1074, exponent(maximum)-1))
        select case(storage)
        case(FP64)
            next%wide = values
        case(FP32)
            next%single = real(values / next%scaling, real32)
        case(FP16, BF16)
            allocate(next%half(size(values)))
            if (storage == FP16) then
                !$omp parallel do simd if(size(values) >= 8192)
                do i=1,size(values)
                    next%half(i) = pack_fp16(values(i) / next%scaling)
                end do
            else
                !$omp parallel do simd if(size(values) >= 8192)
                do i=1,size(values)
                    next%half(i) = pack_bf16(values(i) / next%scaling)
                end do
            end if
        end select
        if (next%scaling > 1.0_real64) then
            do i=1,next%count
                if (abs(next%at_scaled(i)) > huge(1.0_real64)/next%scaling) &
                    return
            end do
        end if
        ! Commit only a successful conversion; failed input preserves old data.
        self%format = next%format
        self%count = next%count
        self%scaling = next%scaling
        call move_alloc(next%wide, self%wide)
        call move_alloc(next%single, self%single)
        call move_alloc(next%half, self%half)
        ok = .true.
    end subroutine

    pure real(real64) function vector_at(self, i) result(x)
        class(precision_vector_t), intent(in) :: self
        integer, intent(in) :: i
        select case(self%format)
        case(FP64); x = self%wide(i)
        case(FP32); x = real(self%single(i), real64) * self%scaling
        case(FP16); x = unpack_fp16(self%half(i)) * self%scaling
        case(BF16); x = unpack_bf16(self%half(i)) * self%scaling
        end select
    end function

    pure real(real64) function vector_at_scaled(self, i) result(x)
        class(precision_vector_t), intent(in) :: self
        integer, intent(in) :: i
        select case(self%format)
        case(FP64); x = self%wide(i)
        case(FP32); x = real(self%single(i), real64)
        case(FP16); x = unpack_fp16(self%half(i))
        case(BF16); x = unpack_bf16(self%half(i))
        end select
    end function

    pure real(real64) function vector_scale(self) result(x)
        class(precision_vector_t), intent(in) :: self
        x = self%scaling
    end function

    pure logical function vector_is_fp64(self) result(wide)
        class(precision_vector_t), intent(in) :: self
        wide = self%format == FP64
    end function

    subroutine vector_load64(self, values, ok)
        class(precision_vector_t), intent(in) :: self
        real(real64), intent(out) :: values(:)
        logical, intent(out) :: ok
        integer :: i
        ok = size(values) == self%count
        if (.not. ok) return
        !$omp parallel do simd if(self%count >= 8192)
        do i=1,self%count
            values(i) = self%at(i)
        end do
    end subroutine

    subroutine vector_load32(self, values, ok)
        class(precision_vector_t), intent(in) :: self
        real(real32), intent(out) :: values(:)
        logical, intent(out) :: ok
        integer :: i
        ok = size(values) == self%count
        if (.not. ok) return
        ! Unscaled fp32 compute requires representable magnitudes.
        do i=1,self%count
            if (abs(self%at(i)) > real(huge(0.0_real32), real64)) then
                ok = .false.
                return
            end if
        end do
        !$omp parallel do simd if(self%count >= 8192)
        do i=1,self%count
            values(i) = real(self%at(i), real32)
        end do
    end subroutine

    pure integer(int64) function vector_bytes(self) result(n)
        class(precision_vector_t), intent(in) :: self
        integer :: width
        width = 2
        if (self%format == FP32) width = 4
        if (self%format == FP64) width = 8
        n = int(width, int64) * self%count
    end function

    pure integer function vector_size(self) result(n)
        class(precision_vector_t), intent(in) :: self
        n = self%count
    end function
end module fortnum_precision_vector
