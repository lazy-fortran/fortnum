module fortnum_tensor_product
    !! Matrix-free tensor-product contractions for structured grids.
    !!
    !! Factor 1 is the innermost (fastest varying) grid dimension.  The
    !! represented matrix is therefore A_d \otimes ... \otimes A_2 \otimes
    !! A_1.  Applying one factor at a time avoids assembling that matrix and
    !! costs sum_i n_i^2 prod_{j /= i} n_j operations.
    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private

    type, public :: tensor_factor_t
        real(dp), allocatable :: values(:, :)
    end type tensor_factor_t

    type, public :: tensor_product_operator_t
        type(tensor_factor_t), allocatable :: factors(:)
        integer, allocatable :: dimensions(:)
        real(dp), pointer :: device_current(:) => null(), device_work(:) => null()
        real(dp), pointer :: device_current_matrix(:, :) => null()
        real(dp), pointer :: device_work_matrix(:, :) => null()
        integer :: total_size = 0
        integer :: device_rhs = 0
        logical :: device_resident = .false.
    contains
        procedure, public :: initialize => tensor_product_initialize
        procedure, public :: matvec => tensor_product_matvec
        procedure, public :: matmat => tensor_product_matmat
        procedure, public :: enter_data => tensor_product_enter_data
        procedure, public :: exit_data => tensor_product_exit_data
        procedure, public :: matvec_device => tensor_product_matvec_device
        procedure, public :: matmat_device => tensor_product_matmat_device
        procedure, public :: diagonal => tensor_product_diagonal
        procedure, public :: element_count => tensor_product_element_count
    end type tensor_product_operator_t

contains

    subroutine tensor_product_initialize(self, factors, status)
        class(tensor_product_operator_t), intent(out) :: self
        type(tensor_factor_t), intent(in) :: factors(:)
        type(fortnum_status_t), intent(out) :: status

        integer :: mode, n_modes, dimension
        integer(int64) :: total_size

        self%total_size = 0
        self%device_resident = .false.
        if (size(factors) < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: at least one factor is required")
            return
        end if

        total_size = 1_int64
        n_modes = size(factors)
        do mode = 1, n_modes
            if (.not. allocated(factors(mode)%values)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "tensor product: factor values are not allocated")
                return
            end if
            dimension = size(factors(mode)%values, 1)
            if (dimension < 1) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "tensor product: factor dimensions must be positive")
                return
            end if
            if (size(factors(mode)%values, 2) /= dimension) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "tensor product: factors must be square")
                return
            end if
            if (total_size > int(huge(0), int64)/int(dimension, int64)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "tensor product: total size overflows the index kind")
                return
            end if
            total_size = total_size*int(dimension, int64)
        end do

        allocate(self%factors(n_modes), self%dimensions(n_modes))
        do mode = 1, n_modes
            dimension = size(factors(mode)%values, 1)
            allocate(self%factors(mode)%values(dimension, dimension))
            self%factors(mode)%values = factors(mode)%values
            self%dimensions(mode) = dimension
        end do
        self%total_size = int(total_size)
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_initialize

    subroutine tensor_product_enter_data(self, status, n_rhs)
        class(tensor_product_operator_t), intent(inout) :: self
        type(fortnum_status_t), intent(out) :: status
        integer, intent(in), optional :: n_rhs

        integer :: mode, requested_rhs

        call validate_ready(self, status)
        if (status%code /= FORTNUM_OK) return
        requested_rhs = 0
        if (present(n_rhs)) requested_rhs = n_rhs
        if (requested_rhs < 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: device RHS count must be nonnegative")
            return
        end if
        if (self%device_resident) then
            if (requested_rhs > self%device_rhs) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "tensor product: device workspace already has fewer RHS")
                return
            end if
            call status_set(status, FORTNUM_OK, "")
            return
        end if
        allocate(self%device_current(self%total_size), &
            self%device_work(self%total_size))
        if (requested_rhs > 0) then
            allocate( &
                self%device_current_matrix(self%total_size, requested_rhs), &
                self%device_work_matrix(self%total_size, requested_rhs))
        end if
        do mode = 1, size(self%factors)
            !$acc enter data copyin(self%factors(mode)%values)
        end do
        !$acc enter data create(self%device_current, self%device_work)
        if (requested_rhs > 0) then
            !$acc enter data create( &
            !$acc& self%device_current_matrix, self%device_work_matrix)
        end if
        self%device_rhs = requested_rhs
        self%device_resident = .true.
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_enter_data

    subroutine tensor_product_exit_data(self, status)
        class(tensor_product_operator_t), intent(inout) :: self
        type(fortnum_status_t), intent(out) :: status

        integer :: mode

        if (.not. self%device_resident) then
            call status_set(status, FORTNUM_OK, "")
            return
        end if
        do mode = size(self%factors), 1, -1
            !$acc exit data delete(self%factors(mode)%values)
        end do
        if (associated(self%device_current_matrix)) then
            !$acc exit data delete( &
            !$acc& self%device_current_matrix, self%device_work_matrix)
            deallocate(self%device_current_matrix, self%device_work_matrix)
        end if
        !$acc exit data delete(self%device_current, self%device_work)
        deallocate(self%device_current, self%device_work)
        self%device_rhs = 0
        self%device_resident = .false.
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_exit_data

    subroutine tensor_product_matvec(self, input, output, status)
        class(tensor_product_operator_t), intent(in) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        type(fortnum_status_t), intent(out) :: status

        real(dp), allocatable :: current(:), work(:)
        integer :: mode, stride, block_count

        call validate_vector_shapes(self, size(input), size(output), status)
        if (status%code /= FORTNUM_OK) return

        allocate(current(self%total_size), work(self%total_size))
        current = input
        do mode = 1, size(self%dimensions)
            stride = mode_stride(self, mode)
            block_count = self%total_size/(stride*self%dimensions(mode))
            call apply_factor_vector( &
                self%factors(mode), self%dimensions(mode), stride, block_count, &
                current, work)
            current = work
        end do
        output = current
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_matvec

    subroutine tensor_product_matmat(self, input, output, status)
        class(tensor_product_operator_t), intent(in) :: self
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)
        type(fortnum_status_t), intent(out) :: status

        real(dp), allocatable :: current(:, :), work(:, :)
        integer :: mode, stride, block_count

        call validate_matrix_shapes( &
            self, size(input, 1), size(input, 2), size(output, 1), &
            size(output, 2), status)
        if (status%code /= FORTNUM_OK) return

        allocate( &
            current(self%total_size, size(input, 2)), &
            work(self%total_size, size(input, 2)))
        current = input
        do mode = 1, size(self%dimensions)
            stride = mode_stride(self, mode)
            block_count = self%total_size/(stride*self%dimensions(mode))
            call apply_factor_matrix( &
                self%factors(mode), self%dimensions(mode), stride, block_count, &
                current, work)
            current = work
        end do
        output = current
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_matmat

    subroutine tensor_product_matvec_device(self, input, output, status)
        class(tensor_product_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        type(fortnum_status_t), intent(out) :: status

        real(dp), pointer :: current(:), work(:)
        integer :: mode, stride, block_count, index, n

        call validate_vector_shapes(self, size(input), size(output), status)
        if (status%code /= FORTNUM_OK) return
        if (.not. self%device_resident) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: enter_data is required for device products")
            return
        end if

        if (.not. associated(self%device_current)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: device workspace is missing")
            return
        end if
        n = self%total_size
        current => self%device_current
        work => self%device_work
        !$acc data present(input, output, current, work)
        !$acc parallel loop
        do index = 1, n
            current(index) = input(index)
        end do
        do mode = 1, size(self%dimensions)
            stride = mode_stride(self, mode)
            block_count = self%total_size/(stride*self%dimensions(mode))
            call apply_factor_vector_acc( &
                self%factors(mode)%values, self%dimensions(mode), stride, &
                block_count, current, work)
            !$acc parallel loop
            do index = 1, n
                current(index) = work(index)
            end do
        end do
        !$acc parallel loop
        do index = 1, n
            output(index) = current(index)
        end do
        !$acc end data
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_matvec_device

    subroutine tensor_product_matmat_device(self, input, output, status)
        class(tensor_product_operator_t), intent(inout) :: self
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)
        type(fortnum_status_t), intent(out) :: status

        real(dp), pointer :: current(:, :), work(:, :)
        integer :: mode, stride, block_count, row, column, n, n_rhs

        call validate_matrix_shapes( &
            self, size(input, 1), size(input, 2), size(output, 1), &
            size(output, 2), status)
        if (status%code /= FORTNUM_OK) return
        if (.not. self%device_resident) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: enter_data is required for device products")
            return
        end if

        if (.not. associated(self%device_current_matrix) .or. &
            size(self%device_current_matrix, 2) < size(input, 2)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: device matrix workspace is too small")
            return
        end if
        n = self%total_size
        n_rhs = size(input, 2)
        current => self%device_current_matrix
        work => self%device_work_matrix
        !$acc data present(input, output, current, work)
        !$acc parallel loop collapse(2)
        do column = 1, n_rhs
            do row = 1, n
                current(row, column) = input(row, column)
            end do
        end do
        do mode = 1, size(self%dimensions)
            stride = mode_stride(self, mode)
            block_count = self%total_size/(stride*self%dimensions(mode))
            call apply_factor_matrix_acc( &
                self%factors(mode)%values, self%dimensions(mode), stride, &
                block_count, current, work)
            !$acc parallel loop collapse(2)
            do column = 1, n_rhs
                do row = 1, n
                    current(row, column) = work(row, column)
                end do
            end do
        end do
        !$acc parallel loop collapse(2)
        do column = 1, n_rhs
            do row = 1, n
                output(row, column) = current(row, column)
            end do
        end do
        !$acc end data
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_matmat_device

    subroutine tensor_product_diagonal(self, output, status)
        class(tensor_product_operator_t), intent(in) :: self
        real(dp), intent(out) :: output(:)
        type(fortnum_status_t), intent(out) :: status

        integer :: flat_index, mode, remainder, factor_index
        real(dp) :: value

        call validate_ready(self, status)
        if (status%code /= FORTNUM_OK) return
        if (size(output) /= self%total_size) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: diagonal has the wrong size")
            return
        end if

        do flat_index = 1, self%total_size
            remainder = flat_index - 1
            value = 1.0_dp
            do mode = 1, size(self%dimensions)
                factor_index = mod(remainder, self%dimensions(mode)) + 1
                remainder = remainder/self%dimensions(mode)
                value = value*self%factors(mode)%values( &
                    factor_index, factor_index)
            end do
            output(flat_index) = value
        end do
        call status_set(status, FORTNUM_OK, "")
    end subroutine tensor_product_diagonal

    integer function tensor_product_element_count(self) result(count)
        class(tensor_product_operator_t), intent(in) :: self

        count = self%total_size
    end function tensor_product_element_count

    subroutine validate_ready(self, status)
        class(tensor_product_operator_t), intent(in) :: self
        type(fortnum_status_t), intent(out) :: status

        if (.not. allocated(self%factors)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: operator is not initialized")
            return
        end if
        if (.not. allocated(self%dimensions)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: operator dimensions are missing")
            return
        end if
        if (self%total_size < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: operator size is invalid")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine validate_ready

    subroutine validate_vector_shapes(self, input_size, output_size, status)
        class(tensor_product_operator_t), intent(in) :: self
        integer, intent(in) :: input_size, output_size
        type(fortnum_status_t), intent(out) :: status

        call validate_ready(self, status)
        if (status%code /= FORTNUM_OK) return
        if (input_size /= self%total_size .or. output_size /= self%total_size) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: vector has the wrong size")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine validate_vector_shapes

    subroutine validate_matrix_shapes( &
            self, input_rows, input_columns, output_rows, output_columns, &
            status)
        class(tensor_product_operator_t), intent(in) :: self
        integer, intent(in) :: input_rows, input_columns, output_rows
        integer, intent(in) :: output_columns
        type(fortnum_status_t), intent(out) :: status

        call validate_ready(self, status)
        if (status%code /= FORTNUM_OK) return
        if (input_rows /= self%total_size .or. output_rows /= self%total_size) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: matrix has the wrong row count")
            return
        end if
        if (input_columns < 1 .or. output_columns /= input_columns) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "tensor product: matrix has the wrong RHS count")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine validate_matrix_shapes

    integer function mode_stride(self, mode) result(stride)
        class(tensor_product_operator_t), intent(in) :: self
        integer, intent(in) :: mode
        integer :: lower_mode

        stride = 1
        do lower_mode = 1, mode - 1
            stride = stride*self%dimensions(lower_mode)
        end do
    end function mode_stride

    subroutine apply_factor_vector( &
            factor, dimension, stride, block_count, input, output)
        type(tensor_factor_t), intent(in) :: factor
        integer, intent(in) :: dimension, stride, block_count
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        integer :: block, column, inner, offset, row, source, target
        real(dp) :: value

        do block = 0, block_count - 1
            offset = block*stride*dimension
            do inner = 1, stride
                do row = 1, dimension
                    value = 0.0_dp
                    do column = 1, dimension
                        source = offset + inner + (column - 1)*stride
                        value = value + factor%values(row, column)*input(source)
                    end do
                    target = offset + inner + (row - 1)*stride
                    output(target) = value
                end do
            end do
        end do
    end subroutine apply_factor_vector

    subroutine apply_factor_matrix( &
            factor, dimension, stride, block_count, input, output)
        type(tensor_factor_t), intent(in) :: factor
        integer, intent(in) :: dimension, stride, block_count
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)

        integer :: block, column, inner, offset, row, source, target

        do block = 0, block_count - 1
            offset = block*stride*dimension
            do inner = 1, stride
                do row = 1, dimension
                    target = offset + inner + (row - 1)*stride
                    output(target, :) = 0.0_dp
                    do column = 1, dimension
                        source = offset + inner + (column - 1)*stride
                        output(target, :) = output(target, :) + &
                            factor%values(row, column)*input(source, :)
                    end do
                end do
            end do
        end do
    end subroutine apply_factor_matrix

    subroutine apply_factor_vector_acc( &
            factor, dimension, stride, block_count, input, output)
        real(dp), intent(in) :: factor(:, :)
        integer, intent(in) :: dimension, stride, block_count
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        integer :: block, column, inner, offset, row, source, target
        real(dp) :: value

        !$acc parallel loop gang collapse(3) present(factor, input, output)
        do block = 0, block_count - 1
            do inner = 1, stride
                do row = 1, dimension
                    offset = block*stride*dimension
                    value = 0.0_dp
                    !$acc loop seq
                    do column = 1, dimension
                        source = offset + inner + (column - 1)*stride
                        value = value + factor(row, column)*input(source)
                    end do
                    target = offset + inner + (row - 1)*stride
                    output(target) = value
                end do
            end do
        end do
    end subroutine apply_factor_vector_acc

    subroutine apply_factor_matrix_acc( &
            factor, dimension, stride, block_count, input, output)
        real(dp), intent(in) :: factor(:, :)
        integer, intent(in) :: dimension, stride, block_count
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)

        integer :: block, column, inner, offset, row, source, target, rhs
        real(dp) :: value

        !$acc parallel loop gang collapse(3) present(factor, input, output)
        do block = 0, block_count - 1
            do inner = 1, stride
                do row = 1, dimension
                    offset = block*stride*dimension
                    do rhs = 1, size(input, 2)
                        value = 0.0_dp
                        !$acc loop seq
                        do column = 1, dimension
                            source = offset + inner + (column - 1)*stride
                            value = value + factor(row, column)* &
                                input(source, rhs)
                        end do
                        target = offset + inner + (row - 1)*stride
                        output(target, rhs) = value
                    end do
                end do
            end do
        end do
    end subroutine apply_factor_matrix_acc

end module fortnum_tensor_product
