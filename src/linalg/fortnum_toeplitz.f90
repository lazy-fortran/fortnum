module fortnum_toeplitz
    !! Matrix-free Toeplitz products through circulant embedding and FFT.
    !!
    !! For column c and row r, the represented matrix is
    !! T(i,j) = c(i-j+1) for i >= j and r(j-i+1) otherwise.  The optional row
    !! defaults to the column, which is the symmetric covariance case.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_fft, only: fft_c2c, fortnum_fft_plan_t, fft_c2c_plan_init
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    implicit none
    private

    type, public :: toeplitz_operator_t
        real(dp), allocatable :: column(:), row(:)
        complex(dp), allocatable :: spectrum(:)
        type(fortnum_fft_plan_t) :: plan
        integer :: n = 0
        integer :: embedding_size = 0
    contains
        procedure, public :: initialize => toeplitz_initialize
        procedure, public :: matvec => toeplitz_matvec
        procedure, public :: matmat => toeplitz_matmat
        procedure, public :: diagonal => toeplitz_diagonal
        procedure, public :: element_count => toeplitz_element_count
    end type toeplitz_operator_t

contains

    subroutine toeplitz_initialize(self, column, status, row)
        class(toeplitz_operator_t), intent(out) :: self
        real(dp), intent(in) :: column(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: row(:)

        complex(dp), allocatable :: embedded(:)
        integer :: index, n

        self%n = 0
        self%embedding_size = 0
        if (size(column) < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "Toeplitz operator: column must be nonempty")
            return
        end if
        if (.not. all(ieee_is_finite(column))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "Toeplitz operator: column contains a nonfinite value")
            return
        end if
        if (present(row)) then
            if (size(row) /= size(column)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "Toeplitz operator: row and column sizes differ")
                return
            end if
            if (.not. all(ieee_is_finite(row))) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "Toeplitz operator: row contains a nonfinite value")
                return
            end if
        end if

        n = size(column)
        allocate(self%column(n), self%row(n), embedded(2*n), &
            self%spectrum(2*n))
        self%column = column
        if (present(row)) then
            self%row = row
        else
            self%row = column
        end if
        embedded = cmplx(0.0_dp, 0.0_dp, dp)
        do index = 1, n
            embedded(index) = cmplx(self%column(index), 0.0_dp, dp)
        end do
        do index = 2, n
            embedded(n + index) = cmplx( &
                self%row(n - index + 2), 0.0_dp, dp)
        end do
        self%spectrum = embedded
        call fft_c2c_plan_init(self%plan, 2*n)
        call fft_c2c(self%spectrum, -1, self%plan)
        self%n = n
        self%embedding_size = 2*n
        call status_set(status, FORTNUM_OK, "")
    end subroutine toeplitz_initialize

    subroutine toeplitz_matvec(self, input, output)
        class(toeplitz_operator_t), intent(in) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        complex(dp), allocatable :: work(:)
        integer :: index

        call validate_vector_shape(self, size(input), size(output))
        allocate(work(self%embedding_size))
        work = cmplx(0.0_dp, 0.0_dp, dp)
        do index = 1, self%n
            work(index) = cmplx(input(index), 0.0_dp, dp)
        end do
        call fft_c2c(work, -1, self%plan)
        work = work*self%spectrum
        call fft_c2c(work, 1, self%plan)
        do index = 1, self%n
            output(index) = real(work(index), dp)/real(self%embedding_size, dp)
        end do
    end subroutine toeplitz_matvec

    subroutine toeplitz_matmat(self, input, output)
        class(toeplitz_operator_t), intent(in) :: self
        real(dp), intent(in) :: input(:, :)
        real(dp), intent(out) :: output(:, :)

        complex(dp), allocatable :: work(:, :)
        integer :: index, column, n_rhs

        call validate_matrix_shape( &
            self, size(input, 1), size(input, 2), size(output, 1), &
            size(output, 2))
        n_rhs = size(input, 2)
        allocate(work(self%embedding_size, n_rhs))
        work = cmplx(0.0_dp, 0.0_dp, dp)
        do column = 1, n_rhs
            do index = 1, self%n
                work(index, column) = cmplx(input(index, column), 0.0_dp, dp)
            end do
            call fft_c2c(work(:, column), -1, self%plan)
            work(:, column) = work(:, column)*self%spectrum
            call fft_c2c(work(:, column), 1, self%plan)
            do index = 1, self%n
                output(index, column) = real(work(index, column), dp)/ &
                    real(self%embedding_size, dp)
            end do
        end do
    end subroutine toeplitz_matmat

    function toeplitz_diagonal(self) result(values)
        class(toeplitz_operator_t), intent(in) :: self
        real(dp), allocatable :: values(:)

        if (self%n < 1) error stop "Toeplitz operator: not initialized"
        allocate(values(self%n))
        values = self%column(1)
    end function toeplitz_diagonal

    integer function toeplitz_element_count(self) result(count)
        class(toeplitz_operator_t), intent(in) :: self

        count = self%n
    end function toeplitz_element_count

    subroutine validate_vector_shape(self, input_size, output_size)
        class(toeplitz_operator_t), intent(in) :: self
        integer, intent(in) :: input_size, output_size

        if (self%n < 1) error stop "Toeplitz operator: not initialized"
        if (input_size /= self%n .or. output_size /= self%n) then
            error stop "Toeplitz operator: vector has the wrong size"
        end if
    end subroutine validate_vector_shape

    subroutine validate_matrix_shape( &
            self, input_rows, input_columns, output_rows, output_columns)
        class(toeplitz_operator_t), intent(in) :: self
        integer, intent(in) :: input_rows, input_columns, output_rows
        integer, intent(in) :: output_columns

        if (self%n < 1) error stop "Toeplitz operator: not initialized"
        if (input_rows /= self%n .or. output_rows /= self%n) then
            error stop "Toeplitz operator: matrix has the wrong row count"
        end if
        if (input_columns < 1 .or. output_columns /= input_columns) then
            error stop "Toeplitz operator: matrix has the wrong RHS count"
        end if
    end subroutine validate_matrix_shape

end module fortnum_toeplitz
