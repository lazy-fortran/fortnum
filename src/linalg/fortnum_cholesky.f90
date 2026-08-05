module fortnum_cholesky
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    interface
        subroutine dpotrf(uplo, n, a, lda, info)
            import :: dp
            character(1), intent(in) :: uplo
            integer, intent(in) :: n, lda
            real(dp), intent(inout) :: a(lda, *)
            integer, intent(out) :: info
        end subroutine dpotrf

        subroutine dpotrs(uplo, n, nrhs, a, lda, b, ldb, info)
            import :: dp
            character(1), intent(in) :: uplo
            integer, intent(in) :: n, nrhs, lda, ldb
            real(dp), intent(in) :: a(lda, *)
            real(dp), intent(inout) :: b(ldb, *)
            integer, intent(out) :: info
        end subroutine dpotrs
    end interface

    type, public :: cholesky_factorization_t
        real(dp), allocatable :: lower(:, :)
        integer :: n = 0
    contains
        procedure, public :: factorize => cholesky_factorize
        procedure, public :: solve_vector => cholesky_solve_vector
        procedure, public :: solve_matrix => cholesky_solve_matrix
        generic, public :: solve => solve_vector, solve_matrix
        procedure, public :: log_determinant => cholesky_log_determinant
    end type cholesky_factorization_t

    public :: cholesky_factorize
    public :: cholesky_solve_vector
    public :: cholesky_solve_matrix
    public :: cholesky_log_determinant

contains

    subroutine cholesky_factorize(self, matrix, status)
        class(cholesky_factorization_t), intent(out) :: self
        real(dp), intent(in) :: matrix(:, :)
        type(fortnum_status_t), intent(out) :: status
        integer :: i, j, n, info

        n = size(matrix, 1)
        if (n < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky factorize: matrix must be nonempty")
            return
        end if
        if (size(matrix, 2) /= n) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky factorize: matrix must be square")
            return
        end if

        allocate(self%lower(n, n))
        self%lower = matrix
        call dpotrf("L", n, self%lower, n, info)
        if (info /= 0) then
            deallocate(self%lower)
            self%n = 0
            if (info > 0) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "cholesky factorize: matrix is not positive definite")
            else
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "cholesky factorize: LAPACK argument error")
            end if
            return
        end if
        do j = 1, n
            do i = 1, j - 1
                self%lower(i, j) = 0.0_dp
            end do
        end do
        self%n = n
        call status_set(status, FORTNUM_OK, "")
    end subroutine cholesky_factorize

    subroutine cholesky_solve_vector(self, rhs, status)
        class(cholesky_factorization_t), intent(in) :: self
        real(dp), intent(inout) :: rhs(:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: work(:, :)
        integer :: info

        if (.not. allocated(self%lower)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky solve: factorization is not initialized")
            return
        end if
        if (self%n < 1 .or. size(rhs) /= self%n) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky solve: right-hand-side shape is invalid")
            return
        end if

        allocate(work(self%n, 1))
        work(:, 1) = rhs
        call dpotrs("L", self%n, 1, self%lower, self%n, work, self%n, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "cholesky solve: LAPACK solve failed")
            return
        end if
        rhs = work(:, 1)
        call status_set(status, FORTNUM_OK, "")
    end subroutine cholesky_solve_vector

    subroutine cholesky_solve_matrix(self, rhs, status)
        class(cholesky_factorization_t), intent(in) :: self
        real(dp), contiguous, intent(inout) :: rhs(:, :)
        type(fortnum_status_t), intent(out) :: status
        integer :: info, nrhs

        if (.not. allocated(self%lower)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky solve: factorization is not initialized")
            return
        end if
        nrhs = size(rhs, 2)
        if (self%n < 1 .or. size(rhs, 1) /= self%n .or. nrhs < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky solve: right-hand-side shape is invalid")
            return
        end if

        call dpotrs("L", self%n, nrhs, self%lower, self%n, rhs, self%n, info)
        if (info /= 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "cholesky solve: LAPACK solve failed")
            return
        end if
        call status_set(status, FORTNUM_OK, "")
    end subroutine cholesky_solve_matrix

    subroutine cholesky_log_determinant(self, value, status)
        class(cholesky_factorization_t), intent(in) :: self
        real(dp), intent(out) :: value
        type(fortnum_status_t), intent(out) :: status
        integer :: i

        value = 0.0_dp
        if (.not. allocated(self%lower)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky log determinant: factorization is not initialized")
            return
        end if
        if (self%n < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cholesky log determinant: factorization size is invalid")
            return
        end if

        do i = 1, self%n
            value = value + 2.0_dp*log(self%lower(i, i))
        end do
        call status_set(status, FORTNUM_OK, "")
    end subroutine cholesky_log_determinant

end module fortnum_cholesky
