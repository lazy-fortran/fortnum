module fortnum_parallel_vector
    !! OpenMP vector kernels for Krylov solvers: axpy-type updates and
    !! reductions.
    !!
    !! Deterministic mode (the default): reductions sum fixed chunks of
    !! chunk_size elements serially inside each chunk and combine the chunk
    !! partials serially in chunk order.  The result depends only on the
    !! input, never on the thread count, scheduling, or whether OpenMP is
    !! enabled.  Fast mode uses an OpenMP reduction(+); its result may differ
    !! by roundoff between thread counts and repeat runs.
    !!
    !! Element-wise updates are independent per element and therefore bitwise
    !! identical in both modes.  No procedure opens a nested parallel region.
    use fortnum_kinds, only: dp
    implicit none

    private

    integer, parameter, public :: PARALLEL_VECTOR_CHUNK = 8192

    logical :: fast_mode = .false.

    public :: parallel_vector_axpy
    public :: parallel_vector_xpay
    public :: parallel_vector_sub
    public :: parallel_vector_scale
    public :: parallel_vector_dot
    public :: parallel_vector_norm2
    public :: parallel_vector_gemv_t
    public :: parallel_vector_gemv_n
    public :: parallel_vector_set_fast
    public :: parallel_vector_is_fast

contains

    subroutine parallel_vector_set_fast(fast)
        !! Select the non-deterministic OpenMP reduction for all reductions.
        logical, intent(in) :: fast

        fast_mode = fast
    end subroutine parallel_vector_set_fast

    function parallel_vector_is_fast() result(fast)
        logical :: fast

        fast = fast_mode
    end function parallel_vector_is_fast

    subroutine parallel_vector_axpy(alpha, x, y)
        !! y = y + alpha*x
        real(dp), intent(in) :: alpha, x(:)
        real(dp), intent(inout) :: y(:)
        integer :: i, n

        n = size(y)
        !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
        do i = 1, n
            y(i) = y(i) + alpha*x(i)
        end do
        !$omp end parallel do simd
    end subroutine parallel_vector_axpy

    subroutine parallel_vector_xpay(x, alpha, y)
        !! y = x + alpha*y
        real(dp), intent(in) :: x(:), alpha
        real(dp), intent(inout) :: y(:)
        integer :: i, n

        n = size(y)
        !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
        do i = 1, n
            y(i) = x(i) + alpha*y(i)
        end do
        !$omp end parallel do simd
    end subroutine parallel_vector_xpay

    subroutine parallel_vector_sub(a, b, c)
        !! c = a - b
        real(dp), intent(in) :: a(:), b(:)
        real(dp), intent(out) :: c(:)
        integer :: i, n

        n = size(c)
        !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
        do i = 1, n
            c(i) = a(i) - b(i)
        end do
        !$omp end parallel do simd
    end subroutine parallel_vector_sub

    subroutine parallel_vector_scale(alpha, x, y)
        !! y = alpha*x
        real(dp), intent(in) :: alpha, x(:)
        real(dp), intent(out) :: y(:)
        integer :: i, n

        n = size(y)
        !$omp parallel do simd schedule(static) if(n >= PARALLEL_VECTOR_CHUNK)
        do i = 1, n
            y(i) = alpha*x(i)
        end do
        !$omp end parallel do simd
    end subroutine parallel_vector_scale

    function parallel_vector_dot(x, y, fast) result(value)
        !! Dot product; deterministic unless fast (or the module mode) is set.
        real(dp), intent(in) :: x(:), y(:)
        logical, intent(in), optional :: fast
        real(dp) :: value

        real(dp), allocatable :: partial(:)
        real(dp) :: chunk_sum
        integer :: chunk, first, last, i, n, n_chunks
        logical :: use_fast

        n = size(x)
        use_fast = fast_mode
        if (present(fast)) use_fast = fast
        value = 0.0_dp
        if (n <= PARALLEL_VECTOR_CHUNK) then
            ! One chunk: identical to the chunked path, no thread overhead.
            do i = 1, n
                value = value + x(i)*y(i)
            end do
            return
        end if

        if (use_fast) then
            !$omp parallel do simd schedule(static) reduction(+:value)
            do i = 1, n
                value = value + x(i)*y(i)
            end do
            !$omp end parallel do simd
            return
        end if

        n_chunks = (n + PARALLEL_VECTOR_CHUNK - 1)/PARALLEL_VECTOR_CHUNK
        allocate (partial(n_chunks))
        !$omp parallel do schedule(static) private(first, last, chunk_sum, i)
        do chunk = 1, n_chunks
            first = (chunk - 1)*PARALLEL_VECTOR_CHUNK + 1
            last = min(n, chunk*PARALLEL_VECTOR_CHUNK)
            chunk_sum = 0.0_dp
            do i = first, last
                chunk_sum = chunk_sum + x(i)*y(i)
            end do
            partial(chunk) = chunk_sum
        end do
        !$omp end parallel do
        do chunk = 1, n_chunks
            value = value + partial(chunk)
        end do
    end function parallel_vector_dot

    function parallel_vector_norm2(x, fast) result(value)
        !! Euclidean norm sqrt(x.x) with the same reduction contract as dot.
        real(dp), intent(in) :: x(:)
        logical, intent(in), optional :: fast
        real(dp) :: value

        value = sqrt(parallel_vector_dot(x, x, fast))
    end function parallel_vector_norm2

    subroutine parallel_vector_gemv_t(basis, columns, w, coefficients)
        !! coefficients(j) = basis(:, j) . w for j = 1..columns.
        !! Always deterministic (fixed chunks, ordered combine), in either mode.
        real(dp), intent(in) :: basis(:, :), w(:)
        integer, intent(in) :: columns
        real(dp), intent(out) :: coefficients(:)

        real(dp), allocatable :: partial(:, :)
        real(dp) :: chunk_sum
        integer :: chunk, first, last, i, j, n, n_chunks

        n = size(w)
        coefficients(1:columns) = 0.0_dp
        if (columns < 1) return
        if (n <= PARALLEL_VECTOR_CHUNK) then
            do j = 1, columns
                chunk_sum = 0.0_dp
                do i = 1, n
                    chunk_sum = chunk_sum + basis(i, j)*w(i)
                end do
                coefficients(j) = chunk_sum
            end do
            return
        end if

        n_chunks = (n + PARALLEL_VECTOR_CHUNK - 1)/PARALLEL_VECTOR_CHUNK
        allocate (partial(columns, n_chunks))
        !$omp parallel do schedule(static) private(first, last, chunk_sum, i, j)
        do chunk = 1, n_chunks
            first = (chunk - 1)*PARALLEL_VECTOR_CHUNK + 1
            last = min(n, chunk*PARALLEL_VECTOR_CHUNK)
            do j = 1, columns
                chunk_sum = 0.0_dp
                do i = first, last
                    chunk_sum = chunk_sum + basis(i, j)*w(i)
                end do
                partial(j, chunk) = chunk_sum
            end do
        end do
        !$omp end parallel do
        do chunk = 1, n_chunks
            do j = 1, columns
                coefficients(j) = coefficients(j) + partial(j, chunk)
            end do
        end do
    end subroutine parallel_vector_gemv_t

    subroutine parallel_vector_gemv_n(alpha, basis, columns, coefficients, w)
        !! w = w + alpha*basis(:, 1:columns)*coefficients(1:columns).
        !! Each row chunk is updated by one thread with a fixed column order.
        real(dp), intent(in) :: alpha, basis(:, :), coefficients(:)
        integer, intent(in) :: columns
        real(dp), intent(inout) :: w(:)

        integer :: chunk, first, last, i, j, n, n_chunks

        n = size(w)
        n_chunks = (n + PARALLEL_VECTOR_CHUNK - 1)/PARALLEL_VECTOR_CHUNK
        !$omp parallel do schedule(static) private(first, last, i, j) &
        !$omp if(n_chunks > 1)
        do chunk = 1, n_chunks
            first = (chunk - 1)*PARALLEL_VECTOR_CHUNK + 1
            last = min(n, chunk*PARALLEL_VECTOR_CHUNK)
            do j = 1, columns
                do i = first, last
                    w(i) = w(i) + (alpha*coefficients(j))*basis(i, j)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine parallel_vector_gemv_n

end module fortnum_parallel_vector
