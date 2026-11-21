module fortnum_blas
    ! Native BLAS on register tiles that fortgen emits in assembly for the
    ! running CPU (src/generated/fortnum_dtile_*.c). GEMM follows the
    ! GotoBLAS/BLIS design: B is packed in KC x NC blocks, A in MC x KC blocks,
    ! both into panels of the tile's height and width, and the tile kernel
    ! multiplies one A panel by one B panel. Complex products use the real tile
    ! on panels that stack real and imaginary parts, so ZGEMM runs at the real
    ! kernel's rate. Arguments and semantics are those of reference BLAS.
    use, intrinsic :: iso_c_binding, only: c_int64_t, c_int32_t, c_double, &
        c_funptr, c_associated, c_f_procpointer
    use fortnum_kinds, only: dp
    !$ use omp_lib, only: omp_in_parallel, omp_get_max_threads, omp_get_thread_num, &
    !$     omp_get_num_threads
    implicit none
    private

    public :: fortnum_dgemm, fortnum_zgemm, fortnum_blas_tile
    public :: fortnum_dgemm_ws, fortnum_zgemm_ws, fortnum_dgemm_lwork, fortnum_zgemm_lwork
    public :: fortnum_dgemm_triangle, thread_count
    public :: fortnum_dtrsm_packed, fortnum_dtrsm_packed_lwork

    abstract interface
        subroutine tile_kernel(k, a, b, acc) bind(c)
            import :: c_int64_t, c_double
            integer(c_int64_t), value :: k
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(out) :: acc(*)
        end subroutine tile_kernel
    end interface

    abstract interface
        subroutine tile_update_kernel(k, a, b, c, ldc) bind(c)
            import :: c_int64_t, c_double
            integer(c_int64_t), value :: k, ldc
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(inout) :: c(*)
        end subroutine tile_update_kernel
    end interface

    abstract interface
        subroutine tile_strided_kernel(k, a, lda, b, c, ldc) bind(c)
            import :: c_int64_t, c_double
            integer(c_int64_t), value :: k, lda, ldc
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(inout) :: c(*)
        end subroutine tile_strided_kernel
        subroutine tile_bstrided_kernel(k, a, b, ldb, c, ldc) bind(c)
            import :: c_int64_t, c_double
            integer(c_int64_t), value :: k, ldb, ldc
            real(c_double), intent(in) :: a(*), b(*)
            real(c_double), intent(inout) :: c(*)
        end subroutine tile_bstrided_kernel
    end interface

    abstract interface
        subroutine tile_trsm_kernel(k, a, b, c, ldc, t, pb) bind(c)
            import :: c_int64_t, c_double
            integer(c_int64_t), value :: k, ldc
            real(c_double), intent(in) :: a(*), b(*), t(*)
            real(c_double), intent(inout) :: c(*)
            real(c_double), intent(out) :: pb(*)
        end subroutine tile_trsm_kernel
    end interface

    interface
        function fortnum_blas_dtrsm_tile(isa) bind(c, name='fortnum_blas_dtrsm_tile') result(kernel)
            import :: c_funptr, c_int32_t
            integer(c_int32_t), value :: isa
            type(c_funptr) :: kernel
        end function fortnum_blas_dtrsm_tile
        function fortnum_blas_dtile(mr, nr, isa, update, strided_update, strided_acc, &
                bstrided_update, bstrided_acc) bind(c, name='fortnum_blas_dtile') result(kernel)
            import :: c_funptr, c_int32_t
            integer(c_int32_t), intent(out) :: mr, nr, isa
            type(c_funptr), intent(out) :: update, strided_update, strided_acc, bstrided_update, bstrided_acc
            type(c_funptr) :: kernel
        end function fortnum_blas_dtile
    end interface

    !> Depth of a packed block, and the targets for the A and B block sizes
    !> (rounded to whole tiles), set with the tile; per-machine values come
    !> from fortnum-bench.
    integer, save :: kc_depth = 512, mc_rows = 120, nc_cols = 4032
    !> Largest tile (rows times columns) any emitted kernel uses.
    integer, parameter :: MAX_TILE = 512
    !> Products with at most this many columns read A unpacked.
    ! fortnum_dgemm multiplies directly up to this many multiply-adds: the
    ! blocked kernel's workspace, packing and threads cost more there.
    real(dp), parameter :: GEMM_SMALL_WORK = 8.0e3_dp
    integer, parameter :: THIN_N = 96
    !> Work per thread below which fewer threads are used (tens of microseconds).
    real(dp), parameter :: THREAD_FLOPS = 2.0e6_dp

    !> FORTNUM_BLAS_PROFILE set: per-thread wall time in packing B, packing A,
    !> kernel with C update, and barriers, printed after each GEMM.
    logical, save :: profile = .false.
    logical, save :: tile_ready = .false.
    integer, save :: tile_mr = 0, tile_nr = 0, tile_isa = 0
    procedure(tile_kernel), pointer, save :: tile => null()
    !> The same tile adding to C in place, used for full tiles.
    procedure(tile_update_kernel), pointer, save :: tile_c => null()
    !> The tile reading A unpacked: added to C, and to a tile buffer.
    procedure(tile_strided_kernel), pointer, save :: tile_sc => null(), tile_sa => null()
    !> The tile reading B unpacked: added to C, and to a tile buffer.
    procedure(tile_bstrided_kernel), pointer, save :: tile_bc => null(), tile_ba => null()
    !> The TRSM tile of the GEMM tile, if there is one.
    procedure(tile_trsm_kernel), pointer, save :: tile_trsm => null()

contains

    !> The tile in use: rows, columns and instruction set (1 AVX2, 2 AVX-512,
    !> 3 NEON, 4 Apple AMX, 5 WebAssembly SIMD128).
    subroutine fortnum_blas_tile(mr, nr, isa)
        integer, intent(out) :: mr, nr, isa
        call select_tile()
        mr = tile_mr
        nr = tile_nr
        isa = tile_isa
    end subroutine fortnum_blas_tile

    subroutine select_tile()
        integer(c_int32_t) :: mr, nr, isa
        type(c_funptr) :: kernel, update, strided_update, strided_acc, bstrided_update, bstrided_acc
        integer :: status
        if (tile_ready) return
        !$omp critical (fortnum_blas_tile_selection)
        if (.not. tile_ready) then
            kernel = fortnum_blas_dtile(mr, nr, isa, update, strided_update, strided_acc, &
                bstrided_update, bstrided_acc)
            if (.not. c_associated(kernel)) &
                error stop 'fortnum_blas: no register-tile kernel for this CPU '// &
                    '(needs AVX2+FMA, AVX-512, NEON or wasm SIMD128)'
            call c_f_procpointer(kernel, tile)
            call c_f_procpointer(update, tile_c)
            call c_f_procpointer(strided_update, tile_sc)
            call c_f_procpointer(strided_acc, tile_sa)
            call c_f_procpointer(bstrided_update, tile_bc)
            call c_f_procpointer(bstrided_acc, tile_ba)
            tile_mr = mr
            tile_nr = nr
            tile_isa = isa
            kernel = fortnum_blas_dtrsm_tile(isa)
            if (c_associated(kernel)) call c_f_procpointer(kernel, tile_trsm)
            if (isa == 4) then
                ! AMX streams its panels from L2; one pass over C for k up to
                ! 2048 (2000^3 on one M3 core: 373 instead of 348 GFLOP/s)
                kc_depth = 2048
                mc_rows = 128
                nc_cols = 2048
            end if
            call get_environment_variable('FORTNUM_BLAS_PROFILE', status=status)
            profile = status == 0
            tile_ready = .true.
        end if
        !$omp end critical (fortnum_blas_tile_selection)
    end subroutine select_tile

    !> C := alpha op(A) op(B) + beta C, op(X) = X or X^T.
    subroutine fortnum_dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        real(dp), allocatable :: work(:)
        integer :: lwork
        if (real(m, dp)*n*k <= GEMM_SMALL_WORK) then
            call dgemm_small(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
            return
        end if
        lwork = 0
        if (m > 0 .and. n > 0 .and. k > 0) &
            lwork = fortnum_dgemm_lwork(m, n, k, thread_count(2.0_dp*m*real(n, dp)*k))
        allocate (work(lwork))
        call fortnum_dgemm_ws(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, work, lwork)
    end subroutine fortnum_dgemm

    !> fortnum_dgemm for small products: the loops of reference DGEMM, with
    !> neither workspace nor packing nor threads.
    subroutine dgemm_small(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        logical :: ta, tb
        integer :: i, j, l
        real(dp) :: t
        ta = transposed(transa, 'fortnum_dgemm', 1, .false.)
        tb = transposed(transb, 'fortnum_dgemm', 2, .false.)
        call check_sizes('fortnum_dgemm', ta, tb, m, n, k, lda, ldb, ldc)
        if (m == 0 .or. n == 0) return
        if ((alpha == 0 .or. k == 0) .and. beta == 1) return
        do j = 1, n
            if (beta == 0) then
                c(1:m, j) = 0
            else if (beta /= 1) then
                c(1:m, j) = beta*c(1:m, j)
            end if
            if (alpha == 0) cycle
            if (.not. ta) then
                do l = 1, k
                    if (tb) then
                        t = alpha*b(j, l)
                    else
                        t = alpha*b(l, j)
                    end if
                    c(1:m, j) = c(1:m, j) + t*a(1:m, l)
                end do
            else
                do i = 1, m
                    t = 0
                    if (tb) then
                        do l = 1, k
                            t = t + a(l, i)*b(j, l)
                        end do
                    else
                        do l = 1, k
                            t = t + a(l, i)*b(l, j)
                        end do
                    end if
                    c(i, j) = c(i, j) + alpha*t
                end do
            end if
        end do
    end subroutine dgemm_small

    !> C := alpha op(A) op(B) + beta C on the triangle uplo ('L' or 'U') of the
    !> square C (m = n) only; the other triangle is not touched. The core of
    !> the SYRK family: tiles outside the triangle are skipped, so the cost is
    !> about half a GEMM at GEMM speed and threading.
    subroutine fortnum_dgemm_triangle(uplo, transa, transb, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: uplo, transa, transb
        integer, intent(in) :: n, k, lda, ldb, ldc
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        real(dp), allocatable :: work(:)
        integer :: lwork
        lwork = 0
        if (n > 0 .and. k > 0) lwork = fortnum_dgemm_lwork(n, n, k, thread_count(1.0_dp*n*real(n, dp)*k))
        allocate (work(lwork))
        call fortnum_dgemm_ws(transa, transb, n, n, k, alpha, a, lda, b, ldb, beta, c, ldc, work, lwork, &
            triangle=uplo)
    end subroutine fortnum_dgemm_triangle

    !> Workspace in real(dp) elements that fortnum_dgemm_ws needs for an
    !> m x n x k product on up to `threads` OpenMP threads (default: as many
    !> as a call outside a parallel region may use). Nondecreasing in each
    !> argument, so the size for the largest dimensions covers every smaller
    !> product.
    integer function fortnum_dgemm_lwork(m, n, k, threads) result(lwork)
        integer, intent(in) :: m, n, k
        integer, intent(in), optional :: threads
        integer :: kd, nt
        call select_tile()
        nt = thread_limit()
        if (present(threads)) nt = max(1, threads)
        kd = min(kc_depth, max(0, k))
        lwork = max(2*b_block(max(0, n), kd, tile_nr) + nt*a_block(max(0, m), kd, tile_mr), &
            tile_nr*kd*((min(max(0, n), THIN_N) + tile_nr - 1)/tile_nr) + tile_mr*kd, &
            tile_mr*kd*((min(max(0, m), THIN_N) + tile_mr - 1)/tile_mr) + tile_nr*kd)
    end function fortnum_dgemm_lwork

    !> fortnum_dgemm with caller-provided packing buffers: work(lwork), lwork
    !> at least fortnum_dgemm_lwork(m, n, k); no memory is allocated. With
    !> less workspace than the threads need, fewer threads are used. With
    !> triangle = 'L' or 'U' (m = n) only that triangle of C is formed.
    subroutine fortnum_dgemm_ws(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, &
            work, lwork, triangle)
        character, intent(in) :: transa, transb
        character, intent(in), optional :: triangle
        integer, intent(in) :: m, n, k, lda, ldb, ldc, lwork
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(in) :: a(lda, *), b(ldb, *)
        real(dp), intent(inout) :: c(ldc, *)
        real(dp), intent(inout), target :: work(lwork)
        logical :: ta, tb
        integer :: jc, pc, nc, kc, mstep, nstep, npb, jr, ir, j, tri, ilo, ihi
        integer :: nthreads, mcb, mt, i0, i1, me, chunk, nchunk, ct, epoch, bsize, asize, boff
        real(dp) :: tp(4), w0, w1
        real(dp), allocatable :: prof(:, :)
        real(dp), pointer, contiguous :: ap(:), bp(:)

        ta = transposed(transa, 'fortnum_dgemm', 1, .false.)
        tb = transposed(transb, 'fortnum_dgemm', 2, .false.)
        call check_sizes('fortnum_dgemm', ta, tb, m, n, k, lda, ldb, ldc)
        ! tri: 0 full C, 1 lower triangle (i >= j), 2 upper triangle (i <= j)
        tri = 0
        if (present(triangle)) then
            select case (triangle)
            case ('L', 'l')
                tri = 1
            case ('U', 'u')
                tri = 2
            case default
                call argument_error('fortnum_dgemm_ws', 16)
            end select
            if (m /= n) call argument_error('fortnum_dgemm_ws', 3)
        end if
        if (m == 0 .or. n == 0) return
        if ((alpha == 0 .or. k == 0) .and. beta == 1) return
        do j = 1, n
            ilo = 1
            ihi = m
            if (tri == 1) ilo = j
            if (tri == 2) ihi = j
            if (beta == 0) then
                c(ilo:ihi, j) = 0
            else if (beta /= 1) then
                c(ilo:ihi, j) = beta*c(ilo:ihi, j)
            end if
        end do
        if (alpha == 0 .or. k == 0) return

        call select_tile()
        ! thin products: A is read unpacked, only B is packed
        if (tri == 0 .and. .not. ta .and. n <= THIN_N .and. m >= tile_mr) then
            call thin_gemm()
            return
        end if
        ! products with few rows: B is read unpacked, only A is packed
        if (tri == 0 .and. .not. tb .and. m <= THIN_N .and. n >= tile_nr .and. associated(tile_bc)) then
            call wide_gemm()
            return
        end if
        mstep = tile_mr
        nstep = tile_nr
        mcb = mc_size(mstep)
        ! buffers for the depth and rows the product has, not the block targets
        bsize = b_block(n, min(kc_depth, k), nstep)
        asize = a_block(m, min(kc_depth, k), mstep)
        if (lwork < 2*bsize + asize) call argument_error('fortnum_dgemm_ws', 15)
        nthreads = min(thread_count(2.0_dp*m*real(n, dp)*k), (lwork - 2*bsize)/asize)
        ! two B buffers: threads may pack the next block while others finish this one
        bp => work(1:2*bsize)
        mt = (m + mstep - 1)/mstep
        ! chunks of A as tall as MC allows but about four per thread, so that
        ! dynamic claiming balances unequal cores while A is packed only once
        ct = max(1, min(mcb/mstep, (mt + 4*nthreads - 1)/(4*nthreads)))
        nchunk = (mt + ct - 1)/ct
        if (profile) allocate (prof(4, 0:nthreads - 1), source=0.0_dp)
        ! One region per call. All threads pack each KC x NC block of B
        ! together; C is then split over a pr x pq grid of threads, each of
        ! which packs its own rows of A in MC chunks and multiplies them with
        ! its range of B panels. Two barriers per block of B, balance to one
        ! tile row or column.
        !$omp parallel num_threads(nthreads) if(nthreads > 1) default(shared) &
        !$omp private(jc, pc, nc, kc, npb, jr, ir, i0, i1, chunk, me, ap, epoch, boff, tp, w0, w1)
        tp = 0
        me = 0
        !$ me = omp_get_thread_num()
        ap => work(2*bsize + me*asize + 1:2*bsize + (me + 1)*asize)
        epoch = 0
        do jc = 1, n, nc_size(nstep)
            nc = min(nc_size(nstep), n - jc + 1)
            npb = (nc + nstep - 1)/nstep
            do pc = 1, k, kc_depth
                kc = min(kc_depth, k - pc + 1)
                boff = mod(epoch, 2)*bsize
                epoch = epoch + 1
                if (profile) w0 = wall()
                !$omp do schedule(static)
                do jr = 1, npb
                    call pack_b_real(tb, b, ldb, pc, kc, jc + (jr - 1)*nstep, &
                        min(nstep, nc - (jr - 1)*nstep), bp(boff + (jr - 1)*tile_nr*kc + 1:))
                end do
                !$omp end do nowait
                if (profile) call lap(tp(1), w0)
                !$omp barrier
                if (profile) call lap(tp(4), w0)
                !$omp do schedule(dynamic, 1)
                do chunk = 0, nchunk - 1
                    i0 = 1 + chunk*ct
                    i1 = min(mt, i0 + ct - 1)
                    do ir = i0, i1
                        call pack_a_real(ta, a, lda, alpha, 1 + (ir - 1)*mstep, &
                                min(mstep, m - (ir - 1)*mstep), pc, kc, ap((ir - i0)*tile_mr*kc + 1:))
                    end do
                    if (profile) call lap(tp(2), w0)
                    do jr = 1, npb
                        do ir = i0, i1
                            call update_real(ir, jr, jc, nc, kc, ap((ir - i0)*tile_mr*kc + 1:), &
                                bp(boff + (jr - 1)*tile_nr*kc + 1:))
                        end do
                    end do
                    if (profile) call lap(tp(3), w0)
                end do
                !$omp end do nowait
            end do
        end do
        if (profile) prof(:, me) = tp
        !$omp end parallel
        if (profile) call report_profile('dgemm', m, n, k, prof)

    contains

        ! C += alpha A op(B) for n <= THIN_N without packing A: B (alpha folded
        ! in) is packed per depth block, full row tiles of A are read in place
        ! by the strided kernels, the last partial row tile is packed.
        subroutine thin_gemm()
            real(dp), pointer, contiguous :: bt(:), ap1(:)
            real(dp) :: acc(MAX_TILE)
            integer :: pcc, kcc, npbt, jrr, irr, mfull, i0, j0, nrr, jj, nth, nb
            npbt = (n + tile_nr - 1)/tile_nr
            mfull = (m/tile_mr)*tile_mr
            nb = tile_nr*min(kc_depth, k)*npbt
            if (lwork < nb + tile_mr*min(kc_depth, k)) call argument_error('fortnum_dgemm_ws', 15)
            bt => work(1:nb)
            ap1 => work(nb + 1:nb + tile_mr*min(kc_depth, k))
            nth = thread_count(2.0_dp*m*real(n, dp)*k)
            do pcc = 1, k, kc_depth
                kcc = min(kc_depth, k - pcc + 1)
                do jrr = 1, npbt
                    call pack_b_real(tb, b, ldb, pcc, kcc, 1 + (jrr - 1)*tile_nr, &
                        min(tile_nr, n - (jrr - 1)*tile_nr), bt((jrr - 1)*tile_nr*kcc + 1:))
                end do
                if (alpha /= 1) bt(1:tile_nr*kcc*npbt) = alpha*bt(1:tile_nr*kcc*npbt)
                !$omp parallel do num_threads(nth) if(nth > 1) schedule(static) private(i0, jrr, j0, nrr, jj, acc)
                do irr = 1, mfull/tile_mr
                    i0 = 1 + (irr - 1)*tile_mr
                    do jrr = 1, npbt
                        j0 = 1 + (jrr - 1)*tile_nr
                        nrr = min(tile_nr, n - j0 + 1)
                        if (nrr == tile_nr) then
                            call tile_sc(int(kcc, c_int64_t), a(i0, pcc), int(lda, c_int64_t), &
                                bt((jrr - 1)*tile_nr*kcc + 1:), c(i0, j0), int(ldc, c_int64_t))
                        else
                            call tile_sa(int(kcc, c_int64_t), a(i0, pcc), int(lda, c_int64_t), &
                                bt((jrr - 1)*tile_nr*kcc + 1:), acc, int(tile_mr, c_int64_t))
                            do jj = 1, nrr
                                c(i0:i0 + tile_mr - 1, j0 + jj - 1) = c(i0:i0 + tile_mr - 1, j0 + jj - 1) + &
                                    acc((jj - 1)*tile_mr + 1:jj*tile_mr)
                            end do
                        end if
                    end do
                end do
                !$omp end parallel do
                if (mfull < m) then
                    ! the last rows through a packed panel (alpha is already in B)
                    call pack_a_real(.false., a, lda, 1.0_dp, mfull + 1, m - mfull, pcc, kcc, ap1)
                    do jrr = 1, npbt
                        j0 = 1 + (jrr - 1)*tile_nr
                        nrr = min(tile_nr, n - j0 + 1)
                        call tile(int(kcc, c_int64_t), ap1, bt((jrr - 1)*tile_nr*kcc + 1:), acc)
                        do jj = 1, nrr
                            c(mfull + 1:m, j0 + jj - 1) = c(mfull + 1:m, j0 + jj - 1) + &
                                acc((jj - 1)*tile_mr + 1:(jj - 1)*tile_mr + m - mfull)
                        end do
                    end do
                end if
            end do
        end subroutine thin_gemm

        ! C += alpha op(A) B for m <= THIN_N without packing B: A (alpha folded
        ! in) is packed per depth block, full column panels of B are read in
        ! place by the B-strided kernels, the last partial panel is packed.
        subroutine wide_gemm()
            real(dp), pointer, contiguous :: at(:), bp1(:)
            real(dp) :: acc(MAX_TILE)
            integer :: pcc, kcc, npat, irr, jrr, nfull, i0, j0, mrr, ii, jj, nth, na
            npat = (m + tile_mr - 1)/tile_mr
            nfull = (n/tile_nr)*tile_nr
            na = tile_mr*min(kc_depth, k)*npat
            if (lwork < na + tile_nr*min(kc_depth, k)) call argument_error('fortnum_dgemm_ws', 15)
            at => work(1:na)
            bp1 => work(na + 1:na + tile_nr*min(kc_depth, k))
            nth = thread_count(2.0_dp*m*real(n, dp)*k)
            do pcc = 1, k, kc_depth
                kcc = min(kc_depth, k - pcc + 1)
                do irr = 1, npat
                    call pack_a_real(ta, a, lda, alpha, 1 + (irr - 1)*tile_mr, &
                        min(tile_mr, m - (irr - 1)*tile_mr), pcc, kcc, at((irr - 1)*tile_mr*kcc + 1:))
                end do
                !$omp parallel do num_threads(nth) if(nth > 1) schedule(static) private(j0, irr, i0, mrr, ii, jj, acc)
                do jrr = 1, nfull/tile_nr
                    j0 = 1 + (jrr - 1)*tile_nr
                    do irr = 1, npat
                        i0 = 1 + (irr - 1)*tile_mr
                        mrr = min(tile_mr, m - i0 + 1)
                        if (mrr == tile_mr) then
                            call tile_bc(int(kcc, c_int64_t), at((irr - 1)*tile_mr*kcc + 1:), b(pcc, j0), &
                                int(ldb, c_int64_t), c(i0, j0), int(ldc, c_int64_t))
                        else
                            call tile_ba(int(kcc, c_int64_t), at((irr - 1)*tile_mr*kcc + 1:), b(pcc, j0), &
                                int(ldb, c_int64_t), acc, int(tile_mr, c_int64_t))
                            do jj = 1, tile_nr
                                do ii = 1, mrr
                                    c(i0 + ii - 1, j0 + jj - 1) = c(i0 + ii - 1, j0 + jj - 1) + acc((jj - 1)*tile_mr + ii)
                                end do
                            end do
                        end if
                    end do
                end do
                !$omp end parallel do
                if (nfull < n) then
                    ! the last columns through a packed panel (alpha is already in A)
                    call pack_b_real(.false., b, ldb, pcc, kcc, nfull + 1, n - nfull, bp1)
                    do irr = 1, npat
                        i0 = 1 + (irr - 1)*tile_mr
                        mrr = min(tile_mr, m - i0 + 1)
                        call tile(int(kcc, c_int64_t), at((irr - 1)*tile_mr*kcc + 1:), bp1, acc)
                        do jj = 1, n - nfull
                            c(i0:i0 + mrr - 1, nfull + jj) = c(i0:i0 + mrr - 1, nfull + jj) + &
                                acc((jj - 1)*tile_mr + 1:(jj - 1)*tile_mr + mrr)
                        end do
                    end do
                end if
            end do
        end subroutine wide_gemm

        subroutine update_real(ir, jr, jc, nc, kc, apanel, bpanel)
            integer, intent(in) :: ir, jr, jc, nc, kc
            real(dp), intent(in) :: apanel(*), bpanel(*)
            real(dp) :: acc(MAX_TILE)
            integer :: i0, j0, mrr, nrr, jj, lo, hi
            logical :: cross
            i0 = 1 + (ir - 1)*mstep
            j0 = jc + (jr - 1)*nstep
            mrr = min(mstep, m - i0 + 1)
            nrr = min(nstep, nc - (jr - 1)*nstep)
            if (tri == 1 .and. i0 + mrr - 1 < j0) return
            if (tri == 2 .and. i0 > j0 + nrr - 1) return
            cross = (tri == 1 .and. i0 < j0 + nrr - 1) .or. (tri == 2 .and. i0 + mrr - 1 > j0)
            if (mrr == mstep .and. nrr == nstep .and. .not. cross) then
                call tile_c(int(kc, c_int64_t), apanel, bpanel, c(i0, j0), &
                    int(ldc, c_int64_t))
                return
            end if
            call tile(int(kc, c_int64_t), apanel, bpanel, acc)
            do jj = 1, nrr
                lo = 1
                hi = mrr
                if (tri == 1) lo = max(1, j0 + jj - i0)
                if (tri == 2) hi = min(mrr, j0 + jj - i0)
                if (lo > hi) cycle
                c(i0 + lo - 1:i0 + hi - 1, j0 + jj - 1) = c(i0 + lo - 1:i0 + hi - 1, j0 + jj - 1) + &
                    acc((jj - 1)*tile_mr + lo:(jj - 1)*tile_mr + hi)
            end do
        end subroutine update_real

    end subroutine fortnum_dgemm_ws

    !> C := alpha op(A) op(B) + beta C, op(X) = X, X^T or X^H.
    subroutine fortnum_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        complex(dp), intent(in) :: alpha, beta
        complex(dp), intent(in) :: a(lda, *), b(ldb, *)
        complex(dp), intent(inout) :: c(ldc, *)
        real(dp), allocatable :: work(:)
        integer :: lwork
        lwork = 0
        if (m > 0 .and. n > 0 .and. k > 0) &
            lwork = fortnum_zgemm_lwork(m, n, k, thread_count(8.0_dp*m*real(n, dp)*k))
        allocate (work(lwork))
        call fortnum_zgemm_ws(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, work, lwork)
    end subroutine fortnum_zgemm

    !> Workspace in real(dp) elements that fortnum_zgemm_ws needs; as
    !> fortnum_dgemm_lwork.
    integer function fortnum_zgemm_lwork(m, n, k, threads) result(lwork)
        integer, intent(in) :: m, n, k
        integer, intent(in), optional :: threads
        integer :: kd, nt
        call select_tile()
        nt = thread_limit()
        if (present(threads)) nt = max(1, threads)
        kd = min(kc_depth, max(0, k))
        lwork = 2*b_block(max(0, n), kd, tile_nr/2) + nt*a_block(max(0, m), kd, tile_mr/2)
    end function fortnum_zgemm_lwork

    !> fortnum_zgemm with caller-provided packing buffers, real(dp) work(lwork)
    !> with lwork at least fortnum_zgemm_lwork(m, n, k); as fortnum_dgemm_ws.
    subroutine fortnum_zgemm_ws(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc, &
            work, lwork)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc, lwork
        complex(dp), intent(in) :: alpha, beta
        complex(dp), intent(in) :: a(lda, *), b(ldb, *)
        complex(dp), intent(inout) :: c(ldc, *)
        real(dp), intent(inout), target :: work(lwork)
        logical :: ta, tb, ca, cb
        integer :: jc, pc, nc, kc, mstep, nstep, npb, jr, ir, j
        integer :: nthreads, mcb, mt, i0, i1, me, chunk, nchunk, ct, epoch, bsize, asize, boff
        real(dp), pointer, contiguous :: ap(:), bp(:)

        ta = transposed(transa, 'fortnum_zgemm', 1, .true., ca)
        tb = transposed(transb, 'fortnum_zgemm', 2, .true., cb)
        call check_sizes('fortnum_zgemm', ta, tb, m, n, k, lda, ldb, ldc)
        if (m == 0 .or. n == 0) return
        if ((alpha == (0.0_dp, 0.0_dp) .or. k == 0) .and. beta == (1.0_dp, 0.0_dp)) return
        do j = 1, n
            if (beta == (0.0_dp, 0.0_dp)) then
                c(1:m, j) = (0.0_dp, 0.0_dp)
            else if (beta /= (1.0_dp, 0.0_dp)) then
                c(1:m, j) = beta*c(1:m, j)
            end if
        end do
        if (alpha == (0.0_dp, 0.0_dp) .or. k == 0) return

        call select_tile()
        ! a real tile mr x nr holds a complex tile mr/2 x nr/2
        mstep = tile_mr/2
        nstep = tile_nr/2
        mcb = mc_size(mstep)
        ! buffers for the depth and rows the product has, not the block targets
        bsize = b_block(n, min(kc_depth, k), nstep)
        asize = a_block(m, min(kc_depth, k), mstep)
        if (lwork < 2*bsize + asize) call argument_error('fortnum_zgemm_ws', 15)
        nthreads = min(thread_count(8.0_dp*m*real(n, dp)*k), (lwork - 2*bsize)/asize)
        ! two B buffers: threads may pack the next block while others finish this one
        bp => work(1:2*bsize)
        mt = (m + mstep - 1)/mstep
        ! chunks of A as tall as MC allows but about four per thread, so that
        ! dynamic claiming balances unequal cores while A is packed only once
        ct = max(1, min(mcb/mstep, (mt + 4*nthreads - 1)/(4*nthreads)))
        nchunk = (mt + ct - 1)/ct
        ! One region per call. All threads pack each KC x NC block of B
        ! together; C is then split over a pr x pq grid of threads, each of
        ! which packs its own rows of A in MC chunks and multiplies them with
        ! its range of B panels. Two barriers per block of B, balance to one
        ! tile row or column.
        !$omp parallel num_threads(nthreads) if(nthreads > 1) default(shared) &
        !$omp private(jc, pc, nc, kc, npb, jr, ir, i0, i1, chunk, me, ap, epoch, boff)
        me = 0
        !$ me = omp_get_thread_num()
        ap => work(2*bsize + me*asize + 1:2*bsize + (me + 1)*asize)
        epoch = 0
        do jc = 1, n, nc_size(nstep)
            nc = min(nc_size(nstep), n - jc + 1)
            npb = (nc + nstep - 1)/nstep
            do pc = 1, k, kc_depth
                kc = min(kc_depth, k - pc + 1)
                boff = mod(epoch, 2)*bsize
                epoch = epoch + 1
                !$omp do schedule(static)
                do jr = 1, npb
                    call pack_b_complex(tb, cb, b, ldb, pc, kc, jc + (jr - 1)*nstep, &
                        min(nstep, nc - (jr - 1)*nstep), nstep, bp(boff + (jr - 1)*tile_nr*kc + 1:))
                end do
                !$omp end do
                !$omp do schedule(dynamic, 1)
                do chunk = 0, nchunk - 1
                    i0 = 1 + chunk*ct
                    i1 = min(mt, i0 + ct - 1)
                    do ir = i0, i1
                        call pack_a_complex(ta, ca, a, lda, alpha, 1 + (ir - 1)*mstep, &
                                min(mstep, m - (ir - 1)*mstep), mstep, pc, kc, ap((ir - i0)*tile_mr*kc + 1:))
                    end do
                    do jr = 1, npb
                        do ir = i0, i1
                            call update_complex(ir, jr, jc, nc, kc, ap((ir - i0)*tile_mr*kc + 1:), &
                                bp(boff + (jr - 1)*tile_nr*kc + 1:))
                        end do
                    end do
                end do
                !$omp end do nowait
            end do
        end do
        !$omp end parallel

    contains

        subroutine update_complex(ir, jr, jc, nc, kc, apanel, bpanel)
            integer, intent(in) :: ir, jr, jc, nc, kc
            real(dp), intent(in) :: apanel(*), bpanel(*)
            real(dp) :: acc(MAX_TILE)
            integer :: i0, j0, mrr, nrr, ii, jj, rr, ri
            i0 = 1 + (ir - 1)*mstep
            j0 = jc + (jr - 1)*nstep
            mrr = min(mstep, m - i0 + 1)
            nrr = min(nstep, nc - (jr - 1)*nstep)
            call tile(int(kc, c_int64_t), apanel, bpanel, acc)
            do jj = 1, nrr
                ! columns jj (real part of B) and nstep+jj (imaginary part)
                rr = (jj - 1)*tile_mr
                ri = (nstep + jj - 1)*tile_mr
                do ii = 1, mrr
                    c(i0 + ii - 1, j0 + jj - 1) = c(i0 + ii - 1, j0 + jj - 1) + &
                        cmplx(acc(rr + ii) - acc(ri + mstep + ii), acc(ri + ii) + acc(rr + mstep + ii), dp)
                end do
            end do
        end subroutine update_complex

    end subroutine fortnum_zgemm_ws

    !> Workspace in real(dp) elements of fortnum_dtrsm_packed for the
    !> triangle of order mt and nt right-hand sides (columns of B for side L,
    !> rows for side R).
    integer function fortnum_dtrsm_packed_lwork(mt, nt) result(lwork)
        integer, intent(in) :: mt, nt
        integer :: kb, np
        call select_tile()
        kb = trsm_block(mt)
        np = kb/tile_mr
        lwork = tile_nr*kb*((min(nt, nc_size(tile_nr)) + tile_nr - 1)/tile_nr) &
            + tile_mr*tile_mr*np*(np + 1)/2 + tile_mr*kb*(mc_size(tile_mr)/tile_mr)
    end function fortnum_dtrsm_packed_lwork

    ! Rows of the diagonal blocks: whole tiles, at most TRSM_KC, and no more
    ! than the triangle needs.
    pure integer function trsm_block(mt) result(kb)
        integer, intent(in) :: mt
        integer, parameter :: TRSM_KC = 256
        kb = tile_mr*((min(mt, TRSM_KC) + tile_mr - 1)/tile_mr)
    end function trsm_block

    !> B := inv(op(A)) B (left) or B inv(op(A)) (right) for real A
    !> triangular, on packed panels and the GEMM register tile; alpha is
    !> applied by the caller. All eight cases are the forward solve L X = B'
    !> of a lower triangle L, read from A with reversed indices (op(A) upper
    !> on the left, lower on the right) and swapped ones (A^T on the left,
    !> A on the right), with B' = B (left) or B^T (right).
    !>
    !> Design after BLIS (BSD-3-Clause, github.com/flame/blis):
    !> frame/3/trsm/bli_trsm_ll_ker_var2.c and the gemmtrsm micro-kernel,
    !> with the diagonal blocks packed with inverted diagonal (invdiag in
    !> frame/1m/packm/bli_packm_struc_cxk.c). For each block of KB
    !> rows, B's rows are packed once into NR-wide panels; each MR-row
    !> micro-panel is updated by the tile kernel against the rows solved
    !> before it in the packed panel, solved in registers-sized tiles and
    !> written back to the panel and to B; the rows below are then updated
    !> from the same packed panel, so the solution is never packed twice.
    subroutine fortnum_dtrsm_packed(left, upper, notrans, unit, m, n, a, lda, b, ldb, work, lwork)
        logical, intent(in) :: left, upper, notrans, unit
        integer, intent(in) :: m, n, lda, ldb, lwork
        real(dp), intent(in) :: a(lda, *)
        real(dp), intent(inout) :: b(ldb, *)
        real(dp), intent(inout), target :: work(lwork)
        real(dp), pointer, contiguous :: bp(:), ad(:), ar(:)
        real(dp) :: acc(MAX_TILE)
        integer :: mt, nt, kb, ncb, jc, ncc, pc, kcc, npb, npd, ic, mcc, npr, nb_, na, nth, &
            jr, ir, q
        logical :: fwd, swap

        call select_tile()
        mt = merge(m, n, left)
        nt = merge(n, m, left)
        if (lwork < fortnum_dtrsm_packed_lwork(mt, nt)) call argument_error('fortnum_dtrsm_ws', 13)
        fwd = left .neqv. (upper .eqv. notrans)
        swap = left .neqv. notrans
        kb = trsm_block(mt)
        ncb = nc_size(tile_nr)
        nb_ = tile_nr*kb*((min(nt, ncb) + tile_nr - 1)/tile_nr)
        npd = kb/tile_mr
        na = tile_mr*tile_mr*npd*(npd + 1)/2
        bp => work(1:nb_)
        ad => work(nb_ + 1:nb_ + na)
        ar => work(nb_ + na + 1:)

        do jc = 1, nt, ncb
            ncc = min(ncb, nt - jc + 1)
            npb = (ncc + tile_nr - 1)/tile_nr
            do pc = 1, mt, kb
                kcc = min(kb, mt - pc + 1)
                npd = (kcc + tile_mr - 1)/tile_mr
                nth = thread_count(real(ncc, dp)*kcc*kcc)
                !$omp parallel num_threads(nth) if(nth > 1) private(jr, ir, q, acc)
                !$omp do schedule(static)
                do jr = 1, npb
                    call pack_rhs(jc + (jr - 1)*tile_nr, min(tile_nr, ncc - (jr - 1)*tile_nr), pc, kcc, &
                        bp((jr - 1)*tile_nr*kcc + 1:jr*tile_nr*kcc))
                end do
                !$omp end do nowait
                !$omp do schedule(static)
                do ir = 1, npd
                    q = tile_mr*tile_mr*(ir - 1)*ir/2
                    call pack_diag(pc, (ir - 1)*tile_mr, min(tile_mr, kcc - (ir - 1)*tile_mr), &
                        ad(q + 1:q + tile_mr*tile_mr*ir))
                end do
                !$omp end do
                !$omp do schedule(static)
                do jr = 1, npb
                    call solve_panel(jc + (jr - 1)*tile_nr, min(tile_nr, ncc - (jr - 1)*tile_nr), pc, kcc, &
                        npd, bp((jr - 1)*tile_nr*kcc + 1:jr*tile_nr*kcc), acc)
                end do
                !$omp end do
                !$omp end parallel
                ! the rows below, from the packed solution
                do ic = pc + kcc, mt, mc_size(tile_mr)
                    mcc = min(mc_size(tile_mr), mt - ic + 1)
                    npr = (mcc + tile_mr - 1)/tile_mr
                    nth = thread_count(2.0_dp*ncc*kcc*real(mcc, dp))
                    !$omp parallel num_threads(nth) if(nth > 1) private(jr, ir, acc)
                    !$omp do schedule(static)
                    do ir = 1, npr
                        call pack_lower(ic + (ir - 1)*tile_mr, min(tile_mr, mcc - (ir - 1)*tile_mr), pc, kcc, &
                            ar((ir - 1)*tile_mr*kcc + 1:ir*tile_mr*kcc))
                    end do
                    !$omp end do
                    !$omp do collapse(2) schedule(static)
                    do jr = 1, npb
                        do ir = 1, npr
                            call update_tile(ic + (ir - 1)*tile_mr, min(tile_mr, mcc - (ir - 1)*tile_mr), &
                                jc + (jr - 1)*tile_nr, min(tile_nr, ncc - (jr - 1)*tile_nr), kcc, &
                                ar((ir - 1)*tile_mr*kcc + 1:ir*tile_mr*kcc), &
                                bp((jr - 1)*tile_nr*kcc + 1:jr*tile_nr*kcc), acc)
                        end do
                    end do
                    !$omp end do
                    !$omp end parallel
                end do
            end do
        end do

    contains

        ! row of A (left) or column (right) holding canonical row i
        pure integer function phys(i)
            integer, intent(in) :: i
            phys = merge(i, mt + 1 - i, fwd)
        end function phys

        pure real(dp) function lval(i, j)
            integer, intent(in) :: i, j
            if (swap) then
                lval = a(phys(j), phys(i))
            else
                lval = a(phys(i), phys(j))
            end if
        end function lval

        ! canonical rows pc..pc+kc-1 of B', columns j0..j0+cols-1, as bp(nr, kc)
        subroutine pack_rhs(j0, cols, pc, kc, p)
            integer, intent(in) :: j0, cols, pc, kc
            real(dp), intent(out) :: p(tile_nr, kc)
            integer :: r, jj
            if (cols < tile_nr) p(cols + 1:, :) = 0
            if (left) then
                ! down the columns of B, which are contiguous
                do jj = 1, cols
                    if (fwd) then
                        p(jj, 1:kc) = b(pc:pc + kc - 1, j0 + jj - 1)
                    else
                        do r = 1, kc
                            p(jj, r) = b(mt + 2 - pc - r, j0 + jj - 1)
                        end do
                    end if
                end do
            else
                do r = 1, kc
                    p(1:cols, r) = b(j0:j0 + cols - 1, phys(pc + r - 1))
                end do
            end if
        end subroutine pack_rhs

        ! micro-panel at offset r0 of the diagonal block at pc: L(rows, pc..)
        ! left of the diagonal as ap(mr, r0), then the mr x mr triangle with
        ! the inverted diagonal (identity padding)
        subroutine pack_diag(pc, r0, rows, p)
            integer, intent(in) :: pc, r0, rows
            real(dp), intent(out) :: p(tile_mr, r0 + tile_mr)
            integer :: ii, c
            if (rows < tile_mr) p(rows + 1:, 1:r0) = 0
            p(:, r0 + 1:) = 0
            if (fwd .and. .not. swap) then
                do c = 1, r0 + rows
                    p(max(1, c - r0):rows, c) = a(pc + r0 + max(1, c - r0) - 1:pc + r0 + rows - 1, pc + c - 1)
                end do
            else if (fwd) then
                do ii = 1, rows
                    p(ii, 1:r0 + ii) = a(pc:pc + r0 + ii - 1, pc + r0 + ii - 1)
                end do
            else
                do c = 1, r0 + rows
                    do ii = max(1, c - r0), rows
                        p(ii, c) = lval(pc + r0 + ii - 1, pc + c - 1)
                    end do
                end do
            end if
            do ii = 1, tile_mr
                if (ii > rows .or. unit) then
                    p(ii, r0 + ii) = 1
                else
                    p(ii, r0 + ii) = 1/p(ii, r0 + ii)
                end if
            end do
        end subroutine pack_diag

        ! -L(rows i0.., columns pc..pc+kc-1) as ap(mr, kc)
        subroutine pack_lower(i0, rows, pc, kc, p)
            integer, intent(in) :: i0, rows, pc, kc
            real(dp), intent(out) :: p(tile_mr, kc)
            integer :: ii, c
            if (rows < tile_mr) p(rows + 1:, :) = 0
            if (fwd .and. .not. swap) then
                do c = 1, kc
                    p(1:rows, c) = -a(i0:i0 + rows - 1, pc + c - 1)
                end do
            else
                do c = 1, kc
                    do ii = 1, rows
                        p(ii, c) = -lval(i0 + ii - 1, pc + c - 1)
                    end do
                end do
            end if
        end subroutine pack_lower

        ! the diagonal block against one NR-wide panel, micro-panel by micro-panel
        subroutine solve_panel(j0, cols, pc, kc, np, p, acc)
            integer, intent(in) :: j0, cols, pc, kc, np
            real(dp), intent(inout) :: p(tile_nr, kc)
            real(dp), intent(inout) :: acc(MAX_TILE)
            integer :: ir, r0, rows, q, ii, jj
            do ir = 1, np
                r0 = (ir - 1)*tile_mr
                rows = min(tile_mr, kc - r0)
                q = tile_mr*tile_mr*(ir - 1)*ir/2
                if (associated(tile_trsm) .and. rows == tile_mr .and. cols == tile_nr) then
                    ! the generated kernel: product, substitution and both stores
                    if (left .and. fwd) then
                        call tile_trsm(int(r0, c_int64_t), ad(q + 1:), p, b(pc + r0, j0), &
                            int(ldb, c_int64_t), ad(q + tile_mr*r0 + 1:), p(:, r0 + 1:))
                    else
                        do jj = 1, tile_nr
                            do ii = 1, tile_mr
                                acc(ii + tile_mr*(jj - 1)) = p(jj, r0 + ii)
                            end do
                        end do
                        call tile_trsm(int(r0, c_int64_t), ad(q + 1:), p, acc, int(tile_mr, c_int64_t), &
                            ad(q + tile_mr*r0 + 1:), p(:, r0 + 1:))
                        call scatter(pc + r0, rows, j0, cols, p(:, r0 + 1:r0 + rows))
                    end if
                    cycle
                end if
                if (r0 > 0) then
                    call tile(int(r0, c_int64_t), ad(q + 1:), p, acc)
                else
                    acc(1:tile_mr*tile_nr) = 0
                end if
                call solve_tile(rows, ad(q + tile_mr*r0 + 1:q + tile_mr*(r0 + tile_mr)), acc, p(:, r0 + 1:r0 + rows))
                call scatter(pc + r0, rows, j0, cols, p(:, r0 + 1:r0 + rows))
            end do
        end subroutine solve_panel

        ! x := inv(T) (x - acc^T) on an mr x mr lower triangle T with
        ! inverted diagonal; x(nr, rows) holds the rows of the panel
        subroutine solve_tile(rows, t, acc, x)
            integer, intent(in) :: rows
            real(dp), intent(in) :: t(tile_mr, tile_mr), acc(tile_mr, tile_nr)
            real(dp), intent(inout) :: x(tile_nr, rows)
            integer, parameter :: W = 8
            real(dp) :: v(W), y(W, MAX_TILE/W)
            integer :: ii, kk
            if (tile_nr > W) then
                do ii = 1, rows
                    x(:, ii) = x(:, ii) - acc(ii, :)
                    do kk = 1, ii - 1
                        x(:, ii) = x(:, ii) - t(ii, kk)*x(:, kk)
                    end do
                    x(:, ii) = x(:, ii)*t(ii, ii)
                end do
                return
            end if
            ! rows padded to a fixed width W, so each step is whole vectors
            y(tile_nr + 1:, 1:rows) = 0
            do ii = 1, rows
                y(1:tile_nr, ii) = x(:, ii) - acc(ii, :)
            end do
            do ii = 1, rows
                v = y(:, ii)
                do kk = 1, ii - 1
                    v = v - t(ii, kk)*y(:, kk)
                end do
                y(:, ii) = v*t(ii, ii)
                x(:, ii) = y(1:tile_nr, ii)
            end do
        end subroutine solve_tile

        subroutine scatter(i0, rows, j0, cols, x)
            integer, intent(in) :: i0, rows, j0, cols
            real(dp), intent(in) :: x(tile_nr, rows)
            integer :: ii, jj
            do ii = 1, rows
                if (left) then
                    do jj = 1, cols
                        b(phys(i0 + ii - 1), j0 + jj - 1) = x(jj, ii)
                    end do
                else
                    b(j0:j0 + cols - 1, phys(i0 + ii - 1)) = x(1:cols, ii)
                end if
            end do
        end subroutine scatter

        subroutine update_tile(i0, rows, j0, cols, kc, ap, p, acc)
            integer, intent(in) :: i0, rows, j0, cols, kc
            real(dp), intent(in) :: ap(*), p(*)
            real(dp), intent(inout) :: acc(MAX_TILE)
            integer :: ii, jj
            if (left .and. fwd .and. rows == tile_mr .and. cols == tile_nr) then
                call tile_c(int(kc, c_int64_t), ap, p, b(i0, j0), int(ldb, c_int64_t))
                return
            end if
            call tile(int(kc, c_int64_t), ap, p, acc)
            if (left) then
                do jj = 1, cols
                    do ii = 1, rows
                        b(phys(i0 + ii - 1), j0 + jj - 1) = b(phys(i0 + ii - 1), j0 + jj - 1) &
                            + acc((jj - 1)*tile_mr + ii)
                    end do
                end do
            else
                ! B^T: the tile's columns run along B's columns
                do ii = 1, rows
                    b(j0:j0 + cols - 1, phys(i0 + ii - 1)) = b(j0:j0 + cols - 1, phys(i0 + ii - 1)) &
                        + acc(ii:ii + (cols - 1)*tile_mr:tile_mr)
                end do
            end if
        end subroutine update_tile

    end subroutine fortnum_dtrsm_packed

    ! A panel: rows i0 .. i0+rows-1 of alpha op(A), depth pc .. pc+kc-1, stored
    ! tile_mr per step of depth and zero-padded to the tile height.
    subroutine pack_a_real(ta, a, lda, alpha, i0, rows, pc, kc, ap)
        logical, intent(in) :: ta
        integer, intent(in) :: lda, i0, rows, pc, kc
        real(dp), intent(in) :: a(lda, *), alpha
        real(dp), intent(out) :: ap(tile_mr, kc)
        integer :: i, p
        if (rows < tile_mr) ap(rows + 1:, :) = 0
        if (.not. ta) then
            do p = 1, kc
                ap(1:rows, p) = alpha*a(i0:i0 + rows - 1, pc + p - 1)
            end do
        else
            do i = 1, rows
                do p = 1, kc
                    ap(i, p) = alpha*a(pc + p - 1, i0 + i - 1)
                end do
            end do
        end if
    end subroutine pack_a_real

    ! B panel: columns j0 .. j0+cols-1 of op(B), depth pc .. pc+kc-1, stored
    ! tile_nr per step of depth and zero-padded to the tile width.
    subroutine pack_b_real(tb, b, ldb, pc, kc, j0, cols, bp)
        logical, intent(in) :: tb
        integer, intent(in) :: ldb, pc, kc, j0, cols
        real(dp), intent(in) :: b(ldb, *)
        real(dp), intent(out) :: bp(tile_nr, kc)
        integer :: j, p
        if (cols < tile_nr) bp(cols + 1:, :) = 0
        if (.not. tb) then
            ! a transpose in 2 x 2 blocks, which compilers turn into paired
            ! loads and zips (one panel 2048 x 16: 3.5 times as fast on M3,
            ! 2.1 on Zen 3, as one column at a time)
            do p = 1, kc - 1, 2
                do j = 1, cols - 1, 2
                    bp(j, p) = b(pc + p - 1, j0 + j - 1)
                    bp(j + 1, p) = b(pc + p - 1, j0 + j)
                    bp(j, p + 1) = b(pc + p, j0 + j - 1)
                    bp(j + 1, p + 1) = b(pc + p, j0 + j)
                end do
                if (mod(cols, 2) == 1) then
                    bp(cols, p) = b(pc + p - 1, j0 + cols - 1)
                    bp(cols, p + 1) = b(pc + p, j0 + cols - 1)
                end if
            end do
            if (mod(kc, 2) == 1) bp(1:cols, kc) = b(pc + kc - 1, j0:j0 + cols - 1)
        else
            do p = 1, kc
                bp(1:cols, p) = b(j0:j0 + cols - 1, pc + p - 1)
            end do
        end if
    end subroutine pack_b_real

    ! Complex A panel of height half = tile_mr/2: real parts in rows 1..half,
    ! imaginary parts in rows half+1..tile_mr.
    subroutine pack_a_complex(ta, ca, a, lda, alpha, i0, rows, half, pc, kc, ap)
        logical, intent(in) :: ta, ca
        integer, intent(in) :: lda, i0, rows, half, pc, kc
        complex(dp), intent(in) :: a(lda, *), alpha
        real(dp), intent(out) :: ap(tile_mr, kc)
        complex(dp) :: v
        integer :: i, p
        ap = 0
        if (.not. ta) then
            do p = 1, kc
                do i = 1, rows
                    v = alpha*a(i0 + i - 1, pc + p - 1)
                    ap(i, p) = real(v, dp)
                    ap(half + i, p) = aimag(v)
                end do
            end do
        else
            do i = 1, rows
                do p = 1, kc
                    if (ca) then
                        v = alpha*conjg(a(pc + p - 1, i0 + i - 1))
                    else
                        v = alpha*a(pc + p - 1, i0 + i - 1)
                    end if
                    ap(i, p) = real(v, dp)
                    ap(half + i, p) = aimag(v)
                end do
            end do
        end if
    end subroutine pack_a_complex

    ! Complex B panel of width half = tile_nr/2: real parts in columns
    ! 1..half, imaginary parts in half+1..tile_nr.
    subroutine pack_b_complex(tb, cb, b, ldb, pc, kc, j0, cols, half, bp)
        logical, intent(in) :: tb, cb
        integer, intent(in) :: ldb, pc, kc, j0, cols, half
        complex(dp), intent(in) :: b(ldb, *)
        real(dp), intent(out) :: bp(tile_nr, kc)
        complex(dp) :: v
        integer :: j, p
        bp = 0
        do j = 1, cols
            do p = 1, kc
                if (.not. tb) then
                    v = b(pc + p - 1, j0 + j - 1)
                else if (cb) then
                    v = conjg(b(j0 + j - 1, pc + p - 1))
                else
                    v = b(j0 + j - 1, pc + p - 1)
                end if
                bp(j, p) = real(v, dp)
                bp(half + j, p) = aimag(v)
            end do
        end do
    end subroutine pack_b_complex

    real(dp) function wall()
        integer(8) :: count, rate
        call system_clock(count, rate)
        wall = real(count, dp)/real(rate, dp)
    end function wall

    ! Adds the time since w0 to t and restarts the clock.
    subroutine lap(t, w0)
        real(dp), intent(inout) :: t, w0
        real(dp) :: w1
        w1 = wall()
        t = t + (w1 - w0)
        w0 = w1
    end subroutine lap

    subroutine report_profile(name, m, n, k, prof)
        character(*), intent(in) :: name
        integer, intent(in) :: m, n, k
        real(dp), intent(in) :: prof(:, 0:)
        character(len=12), parameter :: phase(4) = [character(len=12) :: 'pack B', 'pack A', 'kernel+C', 'barriers']
        integer :: i
        print '(a, " m=", i0, " n=", i0, " k=", i0, ", ", i0, " threads; seconds per thread (min mean max):")', &
            name, m, n, k, size(prof, 2)
        do i = 1, 4
            print '(2x, a12, 3f10.5)', phase(i), minval(prof(i, :)), sum(prof(i, :))/size(prof, 2), &
                maxval(prof(i, :))
        end do
        print '(2x, a, *(f6.1))', 'kernel ms by thread:', 1000*prof(3, :)
    end subroutine report_profile

    pure integer function mc_size(step) result(size)
        integer, intent(in) :: step
        size = step*max(1, mc_rows/step)
    end function mc_size

    pure integer function nc_size(step) result(size)
        integer, intent(in) :: step
        size = step*max(1, nc_cols/step)
    end function nc_size

    !> Thread grid pr x pq = nt over mt tile rows and nb tile columns with the
    !> fewest tiles on the critical path, preferring fewer row groups (each
    !> column group packs its rows of A again).
    pure subroutine thread_grid(nt, mt, nb, pr, pq)
        integer, intent(in) :: nt, mt, nb
        integer, intent(out) :: pr, pq
        integer :: r, q, cost, best
        pr = 1
        pq = nt
        best = huge(1)
        do r = 1, nt
            if (mod(nt, r) /= 0) cycle
            q = nt/r
            cost = ((mt + r - 1)/r)*((nb + q - 1)/q)
            if (cost < best) then
                best = cost
                pr = r
                pq = q
            end if
        end do
    end subroutine thread_grid

    !> Packed B: two KC x NC blocks are in flight (the factor 2 is the caller's).
    pure integer function b_block(n, kd, nstep) result(size)
        integer, intent(in) :: n, kd, nstep
        size = tile_nr*kd*((min(n, nc_size(nstep)) + nstep - 1)/nstep)
    end function b_block

    !> Packed A of one thread: an MC x KC block, or fewer rows if the product has fewer.
    pure integer function a_block(m, kd, mstep) result(size)
        integer, intent(in) :: m, kd, mstep
        size = tile_mr*kd*min(mc_size(mstep)/mstep, (m + mstep - 1)/mstep)
    end function a_block

    !> Threads a GEMM outside a parallel region may use; 1 inside one.
    integer function thread_limit() result(count)
        count = 1
        !$ if (.not. omp_in_parallel()) count = omp_get_max_threads()
    end function thread_limit

    integer function thread_count(flops) result(count)
        real(dp), intent(in) :: flops
        count = 1
        !$ if (.not. omp_in_parallel()) &
        !$     count = max(1, min(omp_get_max_threads(), int(min(flops/THREAD_FLOPS, 1.0e6_dp))))
    end function thread_count

    ! True for a transposed operand. Reference BLAS accepts N, T and C in
    ! either case; for real data C means T.
    logical function transposed(flag, routine, position, is_complex, conjugate) result(t)
        character, intent(in) :: flag
        character(*), intent(in) :: routine
        integer, intent(in) :: position
        logical, intent(in) :: is_complex
        logical, intent(out), optional :: conjugate
        logical :: conj
        select case (flag)
        case ('N', 'n')
            t = .false.
            conj = .false.
        case ('T', 't')
            t = .true.
            conj = .false.
        case ('C', 'c')
            t = .true.
            conj = is_complex
        case default
            call argument_error(routine, position)
        end select
        if (present(conjugate)) conjugate = conj
    end function transposed

    subroutine check_sizes(routine, ta, tb, m, n, k, lda, ldb, ldc)
        character(*), intent(in) :: routine
        logical, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        if (m < 0) call argument_error(routine, 3)
        if (n < 0) call argument_error(routine, 4)
        if (k < 0) call argument_error(routine, 5)
        if (lda < max(1, merge(k, m, ta))) call argument_error(routine, 8)
        if (ldb < max(1, merge(n, k, tb))) call argument_error(routine, 10)
        if (ldc < max(1, m)) call argument_error(routine, 13)
    end subroutine check_sizes

    subroutine argument_error(routine, position)
        character(*), intent(in) :: routine
        integer, intent(in) :: position
        character(len=16) :: number
        write (number, '(i0)') position
        error stop routine//': argument '//trim(number)//' has an illegal value'
    end subroutine argument_error

end module fortnum_blas
