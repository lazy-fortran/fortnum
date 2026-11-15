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

    interface
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
        logical :: ta, tb
        integer :: jc, pc, nc, kc, mstep, nstep, npb, jr, ir, j
        integer :: nthreads, mcb, mt, i0, i1, me, chunk, nchunk, ct, epoch, bsize, boff
        real(dp) :: tp(4), w0, w1
        real(dp), allocatable :: prof(:, :)
        real(dp), allocatable :: ap(:), bp(:)

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
        end do
        if (alpha == 0 .or. k == 0) return

        call select_tile()
        ! thin products: A is read unpacked, only B is packed
        if (.not. ta .and. n <= THIN_N .and. m >= tile_mr) then
            call thin_gemm()
            return
        end if
        ! products with few rows: B is read unpacked, only A is packed
        if (.not. tb .and. m <= THIN_N .and. n >= tile_nr .and. associated(tile_bc)) then
            call wide_gemm()
            return
        end if
        mstep = tile_mr
        nstep = tile_nr
        nthreads = thread_count(2.0_dp*m*real(n, dp)*k)
        mcb = mc_size(mstep)
        ! buffers for the depth and rows the product has, not the block targets
        bsize = tile_nr*min(kc_depth, k)*((min(n, nc_size(nstep)) + nstep - 1)/nstep)
        ! two B buffers: threads may pack the next block while others finish this one
        allocate (bp(2*bsize))
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
        allocate (ap(tile_mr*min(kc_depth, k)*min(mcb/mstep, (m + mstep - 1)/mstep)))
        tp = 0
        me = 0
        !$ me = omp_get_thread_num()
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
        deallocate (ap)
        if (profile) prof(:, me) = tp
        !$omp end parallel
        if (profile) call report_profile('dgemm', m, n, k, prof)

    contains

        ! C += alpha A op(B) for n <= THIN_N without packing A: B (alpha folded
        ! in) is packed per depth block, full row tiles of A are read in place
        ! by the strided kernels, the last partial row tile is packed.
        subroutine thin_gemm()
            real(dp), allocatable :: bt(:), ap1(:)
            real(dp) :: acc(MAX_TILE)
            integer :: pcc, kcc, npbt, jrr, irr, mfull, i0, j0, nrr, jj, nth
            npbt = (n + tile_nr - 1)/tile_nr
            mfull = (m/tile_mr)*tile_mr
            allocate (bt(tile_nr*min(kc_depth, k)*npbt), ap1(tile_mr*min(kc_depth, k)))
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
            real(dp), allocatable :: at(:), bp1(:)
            real(dp) :: acc(MAX_TILE)
            integer :: pcc, kcc, npat, irr, jrr, nfull, i0, j0, mrr, ii, jj, nth
            npat = (m + tile_mr - 1)/tile_mr
            nfull = (n/tile_nr)*tile_nr
            allocate (at(tile_mr*min(kc_depth, k)*npat), bp1(tile_nr*min(kc_depth, k)))
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
            integer :: i0, j0, mrr, nrr, jj
            i0 = 1 + (ir - 1)*mstep
            j0 = jc + (jr - 1)*nstep
            mrr = min(mstep, m - i0 + 1)
            nrr = min(nstep, nc - (jr - 1)*nstep)
            if (mrr == mstep .and. nrr == nstep) then
                call tile_c(int(kc, c_int64_t), apanel, bpanel, c(i0, j0), &
                    int(ldc, c_int64_t))
                return
            end if
            call tile(int(kc, c_int64_t), apanel, bpanel, acc)
            do jj = 1, nrr
                c(i0:i0 + mrr - 1, j0 + jj - 1) = c(i0:i0 + mrr - 1, j0 + jj - 1) + &
                    acc((jj - 1)*tile_mr + 1:(jj - 1)*tile_mr + mrr)
            end do
        end subroutine update_real

    end subroutine fortnum_dgemm

    !> C := alpha op(A) op(B) + beta C, op(X) = X, X^T or X^H.
    subroutine fortnum_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        complex(dp), intent(in) :: alpha, beta
        complex(dp), intent(in) :: a(lda, *), b(ldb, *)
        complex(dp), intent(inout) :: c(ldc, *)
        logical :: ta, tb, ca, cb
        integer :: jc, pc, nc, kc, mstep, nstep, npb, jr, ir, j
        integer :: nthreads, mcb, mt, i0, i1, me, chunk, nchunk, ct, epoch, bsize, boff
        real(dp), allocatable :: ap(:), bp(:)

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
        nthreads = thread_count(8.0_dp*m*real(n, dp)*k)
        mcb = mc_size(mstep)
        ! buffers for the depth and rows the product has, not the block targets
        bsize = tile_nr*min(kc_depth, k)*((min(n, nc_size(nstep)) + nstep - 1)/nstep)
        ! two B buffers: threads may pack the next block while others finish this one
        allocate (bp(2*bsize))
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
        allocate (ap(tile_mr*min(kc_depth, k)*min(mcb/mstep, (m + mstep - 1)/mstep)))
        me = 0
        !$ me = omp_get_thread_num()
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
        deallocate (ap)
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

    end subroutine fortnum_zgemm

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
