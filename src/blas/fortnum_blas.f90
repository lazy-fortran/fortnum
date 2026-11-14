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
    !$ use omp_lib, only: omp_in_parallel
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

    interface
        function fortnum_blas_dtile(mr, nr, isa) bind(c, name='fortnum_blas_dtile') result(kernel)
            import :: c_funptr, c_int32_t
            integer(c_int32_t), intent(out) :: mr, nr, isa
            type(c_funptr) :: kernel
        end function fortnum_blas_dtile
    end interface

    !> Depth of a packed block, and the targets for the A and B block sizes
    !> (rounded to whole tiles); per-machine values come from fortnum-bench.
    integer, parameter :: KC_DEPTH = 256, MC_ROWS = 120, NC_COLS = 4032
    !> Largest tile (rows times columns) any emitted kernel uses.
    integer, parameter :: MAX_TILE = 256
    !> Products below this many flops per thread run on one thread.
    real(dp), parameter :: THREAD_FLOPS = 4.0e6_dp

    logical, save :: tile_ready = .false.
    integer, save :: tile_mr = 0, tile_nr = 0, tile_isa = 0
    procedure(tile_kernel), pointer, save :: tile => null()

contains

    !> The tile in use: rows, columns and instruction set (1 AVX2, 2 AVX-512,
    !> 3 NEON).
    subroutine fortnum_blas_tile(mr, nr, isa)
        integer, intent(out) :: mr, nr, isa
        call select_tile()
        mr = tile_mr
        nr = tile_nr
        isa = tile_isa
    end subroutine fortnum_blas_tile

    subroutine select_tile()
        integer(c_int32_t) :: mr, nr, isa
        type(c_funptr) :: kernel
        if (tile_ready) return
        !$omp critical (fortnum_blas_tile_selection)
        if (.not. tile_ready) then
            kernel = fortnum_blas_dtile(mr, nr, isa)
            if (.not. c_associated(kernel)) &
                error stop 'fortnum_blas: no register-tile kernel for this CPU (needs AVX2+FMA, AVX-512 or NEON)'
            call c_f_procpointer(kernel, tile)
            tile_mr = mr
            tile_nr = nr
            tile_isa = isa
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
        logical :: ta, tb, threaded
        integer :: jc, pc, ic, nc, kc, mc, mstep, nstep, npa, npb, jr, j
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
        mstep = tile_mr
        nstep = tile_nr
        threaded = use_threads(2.0_dp*m*real(n, dp)*k)
        allocate (ap(tile_mr*KC_DEPTH*((min(m, mc_size(mstep)) + mstep - 1)/mstep)))
        allocate (bp(tile_nr*KC_DEPTH*((min(n, nc_size(nstep)) + nstep - 1)/nstep)))
        do jc = 1, n, nc_size(nstep)
            nc = min(nc_size(nstep), n - jc + 1)
            npb = (nc + nstep - 1)/nstep
            do pc = 1, k, KC_DEPTH
                kc = min(KC_DEPTH, k - pc + 1)
                !$omp parallel do if(threaded) schedule(static)
                do jr = 1, npb
                    call pack_b_real(tb, b, ldb, pc, kc, jc + (jr - 1)*nstep, &
                        min(nstep, nc - (jr - 1)*nstep), bp((jr - 1)*tile_nr*kc + 1:))
                end do
                !$omp end parallel do
                do ic = 1, m, mc_size(mstep)
                    mc = min(mc_size(mstep), m - ic + 1)
                    npa = (mc + mstep - 1)/mstep
                    !$omp parallel if(threaded)
                    !$omp do schedule(static)
                    do jr = 1, npa
                        call pack_a_real(ta, a, lda, alpha, ic + (jr - 1)*mstep, &
                            min(mstep, mc - (jr - 1)*mstep), pc, kc, ap((jr - 1)*tile_mr*kc + 1:))
                    end do
                    !$omp end do
                    !$omp do schedule(dynamic, 1)
                    do jr = 1, npb
                        call macro_real(jr)
                    end do
                    !$omp end do
                    !$omp end parallel
                end do
            end do
        end do

    contains

        subroutine macro_real(jr)
            integer, intent(in) :: jr
            real(dp) :: acc(MAX_TILE)
            integer :: ir, i0, j0, mrr, nrr, jj
            j0 = jc + (jr - 1)*nstep
            nrr = min(nstep, nc - (jr - 1)*nstep)
            do ir = 1, npa
                i0 = ic + (ir - 1)*mstep
                mrr = min(mstep, mc - (ir - 1)*mstep)
                call tile(int(kc, c_int64_t), ap((ir - 1)*tile_mr*kc + 1:), &
                    bp((jr - 1)*tile_nr*kc + 1:), acc)
                do jj = 1, nrr
                    c(i0:i0 + mrr - 1, j0 + jj - 1) = c(i0:i0 + mrr - 1, j0 + jj - 1) + &
                        acc((jj - 1)*tile_mr + 1:(jj - 1)*tile_mr + mrr)
                end do
            end do
        end subroutine macro_real

    end subroutine fortnum_dgemm

    !> C := alpha op(A) op(B) + beta C, op(X) = X, X^T or X^H.
    subroutine fortnum_zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: transa, transb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        complex(dp), intent(in) :: alpha, beta
        complex(dp), intent(in) :: a(lda, *), b(ldb, *)
        complex(dp), intent(inout) :: c(ldc, *)
        logical :: ta, tb, ca, cb, threaded
        integer :: jc, pc, ic, nc, kc, mc, mstep, nstep, npa, npb, jr, j
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
        threaded = use_threads(8.0_dp*m*real(n, dp)*k)
        allocate (ap(tile_mr*KC_DEPTH*((min(m, mc_size(mstep)) + mstep - 1)/mstep)))
        allocate (bp(tile_nr*KC_DEPTH*((min(n, nc_size(nstep)) + nstep - 1)/nstep)))
        do jc = 1, n, nc_size(nstep)
            nc = min(nc_size(nstep), n - jc + 1)
            npb = (nc + nstep - 1)/nstep
            do pc = 1, k, KC_DEPTH
                kc = min(KC_DEPTH, k - pc + 1)
                !$omp parallel do if(threaded) schedule(static)
                do jr = 1, npb
                    call pack_b_complex(tb, cb, b, ldb, pc, kc, jc + (jr - 1)*nstep, &
                        min(nstep, nc - (jr - 1)*nstep), nstep, bp((jr - 1)*tile_nr*kc + 1:))
                end do
                !$omp end parallel do
                do ic = 1, m, mc_size(mstep)
                    mc = min(mc_size(mstep), m - ic + 1)
                    npa = (mc + mstep - 1)/mstep
                    !$omp parallel if(threaded)
                    !$omp do schedule(static)
                    do jr = 1, npa
                        call pack_a_complex(ta, ca, a, lda, alpha, ic + (jr - 1)*mstep, &
                            min(mstep, mc - (jr - 1)*mstep), mstep, pc, kc, ap((jr - 1)*tile_mr*kc + 1:))
                    end do
                    !$omp end do
                    !$omp do schedule(dynamic, 1)
                    do jr = 1, npb
                        call macro_complex(jr)
                    end do
                    !$omp end do
                    !$omp end parallel
                end do
            end do
        end do

    contains

        subroutine macro_complex(jr)
            integer, intent(in) :: jr
            real(dp) :: acc(MAX_TILE)
            integer :: ir, i0, j0, mrr, nrr, ii, jj, rr, ri, ir_, ii_
            j0 = jc + (jr - 1)*nstep
            nrr = min(nstep, nc - (jr - 1)*nstep)
            do ir = 1, npa
                i0 = ic + (ir - 1)*mstep
                mrr = min(mstep, mc - (ir - 1)*mstep)
                call tile(int(kc, c_int64_t), ap((ir - 1)*tile_mr*kc + 1:), &
                    bp((jr - 1)*tile_nr*kc + 1:), acc)
                do jj = 1, nrr
                    ! columns jj (real part of B) and nstep+jj (imaginary part)
                    rr = (jj - 1)*tile_mr
                    ri = (nstep + jj - 1)*tile_mr
                    do ii = 1, mrr
                        ir_ = ii
                        ii_ = mstep + ii
                        c(i0 + ii - 1, j0 + jj - 1) = c(i0 + ii - 1, j0 + jj - 1) + &
                            cmplx(acc(rr + ir_) - acc(ri + ii_), acc(ri + ir_) + acc(rr + ii_), dp)
                    end do
                end do
            end do
        end subroutine macro_complex

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
            do j = 1, cols
                do p = 1, kc
                    bp(j, p) = b(pc + p - 1, j0 + j - 1)
                end do
            end do
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

    pure integer function mc_size(step) result(size)
        integer, intent(in) :: step
        size = step*max(1, MC_ROWS/step)
    end function mc_size

    pure integer function nc_size(step) result(size)
        integer, intent(in) :: step
        size = step*max(1, NC_COLS/step)
    end function nc_size

    logical function use_threads(flops) result(threaded)
        real(dp), intent(in) :: flops
        threaded = .false.
        !$ threaded = flops > THREAD_FLOPS .and. .not. omp_in_parallel()
    end function use_threads

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
