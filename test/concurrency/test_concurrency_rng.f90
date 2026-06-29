program test_concurrency_rng
    ! M5.1 concurrency: rng (per-thread rng_t with independent streams).
    ! This is the key reproducibility gate from docs/design/rng.md sections 5
    ! and 7: seed once on the master, split per stream index, and the sequence
    ! thread t draws must be a deterministic function of (seed, t, draw count)
    ! alone -- independent of the number of threads, the schedule, and whether
    ! the draws run serially or under OpenMP.
    !
    ! The test:
    !   1. Build the serial reference: for each stream t in [0, nstream), split
    !      from the master and draw nseq u64 / uniform / normal values.
    !   2. Draw the same ensemble under OpenMP, one rng_t per iteration split
    !      by the loop index. No thread touches another's state.
    !   3. Assert the u64 words are bit-for-bit equal and the reals are exactly
    !      equal (same arithmetic, same order, so exact equality holds).
    !   4. Assert distinct streams differ (independence is not trivially all
    !      streams equal).
    !
    ! PRIMAL concurrency only; the rng module is primal_only (rng.md derivative
    ! classification), so there is no derivative path and none arrives in M6.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_rng, only: rng_t, rng_seed, rng_split, rng_next_u64, &
        rng_uniform, rng_normal
    implicit none

    integer,        parameter :: nstream = 256
    integer,        parameter :: nseq = 64
    integer(int64), parameter :: seed = 20260613_int64
    integer :: nfail, t, k

    integer(int64) :: ser_u64(nseq, nstream), par_u64(nseq, nstream)
    real(dp) :: ser_uni(nseq, nstream), par_uni(nseq, nstream)
    real(dp) :: ser_nrm(nseq, nstream), par_nrm(nseq, nstream)

    nfail = 0

    call draw_serial()
    call draw_parallel()

    ! Exact equality: counter-based generator, identical arithmetic in both
    ! runs, so reals match to the last bit, not just to a tolerance.
    if (any(ser_u64 /= par_u64)) then
        nfail = nfail + 1
        write (error_unit, "(a)") "FAIL: rng u64 parallel != serial"
    end if
    if (any(ser_uni /= par_uni)) then
        nfail = nfail + 1
        write (error_unit, "(a)") "FAIL: rng uniform parallel != serial"
    end if
    if (any(ser_nrm /= par_nrm)) then
        nfail = nfail + 1
        write (error_unit, "(a)") "FAIL: rng normal parallel != serial"
    end if

    ! Independence: stream 0 and stream 1 must not produce the same sequence,
    ! otherwise the per-thread split is not selecting distinct streams.
    if (all(ser_u64(:, 1) == ser_u64(:, 2))) then
        nfail = nfail + 1
        write (error_unit, "(a)") "FAIL: rng streams 0 and 1 identical"
    end if

    ! Order independence: drawing the streams in reverse must match the same
    ! per-stream sequences (rng.md section 5: independent of visit order).
    call check_reverse_order()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_rng: all tests passed"

contains

    subroutine draw_one(g, col_u64, col_uni, col_nrm)
        type(rng_t),    intent(inout) :: g
        integer(int64), intent(out)   :: col_u64(nseq)
        real(dp),       intent(out)   :: col_uni(nseq)
        real(dp),       intent(out)   :: col_nrm(nseq)
        type(rng_t) :: gu, gn
        integer :: i
        ! Draw u64, then uniforms from a fresh copy, then normals from another,
        ! so the three streams are read from the same split state independently
        ! and the counts are fixed per stream.
        gu = g
        gn = g
        do i = 1, nseq
            call rng_next_u64(g, col_u64(i))
            call rng_uniform(gu, col_uni(i))
            call rng_normal(gn, col_nrm(i))
        end do
    end subroutine draw_one

    subroutine draw_serial()
        type(rng_t)            :: master, child
        type(fortnum_status_t) :: st
        call rng_seed(master, seed, st)
        if (.not. status_ok(st)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: rng_seed (serial) not ok"
            return
        end if
        do t = 1, nstream
            call rng_split(master, int(t - 1, int64), child, st)
            if (.not. status_ok(st)) then
                nfail = nfail + 1
                write (error_unit, "(a,i0)") "FAIL: rng_split serial t=", t
                return
            end if
            call draw_one(child, ser_u64(:, t), ser_uni(:, t), ser_nrm(:, t))
        end do
    end subroutine draw_serial

    subroutine draw_parallel()
        type(rng_t)            :: master, child
        type(fortnum_status_t) :: st
        call rng_seed(master, seed, st)
        ! master is shared and read-only (rng_split takes parent intent(in));
        ! child and st are private, split by the loop index, so each iteration
        ! owns an independent stream with no shared mutable state.
        !$omp parallel do default(shared) private(child, st, t) schedule(static)
        do t = 1, nstream
            call rng_split(master, int(t - 1, int64), child, st)
            call draw_one(child, par_u64(:, t), par_uni(:, t), par_nrm(:, t))
        end do
        !$omp end parallel do
    end subroutine draw_parallel

    subroutine check_reverse_order()
        type(rng_t)            :: master, child
        type(fortnum_status_t) :: st
        integer(int64) :: rev_u64(nseq, nstream)
        real(dp)       :: rev_uni(nseq, nstream), rev_nrm(nseq, nstream)
        call rng_seed(master, seed, st)
        !$omp parallel do default(shared) private(child, st, t) schedule(static)
        do t = nstream, 1, -1
            call rng_split(master, int(t - 1, int64), child, st)
            call draw_one(child, rev_u64(:, t), rev_uni(:, t), rev_nrm(:, t))
        end do
        !$omp end parallel do
        if (any(rev_u64 /= ser_u64) .or. any(rev_uni /= ser_uni) &
            .or. any(rev_nrm /= ser_nrm)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: rng reverse-order draw != serial"
        end if
    end subroutine check_reverse_order

end program test_concurrency_rng
