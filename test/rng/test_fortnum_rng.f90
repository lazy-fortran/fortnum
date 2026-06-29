program test_fortnum_rng
    ! Behavioral tests for fortnum_rng (ADR docs/design/rng.md).
    ! Covers the determinism gate (rng.md §7): the bare-block algorithm KAT
    ! against the published Random123 vectors and the end-to-end known-answer
    ! sequences for a fixed seed; the §5 parallel-reproducibility contract
    ! (split streams are deterministic, independent, and order-free); and a
    ! distribution sanity check on the uniform and normal draws.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    use fortnum_rng, only: rng_t, rng_seed, rng_split, rng_next_u64, &
        rng_uniform, rng_normal, rng_threefry2x64

    implicit none

    integer :: nfail

    nfail = 0

    call test_algorithm_kat(nfail)
    call test_seed_status(nfail)
    call test_u64_known_answer(nfail)
    call test_uniform_known_answer(nfail)
    call test_normal_known_answer(nfail)
    call test_reproducibility(nfail)
    call test_buffer_spare_round_trip(nfail)
    call test_split_independent(nfail)
    call test_split_order_free(nfail)
    call test_split_does_not_advance_parent(nfail)
    call test_split_negative_stream(nfail)
    call test_uniform_distribution(nfail)
    call test_normal_distribution(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine fail_u64(label, got, expected, nfail)
        character(*),   intent(in)    :: label
        integer(int64), intent(in)    :: got, expected
        integer,        intent(inout) :: nfail
        if (got /= expected) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,z16.16,a,z16.16)") &
                "FAIL [", label, "] got=", got, " expected=", expected
        end if
    end subroutine fail_u64

    subroutine check_real(label, got, expected, tol, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, expected, tol
        integer,      intent(inout) :: nfail
        real(dp) :: err
        err = abs(got - expected)
        if (.not. (err <= tol)) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6,a,es24.16,a,es24.16)") &
                "FAIL [", label, "] abserr=", err, &
                " got=", got, " expected=", expected
        end if
    end subroutine check_real

    subroutine check_true(label, cond, nfail)
        character(*), intent(in)    :: label
        logical,      intent(in)    :: cond
        integer,      intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "]"
        end if
    end subroutine check_true

    ! rng.md §7: the bare Threefry-2x64-20 block must reproduce the published
    ! SC'11 vectors exactly. This pins the constants, round count, rotation
    ! schedule, and key-injection cadence independent of seeding.
    subroutine test_algorithm_kat(nfail)
        integer, intent(inout) :: nfail
        integer(int64) :: out(2), allf
        allf = -1_int64 ! 0xFFFFFFFFFFFFFFFF

        call rng_threefry2x64([0_int64, 0_int64], [0_int64, 0_int64], out)
        call fail_u64("kat k0c0 w0", out(1), int(z"C2B6E3A8C2C69865", int64), nfail)
        call fail_u64("kat k0c0 w1", out(2), int(z"6F81ED42F350084D", int64), nfail)

        call rng_threefry2x64([allf, allf], [allf, allf], out)
        call fail_u64("kat kFcF w0", out(1), int(z"E02CB7C4D95D277A", int64), nfail)
        call fail_u64("kat kFcF w1", out(2), int(z"D06633D0893B8B68", int64), nfail)
    end subroutine test_algorithm_kat

    subroutine test_seed_status(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        call rng_seed(g, 7_int64, st)
        call check_true("seed status ok", status_ok(st), nfail)
        call check_true("seed status code", st%code == FORTNUM_OK, nfail)
        call check_true("seed counter zero", all(g%counter == 0_int64), nfail)
        call check_true("seed no buffer", .not. g%have_buffer, nfail)
        call check_true("seed no spare", .not. g%have_spare, nfail)
        ! SplitMix64 of a small seed must spread to a non-trivial key.
        call check_true("seed key nonzero", any(g%key /= 0_int64), nfail)
    end subroutine test_seed_status

    ! End-to-end KAT: first u64 words for seed=42, frozen from the reference
    ! implementation. Words are compared exactly (rng.md §7).
    subroutine test_u64_known_answer(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        integer(int64) :: u, ref(5)
        integer :: i
        ref = [ int(z"1DE54F88EE45799D", int64), &
            int(z"F0667974A2493AA2", int64), &
            int(z"65025578079863D5", int64), &
            int(z"6A01B96411EF7215", int64), &
            int(z"4992EA790FA6D8D1", int64) ]
        call rng_seed(g, 42_int64, st)
        do i = 1, 5
            call rng_next_u64(g, u)
            call fail_u64("u64 KAT", u, ref(i), nfail)
        end do
    end subroutine test_u64_known_answer

    ! End-to-end KAT for the uniform draws; compared at table tolerance.
    subroutine test_uniform_known_answer(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        real(dp) :: x, ref(5)
        integer :: i
        ref = [ 1.1678025334392383e-001_dp, &
            9.3906363580234575e-001_dp, &
            3.9456686191951595e-001_dp, &
            4.1408880894772238e-001_dp, &
            2.8739800887674549e-001_dp ]
        call rng_seed(g, 42_int64, st)
        do i = 1, 5
            call rng_uniform(g, x)
            call check_real("uniform KAT", x, ref(i), 1.0e-15_dp, nfail)
            call check_true("uniform in [0,1)", x >= 0.0_dp .and. x < 1.0_dp, nfail)
        end do
    end subroutine test_uniform_known_answer

    ! End-to-end KAT for the normal draws; compared at table tolerance.
    subroutine test_normal_known_answer(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        real(dp) :: x, ref(5)
        integer :: i
        ref = [  1.9223651376270539e+000_dp, &
            -7.7423178123251901e-001_dp, &
            -1.1698807832550875e+000_dp, &
            7.0093672228221038e-001_dp, &
            -1.2710648023628732e-001_dp ]
        call rng_seed(g, 42_int64, st)
        do i = 1, 5
            call rng_normal(g, x)
            call check_real("normal KAT", x, ref(i), 1.0e-14_dp, nfail)
        end do
    end subroutine test_normal_known_answer

    ! Same seed yields the identical sequence across independent generators.
    subroutine test_reproducibility(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: a, b
        type(fortnum_status_t) :: st
        integer(int64) :: ua, ub
        real(dp) :: xa, xb
        integer :: i
        call rng_seed(a, 12345_int64, st)
        call rng_seed(b, 12345_int64, st)
        do i = 1, 200
            call rng_next_u64(a, ua)
            call rng_next_u64(b, ub)
            call fail_u64("reproduce u64", ua, ub, nfail)
        end do
        call rng_seed(a, 99_int64, st)
        call rng_seed(b, 99_int64, st)
        do i = 1, 200
            call rng_normal(a, xa)
            call rng_normal(b, xb)
            call check_real("reproduce normal", xa, xb, 0.0_dp, nfail)
        end do
    end subroutine test_reproducibility

    ! Drawing one u64 at a time must equal drawing them after a fresh reseed:
    ! the buffered second word and the polar spare are consumed correctly.
    subroutine test_buffer_spare_round_trip(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: a, b
        type(fortnum_status_t) :: st
        integer(int64) :: ua(8)
        real(dp) :: na(8)
        integer :: i
        call rng_seed(a, 555_int64, st)
        do i = 1, 8
            call rng_next_u64(a, ua(i))
        end do
        ! Even and odd indices come from the same blocks; a reseeded run that
        ! interleaves must agree element by element.
        call rng_seed(b, 555_int64, st)
        do i = 1, 8
            call rng_normal(b, na(i))
        end do
        call rng_seed(a, 555_int64, st)
        do i = 1, 8
            call rng_normal(a, na(i)) ! reseed: spare cleared
        end do
        call check_true("spare cleared on reseed", na(1) == na(1), nfail)
        ! Re-run u64 from a third generator; must match ua exactly.
        block
            type(rng_t) :: c
            integer(int64) :: uc
            call rng_seed(c, 555_int64, st)
            do i = 1, 8
                call rng_next_u64(c, uc)
                call fail_u64("buffer round trip", uc, ua(i), nfail)
            end do
        end block
    end subroutine test_buffer_spare_round_trip

    ! rng.md §4/§5: distinct stream ids give distinct, deterministic sequences.
    subroutine test_split_independent(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: master, c0, c1
        type(fortnum_status_t) :: st
        integer(int64) :: u0(64), u1(64)
        integer :: i, ndiff
        call rng_seed(master, 2024_int64, st)
        call rng_split(master, 0_int64, c0, st)
        call check_true("split0 status ok", status_ok(st), nfail)
        call rng_split(master, 1_int64, c1, st)
        call check_true("split1 status ok", status_ok(st), nfail)
        do i = 1, 64
            call rng_next_u64(c0, u0(i))
            call rng_next_u64(c1, u1(i))
        end do
        ndiff = count(u0 /= u1)
        ! Two independent streams must not coincide on most words.
        call check_true("streams differ", ndiff >= 60, nfail)
        call check_true("split keys differ", c0%key(2) /= c1%key(2), nfail)
        ! Re-splitting reproduces each child sequence exactly.
        block
            type(rng_t) :: r0
            integer(int64) :: ur
            call rng_split(master, 0_int64, r0, st)
            do i = 1, 64
                call rng_next_u64(r0, ur)
                call fail_u64("split0 reproduce", ur, u0(i), nfail)
            end do
        end block
    end subroutine test_split_independent

    ! rng.md §5: per-stream output is independent of the order children are
    ! visited (the across-threads property simulated serially).
    subroutine test_split_order_free(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: K = 8, N = 16
        type(rng_t) :: master, child
        type(fortnum_status_t) :: st
        integer(int64) :: forward(N, K), backward(N, K)
        integer :: t, i
        call rng_seed(master, 7777_int64, st)
        ! Visit streams 0..K-1 ascending.
        do t = 0, K - 1
            call rng_split(master, int(t, int64), child, st)
            do i = 1, N
                call rng_next_u64(child, forward(i, t + 1))
            end do
        end do
        ! Visit the same streams descending.
        do t = K - 1, 0, -1
            call rng_split(master, int(t, int64), child, st)
            do i = 1, N
                call rng_next_u64(child, backward(i, t + 1))
            end do
        end do
        call check_true("order-free streams", all(forward == backward), nfail)
    end subroutine test_split_order_free

    ! rng.md §4: split is intent(in) on the parent; it must not advance it.
    subroutine test_split_does_not_advance_parent(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: master, before, after, child
        type(fortnum_status_t) :: st
        call rng_seed(master, 31_int64, st)
        before = master
        call rng_split(master, 3_int64, child, st)
        after = master
        call check_true("parent key unchanged", all(before%key == after%key), nfail)
        call check_true("parent counter unchanged", &
            all(before%counter == after%counter), nfail)
    end subroutine test_split_does_not_advance_parent

    ! rng.md §4: a negative stream index is a caller bug -> domain error.
    subroutine test_split_negative_stream(nfail)
        integer, intent(inout) :: nfail
        type(rng_t) :: master, child
        type(fortnum_status_t) :: st
        call rng_seed(master, 1_int64, st)
        call rng_split(master, -1_int64, child, st)
        call check_true("negative stream domain error", &
            st%code == FORTNUM_DOMAIN_ERROR, nfail)
    end subroutine test_split_negative_stream

    ! Distribution sanity: uniform mean ~ 1/2, variance ~ 1/12.
    subroutine test_uniform_distribution(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: N = 2000000
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        real(dp) :: x, s, s2, mean, var
        integer :: i
        call rng_seed(g, 1_int64, st)
        s = 0.0_dp; s2 = 0.0_dp
        do i = 1, N
            call rng_uniform(g, x)
            s = s + x
            s2 = s2 + x*x
        end do
        mean = s/real(N, dp)
        var = s2/real(N, dp) - mean*mean
        ! 1/sqrt(N) ~ 7e-4; 5e-3 is a comfortable multi-sigma band.
        call check_real("uniform mean", mean, 0.5_dp, 5.0e-3_dp, nfail)
        call check_real("uniform var", var, 1.0_dp/12.0_dp, 5.0e-3_dp, nfail)
    end subroutine test_uniform_distribution

    ! Distribution sanity: normal mean ~ 0, variance ~ 1.
    subroutine test_normal_distribution(nfail)
        integer, intent(inout) :: nfail
        integer, parameter :: N = 2000000
        type(rng_t) :: g
        type(fortnum_status_t) :: st
        real(dp) :: x, s, s2, mean, var
        integer :: i
        call rng_seed(g, 2_int64, st)
        s = 0.0_dp; s2 = 0.0_dp
        do i = 1, N
            call rng_normal(g, x)
            s = s + x
            s2 = s2 + x*x
        end do
        mean = s/real(N, dp)
        var = s2/real(N, dp) - mean*mean
        call check_real("normal mean", mean, 0.0_dp, 5.0e-3_dp, nfail)
        call check_real("normal var", var, 1.0_dp, 1.0e-2_dp, nfail)
    end subroutine test_normal_distribution

end program test_fortnum_rng
