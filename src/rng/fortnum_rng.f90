module fortnum_rng
    ! Counter-based pseudorandom generator with caller-owned explicit state
    ! (ADR docs/design/rng.md). The state lives in rng_t and is passed by
    ! reference; there is no module-level mutable state, no global pool, and no
    ! procedure pointer, so two threads that own separate rng_t never race and
    ! the output is a pure function of (key, counter) plus the fixed algorithm.
    !
    ! Algorithm (frozen for issue #20): Threefry-2x64 with 20 rounds (Salmon,
    ! Moraes, Dror, Shaw, "Parallel Random Numbers: As Easy as 1, 2, 3", SC'11),
    ! seeded by SplitMix64 (Steele, Lea, Flood, OOPSLA'14). Threefry needs only
    ! add, xor, and rotate, so the kernel is branch-free and identical across
    ! compilers; that is what the determinism gate checks.
    !
    ! Derivative policy: primal_only (ad.md §4). A pseudorandom draw is not a
    ! differentiable function of its seed, so no derivative entry point exists.
    ! The seed, the stream index, the key, and the counter are inactive
    ! arguments (ad.md §3): they select behavior, they do not carry a derivative.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    private

    public :: rng_t
    public :: rng_seed
    public :: rng_split
    public :: rng_next_u64
    public :: rng_uniform
    public :: rng_normal

    ! Bare Threefry-2x64-20 block, exposed for the determinism gate's algorithm
    ! KAT against the published Random123 vectors (rng.md §7).
    public :: rng_threefry2x64

    ! Caller-owned generator state. The reproducible state is key and counter;
    ! buffer/spare are derived caches that seed and split clear (rng.md §2).
    type :: rng_t
        integer(int64) :: key(2)       = 0_int64
        integer(int64) :: counter(2)   = 0_int64
        integer(int64) :: buffer       = 0_int64
        logical        :: have_buffer  = .false.
        real(dp)       :: spare_normal = 0.0_dp
        logical        :: have_spare   = .false.
    end type rng_t

    ! Threefry-2x64 parameters. Round count is fixed at 20 (Random123 default).
    integer, parameter :: THREEFRY_ROUNDS = 20

    ! Key-schedule parity constant 2^64/phi (SC'11 / Skein).
    integer(int64), parameter :: KS_PARITY = int(z"1BD11BDAA9FC1A22", int64)

    ! Rotation schedule for the 2x64 word size; the cycle has length 8.
    integer, parameter :: ROT(8) = [16, 42, 12, 31, 16, 32, 24, 21]

    ! SplitMix64 constants (OOPSLA'14): golden-ratio increment and the two
    ! finalizer multipliers.
    integer(int64), parameter :: SM_GAMMA = int(z"9E3779B97F4A7C15", int64)
    integer(int64), parameter :: SM_MIX1  = int(z"BF58476D1CE4E5B9", int64)
    integer(int64), parameter :: SM_MIX2  = int(z"94D049BB133111EB", int64)

    ! 2**-53, the scale that maps a 53-bit mantissa to [0, 1).
    real(dp), parameter :: INV_2_53 = 2.0_dp**(-53)

    ! Floor for u1 in Box-Muller so the [0,1) lower endpoint cannot give +inf.
    real(dp), parameter :: U1_FLOOR = INV_2_53

    real(dp), parameter :: TWO_PI = 8.0_dp*atan(1.0_dp)

contains

    ! Turn one 64-bit seed into a 128-bit key by running SplitMix64 twice.
    ! SplitMix64 spreads a low-entropy seed across all key bits, so nearby
    ! seeds give unrelated streams (rng.md §3). The counter and caches reset.
    pure subroutine rng_seed(g, seed, status)
        type(rng_t),            intent(out) :: g
        integer(int64),         intent(in)  :: seed
        type(fortnum_status_t), intent(out) :: status

        integer(int64) :: s

        s = seed
        call splitmix64_next(s, g%key(1))
        call splitmix64_next(s, g%key(2))
        g%counter = 0_int64
        g%buffer = 0_int64
        g%have_buffer = .false.
        g%spare_normal = 0.0_dp
        g%have_spare = .false.

        call status_set(status, FORTNUM_OK, "")
    end subroutine rng_seed

    ! Derive an independent substream. The child keeps the parent key high word
    ! and xors the stream index into the low word, then zeroes the counter; two
    ! distinct indices give two distinct keys and non-overlapping sequences
    ! (rng.md §4). Splitting does not advance the parent (intent(in)).
    pure subroutine rng_split(parent, stream, child, status)
        type(rng_t),            intent(in)  :: parent
        integer(int64),         intent(in)  :: stream
        type(rng_t),            intent(out) :: child
        type(fortnum_status_t), intent(out) :: status

        if (stream < 0_int64) then
            child = parent
            child%counter = 0_int64
            child%buffer = 0_int64
            child%have_buffer = .false.
            child%spare_normal = 0.0_dp
            child%have_spare = .false.
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "rng_split: negative stream index")
            return
        end if

        child%key(2) = ieor(parent%key(2), stream)
        child%key(1) = parent%key(1)
        child%counter = 0_int64
        child%buffer = 0_int64
        child%have_buffer = .false.
        child%spare_normal = 0.0_dp
        child%have_spare = .false.

        call status_set(status, FORTNUM_OK, "")
    end subroutine rng_split

    ! Next uniform 64-bit word. One Threefry call yields two words: the first is
    ! returned, the second cached; the cached word is served on the next call
    ! with no Threefry call (rng.md §6.1). The counter increments one block,
    ! carrying low to high so a stream spans the full 128-bit counter space.
    pure subroutine rng_next_u64(g, value)
        type(rng_t),    intent(inout) :: g
        integer(int64), intent(out)   :: value

        integer(int64) :: out(2)

        if (g%have_buffer) then
            value = g%buffer
            g%have_buffer = .false.
            return
        end if

        call rng_threefry2x64(g%key, g%counter, out)
        value = out(1)
        g%buffer = out(2)
        g%have_buffer = .true.

        ! Advance the 128-bit counter by one block, low word first with carry.
        g%counter(1) = g%counter(1) + 1_int64
        if (g%counter(1) == 0_int64) g%counter(2) = g%counter(2) + 1_int64
    end subroutine rng_next_u64

    ! Uniform real in [0, 1): keep the top 53 bits of a u64 and scale by 2**-53.
    ! The high bits are the best bits of the output; the value is never 1
    ! (rng.md §6.2).
    pure subroutine rng_uniform(g, value)
        type(rng_t), intent(inout) :: g
        real(dp),    intent(out)   :: value

        integer(int64) :: u

        call rng_next_u64(g, u)
        ! Logical (unsigned) shift by 11 leaves the top 53 bits in [0, 2**53).
        value = real(ishft(u, -11), dp)*INV_2_53
    end subroutine rng_uniform

    ! Standard normal by Box-Muller: draw u1, u2; r = sqrt(-2 ln u1),
    ! theta = 2 pi u2; return r cos theta, cache r sin theta (rng.md §6.3). u1 is
    ! floored at a tiny positive value so the [0,1) lower endpoint never gives
    ! +inf. No table and no rejection branch, so the determinism gate extends to
    ! normals.
    pure subroutine rng_normal(g, value)
        type(rng_t), intent(inout) :: g
        real(dp),    intent(out)   :: value

        real(dp) :: u1, u2, r, theta

        if (g%have_spare) then
            value = g%spare_normal
            g%have_spare = .false.
            return
        end if

        call rng_uniform(g, u1)
        call rng_uniform(g, u2)
        u1 = max(u1, U1_FLOOR)
        r = sqrt(-2.0_dp*log(u1))
        theta = TWO_PI*u2
        value = r*cos(theta)
        g%spare_normal = r*sin(theta)
        g%have_spare = .true.
    end subroutine rng_normal

    ! One SplitMix64 step: advance the state by the golden-ratio increment and
    ! return the finalized mix of the new state (OOPSLA'14).
    pure subroutine splitmix64_next(state, z)
        integer(int64), intent(inout) :: state
        integer(int64), intent(out)   :: z

        state = state + SM_GAMMA
        z = state
        z = ieor(z, ishft(z, -30))*SM_MIX1
        z = ieor(z, ishft(z, -27))*SM_MIX2
        z = ieor(z, ishft(z, -31))
    end subroutine splitmix64_next

    ! Threefry-2x64-20 block cipher: map (key, counter) to two output words.
    ! Add/xor/rotate only; key injection every four rounds (SC'11).
    pure subroutine rng_threefry2x64(key, counter, out)
        integer(int64), intent(in)  :: key(2)
        integer(int64), intent(in)  :: counter(2)
        integer(int64), intent(out) :: out(2)

        integer(int64) :: ks(3), x0, x1
        integer        :: round, inject

        ks(1) = key(1)
        ks(2) = key(2)
        ks(3) = ieor(ieor(KS_PARITY, key(1)), key(2))

        x0 = counter(1) + ks(1)
        x1 = counter(2) + ks(2)

        do round = 0, THREEFRY_ROUNDS - 1
            x0 = x0 + x1
            x1 = rotl64(x1, ROT(mod(round, 8) + 1))
            x1 = ieor(x1, x0)

            if (mod(round, 4) == 3) then
                inject = round/4 + 1
                x0 = x0 + ks(mod(inject, 3) + 1)
                x1 = x1 + ks(mod(inject + 1, 3) + 1) + int(inject, int64)
            end if
        end do

        out(1) = x0
        out(2) = x1
    end subroutine rng_threefry2x64

    ! 64-bit left rotation by n in [0, 63]. Built from logical shifts so it is
    ! defined independently of integer sign.
    pure function rotl64(x, n) result(r)
        integer(int64), intent(in) :: x
        integer,        intent(in) :: n
        integer(int64) :: r
        r = ior(ishft(x, n), ishft(x, n - 64))
    end function rotl64

end module fortnum_rng
