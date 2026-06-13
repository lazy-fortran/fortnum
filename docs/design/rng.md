# ADR: fortnum_rng API

Status: accepted (issue #18, M3.1). Normative for the implementation issues #20
(primal generator and the determinism gate) and any later module that draws
random numbers.

This ADR is subordinate to `docs/design/ad.md`. The module is `primal_only`
(ad.md §4): a pseudorandom draw is not a differentiable function of its seed, so
no derivative entry point exists and none is reserved. A caller who needs a
gradient of an estimator built on draws differentiates the estimator, not the
generator. The seed, the stream index, the key, and the counter are all inactive
arguments (ad.md §3: keys and RNG seeds select behavior, they do not carry a
derivative).

An implementer should be able to write the generator from this file without
further design decisions. The algorithm, the round count, and the constants are
fixed here so #20 reproduces this document's known-answer vectors bit for bit.

## 1. Why explicit state, and why a counter-based generator

A generator with module-level state is not thread-safe and not reproducible.
Two threads sharing one hidden stream race on the update; the draw a thread sees
depends on scheduling, so the same seed gives different output on a different
run. fortnum forbids module-level mutable state (ad.md §6, and the house rule on
explicit state), so the generator state lives in a caller-owned `rng_t` passed by
reference. Nothing is global.

Explicit state makes thread safety the caller's to arrange, but it does not by
itself give independent streams. Two threads that seed `rng_t` the same way draw
the same numbers. The generator must produce, from one seed, many streams that
are statistically independent and individually reproducible. Two families of
algorithm do this well:

- A splittable scheme (SplitMix64 seeding xoshiro256**) advances a per-thread
  state and derives child seeds by splitting. Independence between split
  children rests on the seeding mixer; the streams are sequential internally.
- A counter-based generator (Threefry, Philox) is a keyed bijection on a
  counter. The n-th draw of a stream is `cipher(key, n)`, computed directly with
  no sequential dependence. The key selects the stream; the counter indexes
  within it.

fortnum picks Threefry-2x64 with 20 rounds (Salmon, Moraes, Dror, Shaw,
"Parallel Random Numbers: As Easy as 1, 2, 3", SC'11). The counter-based form
fits parallel reproducibility directly: a stream is a key, position n is the
counter value n, and any draw is addressable without replaying the ones before
it. A thread reads its stream by setting the key from its index; no thread
touches another's state, and the result does not depend on thread scheduling.
Threefry over Philox because Threefry needs only add, xor, and rotate, so the
Fortran kernel is short, branch-free, and identical across compilers, which is
what the determinism gate in section 7 checks. The 20-round variant is the
Random123 default and passes BigCrush; fortnum does not ship a reduced-round
variant.

## 2. The rng_t state type

```fortran
type :: rng_t
    integer(int64) :: key(2)       = 0_int64
    integer(int64) :: counter(2)   = 0_int64
    integer(int64) :: buffer       = 0_int64
    logical        :: have_buffer  = .false.
    real(dp)       :: spare_normal = 0.0_dp
    logical        :: have_spare   = .false.
end type rng_t
```

`key` is the 128-bit stream selector; `counter` is the 128-bit position within
the stream. One Threefry call maps `(key, counter)` to two 64-bit outputs.
`buffer` and `have_buffer` hold the second output between `rng_next_u64` calls,
so two draws cost one block. `spare_normal` and `have_spare` hold the second
Box-Muller deviate (section 6). The whole reproducible state is `key` and
`counter`; the buffer and spare fields are derived cache that the seed and split
routines clear.

Every field is inactive. The type carries no `fortnum_status_t`; status rides on
the routine arguments, as in `fortnum_oracle` and `fortnum_ode`. There is no
seed sequence, no global pool, and no procedure pointer; the generator is a pure
function of the explicit state plus the fixed algorithm.

## 3. Seeding

```fortran
subroutine rng_seed(g, seed, status)
    type(rng_t),            intent(out) :: g
    integer(int64),         intent(in)  :: seed
    type(fortnum_status_t), intent(out) :: status
end subroutine rng_seed
```

`rng_seed` turns one 64-bit seed into a full 128-bit key by running SplitMix64
twice (Steele, Lea, Flood, "Fast Splittable Pseudorandom Number Generators",
OOPSLA'14). SplitMix64 is the recommended seeder for any 64-bit-block generator:
it spreads a small or low-entropy seed (`0`, `1`, a date) across all key bits, so
nearby seeds give unrelated streams. The constant is the golden-ratio increment
`0x9E3779B97F4A7C15`; the two finalizer multipliers are `0xBF58476D1CE4E5B9` and
`0x94D049BB133111EB`. After seeding, `counter` is zero and the buffer and spare
caches are empty. `status` is `FORTNUM_OK`; the routine cannot fail on a valid
integer seed but carries status for interface uniformity with the rest of
fortnum.

The seed is inactive (ad.md §3). It selects a stream; it is not a quantity the
output varies smoothly in.

## 4. Stream splitting

```fortran
subroutine rng_split(parent, stream, child, status)
    type(rng_t),            intent(in)  :: parent
    integer(int64),         intent(in)  :: stream
    type(rng_t),            intent(out) :: child
    type(fortnum_status_t), intent(out) :: status
end subroutine rng_split
```

`rng_split` derives an independent substream from a parent. The child copies the
parent key in the high word and xors the stream index into the low word
(`child%key(2) = ieor(parent%key(2), stream)`), then resets the counter to zero.
Two distinct stream indices give two distinct keys, so the Threefry bijection
maps them to non-overlapping output sequences with no statistical relation; this
is the counter-based independence guarantee, not a heuristic about seed spacing.
A negative `stream` reports `FORTNUM_DOMAIN_ERROR`, because the parallel policy
in section 5 maps a thread index `t in [0, nthreads)` onto `stream = t` and a
negative index is a caller bug.

`parent` is `intent(in)`: splitting does not consume or advance the parent, so
the same parent yields the same child for the same index on every call and from
every thread. The stream index is inactive (ad.md §3).

## 5. Parallel reproducibility policy

The contract a caller can rely on: thread `t` gets a deterministic, independent
stream that does not depend on how many threads run or in what order they
execute.

Seed once on the master thread, then split per thread:

```fortran
call rng_seed(master, seed, status)
!$omp parallel private(g, status)
call rng_split(master, int(omp_get_thread_num(), int64), g, status)
! ... draw from g; no other thread touches it ...
!$omp end parallel
```

Each thread owns its `rng_t`. No thread reads or writes another's state, so there
is no race and no lock. Because `rng_split` is a pure function of the parent key
and the index, thread `t` draws the same sequence whether the program runs on one
thread or sixty-four, and whether the OS schedules `t` first or last. The seed
fixes the whole ensemble; the thread index fixes the stream within it. A serial
run that replays the same seeds and indices reproduces every thread's draws.

This is the property a counter-based generator buys: reproducibility is a
function of `(seed, stream index, draw count)` alone, with no shared state to
serialize and no dependence on the schedule.

## 6. Draw routines

### 6.1 Raw 64-bit words

```fortran
subroutine rng_next_u64(g, value)
    type(rng_t),    intent(inout) :: g
    integer(int64), intent(out)   :: value
end subroutine rng_next_u64
```

`rng_next_u64` returns the next uniformly distributed 64-bit word and advances
the state. One Threefry-2x64 call yields two words: the first is returned, the
second is cached in `buffer`, and `counter` is incremented by one block. The next
call returns the cached word without a Threefry call. The counter increment
carries from the low word to the high word, so a stream spans the full 128-bit
counter space before it would repeat.

`value` is `intent(out)` and `g` is `intent(inout)`; the routine is a subroutine
rather than a function because it mutates `g`. It cannot fail and so carries no
status: an advancing draw has no error mode.

### 6.2 Uniform real in [0, 1)

```fortran
subroutine rng_uniform(g, value)
    type(rng_t), intent(inout) :: g
    real(dp),    intent(out)   :: value
end subroutine rng_uniform
```

`rng_uniform` draws one `u64`, keeps the top 53 bits, and scales by `2**-53`. The
top bits are the high-quality bits of the Threefry output, and 53 is the
significand width of `real(dp)`, so every representable multiple of `2**-53` in
`[0, 1)` is reachable and the value is never exactly `1`. This half-open
convention is what samplers expect (a `log(u)` in an exponential draw must not
see `u = 0` from a closed upper end mirrored, and must not see `1`); section 6.3
guards the `log` against the `u = 0` endpoint.

### 6.3 Standard normal

```fortran
subroutine rng_normal(g, value)
    type(rng_t), intent(inout) :: g
    real(dp),    intent(out)   :: value
end subroutine rng_normal
```

`rng_normal` returns one standard normal deviate by the polar-free Box-Muller
transform. It draws two uniforms `u1, u2`, forms `r = sqrt(-2 ln u1)` and
`theta = 2 pi u2`, returns `r cos theta`, and caches `r sin theta` in
`spare_normal` for the next call. `u1` is floored at a tiny positive value before
the `log` so the `[0, 1)` lower endpoint cannot produce a `+inf`. Box-Muller over
the ziggurat for the first implementation: it has no table, no rejection branch,
and no compiler-dependent control flow, so the determinism gate in section 7
extends to normals without a second reference scheme. A ziggurat variant may
arrive later as a separate named routine if profiling justifies it; it does not
change `rng_normal`.

Derivative classification (ad.md §3, §4): the module is `primal_only`. `rng_t`,
`seed`, and `stream` are inactive. No `rng_*_jvp`, `rng_*_vjp`, or `rng_*_grad`
exists, and the oracle CSV `has_derivative` flag (ad.md §6) is `0` for every RNG
table.

## 7. Known-answer-test policy (the determinism gate, #20)

#20 ships a known-answer test (KAT) that fails if the output changes on any
compiler or platform. Two layers:

- Algorithm KAT against the published Random123 vectors. The bare
  Threefry-2x64-20 block must reproduce the reference outputs exactly:
  `(key={0,0}, counter={0,0})` gives `{0xC2B6E3A8C2C69865, 0x6F81ED42F350084D}`,
  and `(key, counter all 0xFFFFFFFFFFFFFFFF)` gives
  `{0xE02CB7C4D95D277A, 0xD06633D0893B8B68}`. These are the canonical SC'11
  vectors and pin the constants, the round count, the rotation schedule, and the
  key-injection cadence. A mismatch means the kernel is wrong, not merely
  differently seeded.
- End-to-end KAT against a fortnum oracle table. For a fixed seed, the first N
  `rng_next_u64`, `rng_uniform`, and `rng_normal` outputs are frozen as a
  reference table under `test/oracle/data/` and checked through `fortnum_oracle`
  with `has_derivative: 0`. The `u64` words are compared exactly; the `real(dp)`
  uniforms and normals are compared at the table tolerance. This catches a
  regression in the seeding mixer, the bit-extraction in `rng_uniform`, or the
  Box-Muller arithmetic, none of which the algorithm KAT alone would see.

A reproducibility test also asserts the section 5 contract: splitting one seeded
parent into K children and drawing from each gives output independent of the
order the children are visited and of whether the draws run serially or under
OpenMP. The gate is exact equality of the per-stream sequences, not a
statistical test; statistical quality rests on the cited BigCrush result for
Threefry-20, which fortnum does not re-run.
