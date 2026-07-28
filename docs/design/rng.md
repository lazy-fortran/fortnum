# Random-number contract

Status: implemented and bitwise tested.

## Algorithm

`fortnum_rng` uses Threefry-2x64 with 20 rounds. SplitMix64 expands one signed
64-bit seed into the 128-bit Threefry key. The algorithm uses fixed-width
integer add, XOR, rotate, and key-injection operations.

The algorithm and round schedule are part of the reproducibility contract.
Changing them creates a new random stream and requires new known-answer data.

## State

`rng_t` is caller-owned:

```fortran
type :: rng_t
    integer(int64) :: key(2)
    integer(int64) :: counter(2)
    integer(int64) :: buffer
    logical :: have_buffer
    real(dp) :: spare_normal
    logical :: have_spare
end type rng_t
```

The key and counter define the stream. The buffer caches the second word of a
Threefry block. The normal cache stores the second Box-Muller value.

Separate `rng_t` objects can be used concurrently. Sharing one mutable object
across threads requires caller synchronization.

## Seeding and splitting

`rng_seed(generator, seed, status)` derives both key words, resets the counter,
and clears both caches.

`rng_split(parent, stream, child, status)` derives a child key from the parent
key and a nonnegative stream index. It resets the child counter and caches.
The parent is not advanced. Negative stream indices are domain errors.

Callers assign stable stream indices from logical work identifiers. Thread
number is unsuitable when scheduling may change.

## Draws

`rng_next_u64` returns one 64-bit word. One Threefry block produces two words,
so every second call uses the cache without advancing the block counter.

`rng_uniform` maps the high 53 bits of a word to a `real(dp)` value in
`[0,1)`.

`rng_normal` uses Box-Muller. It returns one standard normal value and caches
the paired value. Its first uniform input is floored at \(2^{-53}\) to avoid
`log(0)`.

`rng_threefry2x64` exposes one raw block transformation for published
known-answer tests.

## Differentiation

Seeds, stream indices, keys, counters, and discrete draws are inactive. The RNG
module has no derivative products.

Distribution parameters and stochastic estimators belong to a higher layer.
That layer can use reparameterization, score estimators, analytical
expectations, or common random numbers while treating the sampled base noise as
fixed.

## Validation

Tests cover published Threefry block vectors, deterministic seeding, split
reproducibility, parent preservation, buffering, range constraints, normal
cache behavior, and concurrent use of independent states.
