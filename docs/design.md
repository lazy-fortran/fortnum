# fortnum architecture

fortnum is a numerical library for Fortran. The design goals are:
primal-correct values first, no global state, derivative-ready interfaces, and
a test infrastructure that cross-checks against high-precision Python references.

---

## Module map

```
fortnum_kinds          -- kind parameters (dp, sp, i4, i8)
fortnum_status         -- error codes and the fortnum_status_t type

fortnum_special        -- umbrella re-export
  fortnum_special_bessel   -- modified Bessel I_n, K_n
  fortnum_special_dawson   -- Dawson integral
  fortnum_special_gamma    -- lower incomplete gamma and regularized P

fortnum_fft            -- 1D DFT (mixed-radix Stockham + Bluestein)
fortnum_quadrature     -- fixed Gauss-Legendre rule generation

fortnum_integrate_gk   -- single-panel GK pairs (G7K15 through G30K61)
fortnum_integrate      -- globally adaptive driver (QAG/QAGS/QAGP/QAGIU)

fortnum_ode            -- Cash-Karp RK5(4) adaptive integrator
  fortnum_ode_cash_karp  -- stage kernel (internal)
  fortnum_ode_events     -- event detection (internal)
  fortnum_ode_wrapper    -- ode_at: evaluate at prescribed output times

fortnum_roots          -- scalar root-finding (bisect/Newton/Brent)
fortnum_rng            -- Threefry-2x64-20 counter-based PRNG

fortnum_interp         -- binary grid search
fortnum_polynomial     -- Lagrange interpolation weights

fortnum_oracle         -- CSV oracle reader and primal checker (testing)
```

---

## Kinds and status

`fortnum_kinds` re-exports `iso_fortran_env` kind aliases. The primary real kind
is `dp = real64`; `sp = real32` appears only in mixed-precision interfaces.
Integer kinds `i4 = int32` and `i8 = int64` cover index spaces and counter
arithmetic.

`fortnum_status` carries error information without exceptions. Every public
subroutine that can fail takes a `type(fortnum_status_t), intent(out) :: status`
argument. The type is a `(code, msg)` pair; `FORTNUM_OK = 0` means success.
Codes are stable across releases because callers may branch on them. The status
object is inactive in the derivative sense: it is a side channel, not a
differentiable output.

---

## No global state

No module defines a `save` attribute on a mutable variable, no module-level
procedure pointer, and no hidden pool. All state lives in caller-owned derived
types: `fortnum_fft_plan_t`, `integrate_workspace_t`, `ode_workspace_t`,
`rng_t`, and so on. This property is the precondition for thread safety and for
the `transparent` derivative path: Enzyme can differentiate the implementation
only when the implementation has no hidden state.

---

## Derivative-ready design

Every public routine carries a derivative policy. The policy vocabulary and the
naming convention for derivative entry points are defined in
`docs/design/ad.md`. The short summary:

| Policy | Meaning |
|--------|---------|
| `transparent` | Enzyme differentiates the implementation directly. |
| `analytic_rule` | fortnum supplies handwritten JVP/VJP using a closed form or recurrence. |
| `implicit_rule` | Derivative defined by an implicit equation, not the iteration. |
| `trace_rule` | Derivative taken at the frozen schedule the adaptive primal chose. |
| `primal_only` | No derivative is mathematically meaningful; RNG is the canonical case. |

The primal ships first (M1 through M5). Derivative products (`foo_jvp`,
`foo_vjp`, `foo_grad`, `foo_hvp`) ship in issue #40, added beside the primal
as new public names without touching existing signatures. The reserved names
appear in the source comments and in `docs/api.md` but are not yet compiled.

---

## Oracle-table testing

Tests for deterministic functions use Python-generated CSV reference tables.
`test/oracle/gen_oracle.py` calls `mpmath` or `scipy` with high precision and
writes rows of `(index, x, primal, derivative)`. The Fortran test reads the
table with `fortnum_oracle`, calls the Fortran implementation, and asserts that
every entry passes an absolute + relative tolerance check.

The CSV format already reserves the derivative column behind the
`has_derivative` header flag. Tables written before issue #40 set the flag to
`0`; when a derivative product lands, the same file gets `has_derivative: 1`
and the column filled. No format change is needed and the reader is
forward-compatible.

---

## Build

Two build systems are supported:

**CMake + Ninja** (primary):

```
cmake -S . -B build -G Ninja
cmake --build build -j$(nproc)
ctest --test-dir build
```

**fpm** (alternate, for fpm-based consumers):

```
fpm build
fpm test
```

The `fo` tool wraps both paths and adds static analysis and formatting checks.
Run `fo` with no arguments for the full pipeline.

---

## Module ADRs

Detailed design decisions for individual modules:

- `docs/design/ad.md`: derivative contract (normative for all modules)
- `docs/design/integrate.md`: fortnum_integrate API
- `docs/design/ode.md`: fortnum_ode API
- `docs/design/rng.md`: fortnum_rng API

M6 design files (`enzyme_toolchain.md`, `optimizer_api.md`,
`downstream_ad.md`) will be linked here when they are written.
