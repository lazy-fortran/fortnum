# fortnum architecture

fortnum is a numerical library for Fortran. The design goals are:
primal-correct values first, no global state, derivative-plural interfaces, and
a test infrastructure that cross-checks against independent references.

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
reliable autodiff, analytical, and hybrid derivative candidates.

---

## Derivative-plural design

The vocabulary and naming convention are defined in `docs/design/ad.md`.
Each derivative product can have several candidates:

| Term | Meaning |
|---|---|
| `autodiff` | A compiler or source-transformation backend differentiates smooth code. |
| `analytical` | An expression, recurrence, implicit solve, linear operator, sensitivity model, or frozen trace supplies the product. |
| `hybrid` | Autodiff composes code across analytical rules at mathematical operator boundaries. |
| `primal_only` | No derivative is meaningful for the declared active arguments. |

The public products remain `foo_jvp`, `foo_vjp`, `foo_grad`, and `foo_hvp`.
Candidate implementations are independently validated and benchmarked.
Application runtime and peak memory select the winner for each workload class.

`fortnum` will use `../fortsym` for symbolic algebra and code generation.
`fortsym` is unfinished, so no integration API is specified yet.

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
- `docs/design/differentiation_plan.md`: derivative implementation plan
- `docs/design/integrate.md`: fortnum_integrate API
- `docs/design/ode.md`: fortnum_ode API
- `docs/design/rng.md`: fortnum_rng API
- `docs/design/enzyme_toolchain.md`: optional compiler-autodiff toolchain
- `docs/design/optimizer_api.md`: backend-independent product interfaces
- `docs/design/downstream_ad.md`: downstream active-kernel pattern

The rationale for candidate generation and measured selection is
`docs/performance_optimal_differentiation.md`.
