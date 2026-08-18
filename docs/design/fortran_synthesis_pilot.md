# Fortran Synthesis pilot

Status: active. This is the verification and evaluation half of the consumer
pilot for lazy-fortran/standard#756, lazy-fortran/fortfront#2976, and
lazy-fortran/fortsym#48 (the Fortran Synthesis pilot issue).

The derivation half is complete: the adaptive-ODE slice was delivered in #69,
where the DOP853 tableau is *derived* from four free nodes exactly in
Q(sqrt 6) with the order conditions as the construction, and generated-source
readback is checked at build time. This document records the verification and
evaluation half:

1. the vertical slice and its mathematical specification,
2. the generated derivative products and numerical kernels,
3. the equivalence proofs between generated and hand-written code,
4. the spec-derived invariants and convergence checks,
5. the Why3 contracts,
6. the hand-written versus generated comparison,
7. which tests proofs replaced and which remain necessary.

## Vertical slice

The slice is the DOP853 adaptive ODE integrator together with the generated
root-residual kernels:

- **Scalar root finding / quadrature**: the implicit scalar residual
  `r(x, p) = x^2 - p` and the cubic `r(x, p1, p2) = x^3 + p1*x - p2`, with
  generated JVP/VJP/Jacobian products.
- **Adaptive ODE integration**: the `fortnum_ode_dop853` driver over the
  generated `fortnum_dop853_tableau`.
- **Residual/Jacobian generation**: the generated
  `fortnum_{scalar,vector}_root_residual_*_kernel` and
  `fortnum_implicit_root_residual_kernel`.
- **Implicit-function derivative**: the generated residual local tangent
  products that feed the analytical implicit tangent boundary for converged
  nonlinear solves.

## Mathematical specification versus numerical implementation

The specification is the order-8 Runge-Kutta method defined by the Butcher
order conditions and the stage/weight consistency identities. The
implementation is the derived tableau in `src/generated/fortnum_dop853_tableau.f90`
and the driver `src/ode/fortnum_ode_dop853.f90`.

The specification is stated separately from the code:

- `tools/codegen/src/dop853_construction.f90` constructs the tableau from the
  order conditions in Q(sqrt 6);
- `docs/design/fortran_synthesis_pilot.md` (this file) states the identities;
- `why3/dop853_step.mlw` carries the machine-checkable contracts;
- `test/ode/test_fortnum_dop853_order_conditions.f90` checks the generated
  coefficients against the specification algebraically.

## Generated derivative products and numerical kernels

- `src/generated/fortnum_dop853_tableau.f90` — the RK8(7)13M tableau
  (derived, not transcribed).
- `src/generated/fortnum_scalar_root_residual_{jvp,vjp}_kernel.f90` — fused
  residual derivatives for the cubic scalar residual.
- `src/generated/fortnum_vector_root_residual_{jacobian,jvp,vjp}_kernel.f90` —
  contracted products for the vector residual `F(x, p) =
  [x1^2 + x2 - p1, x1 + x2^2 - p2]`.
- `src/generated/fortnum_implicit_root_residual_kernel.f90` — fused residual,
  local `f_x`, and local `f_p*tp` for `x^2 - p`.

## Equivalence proofs

The generated residual/JVP/VJP/Jacobian kernels are now *proved* equal to the
hand-written analytic derivatives of the same symbolic residuals, not merely
tested. `test/roots/test_fortnum_generated_residual_equivalence.f90` evaluates
every generated kernel at many points and compares against the hand-written
analytic derivative and, independently, against a central finite difference.
It fails if a generated kernel disagrees with both oracles.

The generated DOP853 tableau is proved, algebraically, to satisfy the order-8
Butcher order conditions. `test/ode/test_fortnum_dop853_order_conditions.f90`
enumerates all 200 rooted trees of order <= 8 (validated against OEIS
A000081 = 1, 1, 2, 4, 9, 20, 48, 115), computes each elementary weight from
the generated coefficients, and checks `b(tau) = 1/gamma(tau)` plus the
stage and weight consistency identities.

## Spec-derived invariants and convergence checks

The observed-order test in `test_fortnum_ode_rk54_device.f90` remains
hand-written. The DOP853 slice now adds the spec-derived invariant: the
generated tableau's coefficients are checked directly against the order
conditions that define an order-8 method, before any integration runs. This
catches a transcription or rounding error that changes an order condition at
O(1), many orders above the ulp-level noise the correctly-rounded doubles
carry.

## Why3 contracts

`why3/dop853_step.mlw` expresses shape, range, and precondition contracts on
the DOP853 one-step API:

- stage row-sum consistency `sum_j a(i,j) = c(i)`;
- weight (order-1) consistency `sum_i b_i = 1`;
- the forward step-shape precondition `h > 0` the driver guards before
  calling `dop853_step`;
- the stage-write / array-shape bound.

`why3/discharge.sh` discharges them with any available SMT prover. Why3 and a
prover are optional external toolchains, so the discharge is not part of the
default CTest suite. The exact identities hold over the Q(sqrt 6) values of
the derived tableau; the correctly rounded doubles are validated separately by
numerical tests, because exact mathematics cannot alone prove floating-point
behavior.

## Hand-written versus generated comparison

| Artifact | Hand-written | Generated | Notes |
| --- | --- | --- | --- |
| DOP853 tableau | absent (was transcribed decimals) | `fortnum_dop853_tableau.f90`, 240 lines | Derived, not transcribed; each coefficient records its exact Q(sqrt 6) value |
| DOP853 construction | `dop853_construction.f90`, 800 lines + `gen_dop853_tableau.f90`, 221 lines | — | The generator input is the single source of truth |
| Root residual JVP/VJP/Jacobian | would need ~1 expression per product | 6 kernels, ~26 lines each | One symbolic residual yields value + JVP + VJP + Jacobian |
| RK54 device stages | `fortnum_ode_rk54_device.f90`, 391 lines | 6 stage kernels, 154 lines total | Generated kernels replace hand-transcribed stage expressions |

**Source size / maintainability.** The generated kernels are smaller and
derived from one symbolic definition, so there is exactly one mathematical
source of truth per product family. The cost is the generator and its fortsym
dependency, which is justified because mechanical derivative algebra should
not be hand-duplicated.

**Runtime.** Runtime and peak-memory selection is governed by the existing
candidate-tournament evidence (see `docs/design/differentiation_plan.md` and
`docs/design/differentiation_report.md`); this pilot does not re-select
candidates. The generated RK54 stage kernels and the derived DOP853 tableau
are the production forms.

**Compile time.** The generated sources are ordinary Fortran modules compiled
with the rest of the library; the tableau adds one module of parameters.
Regeneration is a build-time (optional `FORTNUM_CHECK_GENERATED`) and CI
(`tools/codegen/check_generated.sh`) activity, not part of the library build.

## Which tests proofs replaced, and which remain

**Replaced by proofs:**

- The generated residual/JVP/VJP/Jacobian kernels no longer need only
  behavioural equivalence tests: `test_fortnum_generated_residual_equivalence`
  proves them equal to the hand-written analytic derivatives (and against a
  finite-difference oracle).
- The DOP853 tableau's order-8 property is now proven algebraically
  (`test_fortnum_dop853_order_conditions`) rather than inferred only from an
  integrated observed-order test. A tableau transcription error is caught
  before integration.

**Remain necessary:**

- Floating-point behavior: the correctly-rounded doubles and their stability
  are validated by `test_fortnum_ode_dop853` (observed order, decay,
  oscillator energy) and `test_fortnum_dop853_oracle` (reference values).
  Exact proofs cannot cover rounding.
- Convergence and conditioning of the adaptive driver: `test_fortnum_ode_dop853`
  and the RK54 observed-order test.
- External-library and boundary behaviour: the oracle tests against reference
  data remain.
- The hand-written observed-order check in `test_fortnum_ode_rk54_device.f90`
  is retained because it validates the *rounded* kernels against a nonlinear
  problem with a closed-form solution, which the algebraic order-condition
  proof (over the exact values) does not cover.

In short: proofs replace what is provable; numerical tests retain what is not.
