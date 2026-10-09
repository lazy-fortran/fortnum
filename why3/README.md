# Why3 contracts for the DOP853 synthesis slice

This directory holds the machine-checkable contracts of the Fortran
Synthesis pilot (lazy-fortran/fortnum issue #62). The derivation half of the
pilot (the DOP853 tableau) is complete; this is the verification half.

`dop853_step.mlw` expresses, in Why3, the shape/range/precondition contracts
of the DOP853 one-step API and the generated tableau:

- stage row-sum consistency `sum_j a(i,j) = c(i)` (stage 1 provable from the
  derived values; the full set is checked numerically in
  `test/ode/test_fortnum_dop853_order_conditions.f90`);
- weight (order-1) consistency `sum_i b_i = 1`;
- the forward step-shape precondition `h > 0` that the driver
  `src/ode/fortnum_ode_dop853.f90` guards before calling `dop853_step`;
- the stage-write / array-shape bound.

## Discharge

The Why3 toolchain and an SMT prover are not part of the default build. To
discharge:

```sh
./discharge.sh
```

or directly, e.g. `why3 prove dop853_step.mlw -P z3`.

## Why exact mathematics is not the whole story

The dischargeable identities hold over the exact Q(sqrt 6) values of the
*derived* tableau. The Fortran code uses correctly-rounded doubles of those
values, and exact mathematics cannot prove floating-point behavior. That
split is intentional: proofs replace what is provable, and numerical tests
(`test_fortnum_dop853_order_conditions`, `test_fortnum_ode_dop853`,
`test_fortnum_dop853_oracle`) retain what is not.
