# Stage precision and mixed refinement

Use `fortnum_precision` for `precision_policy_t`, `precision_pair_t`, `FP64`,
`FP32`, `FP16`, `BF16` and `unit_roundoff(format)`. Storage identities come from
[FortGen](https://github.com/lazy-fortran/fortgen); they are not compiler kinds.

`precision_preset(name,ok)` accepts `cpu-fp64`, `cpu-mixed`, `gpu-consumer`,
and `gpu-hpc`. Presets describe arithmetic policy only: GPU presets do not
select or promise an executable GPU backend. `%get(stage)` returns a pair;
`%set(stage,storage,compute,ok)` permits fp32/fp64 compute and independent
storage, rejecting narrowing fp64 storage into fp32 arithmetic. Stages are
`STAGE_EVALUATE`, `STAGE_ASSEMBLE`, `STAGE_PRECONDITION`, `STAGE_KRYLOV`,
`STAGE_RESIDUAL`, `STAGE_ESTIMATE`, and `STAGE_CERTIFY`. The last three are
locked to fp64 in private policy state. Invalid mutations preserve the policy.

`%set_budget(fraction,ok)` sets the allowed fraction of model/discretization
error (default 0.1). `%assess(stage,rounding_error,model_error,accepted)` checks
a caller-supplied a posteriori rounding estimate **after sensitivity
amplification**. Rejection escalates half storage to fp32, then fp64. The
caller must rebuild that stage's storage and rerun it before accepting its
result. This is an error-budget mechanism, not an estimator or certificate.
At fp64, rejection remains rejection; locked stages never silently lower.

`fortnum_precision_vector:precision_vector_t` provides `%store(values,format,ok)`,
`%load(output,ok)` for real32/real64 outputs, `%at(i)`, `%at_scaled(i)`,
`%scale_factor()`, `%is_fp64()`, `%size()`, and `%bytes()`. Indices passed to
scalar accessors must be valid. Narrow storage uses one exact power-of-two
scale. Nonfinite input or a rounded reconstruction outside fp64 range fails
without replacing old data. A wide dynamic range can still underflow small
entries: the half absolute quantum, not a uniform relative bound, applies.
Loading unscaled magnitudes outside fp32 range fails explicitly. `%bytes()`
counts numeric payload only, excluding object metadata and allocation overhead.
OpenMP parallelizes large loads/stores; setup owns allocations.

`fortnum_refinement:mixed_refine(matvec64,correction,b,x,policy,tolerance,`
`max_steps,info,steps,residual)` accepts an fp64 original-operator callback
and a correction callback `(normalized_rhs,delta,precision_pair,info)`.
The callback builds/applies the requested low-precision solver. Residuals,
updates and acceptance always use fp64. Defects are normalized before the
inner solve; failed/nonfinite solves and measured stagnation escalate the
preconditioner stage. An accepted solution satisfies the infinity-norm relative
residual target. `info=0` means success, `1` step limit, `2` nonfinite original
residual, `3` failed fp64 correction, `-1` invalid arguments. This routine does
not certify a forward error, invent a model error, or choose a physics tolerance.
