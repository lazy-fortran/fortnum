# Performance-Optimal Differentiation in the Age of AI and Computer Algebra

## An evidence-driven strategy for analytical derivatives, automatic differentiation, implicit differentiation, and hybrid numerical software

## Abstract

Automatic differentiation has become the default differentiation strategy in
much modern scientific software because it reduces implementation effort,
tracks changes in primal code, and supplies derivatives to complex optimization
algorithms. That historical argument changes when an AI system can derive
formulas, operate a computer algebra system, generate source code, construct
custom adjoints, and maintain multiple derivative implementations at low human
cost.

In that setting, automatic differentiation should no longer be presumed optimal
merely because it is convenient. Symbolic differentiation should not be
presumed optimal merely because exact formulas are cheap to generate. The
objective is narrower:

$$
\boxed{\text{minimize measured runtime and peak memory for the complete application}}
$$

subject to numerical correctness and the exact derivative contract required by
the caller.

The most efficient implementation may be:

- an explicit analytical formula
- a recurrence relation
- forward-mode automatic differentiation
- reverse-mode automatic differentiation
- a source-transformed or compiler-generated derivative
- an implicit-function derivative
- a tangent or adjoint solve
- a custom Jacobian-vector or vector-Jacobian product
- a frozen-trace derivative of an adaptive algorithm
- a hybrid composition of several of these

No method wins by policy. Every performance-relevant numerical operator should
expose several admissible derivative realizations. AI and computer algebra
infrastructure should generate, simplify, compile, validate, benchmark, and
select among them.

The recommended architecture is:

$$
\boxed{
\text{mathematical decomposition}
+
\text{candidate derivative generation}
+
\text{application-level benchmarking}
+
\text{evidence-based selection}
}
$$

The implicit function theorem is central to this strategy. Whenever an output
is defined by a root, equilibrium, optimization condition, fixed point, linear
system, eigensystem, or other residual equation, differentiation of the
defining equation should be an explicit candidate, and usually the leading
candidate. The iterative algorithm used to obtain the solution is a separate
candidate.

This paper develops that policy generically and applies it to `fortnum`, with
DESC as an instructive comparison.

## 1. The decision problem has changed

### 1.1 The historical compromise

Derivative implementation traditionally involved a three-way tradeoff:

1. **Analytical derivatives**

   - potentially fastest
   - difficult to derive
   - expensive to maintain
   - prone to transcription errors

2. **Automatic differentiation**

   - easier to maintain
   - broadly applicable
   - often sufficiently fast
   - potentially expensive in runtime, memory, compilation, or code complexity

3. **Finite differences**

   - easy to deploy
   - useful for validation
   - generally slower for many inputs
   - less accurate and sensitive to step size

AI-assisted derivation and computer-algebra code generation reduce the
derivation and maintenance costs of analytical derivatives. An AI agent can:

- derive exact formulas
- simplify them with a computer algebra system
- generate generic source code
- regenerate derivatives after primal changes
- construct alternative algebraic factorizations
- implement tangent and adjoint forms
- derive implicit sensitivities
- generate validation tests
- benchmark all variants

Implementation effort should therefore cease to determine the production
method unless generation cost itself becomes operationally significant.
Derivative selection becomes an empirical numerical-kernel optimization
problem.

### 1.2 The new principle

The production implementation is the derivative algorithm that minimizes
measured end-to-end runtime and memory for the required derivative operation.

A fully expanded symbolic derivative can be much worse than reverse
autodifferentiation through a compact recurrence. A custom implicit adjoint can
be orders of magnitude better than reverse autodifferentiation through many
nonlinear-solver iterations. The labels `analytical` and `autodiff` are too
coarse to predict performance.

## 2. Differentiate the mathematical problem, not necessarily the executed algorithm

Suppose the primal computation is represented as

$$
y = \mathcal A(p),
$$

where $\mathcal A$ is an algorithm. The mathematical result may instead be
defined by

$$
R(y,p) = 0.
$$

These descriptions admit different differentiation strategies.

### 2.1 Algorithmic differentiation

Algorithmic differentiation follows every operation performed by $\mathcal A$:

$$
p \longrightarrow \text{iteration 1} \longrightarrow \cdots
\longrightarrow \text{iteration }K \longrightarrow y.
$$

This can require storing the execution trace, differentiating convergence
tests, handling changing iteration counts, propagating derivatives through
preconditioners, and reproducing implementation details irrelevant to the
mathematical solution.

### 2.2 Implicit differentiation

Differentiate the defining equation:

$$
R_y\,dy + R_p\,dp = 0,
$$

which gives

$$
\boxed{\frac{dy}{dp} = -R_y^{-1}R_p}.
$$

The inverse denotes a linear solve, never a requirement to form an explicit
inverse.

For a scalar objective $L(y,p)$,

$$
\frac{dL}{dp}
= L_p - L_yR_y^{-1}R_p.
$$

With an adjoint variable satisfying

$$
R_y^T\lambda = L_y^T,
$$

the gradient becomes

$$
\boxed{\frac{dL}{dp} = L_p - \lambda^T R_p}.
$$

This can reduce one sensitivity solve per parameter to one adjoint solve per
scalar objective.

### 2.3 Systematic use of the implicit-function rule

Implicit candidates should be considered for:

- scalar and vector roots
- nonlinear equilibria and fixed points
- linear systems
- optimization stationarity and KKT systems
- eigenvalues and invariant subspaces
- spline coefficients obtained from a solve
- quadrature nodes defined by orthogonality conditions
- event times defined by zero crossings
- boundary-value problems and differential-algebraic systems
- self-consistency loops
- normalization and projection conditions

The performance question depends on dimensions, sparsity, factorization reuse,
and hardware. Mathematical availability alone does not choose the winner.

## 3. Select derivative objects before derivative methods

For

$$
f:\mathbb R^n\rightarrow\mathbb R^m,
$$

a caller may require

$$
J,\quad Jv,\quad J^Tu,\quad \nabla f,\quad Hv,\quad H,\quad
\operatorname{diag}(H).
$$

The best method can change completely with the requested object.

### 3.1 Full Jacobian

A full Jacobian costs at least $O(nm)$ storage. Produce it only when the
consumer needs it and matrix-free products are not competitive. Candidates
include symbolically generated sparse Jacobians, forward columns, reverse rows,
compressed recovery using coloring, local analytical assembly, and implicit
tangent solves for selected blocks.

### 3.2 Jacobian-vector product

$Jv$ is naturally produced by forward autodifferentiation, a differentiated
recurrence, a tangent linear model, an implicit tangent solve, or direct
computer-algebra generation of the contracted expression. Forming $J$ first is
usually inferior for large problems.

### 3.3 Vector-Jacobian product

$J^Tu$ is naturally produced by reverse autodifferentiation, an analytical or
generated adjoint, an implicit adjoint solve, a transpose recurrence, or direct
generation of the contracted expression.

### 3.4 Hessian-vector product

For scalar $f$,

$$
Hv = \nabla^2 f\,v
$$

can use forward-over-reverse autodifferentiation,
reverse-over-forward autodifferentiation, differentiated adjoint equations,
implicit second-order sensitivities, or a directly generated contraction.

### 3.5 Production API

The numerical API should expose derivative products:

```text
value(...)
jvp(direction, ...)
vjp(cotangent, ...)
gradient(...)
hvp(direction, ...)
jacobian(...)      # only when justified
hessian(...)       # only when justified
```

Each product gets its own candidate-selection process.

## 4. Candidate methods

### 4.1 Explicit analytical formulas

Explicit generated formulas are strong candidates when dimensions are fixed or
small, expressions remain compact, simplification removes substantial work,
common subexpressions can be shared, sparsity is strong, and generated code
preserves vectorization without excessive register pressure.

The computer algebra system should operate on an expression DAG. Value and
derivative should be generated jointly when the caller consumes both.

```text
expr = CAS.parse(primal_expression)

requested = [
    value(expr),
    gradient(expr, active_parameters),
    jvp(expr, active_parameters, direction),
    vjp(expr, active_parameters, cotangent)
]

joint_dag = CAS.combine(requested)
joint_dag = CAS.simplify(joint_dag)
joint_dag = CAS.common_subexpression_elimination(joint_dag)
joint_dag = CAS.factor_for_target(joint_dag, hardware_profile)

codegen(joint_dag, target_language, target_architecture)
```

`fortnum` uses `fortsym` as its build-time symbolic algebra and code-generation
library. Its first pinned integration covers expression DAGs, differentiation,
simplification, fused kernel emission, operation counts, and regeneration
metadata. Further interfaces remain provisional until implemented and tested in
`fortsym`.

### 4.2 Analytical recurrences

For

$$
x_{k+1}=F_k(x_k,p),
$$

the tangent recurrence is

$$
\dot x_{k+1}=F_{k,x}\dot x_k+F_{k,p}\dot p,
$$

and the adjoint recurrence is

$$
\lambda_k=F_{k,x}^T\lambda_{k+1}.
$$

This preserves $O(K)$ complexity and avoids expression growth. It applies to
orthogonal polynomials, special functions, FFT butterflies, time integrators,
continued fractions, filters, and transforms.

### 4.3 Forward autodifferentiation

Forward autodifferentiation is a strong candidate when the active input
direction count is small, many outputs are required, memory should remain near
primal memory, or the computation is naturally streaming. Candidate generation
should test multiple tangent block widths because the best width depends on
register count, cache size, SIMD width, and output dimension.

### 4.4 Reverse autodifferentiation

Reverse autodifferentiation is a strong candidate for many inputs and few
outputs when the trace is moderate, checkpointing or recomputation is
effective, and custom rules exist for expensive primitives. Solver iterations,
adaptive control logic, and opaque library calls require explicit competing
candidates.

### 4.5 Compiler or source-transformation autodifferentiation

Compiler-level autodifferentiation can see optimized intermediate
representation, eliminate dead derivative work, inline across calls, perform
activity analysis, reuse primal values, and exploit alias information. These
advantages make it a required candidate, not an automatic winner.

### 4.6 Implicit tangent and adjoint rules

For $R(y,p)=0$, the tangent implementation evaluates the primal solution,
constructs or applies $R_y$ and $R_p\,dp$, then solves

$$
R_y\,dy=-R_p\,dp.
$$

The adjoint implementation evaluates $L_y$ and $L_p$, solves

$$
R_y^T\lambda=L_y^T,
$$

then forms

$$
\frac{dL}{dp}=L_p-\lambda^TR_p.
$$

Candidates should reuse primal Jacobians, factorizations, preconditioners, and
block structure when available.

### 4.7 Frozen-trace differentiation

Adaptive algorithms make discrete choices about timesteps, subdivision,
polynomial order, pivoting, events, meshes, and termination. A useful
derivative often records the primal schedule and differentiates the resulting
fixed computation.

Candidates include frozen-trace forward and reverse autodifferentiation,
analytical tangent and adjoint recurrences, and continuous sensitivities where
the caller intends the continuous problem.

### 4.8 Finite differences

Finite differences remain validation oracles, emergency fallbacks, benchmark
baselines, and tools for black-box external functions. They should still be
measured.

## 5. Computer algebra does not eliminate autodifferentiation

### 5.1 Expression growth

Repeated substitution can turn an $O(K)$ recurrence into an enormous
straight-line expression. Source size, compilation time, binary size, register
pressure, spills, instruction-cache misses, and GPU occupancy can all worsen.
The symbolic system must preserve algorithmic structure where beneficial.

### 5.2 Numerical stability

Equivalent formulas can differ in cancellation, scaling, overflow, underflow,
branch behavior, and conditioning. A generator should emit numerically stable
algorithms such as `solve(A,b)` and preserve stable recurrences.

### 5.3 Hardware effects

Operation count alone does not determine runtime. Candidate benchmarks should
measure elapsed time, peak memory, allocations, cache misses, branches,
vectorization, register spills, instruction count, GPU occupancy, transfers,
and compilation latency when relevant.

## 6. Hybrid differentiation follows mathematical operator boundaries

### 6.1 Recommended layering

**Layer A: global mathematical sensitivity structure**

Derive implicit equations, tangent and adjoint models, KKT sensitivities,
discrete-time adjoints, checkpointing schedules, Schur complements,
conservation reductions, and symmetry reductions.

**Layer B: local residual and physics kernels**

Generate contracted products such as

$$
R_yv,\quad R_pv,\quad R_y^Tu,\quad R_p^Tu
$$

with analytical code generation or local autodifferentiation.

**Layer C: custom rules for numerical primitives**

Provide candidate rules for roots, linear solves, factorizations, eigensolvers,
interpolation, FFTs, quadrature, ODE integrators, event location, projections,
normalization, and external libraries.

**Layer D: compiler autodifferentiation**

Use autodifferentiation to connect local operators and cover smooth code that
does not justify a specialized rule.

### 6.2 Nonlinear solve inside an objective

For

$$
R(y,p)=0,\qquad L=L(y,p),
$$

a high-performance hybrid can:

1. solve the primal residual
2. reuse the primal factorization or preconditioner
3. generate local $R_yv$ and $R_pv$ kernels
4. solve the adjoint equation
5. use a generated VJP for $\lambda^TR_p$
6. use compiler autodifferentiation for residual regions where symbolic code is
   inefficient

This preserves the mathematical solver boundary while allowing
autodifferentiation inside and around it.

## 7. Lessons from DESC

DESC combines spectral representations, equilibrium solves, optimization,
autodifferentiation, analytical spatial derivatives, and implicit equilibrium
sensitivities.

Its objective interfaces standardize Jacobians, JVPs, VJPs, gradients, and
related products. Its Fourier-Zernike representation supplies analytical
spatial derivatives through the plasma volume. Autodifferentiation supplies
residual and objective derivatives for perturbation and optimization.

For an equilibrium constraint $F(x,c)=0$, DESC uses

$$
\frac{dx}{dc}
=-\left(\frac{\partial F}{\partial x}\right)^{-1}
\frac{\partial F}{\partial c}
$$

and composes this sensitivity with objective products. The architecture is
therefore hybrid: analytical basis derivatives, autodifferentiated residuals
and objectives, and implicit differentiation across the equilibrium solve.

A more aggressive performance policy can also benchmark generated residual
products, forward and reverse modes, explicit sparse and matrix-free products,
tangent block sizes, factorization reuse, and custom rules for transforms and
linear solves. DESC demonstrates the separation of mathematical operators. It
does not establish one universal backend winner.

## 8. Implications for `fortnum`

### 8.1 Current strengths

`fortnum` is primal-first and already exposes additive derivative products. Its
design recognizes autodiff-compatible implementations, analytical rules,
implicit rules, frozen-trace rules, and primal-only operations. It correctly
uses product APIs such as JVP, VJP, gradient, and HVP.

### 8.2 Replace exclusive policy classes

The current requirement that every procedure have exactly one derivative
policy is too restrictive. A procedure can support an autodiff JVP, an
analytical JVP, a recurrence, an implicit JVP, a custom VJP, and
hardware-specific variants.

The revised mapping is:

```text
procedure + derivative product + workload class
    -> admissible implementations
    -> validated candidates
    -> benchmark-selected implementation
```

### 8.3 Candidate mechanisms

Candidate mechanisms include:

```text
autodiff_forward
autodiff_reverse
analytical_expression
analytical_recurrence
analytical_implicit_tangent
analytical_implicit_adjoint
analytical_frozen_trace_forward
analytical_frozen_trace_reverse
hybrid_custom_rule
continuous_sensitivity
discrete_adjoint
sparse_explicit
matrix_free_product
external_custom_rule
finite_difference_reference
primal_only
```

The public terminology remains `autodiff`, `analytical`, and `hybrid`.
Mechanism names refine those categories internally.

### 8.4 Stable public contract

The public API remains:

```text
foo(...)
foo_jvp(...)
foo_vjp(...)
foo_grad(...)
foo_hvp(...)
```

Candidate selection can occur at build time, installation, first use, explicit
configuration, or from a benchmark registry. Dispatch must be resolved before
hot inner loops.

## 9. Module strategy for `fortnum`

### 9.1 Special functions

Candidates should include computer-algebra generated formulas, stable
recurrences, differentiation of the numerical approximation, compiler
autodifferentiation, and high-precision or complex-step references. Winners may
vary by argument region. Dispatch criteria remain inactive.

### 9.2 FFT

At the operator level,

$$
Jv=\operatorname{FFT}(v).
$$

Candidates include the primal FFT for JVP, the adjoint transform for VJP,
compiler autodifferentiation, and custom external-library rules.

### 9.3 Fixed quadrature

For

$$
Q(f)=\sum_i w_i f_i,
$$

the analytical derivative is the same weighted linear map. Parameters affecting
nodes or weights require separate derivative candidates.

### 9.4 Adaptive integration

Candidates include frozen-trace forward and reverse autodifferentiation,
analytical differentiation of the fixed quadrature sum, differentiation under
the integral sign, parameter-dependent boundary terms, and finite-difference
references.

### 9.5 Root solvers

For $F(x,p)=0$, leading candidates are implicit tangent and adjoint solves.
Differentiated iterations and finite differences remain measured candidates.
Implementations should reuse final Jacobians and factorizations, support
multiple right-hand sides, and report conditioning.

### 9.6 Linear algebra

For

$$
A(p)x(p)=b(p),
$$

the tangent equation is

$$
A\,dx=db-dA\,x.
$$

Candidate rules should reuse factorizations, expose transpose solves, and
support matrix-free products.

### 9.7 Interpolation and splines

Cell and knot-span search is inactive. Basis evaluation inside a fixed span can
use analytical formulas, generated expressions, or autodifferentiation.
Coefficient solves use implicit candidates. Parameterized knot crossings need
explicit nonsmooth handling.

### 9.8 ODE solvers

Candidates include forward sensitivities, reverse autodifferentiation of the
frozen discrete solver, generated discrete adjoints, continuous adjoints,
implicit stage rules, checkpointing, recomputation, differentiated dense
output, and implicit event-time derivatives.

Continuous and discrete adjoints solve different finite-step derivative
contracts. Candidate validation must preserve that distinction.

### 9.9 Random number generation

Seeds remain primal-only. Estimator gradients through reparameterization,
score functions, pathwise derivatives, common random numbers, or analytical
expectations belong at the estimator layer.

## 10. The derivative tournament

### 10.1 Candidate generation

For each operator, inspect the primal source, mathematical contract, active
inputs, output dimensions, residual equations, sparsity, recurrence structure,
solver boundaries, and hardware target. Generate every applicable candidate.

### 10.2 Validation

Each candidate must pass:

1. primal consistency
2. directional derivative tests
3. adjoint dot-product tests
4. high-precision or complex-step comparison where applicable
5. branch, singularity, and regime-boundary tests
6. optimization-level consistency
7. deterministic reproducibility

The central identities are:

$$
\frac{f(x+\epsilon v)-f(x-\epsilon v)}{2\epsilon}\approx Jv,
$$

$$
u^T(Jv)\approx(J^Tu)^Tv,
$$

and for implicit candidates,

$$
R_y\,dy+R_p\,dp\approx0.
$$

### 10.3 Workloads

Benchmarks cover realistic dimensions, batches, layouts, direction counts,
cache states, hardware, combined value-and-derivative evaluation, factorization
reuse, solver tolerances, and iteration counts.

### 10.4 Metrics

At minimum, record

$$
T_{\rm value},\quad T_{\rm derivative},\quad
T_{\rm value+derivative},\quad M_{\rm peak}.
$$

Also measure allocations, bytes moved, tape size, checkpoint storage, code
size, compile time, vectorization, cache misses, register spills, factorization
reuse, and solver iterations.

### 10.5 End-to-end selection

The selected candidate minimizes representative application cost, including
compilation, primal work, derivative work, linear algebra, and optimization
iteration count. An isolated kernel winner need not be the application winner.

## 11. Performance registry and dispatch

Benchmark evidence should produce a machine-readable registry:

```text
operator: multiroot
product: vjp
problem_class:
    state_dimension: 128..512
    parameter_dimension: 1000+
    scalar_outputs: 1
    sparse_jacobian: true
hardware:
    architecture: target_A
winner:
    implementation: analytical_implicit_adjoint_matrix_free
    checkpointing: none
    factorization_reuse: true
evidence:
    runtime: ...
    peak_memory: ...
    validation_error: ...
```

Generated selectors may use dimensions, sparsity, direction count, output
shape, hardware, memory, and reusable primal data. Thresholds come from
benchmarks.

## 12. Recommendations for `fortnum`

1. Replace exclusive policy classes with candidate sets.
2. Treat Enzyme as an `autodiff` backend, not the differentiation architecture.
3. Add a common custom-rule layer for `analytical`, `autodiff`, and `hybrid`
   implementations.
4. Make implicit tangent and adjoint products first-class library mechanisms.
5. Generate contracted products directly.
6. Benchmark fused and separate value-and-derivative kernels.
7. Preserve recurrences, transforms, factorizations, solvers, and sparse
   operators unless expansion wins by measurement.
8. Add representative application benchmarks.
9. Use `fortsym` for symbolic algebra and code generation, extending and testing
   it there rather than duplicating symbolic machinery in `fortnum`.

## 13. Development sequence

### Phase 1: derivative contract

Preserve the public product API. Replace exclusive policies with candidate
registration. Define implementation, validation, and performance metadata.

### Phase 2: implicit infrastructure

Implement generic tangent and adjoint interfaces for roots, linear solves,
fixed points, spline coefficient systems, and event times. Include reuse hooks
and conditioning diagnostics.

### Phase 3: symbolic generation

Use the isolated `tools/codegen/` package and `fortsym` for the pipeline:

```text
specification
    -> symbolic model
    -> requested derivative products
    -> simplification variants
    -> code generation
    -> compilation
    -> validation
    -> benchmark
    -> registry
```

This phase states the intended role of `fortsym`. Its implementation waits for
the library to mature.

### Phase 4: autodiff candidates

Generate forward and reverse candidates for smooth kernels. Test activity
layouts, fused products, and custom rules around solver and transform
boundaries.

### Phase 5: module tournaments

Run tournaments for special functions, interpolation, FFT, quadrature, roots,
linear algebra, and ODE solvers. Commit the benchmark evidence.

### Phase 6: downstream applications

Benchmark representative plasma workloads. Application-level selection can
differ from microbenchmark selection because of batching, inlining, cache
behavior, memory pressure, factorization reuse, and optimizer convergence.

## 14. Decision table

| Mathematical structure | Leading candidates | Usually avoid |
|---|---|---|
| Small explicit formula | generated analytical, compiler autodiff | finite differences |
| Long recurrence | analytical tangent/adjoint recurrence, autodiff loop | full symbolic expansion |
| Scalar objective, many parameters | reverse autodiff, analytical adjoint | one forward sweep per parameter |
| Few parameters, many outputs | forward autodiff, tangent model | one reverse sweep per output |
| Root or equilibrium | analytical implicit tangent/adjoint | solver-iteration differentiation as an unmeasured default |
| Linear solve | analytical implicit rule with factorization reuse | Krylov-iteration differentiation as an unmeasured default |
| Eigenproblem | analytical perturbation or implicit formulas | eigensolver-internal differentiation as an unmeasured default |
| Fixed quadrature | analytical linear rule | generic reverse tape |
| Adaptive quadrature | analytical derivative under the integral, frozen trace | refinement-logic differentiation as an unmeasured default |
| ODE integration | forward sensitivity, discrete or continuous adjoint | one universal method |
| Event time | analytical implicit event equation | root-search differentiation as an unmeasured default |
| FFT | analytical transform and adjoint transform | butterfly differentiation as an unmeasured default |
| Spline evaluation | analytical basis derivative | autodiff through span search |
| Spline fitting | analytical implicit coefficient solve | iterative-solver differentiation as an unmeasured default |
| RNG seed to draw | no derivative | artificial seed derivative |
| Random estimator parameters | estimator-level gradient | differentiable RNG kernel |

## 15. Final policy

Derivative implementation is an automated performance search:

1. Identify the mathematical operator.
2. Determine the exact derivative product.
3. Apply the implicit function theorem wherever applicable.
4. Generate `analytical`, `autodiff`, and `hybrid` candidates.
5. Validate every candidate independently.
6. Benchmark runtime and peak memory at realistic scales.
7. Select the measured winner for each workload and target.
8. Repeat when the primal, compiler, hardware, or workload changes.

`fortnum` is primal-first, derivative-plural, and benchmark-selected.

The target architecture combines analytical structure, implicit
differentiation, `fortsym`-generated kernels, compiler autodiff candidates, and
measured selection. Convenience does not choose the production derivative.
Evidence does.
