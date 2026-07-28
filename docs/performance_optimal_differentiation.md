# Performance-optimal differentiation

## Objective

The production objective is:

\[
\boxed{\text{minimize validated complete-application wall clock and peak memory}}
\]

AI-assisted derivation and computer algebra reduce the human cost of explicit
derivatives. They do not make one derivative mechanism universally fastest.
Candidate selection remains a numerical performance problem.

An implementation may use:

- an explicit analytical expression
- a stable recurrence
- forward or reverse autodiff
- an implicit tangent or adjoint solve
- a frozen discrete trace
- a continuous sensitivity equation
- a generated contracted product
- a hybrid composition of these parts

## Product before mechanism

For \(f:\mathbb R^n\rightarrow\mathbb R^m\), the caller may need \(Jv\),
\(J^Tu\), a gradient, an HVP, or the full Jacobian. Their optimal algorithms
depend on different dimensions.

Forward work commonly scales with the number of input directions. Reverse work
commonly scales with the number of output cotangents and may need more retained
or recomputed primal state. Analytical products can follow either scaling
pattern. Finite differences usually scale with perturbed inputs or directions
and add truncation error.

Therefore, a library should expose products directly:

```text
value
jvp(direction)
vjp(cotangent)
gradient
hvp(direction)
```

Full matrices are optional products for consumers that reuse them.

## Mathematical operator boundaries

Differentiate the mathematical definition of a composite operator when that
definition gives a cheaper derivative than its executed algorithm.

For a residual-defined state,

\[
R(y,p)=0,
\]

the tangent and adjoint equations are:

\[
R_y\dot y=-R_p\dot p,
\qquad
R_y^T\lambda=L_y^T,
\qquad
\frac{dL}{dp}=L_p-\lambda^T R_p.
\]

These equations replace differentiation through nonlinear iterations. They
also reuse converged Jacobians, factorizations, and preconditioners.

The same principle applies to:

- linear systems
- fixed points and equilibria
- stationarity and KKT systems
- spline coefficient solves
- event times
- eigenproblems and invariant subspaces

Iteration-level autodiff remains a comparator when it is well defined.

## Candidate methods

### Generated explicit products

Generate contracted JVPs and VJPs for compact local expressions. Generate value
and derivative together when they share expensive intermediates. Preserve an
expression DAG and run common-subexpression elimination before code emission.

### Stable recurrences

Differentiate a recurrence as a recurrence. Expanding a length-\(K\) recurrence
into one symbolic expression can increase source size, register pressure,
compile time, and numerical error.

### Forward autodiff

Forward autodiff is a leading candidate for few input directions, streaming
algorithms, and low-memory tangent propagation. Blocked directions can expose
SIMD or matrix operations.

### Reverse autodiff

Reverse autodiff is a leading candidate for few scalar outputs. Its cost
includes trace storage, checkpointing, or recomputation. Custom operator rules
are necessary around solvers, transforms, and external libraries.

### Frozen traces

Adaptive code makes discrete decisions. A frozen-trace derivative records the
primal schedule and differentiates the resulting fixed map. It is piecewise
valid and differs from the derivative of changing adaptive logic.

### Finite differences

Finite differences provide an independent directional oracle and a black-box
fallback. Step-size error and repeated primal work usually prevent them from
winning exact production tournaments.

## Symbolic generation

Computer algebra can simplify local derivatives, generate several algebraic
forms, and prove symbolic equivalence. Native measurements still decide among
equivalent forms.

Generation must preserve:

- stable evaluation order
- scaling and overflow behavior
- solve structure rather than explicit inverses
- recurrence and sparsity structure
- compiler-visible loops and contiguous access
- bounded live temporaries

Operation count is a filter. Cache traffic, vectorization, register pressure,
instruction cache, and device occupancy can reverse an operation-count ranking.

## Hybrid composition

Useful hybrid boundaries follow the mathematics:

```text
global implicit or adjoint structure
    -> local residual JVP/VJP products
    -> custom rules for transforms and solves
    -> autodiff through remaining smooth composition
```

For example, a nonlinear objective may use a primal root solve, an analytical
implicit adjoint, `fortsym`-generated local products, and Enzyme for an outer
smooth objective. The complete candidate is `hybrid` because autodiff and an
analytical rule both participate.

## Performance tournament

For every product:

1. define the exact returned values
2. generate all admissible candidates
3. validate them independently
4. benchmark representative dimensions and reuse
5. measure wall clock and peak memory
6. inspect CPU/cache or GPU counters where needed
7. select outside hot loops
8. repeat after material source, compiler, hardware, or workload changes

The comparison unit is the complete caller-visible workload. A faster
derivative kernel can lose after value recomputation, data transfer,
factorization, optimizer convergence, or memory pressure.

## `fortnum` policy

`fortnum` combines:

- stable hand-written primal algorithms
- `fortsym` symbolic DAGs and Fortran generation
- Enzyme CPU candidates
- analytical implicit and frozen-trace products
- generated analytical GPU leaves
- machine-readable validation and benchmark selection

Evidence chooses each product and workload. Convenience and mechanism labels do
not.
