# Real S-wave Coulomb functions

`fortnum_special_coulomb:coulomb_fg0(eta,rho,f,g,fp,gp,ok)` returns the
DLMF regular and irregular functions, with derivatives in rho, for l=0.
The Wronskian convention is F'G-FG'=1; eta=0 gives sin(rho), cos(rho).
The supported finite domain is |eta| <= 30, 0 < rho <= 8,
|eta rho| <= 2. Outside it, ok is false and the outputs are zero.
Success means converged floating evaluation, not an enclosure. No general
large-radius or higher-angular-momentum capability is claimed.

For y''+(1-2 eta/rho)y=0, write

    F = C rho sum a_n rho^n,
    G = [sum b_n rho^n + 2 eta rho log(2 rho) sum a_n rho^n]/C,
    C^2 = 2 pi eta/[exp(2 pi eta)-1].

Here a_0=b_0=1, a_1=eta,
b_1=2 eta [Re psi(1+i eta)+2 EulerGamma-1]. Equating powers gives

    a_n = (2 eta a_(n-1)-a_(n-2))/[n(n+1)],
    b_n = [2 eta b_(n-1)-b_(n-2)-2 eta(2n-1)a_(n-1)]/[n(n-1)].

The recurrence is retained in `tools/codegen/app/gen_coulomb_series.f90`.
Its checked output uses binary128 accumulation. Digamma comes from the
existing FortNum implementation; binary64 input and final rounding dominate.
Two consecutive terms must be small: a single zero term at eta=0 is not
convergence. The iteration cap returns failure. Independent mpmath references
cover both signs, the neutral limit, small rho, and barrier parameters to 30.
The oracle tolerance is 2e-10 scaled absolute/relative; it is not a proved bound.

Definitions: [NIST DLMF 33.2](https://dlmf.nist.gov/33.2) and
[33.6](https://dlmf.nist.gov/33.6). Digamma's existing analytical products
remain available; full eta products for this bounded Coulomb evaluator are
currently unavailable. fp/gp provide the rho derivative without differencing.

Regenerate using the [codegen locks](../tools/codegen/fortsym.lock) and
[checkout contract](../tools/codegen/checkout_fortsym.sh):

```sh
cd tools/codegen
fo exec gen_coulomb_series
cd ../..
python test/oracle/gen_coulomb.py  # optional oracle only, mpmath 1.4.1
fo test fortnum_coulomb
```

No CAS or Python runtime is required by the library or its registered test.
