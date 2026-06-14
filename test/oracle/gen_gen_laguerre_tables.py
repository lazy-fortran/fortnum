#!/usr/bin/env python3
"""Regenerate the generalized Gauss-Laguerre oracle table under data/.

fortnum_quadrature.gauss_gen_laguerre produces nodes/weights for the weight
w(x) = x^alpha exp(-x) on [0, inf), alpha > -1, via Golub-Welsch on the Jacobi
matrix (DLMF 18.9 recurrence, DLMF Table 18.3.1 moment).  The reference here is
scipy.special.roots_genlaguerre, an independent Golub-Welsch implementation, so
fortnum links no LAPACK: the verified rule lives in the data file.

Covered: n = 32 at alpha = 5/2 and 7/2, the orders/alphas libneo transport.f90
calc_D_one_over_nu uses.  The generator confirms each rule integrates the
monomials int_0^inf x^(alpha+k) exp(-x) dx = Gamma(alpha+k+1) for k = 0..2n-1
and the zeroth moment mu0 = Gamma(alpha+1) before writing the table, and aborts
if scipy's own rule misses that exactness.

Run: python3 test/oracle/gen_gen_laguerre_tables.py
"""
import math
import os

import numpy as np
import scipy.special as sps

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

# (alpha_id, alpha): the transport.f90 consumer orders.
CASES = [(0, 2.5), (1, 3.5)]
N = 32

# scipy reaches ~3e-15 relative on the worst monomial at n=32; assert a small
# multiple so the float64 Fortran QL rule is held to an honest target.
MONOMIAL_TOL = 1.0e-13


def _verify(alpha, nodes, weights):
    """Confirm scipy's rule integrates monomials and the zeroth moment."""
    mu0 = math.gamma(alpha + 1.0)
    if abs(float(np.sum(weights)) - mu0) > MONOMIAL_TOL * mu0:
        raise SystemExit("scipy zeroth moment off for alpha=%g" % alpha)
    for k in range(2 * N):
        approx = float(np.sum(weights * nodes ** k))
        exact = math.gamma(alpha + k + 1.0)
        if abs(approx - exact) > MONOMIAL_TOL * exact:
            raise SystemExit(
                "scipy monomial k=%d off for alpha=%g (rel %.3e)"
                % (k, alpha, abs(approx - exact) / exact))


def main():
    path = os.path.join(DATA, "gauss_gen_laguerre.csv")
    rows = 0
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: gauss_gen_laguerre\n")
        f.write("# backend: scipy.special.roots_genlaguerre\n")
        f.write("# columns: alpha_id,alpha,n,i,node,weight\n")
        f.write("# note: weight x^alpha exp(-x) on [0,inf); nodes ascending\n")
        for aid, alpha in CASES:
            nodes, weights = sps.roots_genlaguerre(N, alpha)
            _verify(alpha, nodes, weights)
            for i in range(N):
                f.write("%d,%s,%d,%d,%s,%s\n" % (
                    aid, repr(float(alpha)), N, i + 1,
                    repr(float(nodes[i])), repr(float(weights[i]))))
                rows += 1
    print("wrote gauss_gen_laguerre.csv (%d rows)" % rows)


if __name__ == "__main__":
    main()
