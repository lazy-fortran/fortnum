#!/usr/bin/env python3
"""Regenerate the CQUAD oracle table under test/oracle/data/cquad.csv.

scipy.integrate.quad IS QUADPACK; run at epsrel=1e-13 it delivers the true
integral to full double precision. The doubly-adaptive Clenshaw-Curtis rule in
fortnum_cquad converges to that same true value, so scipy is the reference.

Cases are finite-valued at the interval endpoints (Clenshaw-Curtis samples the
endpoints) and are representative of the collision-operator moments where NEO-2
uses CQUAD: Maxwellian-weighted polynomials and decaying oscillatory integrands.

Run: python3 test/oracle/gen_cquad_tables.py
"""
import math
import os

import scipy.integrate as si

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

# case_id, label, integrand, a, b. case_id matches the Fortran dispatch.
CASES = [
    (0, "3*x**2 + 2*x + 1", lambda x: 3.0 * x**2 + 2.0 * x + 1.0, 0.0, 2.0),
    (1, "exp(x)", lambda x: math.exp(x), 0.0, 1.0),
    (2, "1/(1+x**2)", lambda x: 1.0 / (1.0 + x * x), 0.0, 1.0),
    (3, "exp(-x**2)", lambda x: math.exp(-x * x), -2.0, 2.0),
    (4, "sqrt(x)", lambda x: math.sqrt(x), 0.0, 1.0),
    (5, "x**2*exp(-x**2)", lambda x: x * x * math.exp(-x * x), 0.0, 5.0),
    (6, "x**4*exp(-x**2)", lambda x: x**4 * math.exp(-x * x), 0.0, 6.0),
    (7, "x*exp(-x)", lambda x: x * math.exp(-x), 0.0, 10.0),
    (8, "cos(5*x)*exp(-x)", lambda x: math.cos(5.0 * x) * math.exp(-x), 0.0, 8.0),
]


def main():
    rows = []
    for cid, label, fn, a, b in CASES:
        val, _ = si.quad(fn, a, b, epsabs=0.0, epsrel=1.0e-13, limit=500)
        rows.append((cid, a, b, val, label))
    path = os.path.join(DATA, "cquad.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: integrate_cquad\n")
        f.write("# backend: scipy.integrate.quad (QUADPACK), epsrel=1e-13\n")
        f.write("# columns: case_id,a,b,expected\n")
        f.write("# has_derivative: 0\n")
        for cid, a, b, val, label in rows:
            f.write("# case %d: %s\n" % (cid, label))
        for cid, a, b, val, _label in rows:
            f.write("%d,%.17e,%.17e,%.17e\n" % (cid, a, b, val))
    print("wrote cquad.csv")


if __name__ == "__main__":
    main()
