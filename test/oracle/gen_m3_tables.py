#!/usr/bin/env python3
"""Regenerate the M3 root-finding oracle tables under test/oracle/data/.

High-precision reference roots from mpmath (50 digits); the Fortran oracle
tests assert agreement within a documented tolerance, so fortnum carries the
reference in the data file rather than linking an external root finder.

Run: python3 test/oracle/gen_m3_tables.py
"""
import os

import mpmath

mpmath.mp.dps = 50

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

# bisection / Newton cases
ROOTS = [
    (0, 1.0, 2.0, lambda x: x**2 - 2, "x^2 - 2, root = sqrt(2)"),
    (1, 0.5, 1.0, lambda x: mpmath.cos(x) - x, "cos(x) - x, Dottie number"),
    (2, 1.0, 2.0, lambda x: x**3 - x - 2, "x^3 - x - 2"),
    (3, 0.0, 2.0, lambda x: mpmath.exp(x) - 3, "exp(x) - 3, root = ln(3)"),
    (4, 0.5, 1.5, lambda x: x**3 - 7 * x + 6, "x^3 - 7*x + 6, root = 1"),
]

# Brent cases (wider brackets, near-flat regions)
BRENT = [
    (0, 0.0, 1.0, lambda x: x * mpmath.exp(x) - 1),
    (1, 0.0, 2.0, lambda x: mpmath.atan(x) - mpmath.pi / 6),
    (2, 1.0, 2.0, lambda x: x**5 - x - 1),
    (3, 0.0, mpmath.pi / 2, lambda x: mpmath.sin(x) - mpmath.mpf("0.8")),
    (4, 1.0, 3.0, lambda x: x**2 - 3),
    (5, 1.0, 2.0, lambda x: x**3 + 4 * x**2 - 10),
]


def _write(name, header, cases):
    path = os.path.join(DATA, name)
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write(f"# function: {header}\n")
        f.write("# backend: mpmath.findroot (dps=50)\n")
        f.write("# columns: index,a,b,root\n")
        f.write("# has_derivative: 0\n")
        for idx, a, b, fn in cases:
            mid = (mpmath.mpf(a) + mpmath.mpf(b)) / 2
            root = mpmath.findroot(fn, mid)
            f.write(f"{idx},{float(a)},{float(b)},{mpmath.nstr(root, 22)}\n")
    return len(cases)


if __name__ == "__main__":
    n1 = _write("roots.csv", "scalar roots (bisection / Newton)",
                [(i, a, b, fn) for i, a, b, fn, _ in ROOTS])
    n2 = _write("roots_brent.csv", "scalar roots (Brent)", BRENT)
    print(f"roots.csv: {n1} cases")
    print(f"roots_brent.csv: {n2} cases")
