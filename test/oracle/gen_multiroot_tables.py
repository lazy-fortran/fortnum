#!/usr/bin/env python3
"""Regenerate the multiroot / deriv_central / argsort oracle tables.

Three reference files under test/oracle/data/:

  multiroot.csv      converged roots of standard n-dim test systems from
                     scipy.optimize.root(method='hybr') (MINPACK hybrj).
  deriv_central.csv  analytic first derivatives of smooth scalar functions,
                     the reference for the central finite-difference helper.
  argsort.csv        ascending index permutations from numpy.argsort.

Run: python3 test/oracle/gen_multiroot_tables.py
"""
import os

import numpy as np
from scipy.optimize import root

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


# --- n-dim root systems (index matches the Fortran case dispatch) ----------
#
# 0: Rosenbrock gradient = 0.
#    F = grad of 100(x2 - x1^2)^2 + (1 - x1)^2; root (1, 1).
def rosenbrock_grad(v):
    x1, x2 = v
    return [
        -400.0 * x1 * (x2 - x1 * x1) - 2.0 * (1.0 - x1),
        200.0 * (x2 - x1 * x1),
    ]


# 1: Powell singular function (More-Garbow-Hillstrom #13), root at origin.
def powell_singular(v):
    x1, x2, x3, x4 = v
    return [
        x1 + 10.0 * x2,
        np.sqrt(5.0) * (x3 - x4),
        (x2 - 2.0 * x3) ** 2,
        np.sqrt(10.0) * (x1 - x4) ** 2,
    ]


# 2: 2x2 linear-ish smooth system with a unique root (KiLCA n=2 shape).
#    F1 = x1^2 + x2^2 - 2, F2 = x1 - x2; root (1, 1).
def circle_line(v):
    x1, x2 = v
    return [x1 * x1 + x2 * x2 - 2.0, x1 - x2]


SYSTEMS = [
    (0, "rosenbrock_grad", rosenbrock_grad, [0.5, 0.5]),
    (1, "powell_singular", powell_singular, [3.0, -1.0, 0.0, 1.0]),
    (2, "circle_line", circle_line, [2.0, 0.5]),
]


def _write_multiroot():
    path = os.path.join(DATA, "multiroot.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: multidimensional roots F(x)=0\n")
        f.write("# backend: scipy.optimize.root(method='hybr')\n")
        f.write("# columns: index,n,x0...,root...\n")
        f.write("# has_derivative: 0\n")
        for idx, name, fn, x0 in SYSTEMS:
            sol = root(fn, x0, method="hybr", tol=1e-12)
            # Powell singular has a singular Jacobian at the root, so MINPACK
            # clears the residual but flags success=False; gate on the actual
            # residual instead of the convergence flag.
            res = float(np.max(np.abs(fn(sol.x))))
            assert res < 1e-9, f"{name}: residual {res} too large"
            n = len(x0)
            cols = [str(idx), str(n)]
            cols += [repr(float(v)) for v in x0]
            cols += [repr(float(v)) for v in sol.x]
            f.write(",".join(cols) + "\n")
    return len(SYSTEMS)


# --- central finite difference: analytic derivative reference --------------
#
# index, x, h, f'(x) analytic.  h chosen as a small fixed fraction (~1e-3)
# of a unit grid spacing, matching the KiLCA r-grid usage.
DERIV = [
    (0, 1.3, 1.0e-3, np.cos(1.3)),                 # f = sin x
    (1, 0.7, 1.0e-3, np.exp(0.7)),                 # f = exp x
    (2, 2.0, 1.0e-3, 3.0 * 2.0 * 2.0),             # f = x^3 -> 3x^2
    (3, 0.4, 1.0e-3, 1.0 / 0.4),                   # f = ln x -> 1/x
    (4, 1.1, 1.0e-3, -np.sin(1.1)),                # f = cos x
]


def _write_deriv():
    path = os.path.join(DATA, "deriv_central.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: central finite difference (deriv_central)\n")
        f.write("# backend: analytic derivatives (numpy)\n")
        f.write("# columns: index,x,h,deriv\n")
        f.write("# has_derivative: 1\n")
        for idx, x, h, d in DERIV:
            f.write(f"{idx},{float(x)!r},{float(h)!r},{float(d)!r}\n")
    return len(DERIV)


# --- argsort: ascending index permutation reference ------------------------
ARGSORT = [
    [3.0, 1.0, 2.0, 5.0, 4.0],
    [-1.5, -3.2, 0.0, 2.7, 1.1, -0.4],
    [10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0],
]


def _write_argsort():
    path = os.path.join(DATA, "argsort.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: ascending index sort (argsort)\n")
        f.write("# backend: numpy.argsort\n")
        f.write("# columns: index,n,x...,perm... (perm 1-based)\n")
        f.write("# has_derivative: 0\n")
        for idx, x in enumerate(ARGSORT):
            arr = np.asarray(x, dtype=float)
            perm = np.argsort(arr, kind="stable") + 1  # 1-based for Fortran
            cols = [str(idx), str(len(arr))]
            cols += [repr(float(v)) for v in arr]
            cols += [str(int(p)) for p in perm]
            f.write(",".join(cols) + "\n")
    return len(ARGSORT)


if __name__ == "__main__":
    n1 = _write_multiroot()
    n2 = _write_deriv()
    n3 = _write_argsort()
    print(f"multiroot.csv: {n1} systems")
    print(f"deriv_central.csv: {n2} cases")
    print(f"argsort.csv: {n3} cases")
