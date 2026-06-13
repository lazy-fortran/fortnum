#!/usr/bin/env python3
"""Regenerate adaptive-integration oracle tables under test/oracle/data/.

scipy.integrate.quad IS QUADPACK (the same dqag/dqags this module reimplements),
so it is the reference for the QAG/QAGS primal. Each case has an analytic value;
the CSV carries the scipy value, and the Fortran oracle test asserts agreement
within the requested relative tolerance.

Run: python3 test/oracle/gen_integrate_tables.py
"""
import math
import os

import numpy as np
import scipy.integrate as si

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


def _path(name):
    return os.path.join(DATA, name)


# case_id, mode (0=QAG, 1=QAGS), a, b, integrand, label.
# Integrand ids match the Fortran dispatch in the oracle test.
QAG_CASES = [
    (0, "3*x**2 + 2*x + 1", lambda x: 3.0 * x**2 + 2.0 * x + 1.0, 0.0, 2.0),
    (1, "exp(x)", lambda x: math.exp(x), 0.0, 1.0),
    (2, "cos(x)", lambda x: math.cos(x), 0.0, 0.5 * math.pi),
    (3, "sin(10*x) + 2", lambda x: math.sin(10.0 * x) + 2.0, 0.0, math.pi),
    (4, "1/(1+x**2)", lambda x: 1.0 / (1.0 + x * x), 0.0, 1.0),
    (5, "exp(-x**2)", lambda x: math.exp(-x * x), -2.0, 2.0),
]

# QAGS: integrable endpoint singularities and sharply peaked integrands. The
# lower bound is the true singular endpoint a=0; Gauss-Kronrod nodes are
# strictly interior, so neither scipy nor fortnum ever evaluates f at 0, and
# both extrapolate the same singularity. This is the genuine QAGS regime.
QAGS_CASES = [
    (0, "1/sqrt(x)", lambda x: 1.0 / math.sqrt(x), 0.0, 1.0),
    (1, "log(x)", lambda x: math.log(x), 0.0, 1.0),
    (2, "lorentzian w=1e-2 c=0.5",
     lambda x: 1.0 / (1.0 + ((x - 0.5) / 1.0e-2) ** 2), 0.0, 1.0),
    (3, "lorentzian w=1e-3 c=0.3",
     lambda x: 1.0 / (1.0 + ((x - 0.3) / 1.0e-3) ** 2), 0.0, 1.0),
    (4, "x**(-0.5)*exp(-x)", lambda x: x ** (-0.5) * math.exp(-x), 0.0, 1.0),
]


def _emit(name, header_func, cases):
    rows = []
    for cid, label, fn, a, b in cases:
        val, _ = si.quad(fn, a, b, epsabs=0.0, epsrel=1.0e-12, limit=500)
        rows.append((cid, a, b, val, label))
    with open(_path(name), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: %s\n" % header_func)
        f.write("# backend: scipy.integrate.quad (QUADPACK)\n")
        f.write("# columns: case_id,a,b,expected\n")
        f.write("# has_derivative: 0\n")
        for cid, a, b, val, label in rows:
            f.write("# case %d: %s\n" % (cid, label))
        for cid, a, b, val, _label in rows:
            f.write("%d,%.17e,%.17e,%.17e\n" % (cid, a, b, val))


def main():
    _emit("integrate_qag_smooth.csv", "integrate_qag", QAG_CASES)
    _emit("integrate_qags_singular.csv", "integrate_qags", QAGS_CASES)
    print("wrote integrate_qag_smooth.csv, integrate_qags_singular.csv")


if __name__ == "__main__":
    main()
