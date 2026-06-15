#!/usr/bin/env python3
"""Scipy DOP853 oracle generator for the fortnum RK8(7)13M integrator.

Produces reference CSV files under test/oracle/data/dop853_<case>.csv. Each
file records (t, y1, y2, ...) at fixed output points from
scipy.integrate.solve_ivp(method='DOP853') at very tight tolerance, the same
method family the fortnum integrator implements. The fortnum
ode_integrate_dop run in the oracle test integrates the same IVPs and must
match these points.

Cases mirror the consumers KiLCA drives through its RK8(7) stepper: a tight
exponential decay (background-style stiffly tight tolerance), exponential
growth, and the harmonic oscillator (flow-style oscillatory transport). All
have closed-form solutions, so the test also checks the analytic value.

decay      y' = -y,        y(0) = 1,     t in [0, 5]
growth     y' = y,         y(0) = 1,     t in [0, 3]
harmonic   y1' = y2,       y(0) = [1,0], t in [0, 6*pi]
           y2' = -y1
"""

from __future__ import annotations

import math
import os

try:
    from scipy.integrate import solve_ivp
    import numpy as np
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# Very tight ground truth. The KiLCA background consumer runs rk8pd at
# eps_abs=eps_rel=1e-16; scipy clamps atol/rtol above machine epsilon, so the
# tightest reproducible ground truth is a few ulp above eps.
SCIPY_RTOL = 1e-14
SCIPY_ATOL = 1e-16
SCIPY_METHOD = "DOP853"


def _write_csv(name: str, t_eval, y_ref) -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, f"dop853_{name}.csv")
    neq = y_ref.shape[0]
    cols = ["t"] + [f"y{k+1}" for k in range(neq)]
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("# fortnum RK8(7)13M oracle table\n")
        fh.write(f"# case: {name}\n")
        fh.write(f"# backend: scipy.integrate.solve_ivp ({SCIPY_METHOD})\n")
        fh.write(f"# scipy_rtol: {SCIPY_RTOL!r}\n")
        fh.write(f"# scipy_atol: {SCIPY_ATOL!r}\n")
        fh.write(f"# columns: {','.join(cols)}\n")
        for i, ti in enumerate(t_eval):
            row = [repr(float(ti))] + [repr(float(y_ref[k, i])) for k in range(neq)]
            fh.write(",".join(row) + "\n")
    return path


def _solve(rhs, t_span, y0, t_eval):
    sol = solve_ivp(
        rhs,
        t_span,
        y0,
        method=SCIPY_METHOD,
        t_eval=t_eval,
        rtol=SCIPY_RTOL,
        atol=SCIPY_ATOL,
        dense_output=False,
    )
    if not sol.success:
        raise RuntimeError(f"scipy solve_ivp failed: {sol.message}")
    return sol.y


def generate_decay():
    t_eval = np.linspace(0.0, 5.0, 21)
    y_ref = _solve(lambda t, y: [-y[0]], (0.0, 5.0), [1.0], t_eval)
    return _write_csv("decay", t_eval, y_ref)


def generate_growth():
    t_eval = np.linspace(0.0, 3.0, 16)
    y_ref = _solve(lambda t, y: [y[0]], (0.0, 3.0), [1.0], t_eval)
    return _write_csv("growth", t_eval, y_ref)


def generate_harmonic():
    t_eval = np.linspace(0.0, 6.0 * math.pi, 31)
    y_ref = _solve(
        lambda t, y: [y[1], -y[0]],
        (0.0, 6.0 * math.pi),
        [1.0, 0.0],
        t_eval,
    )
    return _write_csv("harmonic", t_eval, y_ref)


def generate_all() -> list[str]:
    if not _HAS_SCIPY:
        raise ImportError("scipy is required to generate DOP853 oracle tables")
    return [
        generate_decay(),
        generate_growth(),
        generate_harmonic(),
    ]


if __name__ == "__main__":
    paths = generate_all()
    for p in paths:
        print(f"wrote {p}")
