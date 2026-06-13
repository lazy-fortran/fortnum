#!/usr/bin/env python3
"""Scipy oracle generator for fortnum ODE tests.

Produces reference CSV files under test/oracle/data/ode_<case>.csv.
Each file records (t, y1, y2, ...) at fixed output points, computed by
scipy.integrate.solve_ivp at tight tolerance (rtol=1e-12, atol=1e-12).

Cases
-----
decay       y' = -y,           y(0) = 1,     t in [0, 5]
growth      y' = y,            y(0) = 1,     t in [0, 3]
harmonic    y1' = y2,          y(0) = [1,0], t in [0, 6*pi]
            y2' = -y1
linear      y1' = -2*y1+y2,   y(0) = [1,0], t in [0, 2]
            y2' =  y1 - 3*y2

Derivative extension note
-------------------------
Derivative reference tables for ode_at_jvp / ode_at_vjp belong under
issue #40. Do not add derivative columns here by mechanically differentiating
the solution: the sensitivity computation requires differentiating through the
integrator, not through the output values. Add a separate registered case in
a new function (e.g. ode_decay_jvp) when #40 is implemented.
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

# Tight tolerance for scipy ground truth. Cash-Karp and DOP853 use different
# tableaux; comparing fortnum (loose tol) against scipy (tight tol) is the
# correct oracle strategy -- see docs/design/ode.md §5.
SCIPY_RTOL = 1e-12
SCIPY_ATOL = 1e-12
SCIPY_METHOD = "DOP853"


def _write_csv(name: str, t_eval, y_ref) -> str:
    """Write (t, y1, ...) rows to data/ode_<name>.csv."""
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, f"ode_{name}.csv")
    neq = y_ref.shape[0]
    cols = ["t"] + [f"y{k+1}" for k in range(neq)]
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("# fortnum ODE oracle table\n")
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
    return sol.y   # shape (neq, n_eval)


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


def generate_linear():
    # y1' = -2*y1 + y2,  y2' = y1 - 3*y2;  exact: exp of a 2x2 matrix
    t_eval = np.linspace(0.0, 2.0, 21)
    y_ref = _solve(
        lambda t, y: [-2.0 * y[0] + y[1], y[0] - 3.0 * y[1]],
        (0.0, 2.0),
        [1.0, 0.0],
        t_eval,
    )
    return _write_csv("linear", t_eval, y_ref)


def generate_all() -> list[str]:
    if not _HAS_SCIPY:
        raise ImportError("scipy is required to generate ODE oracle tables")
    return [
        generate_decay(),
        generate_growth(),
        generate_harmonic(),
        generate_linear(),
    ]


if __name__ == "__main__":
    paths = generate_all()
    for p in paths:
        print(f"wrote {p}")
