#!/usr/bin/env python3
"""Reference ("oracle") table generator for fortnum.

Emits high-precision reference tables that the Fortran test suite reads back
to verify ported functions against an independent implementation. NumPy and
SciPy are the preferred backends; the generator falls back to mpmath or the
standard-library ``math`` module so the seed tables can be produced on a
machine that has neither installed.

Adding a function is a registry entry: name a callable for the primal value,
optionally one for the derivative, and an input grid. See ``_register_seed``
for the worked gamma example.

Output format (CSV, one file per registered function under data/):

    # fortnum oracle table
    # function: <name>
    # backend: <which library produced the values>
    # columns: index,x,primal,derivative
    # has_derivative: <0|1>
    0,<x0>,<primal0>,<deriv0>
    ...

Values are written with full float64 round-trip precision (repr). The
``derivative`` column is always present so derivative reference tables
(issue #40) can be added without changing the file format or the Fortran
reader. When no derivative is registered the column holds 0 and the
``has_derivative`` header flag is 0, telling consumers to ignore it.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from typing import Callable, Optional, Sequence

# Optional backends. The generator works with whatever is available and
# records which one produced each table in the file header.
try:
    import numpy as _np
except ImportError:  # pragma: no cover - exercised only on minimal hosts
    _np = None

try:
    import scipy.special as _sps
except ImportError:  # pragma: no cover
    _sps = None

try:
    import mpmath as _mp
except ImportError:  # pragma: no cover
    _mp = None


DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@dataclass(frozen=True)
class OracleEntry:
    """One registered reference function.

    primal/derivative take a single float and return a float. derivative may
    be None when no analytic derivative is registered yet; the data column is
    then filled with zeros and flagged unused.
    """

    name: str
    grid: Sequence[float]
    primal: Callable[[float], float]
    backend: str
    derivative: Optional[Callable[[float], float]] = None


_REGISTRY: list[OracleEntry] = []


def register(entry: OracleEntry) -> None:
    """Add a reference function. Names must be unique."""
    if any(e.name == entry.name for e in _REGISTRY):
        raise ValueError(f"duplicate oracle entry: {entry.name}")
    _REGISTRY.append(entry)


def _write_table(entry: OracleEntry) -> str:
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, f"{entry.name}.csv")
    has_deriv = 1 if entry.derivative is not None else 0
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("# fortnum oracle table\n")
        fh.write(f"# function: {entry.name}\n")
        fh.write(f"# backend: {entry.backend}\n")
        fh.write("# columns: index,x,primal,derivative\n")
        fh.write(f"# has_derivative: {has_deriv}\n")
        for i, x in enumerate(entry.grid):
            xf = float(x)
            p = float(entry.primal(xf))
            d = float(entry.derivative(xf)) if entry.derivative is not None else 0.0
            # repr() round-trips float64 exactly; the Fortran reader parses
            # these with list-directed input into real64.
            fh.write(f"{i},{xf!r},{p!r},{d!r}\n")
    return path


def generate_all() -> list[str]:
    """Produce every registered table; return the written file paths."""
    return [_write_table(e) for e in _REGISTRY]


# --- Seed registration -------------------------------------------------------
#
# Complete gamma function. fortnum's special-function port (issue family
# #40) will expose a primal gamma; this seed lets the oracle harness be
# validated end-to-end today against gfortran's intrinsic gamma().


def _gamma_backend():
    """Pick the most authoritative gamma available and return (fn, label)."""
    if _sps is not None:
        return (lambda x: float(_sps.gamma(x)), "scipy.special.gamma")
    if _mp is not None:
        return (lambda x: float(_mp.gamma(x)), "mpmath.gamma")
    return (math.gamma, "math.gamma")


def _register_seed() -> None:
    gamma_fn, gamma_label = _gamma_backend()
    # Positive arguments only: the gamma poles at non-positive integers make
    # the negative axis a poor first verification grid, and the seeded primal
    # (gfortran gamma) agrees most cleanly on x > 0.
    grid = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 10.0]
    register(
        OracleEntry(
            name="gamma",
            grid=grid,
            primal=gamma_fn,
            backend=gamma_label,
            # Derivative gamma'(x) = gamma(x) * digamma(x); registered when a
            # digamma backend exists so #40 inherits a populated column.
            derivative=_gamma_derivative(gamma_fn),
        )
    )


def _gamma_derivative(gamma_fn) -> Optional[Callable[[float], float]]:
    if _sps is not None:
        return lambda x: float(_sps.gamma(x) * _sps.digamma(x))
    if _mp is not None:
        return lambda x: float(_mp.gamma(x) * _mp.digamma(x))
    return None  # no stdlib digamma; leave the derivative column unused


_register_seed()


if __name__ == "__main__":
    written = generate_all()
    for p in written:
        print(f"wrote {p}")
