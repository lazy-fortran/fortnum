#!/usr/bin/env python3
"""B-spline basis oracle generator for fortnum (fortnum_bspline).

Emits reference basis-function values and derivatives that the Fortran oracle
test reads back. The reference is scipy.interpolate.BSpline: each basis
function B_{i,k} is recovered by evaluating a spline whose coefficient vector is
the i-th unit vector, which is exactly splev/BSpline with unit coefficients.

This reproduces the clamped-knot convention NEO-2 uses (collop_bspline):
spline order k = degree + 1, nbreak breakpoints, an augmented
clamped knot vector with end multiplicity k, and ncoef = nbreak + k - 2 basis
functions.

The breakpoint sets mirror the NEO-2 collop domain: a geometrically stretched
grid on [0, phi_x_max] (collop_bspline_dist) plus a uniform control case, over
the orders collop_bspline uses (k = 3, 4, 5).

Output CSV (one file, data/bspline.csv). Header lines start with '#':

    # fortnum bspline oracle table
    # backend: scipy.interpolate.BSpline
    # columns: case,order,nbreak,nx,nbreak_values...,
    # then per case a knots block and per (x,deriv) rows
    # row format: case,order,nbreak,x,deriv,ncoef,v_0,...,v_{ncoef-1}

A self-describing flat layout keeps the Fortran reader simple: every data row is
case,order,nbreak,x,deriv,ncoef followed by ncoef basis values. The breakpoints
are recoverable in Fortran from order+nbreak+domain via the same geometric rule,
so the test rebuilds the workspace independently and only the basis values are
compared. To keep the test self-contained, breakpoints are also emitted in
dedicated 'brk' rows: brk,case,order,nbreak,b_0,...,b_{nbreak-1}.
"""

from __future__ import annotations

import os

import numpy as np
from scipy.interpolate import BSpline

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def clamped_knots(breakpts: np.ndarray, order: int) -> np.ndarray:
    """Augmented clamped knot vector: end breakpoints with multiplicity order."""
    k = order
    interior = breakpts[1:-1]
    left = np.repeat(breakpts[0], k)
    right = np.repeat(breakpts[-1], k)
    return np.concatenate([left, interior, right])


def geom_breakpoints(nbreak: int, xmax: float, dist: float) -> np.ndarray:
    """NEO-2 collop geometric breakpoint stretch on [0, xmax] (init_phi_bspline)."""
    gam_all = sum(dist**kk for kk in range(1, nbreak))
    b = np.zeros(nbreak)
    x_del = xmax / gam_all
    for kk in range(1, nbreak):
        b[kk] = b[kk - 1] + x_del * dist**kk
    b[-1] = xmax
    return b


def basis_table(breakpts: np.ndarray, order: int, xs, max_deriv: int):
    """Return list of (x, deriv, values[ncoef]) for the basis over breakpts."""
    k = order
    degree = k - 1
    t = clamped_knots(breakpts, order)
    ncoef = len(t) - k
    rows = []
    for x in xs:
        for d in range(max_deriv + 1):
            vals = np.zeros(ncoef)
            if d > degree:
                # A degree-p spline has identically zero derivatives above p.
                rows.append((x, d, vals))
                continue
            for i in range(ncoef):
                c = np.zeros(ncoef)
                c[i] = 1.0
                spl = BSpline(t, c, degree, extrapolate=False)
                if d == 0:
                    vals[i] = spl(x)
                else:
                    vals[i] = spl.derivative(d)(x)
            # BSpline returns nan just outside support with extrapolate=False;
            # the test x grid stays strictly inside [b0, b_last].
            rows.append((x, d, vals))
    return ncoef, rows


def make_cases():
    cases = []
    # Case 0: NEO-2-like geometric stretch, order 3 (collop_bspline_order=3).
    # phi_x_max and collop_bspline_dist representative values; lagmax-order+3
    # breakpoints. Use lagmax=8 -> nbreak = 8 - 3 + 3 = 8.
    b0 = geom_breakpoints(8, 5.0, 1.3)
    cases.append((3, b0))
    # Case 1: geometric stretch, order 4, nbreak = lagmax-order+3 with lagmax=9.
    b1 = geom_breakpoints(8, 5.0, 1.2)
    cases.append((4, b1))
    # Case 2: geometric stretch, order 5, nbreak = 7.
    b2 = geom_breakpoints(7, 4.0, 1.25)
    cases.append((5, b2))
    # Case 3: uniform breakpoints, order 4 (control / interior-multiplicity-1).
    b3 = np.linspace(0.0, 1.0, 6)
    cases.append((4, b3))
    return cases


def eval_grid(breakpts: np.ndarray, n: int) -> np.ndarray:
    """Interior evaluation points, strictly inside [b0, b_last], avoiding knots."""
    b0, bn = breakpts[0], breakpts[-1]
    # Offset off the endpoints and off interior knots by an irrational fraction.
    frac = np.linspace(0.0, 1.0, n + 2)[1:-1]
    return b0 + (bn - b0) * frac


def write_table(path: str, cases):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("# fortnum bspline oracle table\n")
        fh.write("# backend: scipy.interpolate.BSpline\n")
        fh.write("# brk row:  brk,case,order,nbreak,b_0,...,b_{nbreak-1}\n")
        fh.write("# val row:  val,case,order,nbreak,x,deriv,ncoef,v_0,...,v_{ncoef-1}\n")
        for ci, (order, breakpts) in enumerate(cases):
            nbreak = len(breakpts)
            brk = ",".join(repr(float(b)) for b in breakpts)
            fh.write(f"brk,{ci},{order},{nbreak},{brk}\n")
            xs = eval_grid(breakpts, 5)
            ncoef, rows = basis_table(breakpts, order, xs, max_deriv=order)
            for (x, d, vals) in rows:
                vv = ",".join(repr(float(v)) for v in vals)
                fh.write(
                    f"val,{ci},{order},{nbreak},{float(x)!r},{d},{ncoef},{vv}\n"
                )


if __name__ == "__main__":
    out = os.path.join(DATA_DIR, "bspline.csv")
    write_table(out, make_cases())
    print(f"wrote {out}")
