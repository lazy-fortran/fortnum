#!/usr/bin/env python3
"""Regenerate the complex analytic-zero-finder oracle table.

One reference file under test/oracle/data/:

  complex_roots.csv  distinct zeros and multiplicities of analytic functions
                     inside a rectangular box, plus the total winding number.

The Fortran finder (fortnum_roots_complex.complex_region_roots) evaluates the
same functions; the references here are computed independently:

  * polynomial cases: roots are assigned by construction, multiplicities are
    the assigned exponents; cross-checked against numpy.roots of the expanded
    polynomial (clustered to recover the multiple root).
  * sin(z): zeros are k*pi for integer k inside the box; cross-checked against
    mpmath.findroot from a nearby seed.

Each row:
  case_id, kind, ll_re, ll_im, ur_re, ur_im, ntotal, ndistinct,
  then ndistinct triples (root_re, root_im, mult), each root sorted by
  (real, imag).

Run: python3 test/oracle/gen_complex_roots_tables.py
"""
import os

import mpmath as mp
import numpy as np

mp.mp.dps = 40

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


def _verify_poly(roots_with_mult):
    """Expand the polynomial, recompute its roots, and confirm clustering."""
    coeffs = np.array([1.0 + 0.0j])
    flat = []
    for z, m in roots_with_mult:
        for _ in range(m):
            coeffs = np.convolve(coeffs, np.array([1.0 + 0.0j, -z]))
            flat.append(z)
    found = np.roots(coeffs)
    # Each assigned root must have m numpy roots within a cluster radius that
    # scales with the multiplicity: a root of multiplicity m perturbs by
    # O(eps^(1/m)) under the unavoidable rounding of np.roots.
    for z, m in roots_with_mult:
        radius = 1e-6 if m == 1 else 5e-4
        close = np.sum(np.abs(found - z) < radius)
        assert close == m, f"poly root {z} mult {m}: numpy found {close}"


def _verify_sin(zeros):
    for z in zeros:
        r = mp.findroot(mp.sin, mp.mpc(z.real + 0.3, z.imag + 0.2))
        assert abs(complex(r) - z) < 1e-12, f"sin zero {z}: mpmath {r}"


def _sorted_rm(rm):
    return sorted(rm, key=lambda t: (round(t[0].real, 9), round(t[0].imag, 9)))


def _row(case_id, kind, ll, ur, rm):
    rm = _sorted_rm(rm)
    ntot = sum(m for _, m in rm)
    nd = len(rm)
    cols = [str(case_id), kind,
            repr(ll.real), repr(ll.imag), repr(ur.real), repr(ur.imag),
            str(ntot), str(nd)]
    for z, m in rm:
        cols += [repr(float(z.real)), repr(float(z.imag)), str(m)]
    return ",".join(cols)


def main():
    rows = []

    # Case 0: simple polynomial, two real-ish simple roots well separated.
    rm = [(complex(0.5, 0.5), 1), (complex(-0.6, -0.3), 1)]
    _verify_poly(rm)
    rows.append(_row(0, "poly_two_simple",
                     complex(-2.0, -2.0), complex(2.0, 2.0), rm))

    # Case 1: a double root and a simple root (multiplicity recovery).
    rm = [(complex(0.4, 0.7), 2), (complex(-0.8, 0.1), 1)]
    _verify_poly(rm)
    rows.append(_row(1, "poly_double_simple",
                     complex(-2.0, -2.0), complex(2.0, 2.0), rm))

    # Case 2: a triple root alone (multiplicity 3).
    rm = [(complex(-0.3, 0.6), 3)]
    _verify_poly(rm)
    rows.append(_row(2, "poly_triple",
                     complex(-2.0, -2.0), complex(2.0, 2.0), rm))

    # Case 3: TWO NEARBY SIMPLE ROOTS plus a double root in the same box.
    # The pair sits 0.3 apart; the double root is elsewhere in the box.
    rm = [(complex(0.20, 0.10), 1), (complex(0.50, 0.10), 1),
          (complex(-0.9, -0.6), 2)]
    _verify_poly(rm)
    rows.append(_row(3, "poly_nearby_pair_and_double",
                     complex(-2.0, -2.0), complex(2.0, 2.0), rm))

    # Case 4: entire function sin(z) on a strip enclosing 0, pi, -pi.
    zeros = [complex(-float(mp.pi), 0.0), complex(0.0, 0.0),
             complex(float(mp.pi), 0.0)]
    _verify_sin(zeros)
    rm = [(z, 1) for z in zeros]
    rows.append(_row(4, "sin_strip",
                     complex(-4.0, -1.0), complex(4.0, 1.0), rm))

    # Case 5: sin(z) on a box around a single interior zero (pi), exercising
    # the independent winding-number count = 1.
    rm = [(complex(float(mp.pi), 0.0), 1)]
    _verify_sin([complex(float(mp.pi), 0.0)])
    rows.append(_row(5, "sin_single",
                     complex(2.0, -1.0), complex(4.5, 1.0), rm))

    path = os.path.join(DATA, "complex_roots.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: complex analytic-zero finder over a rectangle\n")
        f.write("# backend: numpy.roots (polynomials) / mpmath.findroot (sin)\n")
        f.write("# columns: case_id,kind,ll_re,ll_im,ur_re,ur_im,"
                "ntotal,ndistinct,(root_re,root_im,mult)*ndistinct\n")
        f.write("# has_derivative: 0\n")
        for r in rows:
            f.write(r + "\n")
    print(f"wrote {path}: {len(rows)} cases")


if __name__ == "__main__":
    main()
