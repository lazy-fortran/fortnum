#!/usr/bin/env python3
"""Regenerate the confluent hypergeometric oracle table under test/oracle/data/.

The reference for 1F1(a;b;z) (Kummer M) at complex b and z is mpmath.hyp1f1
evaluated at high working precision, cross-checked against scipy.special.hyp1f1
on the real axis. fortnum never links GSL/scipy at runtime: the Fortran oracle
test reads this CSV and asserts agreement within a documented tolerance.

The grid covers the (a=1, b=1+t2, z=t1) domain MEPHIT and KiLCA use, where
t1 = x1^2 and t2 = -i*x2 + t1 from the FLR / plasma-dispersion sums, plus
general a != 1 points and a negative-Re(z) block to exercise the Kummer
transformation branch and a large-|z| block to exercise the asymptotic branch.

Run: python3 test/oracle/gen_hyperg_tables.py
"""
import os

import mpmath as mp
import scipy.special as sps

mp.mp.dps = 50

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


def _rows():
    rows = []
    # Consumer domain: a = 1, b = 1 + t2, z = t1, with t1 = x1^2 and
    # t2 = -1j*x2 + t1. x1, x2 sweep small-to-moderate FLR arguments.
    x1_vals = [0.1, 0.3, 0.7, 1.0, 1.5, 2.0, 3.0]
    x2_vals = [0.0, 0.2, 0.5, 1.0, 2.0, 4.0]
    for x1 in x1_vals:
        for x2 in x2_vals:
            t1 = x1 * x1
            t2 = -1j * x2 + t1
            rows.append((1.0 + 0.0j, 1.0 + t2, t1 + 0.0j))

    # General a != 1, real and complex b, moderate z (Taylor branch).
    general = [
        (0.5 + 0.0j, 1.5 + 0.0j, 0.4 + 0.0j),
        (2.0 + 0.0j, 3.0 + 0.0j, 1.0 + 0.5j),
        (1.5 + 0.3j, 2.5 - 0.4j, 0.8 + 0.2j),
        (-0.5 + 0.0j, 2.0 + 0.0j, 0.6 + 0.0j),
        (1.0 + 0.0j, 0.5 + 0.0j, 1.2 + 0.0j),
    ]
    rows.extend(general)

    # Negative Re(z): exercises the Kummer transformation M=e^z M(b-a,b,-z).
    negz = [
        (1.0 + 0.0j, 2.0 + 0.0j, -1.0 + 0.0j),
        (1.0 + 0.0j, 1.5 + 0.5j, -2.0 + 1.0j),
        (0.5 + 0.0j, 2.0 + 0.0j, -3.0 + 0.0j),
        (1.0 + 0.0j, 2.0 + 1.0j, -5.0 + 2.0j),
    ]
    rows.extend(negz)

    # Large |z|, Re(z) > 0: exercises the asymptotic branch (|z| > 60).
    bigz = [
        (1.0 + 0.0j, 2.0 + 0.0j, 80.0 + 0.0j),
        (1.0 + 0.0j, 2.5 + 0.5j, 100.0 + 10.0j),
        (1.0 + 0.0j, 3.0 + 0.0j, 120.0 + 30.0j),
    ]
    rows.extend(bigz)
    return rows


def gen_hyperg():
    rows = _rows()
    path = os.path.join(DATA, "hyperg_1f1.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: hyperg_1f1 (Kummer M, complex a,b,z)\n")
        f.write("# backend: mpmath.hyp1f1 (dps=50), real-axis cross-check scipy.special.hyp1f1\n")
        f.write("# convention: M(a,b,z) = sum_k (a)_k/(b)_k z^k/k!  (DLMF 13.2.2)\n")
        f.write("# columns: index,a_re,a_im,b_re,b_im,z_re,z_im,m_re,m_im\n")
        for i, (a, b, z) in enumerate(rows):
            m = mp.hyp1f1(
                mp.mpc(a.real, a.imag),
                mp.mpc(b.real, b.imag),
                mp.mpc(z.real, z.imag),
            )
            m_re = float(m.real)
            m_im = float(m.imag)
            # Real-axis sanity cross-check against scipy where applicable.
            if a.imag == 0.0 and b.imag == 0.0 and z.imag == 0.0:
                ref = float(sps.hyp1f1(a.real, b.real, z.real))
                if ref != 0.0 and abs(ref - m_re) > 1e-9 * abs(ref):
                    raise SystemExit(
                        f"scipy/mpmath disagree at row {i}: {ref} vs {m_re}"
                    )
            f.write(
                f"{i},{a.real!r},{a.imag!r},{b.real!r},{b.imag!r},"
                f"{z.real!r},{z.imag!r},{m_re!r},{m_im!r}\n"
            )
    return len(rows)


if __name__ == "__main__":
    print(f"hyperg_1f1: {gen_hyperg()} rows")
