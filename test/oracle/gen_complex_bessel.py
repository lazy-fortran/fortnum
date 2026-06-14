#!/usr/bin/env python3
"""Regenerate the complex-Bessel oracle table under test/oracle/data/.

Reference for fortnum_special_complex_bessel: J_n(z), I_n(z), K_n(z) of
complex z, integer order. scipy.special.jv / iv / ive / kv / kve supply the
values; the Fortran oracle test reads this CSV and asserts agreement within a
documented tolerance, so fortnum never links AMOS.

The grid spans the domain KiLCA (KAMEL) drives the replaced zbesj/zbesi/zbesk
over: complex z = gamma*r at moderate |z| (hom_medium), z on the imaginary
axis (besseli via J(i z)), and real z >= 0 with KODE=2 scaling out to large
Re z (flre conductivity).

Run: python3 test/oracle/gen_complex_bessel.py
"""
import os

import numpy as np
import scipy.special as sps

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

# func codes in the CSV
F_J = 0       # J_n(z)               unscaled
F_I = 1       # I_n(z)               unscaled
F_I_SC = 2    # e^{-Re z} I_n(z)     scaled  (AMOS KODE=2)
F_K = 3       # K_n(z)               unscaled
F_K_SC = 4    # e^{z} K_n(z)         scaled  (AMOS KODE=2)


def _grid_moderate():
    # complex z = r e^{i theta}, moderate |z| as in hom_medium (Re z may be
    # either sign for J/I; K requires Re z > 0 and is filtered below).
    radii = [0.1, 0.5, 1.0, 2.0, 3.5, 5.0, 7.0, 9.0, 11.0, 13.0]
    angles = [0.0, 0.4, -0.4, 0.9, -0.9, 1.3, -1.3]
    pts = []
    for r in radii:
        for a in angles:
            pts.append(complex(r * np.cos(a), r * np.sin(a)))
    # imaginary axis (besseli wrapper feeds z = i*zarg)
    for t in [0.2, 0.5, 1.0, 2.0, 4.0, 7.0, 10.0]:
        pts.append(complex(0.0, t))
        pts.append(complex(0.0, -t))
    return pts


def _grid_large_j():
    # Hankel-asymptotic regime for J, away from the negative real axis.
    radii = [16.0, 20.0, 30.0, 50.0]
    angles = [0.0, 0.5, -0.5, 1.0, -1.0]
    return [complex(r * np.cos(a), r * np.sin(a)) for r in radii for a in angles]


def _grid_real_scaled():
    # flre conductivity: x2 = (ks vT/omc)^2 real >= 0, KODE=2, high order.
    return [complex(x, 0.0) for x in
            [0.05, 0.5, 2.0, 8.0, 25.0, 80.0, 200.0, 500.0]]


def gen_complex_bessel():
    orders = [0, 1, 2, 3, 5]
    rows = []
    idx = 0

    def emit(func, n, z, val):
        nonlocal idx
        if not (np.isfinite(val.real) and np.isfinite(val.imag)):
            return
        rows.append((idx, func, n, z.real, z.imag, val.real, val.imag))
        idx += 1

    mod = _grid_moderate()
    for n in orders:
        for z in mod:
            emit(F_J, n, z, complex(sps.jv(n, z)))
            emit(F_I, n, z, complex(sps.iv(n, z)))
            if z.real > 0.0:
                emit(F_K, n, z, complex(sps.kv(n, z)))
                emit(F_K_SC, n, z, complex(sps.kve(n, z)))
    for n in orders:
        for z in _grid_large_j():
            emit(F_J, n, z, complex(sps.jv(n, z)))
    for n in orders:
        for x in _grid_real_scaled():
            emit(F_I_SC, n, x, complex(sps.ive(n, x)))

    path = os.path.join(DATA, "complex_bessel.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: bessel_j_complex / bessel_i_complex / "
                "bessel_k_complex\n")
        f.write("# backend: scipy.special jv / iv / ive / kv / kve\n")
        f.write("# columns: index,func,n,re_z,im_z,re_expected,im_expected\n")
        f.write("# func: 0=J 1=I 2=I_scaled(KODE2) 3=K 4=K_scaled(KODE2)\n")
        for row in rows:
            f.write("%d,%d,%d,%.17e,%.17e,%.17e,%.17e\n" % row)
    return len(rows)


if __name__ == "__main__":
    print(f"complex_bessel: {gen_complex_bessel()} rows")
