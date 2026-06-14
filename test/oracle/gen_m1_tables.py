#!/usr/bin/env python3
"""Regenerate the M1 oracle reference tables under test/oracle/data/.

Each table is the scipy/numpy reference for one fortnum module. The Fortran
oracle tests read these CSVs and assert agreement within a documented
tolerance, so fortnum never links GSL: the reference lives in the data file.

Run: python3 test/oracle/gen_m1_tables.py
"""
import math
import os

import numpy as np
import scipy.integrate as si
import scipy.special as sps

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


def _path(name):
    return os.path.join(DATA, name)


def gen_bessel():
    # I_n = scipy.special.iv, K_n = scipy.special.kn. Underflow below 1e-280
    # is written as 0.0; overflow/non-finite rows are dropped.
    underflow = 1.0e-280
    orders = [0, 1, 2, 5, 10, 50, 100, 200]
    x_i = np.logspace(-10, np.log10(700.0), 200)
    x_k = np.logspace(-4, np.log10(700.0), 200)
    rows = []
    idx = 0
    for n in orders:
        for x in np.concatenate([x_i, -x_i]):
            val = float(sps.iv(n, x))
            if not np.isfinite(val):
                continue
            if abs(val) < underflow:
                val = 0.0
            rows.append((idx, 0, n, x, val))
            idx += 1
    for n in orders:
        for x in x_k:
            val = float(sps.kn(n, x))
            if not np.isfinite(val) or val < 0.0:
                continue
            if abs(val) < underflow:
                val = 0.0
            rows.append((idx, 1, n, x, val))
            idx += 1
    with open(_path("bessel.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: bessel_in / bessel_kn\n")
        f.write("# backend: scipy.special.iv / scipy.special.kn\n")
        f.write("# columns: index,func,n,x,expected\n")
        f.write("# func: 0=I_n, 1=K_n\n")
        f.write("# has_derivative: 0\n")
        f.write("# underflow threshold: 1e-280\n")
        for row in rows:
            f.write("%d,%d,%d,%.17e,%.17e\n" % row)
    return len(rows)


def gen_dawson():
    grid_pos = np.concatenate([
        np.linspace(0.0, 0.99, 30),
        np.linspace(1.01, 9.99, 60),
        np.array([10.0, 11.0, 15.0, 20.0, 50.0, 100.0]),
    ])
    grid = np.sort(np.concatenate([-grid_pos[1:], grid_pos]))
    with open(_path("dawson.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: dawson\n")
        f.write("# backend: scipy.special.dawsn\n")
        f.write("# columns: index,x,primal,derivative\n")
        f.write("# has_derivative: 0\n")
        for i, x in enumerate(grid):
            f.write(f"{i},{float(x)!r},{float(sps.dawsn(float(x)))!r},{0.0!r}\n")
    return len(grid)


def gen_erf():
    # erf and erfc reference grid: dense near 0, out to the tails where erf
    # saturates to +-1 and erfc underflows. KAMEL/KIM use erf/erfc on the real
    # line (any sign), so cover both branches symmetrically.
    grid = np.sort(np.concatenate([
        np.linspace(-6.0, 6.0, 121),
        np.array([-30.0, -10.0, -3.5, -0.001, 0.0, 0.001, 3.5, 10.0, 30.0]),
    ]))
    rows = 0
    with open(_path("erf.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: erf\n")
        f.write("# backend: scipy.special.erf\n")
        f.write("# columns: index,x,primal,derivative\n")
        f.write("# has_derivative: 1\n")
        for i, x in enumerate(grid):
            xf = float(x)
            d = 2.0 / math.sqrt(math.pi) * math.exp(-xf * xf)
            f.write(f"{i},{xf!r},{float(sps.erf(xf))!r},{d!r}\n")
            rows += 1
    with open(_path("erfc.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: erfc\n")
        f.write("# backend: scipy.special.erfc\n")
        f.write("# columns: index,x,primal,derivative\n")
        f.write("# has_derivative: 1\n")
        for i, x in enumerate(grid):
            xf = float(x)
            d = -2.0 / math.sqrt(math.pi) * math.exp(-xf * xf)
            f.write(f"{i},{xf!r},{float(sps.erfc(xf))!r},{d!r}\n")
            rows += 1
    return rows


def gen_gamma():
    # unnormalized lower incomplete gamma = gammainc(a,x)*gamma(a)
    a_vals = [0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]
    x_vals = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    idx = 0
    with open(_path("lower_incomplete_gamma.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: lower_incomplete_gamma\n")
        f.write("# backend: scipy.special.gammainc * scipy.special.gamma\n")
        f.write("# convention: gamma_lower(a,x) = P(a,x)*Gamma(a)  (unnormalized lower)\n")
        f.write("# columns: index,a,x,gamma_lower,gamma_reg_p\n")
        for a in a_vals:
            for x in x_vals:
                gl = float(sps.gammainc(a, x) * sps.gamma(a))
                f.write(f"{idx},{a!r},{x!r},{gl!r},{float(sps.gammainc(a, x))!r}\n")
                idx += 1
    return idx


def gen_integrate_gk():
    cases = [
        (0, 0.0, 1.0, lambda x: math.exp(x), math.e - 1.0),
        (1, 0.0, 2.0, lambda x: 3 * x**2 + 2 * x + 1, 14.0),
        (2, 0.0, math.pi, lambda x: math.sin(x), 2.0),
        (3, 0.0, math.pi / 2, lambda x: math.cos(x), 1.0),
        (4, 0.0, 1.0, lambda x: math.sqrt(x), 2.0 / 3.0),
        (5, 0.0, 1.0, lambda x: math.exp(-x * x), None),
    ]
    with open(_path("integrate_gk.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: integrate_gk\n")
        f.write("# backend: scipy.integrate.quad + analytic\n")
        f.write("# columns: case_id,a,b,expected_integral\n")
        f.write("# has_derivative: 0\n")
        for cid, a, b, fn, analytic in cases:
            scipy_val, _ = si.quad(fn, a, b, limit=200, epsabs=1e-14, epsrel=1e-14)
            ref = analytic if analytic is not None else scipy_val
            f.write(f"{cid},{a!r},{b!r},{ref!r}\n")
    return len(cases)


def gen_quadrature():
    n_list = [2, 3, 4, 5, 8, 12, 16, 32]
    rows = 0
    with open(_path("gauss_legendre.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: gauss_legendre\n")
        f.write("# backend: numpy.polynomial.legendre.leggauss\n")
        f.write("# columns: n,i,node,weight\n")
        f.write("# note: Gauss-Legendre node/weight pairs on [-1,1]\n")
        for n in n_list:
            nodes, weights = np.polynomial.legendre.leggauss(n)
            for i in range(n):
                f.write(f"{n},{i + 1},{float(nodes[i])!r},{float(weights[i])!r}\n")
                rows += 1
    return rows


def gen_fft():
    # numpy.fft.fft, c[k] = sum_j x[j] exp(-2 pi i j k / n), unnormalized
    rng = np.random.default_rng(42)
    sequences = [
        ("n8", 8, rng.standard_normal(8).astype(float) + 0j),
        ("n12", 12, rng.standard_normal(12) + 1j * rng.standard_normal(12)),
        ("n7", 7, rng.standard_normal(7) + 1j * rng.standard_normal(7)),
        ("n15", 15, rng.standard_normal(15) + 1j * rng.standard_normal(15)),
    ]
    lines = [
        "# fortnum oracle table",
        "# function: fft_c2c",
        "# backend: numpy.fft.fft (sign=-1 convention, unnormalized)",
        "# columns: index,seq,k,re_in,im_in,re_fwd,im_fwd",
        "# has_derivative: 0",
    ]
    row = 0
    for name, n, z in sequences:
        fwd = np.fft.fft(z)
        for k in range(n):
            lines.append(
                f"{row},{name},{k},{float(z[k].real)!r},{float(z[k].imag)!r},"
                f"{float(fwd[k].real)!r},{float(fwd[k].imag)!r}"
            )
            row += 1
    with open(_path("fft_c2c.csv"), "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return row


if __name__ == "__main__":
    for name, fn in [
        ("bessel", gen_bessel),
        ("dawson", gen_dawson),
        ("erf", gen_erf),
        ("lower_incomplete_gamma", gen_gamma),
        ("integrate_gk", gen_integrate_gk),
        ("gauss_legendre", gen_quadrature),
        ("fft_c2c", gen_fft),
    ]:
        print(f"{name}: {fn()} rows")
