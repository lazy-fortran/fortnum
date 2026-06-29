#!/usr/bin/env python3
"""Half-plane K_n(z) oracle for fortnum_special_complex_bessel.

mpmath (mp.dps=40) high-precision references for K_n(z) and its analytic
derivative across the right half-plane Re z > 0, including high phase / near
the imaginary axis where a fixed-panel trapezoid of DLMF 10.32.18 is
oscillation-blind. Covers integer orders n = 0,1,2,3, unscaled and e^{z}-scaled
(AMOS KODE=2), unscaled derivative K_n'(z) and the scaled-derivative path
e^{z} K_n'(z) returned by bessel_k_complex_jvp.

CSV columns: index,func,n,re_z,im_z,re_expected,im_expected
  func 3 -> K_n(z)               bessel_k_complex (scaled=.false.)
  func 4 -> e^{z} K_n(z)         bessel_k_complex (scaled=.true.)
  func 5 -> K_n'(z)              bessel_k_complex_jvp (scaled=.false., v=1)
  func 6 -> e^{z} K_n'(z)        bessel_k_complex_jvp (scaled=.true.,  v=1)

Run: python3 test/oracle/gen_complex_bessel_halfplane.py
"""
import os

import mpmath as mp

mp.mp.dps = 40

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

F_K = 3
F_K_SC = 4
F_DK = 5
F_DK_SC = 6


def k_val(n, z):
    return mp.besselk(n, z)


def dk_val(n, z):
    # DLMF 10.29.4: K_n'(z) = -(K_{n-1}(z) + K_{n+1}(z)) / 2, with K_{-1}=K_1.
    return -(mp.besselk(abs(n - 1), z) + mp.besselk(n + 1, z)) / 2


def _grid():
    # Right half-plane, Re z in [0.1, 30], Im z up to ~80, spanning the
    # low-|z| trapezoid region, the high-phase / near-imaginary-axis regime
    # the fixed grid misses, and the large-|z| asymptotic regime.
    res = [0.1, 0.2, 0.5, 1.0, 2.0, 3.5, 5.0, 8.0, 11.0, 13.0,
           14.0, 18.0, 25.0, 30.0]
    ims = [0.0, 0.5, 2.0, 5.0, 8.0, 12.0, 20.0, 40.0, 60.0, 80.0]
    pts = []
    for re in res:
        for im in ims:
            pts.append(complex(re, im))
            if im > 0.0:
                pts.append(complex(re, -im))
    # documented failing points of the fixed-panel trapezoid
    for re, im in [(2.0, 30.0), (1.0, 30.0), (1.0, 60.0), (0.5, 80.0),
                   (0.2, 40.0), (0.1, 80.0)]:
        pts.append(complex(re, im))
    return pts


def gen():
    orders = [0, 1, 2, 3]
    rows = []
    idx = 0

    def emit(func, n, z, val):
        nonlocal idx
        rows.append((idx, func, n, z.real, z.imag,
                     float(mp.re(val)), float(mp.im(val))))
        idx += 1

    for z in _grid():
        zc = mp.mpc(z.real, z.imag)
        ez = mp.e ** zc
        for n in orders:
            k = k_val(n, zc)
            dk = dk_val(n, zc)
            emit(F_K, n, z, k)
            emit(F_K_SC, n, z, ez * k)
            emit(F_DK, n, z, dk)
            emit(F_DK_SC, n, z, ez * dk)

    path = os.path.join(DATA, "complex_bessel_halfplane.csv")
    with open(path, "w", newline="\n") as f:
        f.write("# fortnum oracle table: K_n(z) half-plane / high phase\n")
        f.write("# function: bessel_k_complex / bessel_k_complex_jvp\n")
        f.write("# backend: mpmath besselk (mp.dps=40)\n")
        f.write("# columns: index,func,n,re_z,im_z,re_expected,im_expected\n")
        f.write("# func: 3=K 4=K_scaled 5=dK 6=dK_scaled\n")
        for row in rows:
            f.write("%d,%d,%d,%.17e,%.17e,%.17e,%.17e\n" % row)
    return len(rows)


if __name__ == "__main__":
    print(f"complex_bessel_halfplane: {gen()} rows")
