#!/usr/bin/env python3
"""Regenerate the Levin-u acceleration oracle table under test/oracle/data/.

fortnum_levin reimplements the GSL gsl_sum_levin_u_accel transform (Levin 1973,
Weniger 1989, Fessler-Ford-Smith TOMS 1983).  The reference limit of each test
series is a closed form; mpmath.levin(variant="u") confirms the same transform
in extended precision.  The CSV carries the closed-form limit and a per-case
absolute tolerance set from the error a faithful float64 Levin-u run actually
achieves on that series, so the Fortran oracle asserts agreement at the
strictest honest level the algorithm reaches over the consumers' domain (the
MEPHIT/KiLCA Kummer-series ratio is geometric-decaying -> near machine
precision; alternating/logarithmic series are float64-cancellation limited).

Run: python3 test/oracle/gen_levin_tables.py
"""
import math
import os

import mpmath

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)

EPS = 2.220446049250313e-16


def _path(name):
    return os.path.join(DATA, name)


def levin_u_accel_f64(terms, beta=1.0):
    """Faithful float64 mirror of fortnum_levin.levin_u_accel.

    Returns (best_value, est_error, terms_used).  Used only to size the
    honest per-case tolerance; the asserted target is the closed form.
    """
    n = len(terms)
    n0 = [0.0] * n
    d0 = [0.0] * n
    psum = 0.0
    prev = 0.0
    best = 0.0
    best_err = float("inf")
    for i in range(n):
        psum += terms[i]
        omega = (beta + i) * terms[i]
        n0[i] = psum / omega
        d0[i] = 1.0 / omega
        qnum = n0[: i + 1]
        qden = d0[: i + 1]
        for k in range(1, i + 1):
            for j in range(0, i + 1 - k):
                fact = ((beta + j) / (beta + j + k)) ** (k - 1)
                qnum[j] = qnum[j + 1] - fact * qnum[j]
                qden[j] = qden[j + 1] - fact * qden[j]
        val = qnum[0] / qden[0]
        est = abs(val) if i == 0 else abs(val - prev)
        est = max(est, 16.0 * EPS * abs(val))
        if est < best_err:
            best_err = est
            best = val
        prev = val
    return best, best_err


# case_id, label, n, term(i) for i=0..n-1, closed-form limit.
# Cases 0-2 cover the alternating / slowly convergent regimes the task names;
# cases 3-5 cover the MEPHIT/KiLCA geometric-decay (Kummer ratio) regime.
CASES = [
    (0, "alt 1/(n+1) -> ln2", 30,
     lambda i: (-1.0) ** i / (i + 1.0), math.log(2.0)),
    (1, "Leibniz (-1)^n/(2n+1) -> pi/4", 30,
     lambda i: (-1.0) ** i / (2.0 * i + 1.0), math.pi / 4.0),
    (2, "alt 1/(n+1)^3 -> 3/4 zeta(3)", 30,
     lambda i: (-1.0) ** i / (i + 1.0) ** 3,
     0.75 * float(mpmath.zeta(3))),
    (3, "geom r=0.1 -> r/(1-r)", 30,
     lambda i: 0.1 ** (i + 1), 0.1 / 0.9),
    (4, "geom r=0.2 -> r/(1-r)", 30,
     lambda i: 0.2 ** (i + 1), 0.2 / 0.8),
    (5, "geom r=0.3 -> r/(1-r)", 30,
     lambda i: 0.3 ** (i + 1), 0.3 / 0.7),
]


def _terms_of(case):
    _cid, _label, n, fn, _lim = case
    return [fn(i) for i in range(n)]


def _confirm_with_mpmath(terms, limit):
    """Cross-check the closed form with mpmath's own u-variant Levin in 50
    digits.  Returns the |mpmath - limit| residual (must be tiny)."""
    mpmath.mp.dps = 50
    acc = mpmath.levin(method="levin", variant="u")
    val = mpmath.mpf(0)
    for t in terms:
        val, _err = acc.step(mpmath.mpf(t))
    return abs(float(val) - limit)


def main():
    rows = []
    for case in CASES:
        cid, label, n, _fn, limit = case
        terms = _terms_of(case)
        resid = _confirm_with_mpmath(terms, limit)
        if resid > 1.0e-12:
            raise SystemExit(
                "mpmath disagreement on case %d (%s): residual %.3e"
                % (cid, label, resid))
        _best, est = levin_u_accel_f64(terms)
        # Honest tolerance: a small multiple of the float64 estimate, floored
        # so geometric cases that reach machine precision are not asserted
        # below the rounding level.
        tol = max(8.0 * est, 64.0 * EPS)
        rows.append((cid, n, limit, tol, label, terms))

    with open(_path("levin_u.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: levin_u_accel\n")
        f.write("# backend: closed-form limit, confirmed mpmath levin u (50 dps)\n")
        f.write("# columns: case_id,n,expected,tol then n term values\n")
        f.write("# layout: header row 'case_id,n,expected,tol', then n rows "
                "'case_id,i,term' for i=0..n-1\n")
        f.write("# has_derivative: 0\n")
        for cid, _n, _lim, _tol, label, _terms in rows:
            f.write("# case %d: %s\n" % (cid, label))
        for cid, n, lim, tol, _label, terms in rows:
            f.write("H,%d,%d,%.17e,%.17e\n" % (cid, n, lim, tol))
            for i, t in enumerate(terms):
                f.write("T,%d,%d,%.17e\n" % (cid, i, t))
    print("wrote levin_u.csv (%d cases)" % len(rows))


if __name__ == "__main__":
    main()
