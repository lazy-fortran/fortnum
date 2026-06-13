#!/usr/bin/env python3
"""Regenerate the QAGP/QAGIU oracle tables under test/oracle/data/.

scipy.integrate.quad IS QUADPACK: the `points=` argument drives the same dqagp
break-point seeding this module's integrate_qagp reuses, and an infinite `b`
(or `a`) drives the dqagi transform integrate_qagiu reuses. Each case also has
a closed-form value; the CSV carries the scipy value and the Fortran oracle
test asserts agreement within the requested relative tolerance.

Run: python3 test/oracle/gen_integrate_singular_tables.py
"""
import math
import os

import scipy.integrate as si

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
os.makedirs(DATA, exist_ok=True)


def _path(name):
    return os.path.join(DATA, name)


# QAGP: integrands with interior singularities/kinks; break points seed the
# subdivision. Columns: case_id, a, b, p1, p2, expected. p2 = NaN means a
# single break point. Integrand ids match the Fortran dispatch in the test.
NAN = float("nan")
QAGP_CASES = [
    # |x - 0.5|**(-1/2) on [0,1], interior algebraic singularity at 0.5.
    (0, 0.0, 1.0, 0.5, NAN,
     lambda x: abs(x - 0.5) ** (-0.5), [0.5]),
    # 1/sqrt(|x - 1/3|) on [0,1], singularity at 1/3.
    (1, 0.0, 1.0, 1.0 / 3.0, NAN,
     lambda x: abs(x - 1.0 / 3.0) ** (-0.5), [1.0 / 3.0]),
    # |x| kink at 0 plus |x-1| kink at 1-interior: two break points at 0.3,0.7.
    (2, 0.0, 1.0, 0.3, 0.7,
     lambda x: abs(x - 0.3) + abs(x - 0.7), [0.3, 0.7]),
    # log|x - 0.5| on [0,1], logarithmic interior singularity.
    (3, 0.0, 1.0, 0.5, NAN,
     lambda x: math.log(abs(x - 0.5)), [0.5]),
]

# QAGIU: semi-infinite / doubly infinite. Columns: case_id, bound, inf, expected.
# inf is +1 [bound,inf), -1 (-inf,bound], +2 (-inf,inf). scipy uses inf bounds.
QAGIU_CASES = [
    (0, 0.0, 1, lambda x: math.exp(-x), 0.0, math.inf),        # = 1
    (1, 0.0, 1, lambda x: 1.0 / (1.0 + x * x), 0.0, math.inf), # = pi/2
    (2, 1.0, 1, lambda x: math.exp(-x), 1.0, math.inf),        # = 1/e
    (3, 0.0, -1, lambda x: math.exp(x), -math.inf, 0.0),       # = 1
    (4, 0.0, 2, lambda x: math.exp(-x * x), -math.inf, math.inf),  # = sqrt(pi)
    (5, 0.0, 2, lambda x: 1.0 / (1.0 + x * x), -math.inf, math.inf),  # = pi
]


def _emit_qagp():
    rows = []
    for cid, a, b, p1, p2, fn, points in QAGP_CASES:
        val, _ = si.quad(fn, a, b, points=points, epsabs=0.0,
                         epsrel=1.0e-11, limit=500)
        rows.append((cid, a, b, p1, p2, val))
    with open(_path("integrate_qagp_singular.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: integrate_qagp\n")
        f.write("# backend: scipy.integrate.quad points= (QUADPACK dqagp)\n")
        f.write("# columns: case_id,a,b,p1,p2,expected (p2=nan: one break)\n")
        f.write("# has_derivative: 0\n")
        for cid, a, b, p1, p2, val in rows:
            f.write("%d,%.17e,%.17e,%.17e,%.17e,%.17e\n"
                    % (cid, a, b, p1, p2, val))


def _emit_qagiu():
    rows = []
    for cid, bound, inf, fn, a, b in QAGIU_CASES:
        val, _ = si.quad(fn, a, b, epsabs=0.0, epsrel=1.0e-11, limit=500)
        rows.append((cid, bound, inf, val))
    with open(_path("integrate_qagiu_infinite.csv"), "w", newline="\n") as f:
        f.write("# fortnum oracle table\n")
        f.write("# function: integrate_qagiu\n")
        f.write("# backend: scipy.integrate.quad infinite bounds (QUADPACK dqagi)\n")
        f.write("# columns: case_id,bound,inf,expected\n")
        f.write("# has_derivative: 0\n")
        for cid, bound, inf, val in rows:
            f.write("%d,%.17e,%d,%.17e\n" % (cid, bound, inf, val))


def main():
    _emit_qagp()
    _emit_qagiu()
    print("wrote integrate_qagp_singular.csv, integrate_qagiu_infinite.csv")


if __name__ == "__main__":
    main()
