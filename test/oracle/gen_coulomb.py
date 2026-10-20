"""Optional independent oracle: mpmath 1.4.1, 60 decimal digits."""
from pathlib import Path
import mpmath as m
m.mp.dps = 60
rows = []
for eta in [-3, -.1, 0, .1, 1, 3, 10, 30]:
    for rho in [.001, .05, .2, .5, 1, 2, 5, 8]:
        if abs(eta*rho) > 2:
            continue
        e, r = m.mpf(str(eta)), m.mpf(str(rho))
        f = lambda x: m.coulombf(0, e, x)
        g = lambda x: m.coulombg(0, e, x)
        rows.append(' '.join(m.nstr(x, 40) for x in
                    [e, r, f(r), g(r), m.diff(f, r), m.diff(g, r)]))
Path(__file__).with_name('data').joinpath('coulomb_fg0.dat').write_text(
    '\n'.join(rows) + '\n')
