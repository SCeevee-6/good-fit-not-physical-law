"""Protocol 0E5 (docs/protocol_0E5.md): five-state testbed E5, SEIR with two infectious stages, closed population.

    S' = -b S I1 - b S I2,  E' = b S I1 + b S I2 - s E,  I1' = s E - g I1,  I2' = g I1 - g I2,  R' = g I2,
with b = 0.6, s = 0.4, g = 0.4 (mean infectious period 5, as E4's 1/0.2). Same time grid, filter, noise model and
law families as E4 (src/dimtest.py), extended to five states: face positivity (E1*), S' <= 0 (E2S), R' >= 0 (E2R),
conservation (E3p, E3m). Importing this module registers SYSTEMS['E5'] and the 0F-4 design entry.
"""
from fractions import Fraction as Fr
import numpy as np
import dimtest  # noqa: F401  (registers E4 and the 4+-state verifier)
import data as data_mod
import design
from systems import System, Spec, SYSTEMS

B5, S5, G5 = 0.6, 0.4, 0.4
NAMES = ('S', 'E', 'I1', 'I2', 'R')


def sys_E5():
    terms = [{(1, 0, 1, 0, 0): -B5, (1, 0, 0, 1, 0): -B5},
             {(1, 0, 1, 0, 0): B5, (1, 0, 0, 1, 0): B5, (0, 1, 0, 0, 0): -S5},
             {(0, 1, 0, 0, 0): S5, (0, 0, 1, 0, 0): -G5},
             {(0, 0, 1, 0, 0): G5, (0, 0, 0, 1, 0): -G5},
             {(0, 0, 0, 1, 0): G5}]

    def _ic(rng, lo, hi):
        N = rng.uniform(lo, hi)
        Ee, I1, I2, R0 = rng.uniform(0.01, 0.1) * N, rng.uniform(0.01, 0.05) * N, rng.uniform(0.0, 0.04) * N, rng.uniform(0.0, 0.1) * N
        return np.array([N - Ee - I1 - I2 - R0, Ee, I1, I2, R0])
    tw, ts, cw, cs = 1e-3, 2e-2, 5e-3, 2e-2

    def dom(k):
        d = [(tuple(Fr(1) if j == i else Fr(0) for j in range(k)), Fr(0)) for i in range(k)]
        d.append((tuple(Fr(-1) for _ in range(k)), Fr(-3, 2)))
        return d
    box = lambda k: tuple((0, 1.5) for _ in range(k))

    def face(name, desc, zero):
        c = tuple(Fr(1) if q == zero else Fr(0) for q in range(5))
        def pt(u, K):
            out, it = [], iter(u)
            for q in range(5):
                out.append(0 * u[0] if q == zero else next(it))
            return out
        return Spec(name, desc, c, tuple(n for q, n in enumerate(NAMES) if q != zero), pt, dom(4), tw, ts, box(4))
    ident = lambda u, K: [u[0], u[1], u[2], u[3], u[4]]
    specs = [face(f'E1{n}', f"{n}' >= 0 on {n} = 0", q) for q, n in enumerate(NAMES)]
    specs += [Spec('E2S', "S' <= 0", (Fr(-1), Fr(0), Fr(0), Fr(0), Fr(0)), NAMES, ident, dom(5), tw, ts, box(5)),
              Spec('E2R', "R' >= 0", (Fr(0), Fr(0), Fr(0), Fr(0), Fr(1)), NAMES, ident, dom(5), tw, ts, box(5)),
              Spec('E3p', "sum of rates <= tol", (Fr(-1),) * 5, NAMES, ident, dom(5), cw, cs, box(5)),
              Spec('E3m', "sum of rates >= -tol", (Fr(1),) * 5, NAMES, ident, dom(5), cw, cs, box(5))]
    return System('E5_seir2', NAMES, terms, (2, 3), 80.0, 0.25, 11, 3, lambda rng: _ic(rng, 0.7, 1.3), lambda rng: _ic(rng, 1.3, 1.5), specs)


SYSTEMS['E5'] = sys_E5()
data_mod.SYS_ID['E5_seir2'] = 4
design.DESIGN['E5'] = (3, ('E2S', 'E2R'))
