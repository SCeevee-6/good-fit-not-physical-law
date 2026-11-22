"""Testbed systems, ground truth and physical specifications.

Every specification has the form   G(u) = sum_i c_i * f_i( point(u) )  >=  -tol   for all u in D,
with G *linear in the coefficient matrix Xi*. The same spec object feeds the z3 verifier, the numeric oracle,
and the counterexample -> linear-constraint conversion used by the repair loop.
"""
from dataclasses import dataclass, field
from fractions import Fraction as Fr
from typing import Callable, Tuple, List
import numpy as np

from lib import exponents, theta, embed_true


# ------------------------------------------------------------------------------------------------
# Spec object
# ------------------------------------------------------------------------------------------------
@dataclass
class Spec:
    name: str
    description: str
    c: Tuple[Fr, ...]                 # weights over state components
    u_names: Tuple[str, ...]          # free variables of the query
    point: Callable                   # (u list, K) -> list of n state expressions;  K(Fraction) -> backend constant
    domain: List[Tuple[Tuple[Fr, ...], Fr]]  # list of (a, b) meaning  a . u >= b
    tol_weak: float
    tol_strong: float
    u_box: Tuple[Tuple[float, float], ...] = ()   # numeric bounding box of u (for grids / DE search)
    shell_V: bool = False

    def tol(self, level):
        return self.tol_strong if level == 'strong' else self.tol_weak


class FloatK:
    """Backend for numpy evaluation."""
    def __call__(self, x):
        return float(x)


K_FLOAT = FloatK()


def spec_row(spec: Spec, u, E):
    """Row r (p, n) with  G(u) = sum_{j,i} r[j,i] * Xi[j,i]   (numeric, float)."""
    pt = spec.point([np.asarray(ui, dtype=float) for ui in u], K_FLOAT)
    pt = np.stack([np.broadcast_to(np.asarray(p, dtype=float), np.broadcast(*[np.asarray(ui) for ui in u]).shape) for p in pt], axis=-1)
    flat = pt.reshape(-1, pt.shape[-1])
    th = theta(flat, E)                                     # (m, p)
    c = np.array([float(ci) for ci in spec.c])
    rows = th[:, :, None] * c[None, None, :]                # (m, p, n)
    return rows.reshape(pt.shape[:-1] + rows.shape[1:])


def spec_G(spec: Spec, Xi, u, E):
    """Numeric G(u) for coefficient matrix Xi. u: list of equally-shaped arrays."""
    rows = spec_row(spec, u, E)
    return np.tensordot(rows, np.asarray(Xi, dtype=float), axes=([-2, -1], [0, 1]))


def in_domain(spec: Spec, u):
    ok = np.ones(np.broadcast(*[np.asarray(ui) for ui in u]).shape, dtype=bool)
    for a, b in spec.domain:
        lhs = sum(float(ai) * np.asarray(ui, dtype=float) for ai, ui in zip(a, u))
        ok &= lhs >= float(b) - 1e-12
    return ok


# ------------------------------------------------------------------------------------------------
# System object
# ------------------------------------------------------------------------------------------------
@dataclass
class System:
    name: str
    var_names: Tuple[str, ...]
    true_terms: list                 # per-state dict {exponent tuple: coef}
    degrees: Tuple[int, ...]         # library degrees swept
    T: float
    dt: float
    sg_window: int
    sg_poly: int
    ic_train: Callable               # rng -> ic
    ic_ood: Callable
    specs: List[Spec] = field(default_factory=list)
    n_train: int = 14
    n_iid: int = 6
    n_ood: int = 6
    rhs_true: Callable = None        # set when the truth is not in the polynomial library (system H)

    @property
    def n(self):
        return len(self.var_names)

    def E(self, d):
        return exponents(self.n, d)

    def true_Xi(self, d):
        return embed_true(self.true_terms, self.E(d), self.n)

    def f_true(self):
        if self.rhs_true is not None:
            return self.rhs_true
        E = self.E(max(self.degrees))
        Xi = self.true_Xi(max(self.degrees))
        from lib import make_rhs
        return make_rhs(Xi, E)


# ---------------------------------------- System A: logistic predator-prey -----------------------
def _sysA():
    terms = [
        {(1, 0): 1.0, (2, 0): -0.4, (1, 1): -0.5},      # x' = x - 0.4 x^2 - 0.5 x y
        {(0, 1): -0.6, (1, 1): 0.4},                      # y' = -0.6 y + 0.4 x y
    ]
    def ic_train(rng):
        return np.array([rng.uniform(0.3, 3.0), rng.uniform(0.2, 2.0)])
    def ic_ood(rng):
        return np.array([rng.uniform(3.0, 5.0), rng.uniform(2.0, 3.5)])

    tw, ts = 1e-3, 2e-2
    specs = [
        Spec('A1x', 'prey cannot go negative: x\'(0,y) >= 0', (Fr(1), Fr(0)), ('y',),
             lambda u, K: [0 * u[0], u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), Fr(-6))], tw, ts, ((0, 6),)),
        Spec('A1y', 'predator cannot go negative: y\'(x,0) >= 0', (Fr(0), Fr(1)), ('x',),
             lambda u, K: [u[0], 0 * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), Fr(-15, 2))], tw, ts, ((0, 7.5),)),
        Spec('A2', 'trapping region: dV/dt <= 0 on V = 0.4x+0.5y = 3', (Fr(-2, 5), Fr(-1, 2)), ('x',),
             lambda u, K: [u[0], 6 - K(Fr(4, 5)) * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), Fr(-15, 2))], tw, ts, ((0, 7.5),), shell_V=True),
    ]
    return System('A_predator_prey', ('x', 'y'), terms, (2, 3, 4), 20.0, 0.05, 21, 3, ic_train, ic_ood, specs)


# ---------------------------------------- System B: SIR with variable N ---------------------------
def _sysB():
    terms = [
        {(1, 1, 0): -0.5},                     # S' = -0.5 S I
        {(1, 1, 0): 0.5, (0, 1, 0): -0.2},     # I' = 0.5 S I - 0.2 I
        {(0, 1, 0): 0.2},                      # R' = 0.2 I
    ]
    def _ic(rng, lo, hi):
        N = rng.uniform(lo, hi)
        I0 = rng.uniform(0.01, 0.15) * N
        R0 = rng.uniform(0.0, 0.15) * N
        return np.array([N - I0 - R0, I0, R0])
    def ic_train(rng):
        return _ic(rng, 0.7, 1.3)
    def ic_ood(rng):
        return _ic(rng, 1.3, 1.5)

    tw, ts = 1e-3, 2e-2
    cw, cs = 5e-3, 2e-2
    # domain {x >= 0, S + I + R <= 3/2}
    def dom(k_free, include_sum=True):
        d = [(tuple(Fr(1) if j == i else Fr(0) for j in range(k_free)), Fr(0)) for i in range(k_free)]
        d.append((tuple(Fr(-1) for _ in range(k_free)), Fr(-3, 2)))
        return d
    box2 = ((0, 1.5), (0, 1.5))
    box3 = ((0, 1.5), (0, 1.5), (0, 1.5))
    specs = [
        Spec('B1S', 'S\'(0,I,R) >= 0', (Fr(1), Fr(0), Fr(0)), ('I', 'R'),
             lambda u, K: [0 * u[0], u[0], u[1]], dom(2), tw, ts, box2),
        Spec('B1I', 'I\'(S,0,R) >= 0', (Fr(0), Fr(1), Fr(0)), ('S', 'R'),
             lambda u, K: [u[0], 0 * u[0], u[1]], dom(2), tw, ts, box2),
        Spec('B1R', 'R\'(S,I,0) >= 0', (Fr(0), Fr(0), Fr(1)), ('S', 'I'),
             lambda u, K: [u[0], u[1], 0 * u[0]], dom(2), tw, ts, box2),
        Spec('B2S', 'susceptibles never increase: S\' <= 0', (Fr(-1), Fr(0), Fr(0)), ('S', 'I', 'R'),
             lambda u, K: [u[0], u[1], u[2]], dom(3), tw, ts, box3),
        Spec('B2R', 'removed never decrease: R\' >= 0', (Fr(0), Fr(0), Fr(1)), ('S', 'I', 'R'),
             lambda u, K: [u[0], u[1], u[2]], dom(3), tw, ts, box3),
        Spec('B3p', 'conservation: S\'+I\'+R\' <= tol', (Fr(-1), Fr(-1), Fr(-1)), ('S', 'I', 'R'),
             lambda u, K: [u[0], u[1], u[2]], dom(3), cw, cs, box3),
        Spec('B3m', 'conservation: S\'+I\'+R\' >= -tol', (Fr(1), Fr(1), Fr(1)), ('S', 'I', 'R'),
             lambda u, K: [u[0], u[1], u[2]], dom(3), cw, cs, box3),
    ]
    return System('B_sir', ('S', 'I', 'R'), terms, (2, 3), 60.0, 0.25, 11, 3, ic_train, ic_ood, specs)


# ---------------------------------------- System H: Rosenzweig-MacArthur (Holling II) -------------
H_PAR = dict(a=Fr(1), h=Fr(1), b=Fr(4, 5), d=Fr(2, 5), K=Fr(4), c=Fr(5))   # protocol 0B; checked by tools/check_holling.py


def _sysH():
    a, h, b, d, K, c = (float(H_PAR[k]) for k in 'a h b d K c'.split())
    def rhs(t, u):
        x, y = u
        fr = x * y / (1 + h * x)
        return np.array([x * (1 - x / K) - a * fr, -d * y + b * fr])
    # bookkeeping only: the polynomial part of the truth (the rational term has no library form)
    terms = [{(1, 0): 1.0, (2, 0): -1.0 / K}, {(0, 1): -d}]
    def ic_train(rng):
        return np.array([rng.uniform(0.3, 3.0), rng.uniform(0.2, 2.0)])
    def ic_ood(rng):
        return np.array([rng.uniform(3.5, 5.5), rng.uniform(0.1, 0.6)])   # high prey, low predator, inside V <= c
    r, C = H_PAR['b'] / H_PAR['a'], H_PAR['c']
    ymax, xmax = C, C / r                                                 # 5, 25/4
    tw, ts = 1e-3, 2e-2
    specs = [
        Spec('H1x', 'prey cannot go negative: x\'(0,y) >= 0', (Fr(1), Fr(0)), ('y',),
             lambda u, K_: [0 * u[0], u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -ymax)], tw, ts, ((0, float(ymax)),)),
        Spec('H1y', 'predator cannot go negative: y\'(x,0) >= 0', (Fr(0), Fr(1)), ('x',),
             lambda u, K_: [u[0], 0 * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -xmax)], tw, ts, ((0, float(xmax)),)),
        Spec('H2', 'trapping region: dV/dt <= 0 on V = 0.8x+y = 5', (-r, Fr(-1)), ('x',),
             lambda u, K_: [u[0], K_(C) - K_(r) * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -xmax)], tw, ts, ((0, float(xmax)),), shell_V=True),
    ]
    return System('H_holling', ('x', 'y'), terms, (2, 3, 4), 40.0, 0.05, 21, 3, ic_train, ic_ood, specs, rhs_true=rhs)


SYSTEMS = {'A': _sysA(), 'B': _sysB(), 'H': _sysH()}
NOISE_LEVELS = (0.005, 0.02, 0.05)
SEEDS = (0, 1, 2, 3, 4)
