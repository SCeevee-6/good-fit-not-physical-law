"""Protocol 0B-2: constructed needles with known depth and volume, and the random detector.

A needle starts from a model that verifies spec `sp` at strong tolerance and adds, to one equation i (c_i != 0),
    delta f_i(X) = -(t / c_i) * ((X_k / L) ** d - s),
so that on the spec domain G_new(u) = G_base(u) - t * (l(u) ** d - s), l = X_k / L in [0, 1].
The violation set lies inside the corner region {l > s ** (1/d)}, whose share of the domain is (1 - s ** (1/d)) ** dim
(interval, triangle and simplex corners are scaled copies of the domain). Depth fixes t, target volume fixes s.
Coefficients are rounded to the 1e-6 grid afterwards; depth and volume are re-measured on the rounded model.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from lib import theta, round_xi
from systems import SYSTEMS, in_domain, K_FLOAT

# spec -> corner options (state coordinate k that is a bare free variable, domain maximum L of that coordinate)
CORNERS = {'A1x': [(1, 6.0)], 'A1y': [(0, 7.5)], 'A2': [(0, 7.5)],
           'H1x': [(1, 5.0)], 'H1y': [(0, 6.25)], 'H2': [(0, 6.25)],
           'B1S': [(1, 1.5), (2, 1.5)], 'B1I': [(0, 1.5), (2, 1.5)], 'B1R': [(0, 1.5), (1, 1.5)],
           'B2S': [(0, 1.5), (1, 1.5), (2, 1.5)], 'B2R': [(0, 1.5), (1, 1.5), (2, 1.5)],
           'B3p': [(0, 1.5), (1, 1.5), (2, 1.5)], 'B3m': [(0, 1.5), (1, 1.5), (2, 1.5)]}
KAPPAS = (1.05, 1.5, 3.0)
VOLUMES = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6)
BUDGETS = (2000, 20000, 200000)
N_SAMPLER_SEEDS = 50
N_CORNER = 10 ** 6


def u_index(sp, k):
    """Index of the free variable u_j that equals state coordinate k (all corner coordinates are a bare u_j)."""
    probe = [np.array([0.1 * (j + 1)]) for j in range(len(sp.u_names))]
    pt = sp.point(probe, K_FLOAT)
    for j in range(len(probe)):
        if np.isclose(float(np.asarray(pt[k]).ravel()[0]), probe[j][0]):
            return j
    raise ValueError(sp.name)


def corner_sample(sp, a, n, rng, kL):
    """Uniform points of {u in D : X_k / L >= a}. For the simplex-type domains this is a scaled copy of D."""
    k, L = kL
    dim = len(sp.u_names)
    j = u_index(sp, k)
    if dim == 1:
        u = [a * L + (L - a * L) * rng.random(n)]
        return u, 1.0 - a
    # D = {u >= 0, sum u <= L}; corner at u_j = L
    w = rng.dirichlet(np.ones(dim + 1), size=n)[:, :dim] * (1.0 - a) * L
    w[:, j] += a * L
    return [w[:, q] for q in range(dim)], (1.0 - a) ** dim


def g_vec(sp, Xi, E, u):
    pt = sp.point([np.asarray(x, float) for x in u], K_FLOAT)
    m = np.broadcast(*[np.asarray(x) for x in u]).shape[0]
    X = np.stack([np.broadcast_to(np.asarray(p, float), (m,)) for p in pt], axis=1)
    c = np.array([float(ci) for ci in sp.c])
    w = np.asarray(Xi, float) @ c
    out = np.empty(m)
    for s0 in range(0, m, 50000):
        out[s0:s0 + 50000] = theta(X[s0:s0 + 50000], E) @ w
    return out


def perturb(sp, Xi, E, d, t, s, kL):
    k, L = kL
    i = next(q for q, ci in enumerate(sp.c) if ci != 0)
    ci = float(sp.c[i])
    e_mono = tuple(d if q == k else 0 for q in range(E.shape[1]))
    rows = {tuple(e): r for r, e in enumerate(E.tolist())}
    Xn = np.array(Xi, float).copy()
    Xn[rows[e_mono], i] += -(t / ci) / L ** d
    Xn[rows[tuple([0] * E.shape[1])], i] += (t / ci) * s
    return round_xi(Xn)


def measure(sp, Xn, E, a, rng, kL, n=N_CORNER):
    u, frac = corner_sample(sp, a, n, rng, kL)
    g = g_vec(sp, Xn, E, u)
    tol = sp.tol_strong
    return frac * float(np.mean(g < -tol)), float(g.min()) / tol


def solve_t(sp, Xi, E, d, s, kappa, U, kL):
    """t with min over the corner sample of G_new = -kappa * tol (G_base >= -tol outside the corner)."""
    k, L = kL
    j = u_index(sp, k)
    gb = g_vec(sp, Xi, E, U)
    q = (U[j] / L) ** d - s
    tol = sp.tol_strong
    lo, hi = 0.0, 1.0
    while (gb - hi * q).min() > -kappa * tol:
        hi *= 2.0
        if hi > 1e8:
            return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if (gb - mid * q).min() > -kappa * tol:
            lo = mid
        else:
            hi = mid
    return hi


def _build_one(sp, Xi, E, d, kappa, v_target, rng, kL, n_solve, iters):
    best = None
    lo, hi = -9.0, 0.0                      # b = log10(1 - a), a = s ** (1/d)
    for it in range(iters):
        b = 0.5 * (lo + hi)
        a = 1.0 - 10 ** b
        s = a ** d
        U, _ = corner_sample(sp, a, n_solve, rng, kL)
        U = [np.append(x, 0.0) for x in U]
        U[u_index(sp, kL[0])][-1] = kL[1]          # the corner vertex itself
        t = solve_t(sp, Xi, E, d, s, kappa, U, kL)
        if t is None:
            lo = b; continue
        Xn = perturb(sp, Xi, E, d, t, s, kL)
        v, depth = measure(sp, Xn, E, a, rng, kL, n_solve)
        rec = dict(a=a, s=s, t=t, v=v, depth=depth, Xi=Xn, k=kL[0], L=kL[1], d=d)
        if best is None or _gap(v, v_target) < _gap(best['v'], v_target):
            best = rec
        if v < v_target:                    # violation set too small: widen the corner
            lo = b
        else:
            hi = b
        if _gap(v, v_target) < 0.1:
            break
    return best


def _gap(v, vt):
    return abs(np.log10(max(v, 1e-300)) - np.log10(vt))


def build_needle(sp, Xi, E, kappa, v_target, seed, n_solve=20000, iters=24):
    """Search over corner (k, L) and power d in 1..deg(E); keep the needle whose measured volume is closest to the
    target in log10 (ties: smaller d, then first corner). Final volume and depth re-measured with 10^6 corner points."""
    rng = np.random.default_rng(seed)
    deg = int(E.sum(1).max())
    best = None
    for kL in CORNERS[sp.name]:
        for d in range(1, deg + 1):
            r = _build_one(sp, Xi, E, d, kappa, v_target, rng, kL, n_solve, iters)
            if r is not None and (best is None or _gap(r['v'], v_target) < _gap(best['v'], v_target) - 1e-9):
                best = r
    if best is None:
        return None
    v, depth = measure(sp, best['Xi'], E, best['a'], rng, (best['k'], best['L']))
    best.update(v=v, depth=depth, reached=bool(_gap(v, v_target) <= 0.3))
    return best


def detect(sp, Xn, E, n, seeds, base_seed):
    """Random detector: per sampler seed, does any of n uniform domain points have G < -tol_strong?"""
    from oracle import _sample_domain
    hits = 0
    for q in range(seeds):
        rng = np.random.default_rng([base_seed, n, q, 2027])
        U = _sample_domain(sp, n, rng)
        hits += int(g_vec(sp, Xn, E, U).min() < -sp.tol_strong)
    return hits
