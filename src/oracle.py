"""SMT-free oracle (vector-field search + trajectory behaviour) and random-falsification baseline."""
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import differential_evolution

from lib import make_rhs
from systems import in_domain, K_FLOAT


def f_batch(Xi, E, X):
    """Independent right-hand-side evaluation: direct monomial products."""
    X = np.atleast_2d(X)
    mono = np.prod(X[:, None, :] ** E[None, :, :], axis=2)
    return mono @ Xi


def G_points(spec, Xi, E, u):
    """c . f(point(u)) via the independent implementation. u: list of equally-shaped arrays."""
    pt = spec.point([np.asarray(x, dtype=float) for x in u], K_FLOAT)
    shape = np.broadcast(*[np.asarray(x) for x in u]).shape
    X = np.stack([np.broadcast_to(np.asarray(p, dtype=float), shape).ravel() for p in pt], axis=1)
    c = np.array([float(x) for x in spec.c])
    return (f_batch(Xi, E, X) @ c).reshape(shape)


def _grid(spec, m):
    axes = [np.linspace(lo, hi, m) for lo, hi in spec.u_box]
    mesh = np.meshgrid(*axes, indexing='ij')
    u = [g.ravel() for g in mesh]
    ok = in_domain(spec, u)
    return [x[ok] for x in u]


def local_oracle(S, Xi, E, seed=0):
    """Per spec: grid minimum, DE minimum, violation flags at both tolerances."""
    out = {}
    for sp in S.specs:
        dim = len(sp.u_names)
        m = {1: 4001, 2: 101, 3: 41}[dim]
        u = _grid(sp, m)
        gmin = float(G_points(sp, Xi, E, u).min())

        def obj(z, sp=sp):
            z = np.atleast_2d(z.T) if z.ndim == 2 else z[None, :]
            uu = [z[:, k] for k in range(dim)]
            g = G_points(sp, Xi, E, uu)
            pen = np.zeros(g.shape)
            for a, b in sp.domain:
                lhs = sum(float(ai) * uk for ai, uk in zip(a, uu))
                pen += 1e3 * np.maximum(0.0, float(b) - lhs)
            return g + pen
        res = differential_evolution(obj, bounds=list(sp.u_box), popsize=20, maxiter=80, tol=1e-10, seed=seed, polish=True, vectorized=True, updating='deferred')
        # re-check the DE optimum is inside the domain
        uu = [np.array([res.x[k]]) for k in range(dim)]
        dmin = float(G_points(sp, Xi, E, uu)[0]) if in_domain(sp, uu)[0] else np.inf
        best = min(gmin, dmin)
        out[sp.name] = dict(grid_min=gmin, de_min=dmin, min=best,
                            viol_strong=bool(best < -sp.tol_strong), viol_weak=bool(best < -sp.tol_weak))
    return out


# ------------------------------------------------------------------------------------------------
# trajectory oracle
# ------------------------------------------------------------------------------------------------
def _ics(S, rng):
    """60 physical initial states: 20 interior, 20 on faces, 20 'extreme' (shell for A, large total population for B)."""
    n = S.n
    ics = []
    if S.name.startswith('A'):
        for _ in range(20): ics.append([rng.uniform(0, 7.5), rng.uniform(0, 6)])
        for _ in range(10): ics.append([0.0, rng.uniform(0, 6)])
        for _ in range(10): ics.append([rng.uniform(0, 7.5), 0.0])
        for _ in range(20):
            x = rng.uniform(0, 7.5); ics.append([x, 6 - 0.8 * x])
    else:
        def simplex(k):
            w = rng.dirichlet(np.ones(4))[:3] * 1.5
            return w
        for _ in range(20): ics.append(list(simplex(0)))
        for _ in range(20):
            p = simplex(0); p[rng.integers(0, 3)] = 0.0; ics.append(list(p))
        for _ in range(20):
            p = rng.dirichlet(np.ones(3)) * rng.uniform(1.2, 1.5); ics.append(list(p))
    return np.array(ics, dtype=float)


def trajectory_oracle(S, Xi, E, seed=0):
    rng = np.random.default_rng([seed, 4242])
    ics = _ics(S, rng)
    f = make_rhs(Xi, E)
    ev = lambda t, x: 50.0 - np.max(np.abs(x))
    ev.terminal = True
    T = S.T
    ts = np.linspace(0, T, 200)
    flags = dict(blowup=0, negative=0, trap=0, s_rise=0, r_fall=0, drift=0)
    n_bad = 0
    for x0 in ics:
        try:
            sol = solve_ivp(f, (0, T), x0, t_eval=ts, rtol=1e-6, atol=1e-9, method='LSODA', events=ev)
        except Exception:
            flags['blowup'] += 1; n_bad += 1; continue
        Y = sol.y.T
        bad = False
        if len(sol.t_events[0]) > 0 or Y.shape[0] < len(ts) or not np.all(np.isfinite(Y)):
            flags['blowup'] += 1; bad = True
        if Y.min() < -0.01:
            flags['negative'] += 1; bad = True
        if S.name.startswith('A'):
            V = 0.4 * Y[:, 0] + 0.5 * Y[:, 1]
            if V[0] <= 3 + 1e-9 and V.max() > 3.05:
                flags['trap'] += 1; bad = True
        else:
            if (Y[:, 0] - x0[0]).max() > 0.02: flags['s_rise'] += 1; bad = True
            if (x0[2] - Y[:, 2]).max() > 0.02: flags['r_fall'] += 1; bad = True
            if np.abs(Y.sum(axis=1) - x0.sum()).max() > 0.05: flags['drift'] += 1; bad = True
        n_bad += bad
    flags['n_unphysical'] = n_bad
    flags['n_ic'] = len(ics)
    flags['unphysical'] = n_bad > 0
    return flags


# ------------------------------------------------------------------------------------------------
# random-falsification baseline
# ------------------------------------------------------------------------------------------------
def _sample_domain(spec, n, rng):
    """Uniform samples of the spec domain by rejection from its bounding box."""
    out = []
    got = 0
    while got < n:
        cand = [rng.uniform(lo, hi, size=4 * n) for lo, hi in spec.u_box]
        ok = in_domain(spec, cand)
        sel = [c[ok] for c in cand]
        out.append(sel); got += sel[0].shape[0]
    return [np.concatenate([o[k] for o in out])[:n] for k in range(len(spec.u_box))]


def random_falsify(S, Xi, E, n, rng):
    """Returns (flag at strong tol, continuous score = max over specs of -minG/tol_strong)."""
    score = -np.inf
    for sp in S.specs:
        u = _sample_domain(sp, n, rng)
        g = G_points(sp, Xi, E, u).min()
        score = max(score, -g / sp.tol_strong)
    return bool(score > 1.0), float(score)
