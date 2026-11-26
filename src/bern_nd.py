"""Bernstein-form branch-and-bound on the d-simplex {x >= 0, sum x <= L}, any d (the 3-D version in bern.py is the d = 3 case).

A polynomial G of degree <= n in d variables is homogenised to a form of degree n in d + 1 barycentric coordinates; its Bernstein
coefficients bound it from below on the simplex. Subdivision is longest-edge bisection as a precomputed linear map on coefficient
vectors. 'verified' means every leaf has min Bernstein coefficient >= -tol + slack (float64). 'counterexample' is returned only for a
witness (vertex or centroid of a leaf) that the caller re-checks exactly, so a float error cannot make a false counterexample.
Input is a dict {exponent tuple (length d): coefficient}.
"""
import itertools, math, time
import numpy as np
import sympy as sp

_TAB = {}


def _tables(n, d):
    if (n, d) in _TAB:
        return _TAB[(n, d)]
    k = d + 1
    lam = sp.symbols(f'l0:{k}'); mu = sp.symbols(f'm0:{k}')
    alphas = [a for a in itertools.product(range(n + 1), repeat=k) if sum(a) == n]
    idx = {a: i for i, a in enumerate(alphas)}
    hom = {}
    for e in itertools.product(range(n + 1), repeat=d):
        m = sum(e)
        if m > n:
            continue
        expr = (sp.Integer(1) * sp.prod([lam[i + 1] ** e[i] for i in range(d)])) * (sum(lam)) ** (n - m)
        v = np.zeros(len(alphas))
        for mon, coef in sp.Poly(sp.expand(expr), *lam).terms():
            v[idx[mon]] = float(coef)
        hom[e] = v
    maps = {}
    for i, j in itertools.combinations(range(k), 2):
        for child in (0, 1):
            sub = {l: l for l in lam}
            if child == 0:
                sub[lam[i]] = mu[i] + mu[j] / 2; sub[lam[j]] = mu[j] / 2
            else:
                sub[lam[i]] = mu[i] / 2; sub[lam[j]] = mu[j] + mu[i] / 2
            for q in range(k):
                if q not in (i, j):
                    sub[lam[q]] = mu[q]
            M = np.zeros((len(alphas), len(alphas)))
            for col, alpha in enumerate(alphas):
                expr = sp.Integer(1)
                for q in range(k):
                    expr *= sub[lam[q]] ** alpha[q]
                for mon, coef in sp.Poly(sp.expand(expr), *mu).terms():
                    M[idx[mon], col] = float(coef)
            maps[(i, j, child)] = M
    wfac = np.array([np.prod([math.factorial(x) for x in a]) / math.factorial(n) for a in alphas])
    vert_idx = [idx[tuple(n if q == v else 0 for q in range(k))] for v in range(k)]
    _TAB[(n, d)] = dict(alphas=alphas, hom=hom, maps=maps, wfac=wfac, vert_idx=vert_idx)
    return _TAB[(n, d)]


def bernstein_nd(poly, d, L, tol, check, time_budget=60.0, max_nodes=3_000_000, slack=1e-9):
    """poly: {exponent tuple: coef}. check(u_float_list) -> True if u is an exact counterexample (G(u) < -tol, u in domain).
    Returns dict(status in verified/counterexample/timeout, u, nodes, time)."""
    t0 = time.time()
    k = d + 1
    nz = {e: c for e, c in poly.items() if c != 0.0}
    n = max([sum(e) for e in nz], default=0)
    out = dict(status=None, u=None, nodes=0, time=0.0)
    if n == 0:
        c0 = nz.get(tuple([0] * d), 0.0)
        if c0 >= -tol:
            out['status'] = 'verified'
        else:
            u0 = [0.0] * d
            out['status'] = 'counterexample' if check(u0) else 'timeout'; out['u'] = u0
        return out
    T = _tables(n, d)
    f = np.zeros(len(T['alphas']))
    for e, c in nz.items():
        f += c * (L ** sum(e)) * T['hom'][e]
    verts = np.zeros((k, d))
    for q in range(1, k):
        verts[q, q - 1] = L
    stack = [(f, verts)]
    nodes = 0
    thr = -tol + slack
    pairs = list(itertools.combinations(range(k), 2))
    while stack:
        f, V = stack.pop()
        nodes += 1
        if nodes > max_nodes or time.time() - t0 > time_budget:
            out.update(status='timeout', nodes=nodes, time=time.time() - t0); return out
        if (f * T['wfac']).min() >= thr:
            continue
        cands = [(f[T['vert_idx'][q]], V[q]) for q in range(k)] + [(f.sum() / k ** n, V.mean(axis=0))]
        for val, pt in cands:
            if val < -tol - 1e-12 and check([float(x) for x in pt]):
                out.update(status='counterexample', u=[float(x) for x in pt], nodes=nodes, time=time.time() - t0); return out
        best, bi = -1.0, None
        for i, j in pairs:
            dd = float(np.sum((V[i] - V[j]) ** 2))
            if dd > best:
                best, bi = dd, (i, j)
        i, j = bi
        mid = 0.5 * (V[i] + V[j])
        V1 = V.copy(); V1[j] = mid
        V2 = V.copy(); V2[i] = mid
        stack.append((T['maps'][(i, j, 0)] @ f, V1)); stack.append((T['maps'][(i, j, 1)] @ f, V2))
    out.update(status='verified', nodes=nodes, time=time.time() - t0)
    return out
