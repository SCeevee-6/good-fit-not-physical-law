"""Bernstein-form branch-and-bound on the 3-simplex  {x >= 0, sum x <= L}.

Independent decision procedure used (a) as a fallback when z3 times out, and (b) as a cross-check of z3 on the 3-D queries.
Polynomial G of degree <= n is homogenised to a form F of degree n in 4 barycentric coordinates; its Bernstein coefficients
b_alpha = f_alpha * alpha!/n!  bound F from below on the simplex.  Subdivision = longest-edge bisection, implemented as a precomputed
linear map on coefficient vectors.  'Verified' means every leaf has min Bernstein coefficient >= -tol + slack (float64, slack 1e-9);
'counterexample' is only returned after exact rational re-evaluation of the witness.
"""
import itertools, math, time
from fractions import Fraction
import numpy as np
import sympy as sp

from lib import to_fraction

_TAB = {}


def _tables(n):
    if n in _TAB:
        return _TAB[n]
    lam = sp.symbols('l0:4'); mu = sp.symbols('m0:4')
    alphas = [a for a in itertools.product(range(n + 1), repeat=4) if sum(a) == n]
    idx = {a: i for i, a in enumerate(alphas)}
    # homogenisation vectors for monomials S^a I^b R^c (unscaled by L^m):  l1^a l2^b l3^c (sum l)^(n-m)
    hom = {}
    for a in range(n + 1):
        for b in range(n + 1 - a):
            for c in range(n + 1 - a - b):
                m = a + b + c
                expr = sp.expand(lam[1] ** a * lam[2] ** b * lam[3] ** c * (lam[0] + lam[1] + lam[2] + lam[3]) ** (n - m))
                v = np.zeros(len(alphas))
                for mon, coef in sp.Poly(expr, *lam).terms():
                    v[idx[mon]] = float(coef)
                hom[(a, b, c)] = v
    # subdivision maps
    maps = {}
    for i, j in itertools.combinations(range(4), 2):
        for child in (0, 1):
            sub = {l: l for l in lam}
            if child == 0:   # replace vertex j by midpoint(i, j)
                sub = {lam[i]: mu[i] + mu[j] / 2, lam[j]: mu[j] / 2}
            else:            # replace vertex i by midpoint(i, j)
                sub = {lam[i]: mu[i] / 2, lam[j]: mu[j] + mu[i] / 2}
            for k in range(4):
                if k not in (i, j):
                    sub[lam[k]] = mu[k]
            M = np.zeros((len(alphas), len(alphas)))
            for col, alpha in enumerate(alphas):
                expr = sp.Integer(1)
                for k in range(4):
                    expr *= sub[lam[k]] ** alpha[k]
                for mon, coef in sp.Poly(sp.expand(expr), *mu).terms():
                    M[idx[mon], col] = float(coef)
            maps[(i, j, child)] = M
    wfac = np.array([np.prod([math.factorial(x) for x in a]) / math.factorial(n) for a in alphas])
    vert_idx = [idx[tuple(n if k == v else 0 for k in range(4))] for v in range(4)]
    _TAB[n] = dict(alphas=alphas, idx=idx, hom=hom, maps=maps, wfac=wfac, vert_idx=vert_idx, n=n)
    return _TAB[n]


def _poly_coeffs(spec, Xi, E):
    """Monomial coefficients g_j of G in the (S,I,R) library basis (identity point map assumed)."""
    c = np.array([float(ci) for ci in spec.c])
    return Xi @ c


def bernstein_verify(spec, Xi, E, level='strong', time_budget=60.0, max_nodes=2_000_000, slack=1e-9):
    from verify import _exact_G, _in_domain_exact
    t0 = time.time()
    assert len(spec.u_names) == 3
    L = -float(spec.domain[-1][1]) / float(spec.domain[-1][0][0]) * -1  # domain row: -(u1+u2+u3) >= -L  ->  b=-L, a=-1
    L = float(-spec.domain[-1][1])
    tol = spec.tol(level)
    g = _poly_coeffs(spec, Xi, E)
    nz = [j for j in range(E.shape[0]) if g[j] != 0]
    n = max([int(E[j].sum()) for j in nz], default=0)
    out = dict(spec=spec.name, level=level, u=None, G=None, status=None, witness_exact=None, nodes=0, time=0.0)
    if n == 0:
        const = float(g[0]) if E[0].sum() == 0 else 0.0
        out['status'] = 'verified' if const >= -tol else 'counterexample'
        if out['status'] == 'counterexample':
            out.update(u=[0.0, 0.0, 0.0], G=const, witness_exact=True)
        return out
    T = _tables(n)
    f = np.zeros(len(T['alphas']))
    for j in nz:
        a, b, c_ = (int(x) for x in E[j])
        f += g[j] * (L ** (a + b + c_)) * T['hom'][(a, b, c_)]
    verts = np.array([[0, 0, 0], [L, 0, 0], [0, L, 0], [0, 0, L]], dtype=float)
    stack = [(f, verts)]
    nodes = 0
    thr = -tol + slack
    while stack:
        f, V = stack.pop()
        nodes += 1
        if nodes > max_nodes or time.time() - t0 > time_budget:
            out.update(status='timeout', nodes=nodes, time=time.time() - t0)
            return out
        b = f * T['wfac']
        if b.min() >= thr:
            continue
        # candidate witnesses: vertices and centroid
        cands = [(f[T['vert_idx'][k]], V[k]) for k in range(4)] + [(f.sum() / 4 ** n, V.mean(axis=0))]
        for val, pt in cands:
            if val < -tol - 1e-12:
                uF = [Fraction(float(x)) for x in pt]
                Gx = _exact_G(spec, Xi, E, uF)
                if _in_domain_exact(spec, uF) and Gx < -Fraction(tol).limit_denominator(10 ** 6):
                    out.update(status='counterexample', u=[float(x) for x in uF], G=float(Gx), witness_exact=True, nodes=nodes, time=time.time() - t0)
                    return out
        # bisect longest edge
        best, bi = -1.0, None
        for i, j in itertools.combinations(range(4), 2):
            d = float(np.sum((V[i] - V[j]) ** 2))
            if d > best:
                best, bi = d, (i, j)
        i, j = bi
        mid = 0.5 * (V[i] + V[j])
        V1 = V.copy(); V1[j] = mid
        V2 = V.copy(); V2[i] = mid
        stack.append((T['maps'][(i, j, 0)] @ f, V1))
        stack.append((T['maps'][(i, j, 1)] @ f, V2))
    out.update(status='verified', nodes=nodes, time=time.time() - t0)
    return out
