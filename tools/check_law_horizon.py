"""Numerical checks for docs/law_horizon.md (unit T). Pass rules were fixed in that doc before the first run.

H1 Lemma 1: K(u) = max_q q(u)^2 / ||q||^2 over the library, attained at c = G^-1 phi(u) (correlated noise).
H2 Lemma 2: max |q(u)| over deg <= d, |q| <= 1 on a symmetric convex body, equals T_d(g(u)) (disk, square; LP).
H3 Theorems 3, 4: T_d(g)^2 <= T_d test <= projected bound <= Kbar <= Nbar^2 T_d(g)^2; equality of the
   projected bound in 1-D.
H4 Theorem 5: zone frequencies under exact regressors and Gaussian noise match Phi(z - sqrt q), Phi(-z - sqrt q).
H5 Corollary 6: Chebyshev-node data give Kbar = 1 + 2 sum T_k^2 and the median horizon lies in
   [T_d^-1(R / Nbar), T_d^-1(R)].
H6 Proposition 7: d > d_m gives |z| < sqrt q beyond cosh(log R_m / (d - d_m)); d = d_m has the stated limit of z.
H7 Proposition 8: restricted library (degree 2 without x^2): monotone in the library, line degree bound on the
   x axis, lower bound from T_{d'} and monomial tests.

Usage: python tools/check_law_horizon.py [out_file]
"""
import sys
import itertools
import numpy as np
from numpy.polynomial import chebyshev
from scipy.optimize import linprog, minimize, brentq
from scipy.stats import chi2, norm

TOL = 1e-9


def expo(s, d, drop=()):
    E = [e for k in range(d + 1) for e in itertools.product(range(k + 1), repeat=s) if sum(e) == k]
    E = sorted(set(E), key=lambda e: (sum(e), e))
    return np.array([e for e in E if e not in drop])


def phi(X, E):
    X = np.atleast_2d(X)
    return np.prod(X[:, None, :] ** E[None, :, :], axis=2)


def T(d, x):
    return chebyshev.chebval(x, [0] * d + [1])


def Tinv(d, y):
    return np.cosh(np.arccosh(y) / d)


def kbar(Xd, E, U, Sinv=None):
    """Kbar(u) = phi(u)' (Theta' S^-1 Theta / n)^-1 phi(u)."""
    Th = phi(Xd, E); n = len(Xd)
    M = (Th.T @ (Sinv if Sinv is not None else np.eye(n)) @ Th) / n
    P = phi(U, E)
    return np.einsum('ij,ij->i', P, np.linalg.solve(M, P.T).T)


def disk(rng, n):
    r = np.sqrt(rng.uniform(size=n)); a = rng.uniform(0, 2 * np.pi, n)
    return np.c_[r * np.cos(a), r * np.sin(a)]


def sup_on_disk(f, rng, n=20000):
    V = np.vstack([disk(rng, n), np.c_[np.cos(np.linspace(0, 2 * np.pi, 4000)), np.sin(np.linspace(0, 2 * np.pi, 4000))]])
    k = f(V); best = float(k.max())
    for i in np.argsort(-k)[:10]:
        def neg(v):
            v = np.asarray(v); nv = np.linalg.norm(v)
            return -float(f((v / nv if nv > 1 else v)[None])[0])
        best = max(best, -minimize(neg, V[i], method='Nelder-Mead', options=dict(xatol=1e-10, fatol=1e-14, maxiter=400)).fun)
    return best


def h1(rng):
    n, d = 120, 3; E = expo(2, d); Xd = disk(rng, n)
    S = 0.6 ** np.abs(np.subtract.outer(np.arange(n), np.arange(n))); Si = np.linalg.inv(S)
    Th = phi(Xd, E); G = Th.T @ Si @ Th
    U = disk(rng, 200) * rng.uniform(1, 4, (200, 1)); P = phi(U, E)
    K = np.einsum('ij,ij->i', P, np.linalg.solve(G, P.T).T)
    C = rng.normal(size=(1000, len(E)))
    ratio = (P @ C.T) ** 2 / np.einsum('ij,jk,ik->i', C, G, C)[None, :] / K[:, None]
    cs = np.linalg.solve(G, P.T).T
    att = (np.einsum('ij,ij->i', P, cs) ** 2 / np.einsum('ij,jk,ik->i', cs, G, cs)) / K
    return dict(max_random_ratio=float(ratio.max()), attain_err=float(np.abs(att - 1).max()),
                ok=bool(ratio.max() <= 1 + TOL and np.abs(att - 1).max() <= TOL))


def h2(rng):
    worst_lo, worst_hi = np.inf, 0.0
    tb = np.linspace(0, 2 * np.pi, 2000, endpoint=False)
    bodies = {
        'disk': (np.vstack([np.c_[np.cos(tb), np.sin(tb)], disk(rng, 3000)]), lambda u: np.linalg.norm(u)),
        'square': (np.vstack([np.c_[np.ones(500), np.linspace(-1, 1, 500)], np.c_[-np.ones(500), np.linspace(-1, 1, 500)],
                              np.c_[np.linspace(-1, 1, 500), np.ones(500)], np.c_[np.linspace(-1, 1, 500), -np.ones(500)],
                              rng.uniform(-1, 1, (3000, 2))]), lambda u: np.abs(u).max())}
    for name, (Grid, gauge) in bodies.items():
        for d in range(1, 6):
            E = expo(2, d); A = phi(Grid, E)
            for _ in range(3):
                e = rng.normal(size=2); e /= gauge(e)
                for g in (1.2, 2.0, 4.0):
                    u = g * e; c = -phi(u, E)[0]
                    r = linprog(c, A_ub=np.vstack([A, -A]), b_ub=np.ones(2 * len(A)), bounds=[(None, None)] * len(E), method='highs')
                    ratio = -r.fun / T(d, g)
                    worst_lo = min(worst_lo, ratio); worst_hi = max(worst_hi, ratio)
    return dict(min_LP_over_Td=float(worst_lo), max_LP_over_Td=float(worst_hi),
                ok=bool(worst_lo >= 1 - 1e-7 and worst_hi <= 1.02))


def proj_bound(Xd, u, ell, d):
    t = Xd @ ell; tu = float(u @ ell)
    V = np.array([T(k, t) for k in range(d + 1)]).T; M = V.T @ V / len(Xd)
    v = np.array([T(k, tu) for k in range(d + 1)])
    return float(v @ np.linalg.solve(M, v)), float(T(d, tu) ** 2 / np.mean(T(d, t) ** 2))


def h3(rng):
    out = {}; ok = True
    Xd = disk(rng, 300)
    for d in (2, 3, 4):
        E = expo(2, d); f = lambda V: kbar(Xd, E, V)
        N2 = sup_on_disk(f, rng)
        bad = 0; r_proj, r_test = [], []
        for _ in range(400):
            e = rng.normal(size=2); e /= np.linalg.norm(e); g = rng.uniform(1.25, 8)
            u = g * e; K = f(u[None])[0]; pb, tb = proj_bound(Xd, u, e, d)
            pr = max([pb] + [proj_bound(Xd, u, l / np.linalg.norm(l), d)[0] for l in rng.normal(size=(3, 2))])
            chain = [T(d, g) ** 2, tb, pb, K, N2 * T(d, g) ** 2]
            bad += any(chain[i] > chain[i + 1] * (1 + TOL) for i in range(4)) or pr > K * (1 + TOL)
            r_proj.append(K / pb); r_test.append(K / tb)
        out[f'd{d}'] = dict(N2=round(N2, 3), p=len(E), violations=bad, med_K_over_proj=round(float(np.median(r_proj)), 3),
                            med_K_over_Tdtest=round(float(np.median(r_test)), 3))
        ok &= bad == 0 and N2 >= len(E) * (1 - TOL)
    x1 = rng.uniform(-1, 1, (200, 1)); errs = []
    for d in (2, 3, 4):
        E = expo(1, d)
        for g in (1.5, 3.0, 6.0):
            errs.append(abs(proj_bound(x1, np.array([g]), np.array([1.0]), d)[0] / kbar(x1, E, np.array([[g]]))[0] - 1))
    out['oneD_max_rel_err'] = float(max(errs)); ok &= max(errs) <= 1e-8
    out['ok'] = bool(ok); return out


def h4(rng, reps=40000, n=200, sigma=0.05):
    E = expo(2, 2); Xd = disk(rng, n); Th = phi(Xd, E); p = len(E); q = chi2.ppf(0.95, p)
    xi = np.array([0.3, -0.2, 0.1, 0.05, -0.08, 0.1])
    G = Th.T @ Th / sigma ** 2; A = np.linalg.solve(Th.T @ Th, Th.T)
    cand = disk(rng, 4000) * rng.uniform(0.5, 4, (4000, 1)); P = phi(cand, E)
    K = np.einsum('ij,ij->i', P, np.linalg.solve(G, P.T).T); z = P @ xi / np.sqrt(K)
    pick = []
    for lo in np.linspace(-4, 3.6, 20):
        idx = np.nonzero((z >= lo) & (z < lo + 0.4))[0]
        if len(idx):
            pick.append(idx[0])
    P, K, z = P[pick], K[pick], z[pick]
    worst = 0.0
    for b in range(0, reps, 10000):
        Y = Th @ xi + sigma * rng.normal(size=(10000, n))
        M = Y @ (A.T @ P.T)
        cert = (M >= np.sqrt(q * K)).mean(0); ref = (M < -np.sqrt(q * K)).mean(0)
        if b == 0:
            C, Rf = cert * 0, ref * 0
        C += cert / (reps / 10000); Rf += ref / (reps / 10000)
    for f, pth in ((C, norm.cdf(z - np.sqrt(q))), (Rf, norm.cdf(-z - np.sqrt(q)))):
        se = np.sqrt(pth * (1 - pth) / reps)
        worst = max(worst, float(np.max(np.abs(f - pth) / (4 * se + 1e-3))))
    return dict(points=len(pick), z_range=[round(float(z.min()), 2), round(float(z.max()), 2)],
                worst_scaled_dev=round(worst, 3), ok=bool(worst <= 1))


def h5(rng):
    out = []; ok = True
    for n in (50, 200, 800):
        x = np.cos((2 * np.arange(1, n + 1) - 1) * np.pi / (2 * n))[:, None]
        for d in (2, 3, 4):
            E = expo(1, d); q = chi2.ppf(0.95, d + 1); Nb = np.sqrt(2 * d + 1)
            gg = np.linspace(1, 20, 400)
            Kc = 1 + 2 * sum(T(k, gg) ** 2 for k in range(1, d + 1)); Kn = kbar(x, E, gg[:, None])
            ok &= float(np.max(np.abs(Kn / Kc - 1))) <= 1e-8
            for mu0 in (0.1, 0.3):
                sigma = 0.05; Rr = mu0 * np.sqrt(n) / (sigma * np.sqrt(q))
                zf = lambda g: mu0 * np.sqrt(n) / (sigma * np.sqrt(kbar(x, E, np.array([[g]]))[0])) - np.sqrt(q)
                gmed = brentq(zf, 1.0, 1e6) if zf(1.0) > 0 else 1.0
                glo = Tinv(d, Rr / Nb) if Rr / Nb >= 1 else 1.0; ghi = Tinv(d, Rr)
                ok &= glo * (1 - TOL) <= gmed <= ghi * (1 + TOL)
                out.append((n, d, mu0, round(glo, 4), round(gmed, 4), round(ghi, 4)))
    return dict(rows_n_d_mu0_glo_gmed_ghi=out, ok=bool(ok))


def h6(rng, n=200, sigma=0.05):
    x = np.cos((2 * np.arange(1, n + 1) - 1) * np.pi / (2 * n))[:, None]
    mu = lambda g: 0.2 + 0.5 * g; C = 0.7; res = {}; ok = True
    for d in (1, 2, 3):
        E = expo(1, d); q = chi2.ppf(0.95, d + 1)
        z = lambda g: np.sqrt(n) * mu(g) / (sigma * np.sqrt(kbar(x, E, np.atleast_1d(g)[:, None])))
        if d == 1:
            Th = phi(x, E); Mi = np.linalg.inv(Th.T @ Th / n); zinf = np.sqrt(n) * 0.5 / (sigma * np.sqrt(Mi[1, 1]))
            rel = abs(z(1e4)[0] / zinf - 1); gs = np.geomspace(1, 1e4, 2000)
            res['d1'] = dict(z_inf=round(float(zinf), 3), rel_err_at_1e4=float(rel), min_z=round(float(z(gs).min()), 3))
            ok &= rel <= 1e-3 and z(gs).min() >= np.sqrt(q)
        else:
            Rm = 2 * C * np.sqrt(n) / (sigma * np.sqrt(q)); gb = np.cosh(np.log(Rm) / (d - 1))
            gs = np.geomspace(gb, 1e4, 2000); zmax = float(z(gs).max())
            cross = brentq(lambda g: z(g)[0] - np.sqrt(q), 1.0, gb * 1.0000001) if z(1.0)[0] > np.sqrt(q) else 1.0
            res[f'd{d}'] = dict(bound_radius=round(float(gb), 3), actual_crossing=round(float(cross), 3),
                                max_z_beyond=round(zmax, 3), sqrt_q=round(float(np.sqrt(q)), 3))
            ok &= zmax < np.sqrt(q)
    res['ok'] = bool(ok); return res


def h7(rng, n=300):
    Xd = disk(rng, n); Ef = expo(2, 2); Em = expo(2, 2, drop=((2, 0),))
    U = disk(rng, 2000) * rng.uniform(0, 8, (2000, 1))
    Kf, Km = kbar(Xd, Ef, U), kbar(Xd, Em, U)
    a_ok = bool(np.all(Km <= Kf * (1 + TOL)))
    N2m = sup_on_disk(lambda V: kbar(Xd, Em, V), rng); N2f = sup_on_disk(lambda V: kbar(Xd, Ef, V), rng)
    gs = np.linspace(1.25, 8, 200); Ux = np.vstack([np.c_[gs, 0 * gs], np.c_[-gs, 0 * gs]]); g2 = np.r_[gs, gs]
    b_ok = bool(np.all(kbar(Xd, Em, Ux) <= N2m * T(1, g2) ** 2 * (1 + TOL)))
    full_exceeds = int(np.sum(kbar(Xd, Ef, Ux) > N2m * T(1, g2) ** 2))
    gU = np.linalg.norm(U, axis=1); sel = gU >= 1
    mono = np.max(phi(U[sel], Em[1:]) ** 2 / np.array([1, 1, 1, 0.25])[None, :], axis=1)  # Em[1:] = y, x, y^2, xy; max over disk of squares
    lo = np.maximum(T(1, gU[sel]) ** 2, mono)
    c_ok = bool(np.all(lo <= Km[sel] * (1 + TOL)))
    return dict(monotone=a_ok, line_degree_bound=b_ok, full_library_points_above_line_bound=full_exceeds,
                lower_bound=c_ok, N2_restricted=round(N2m, 3), N2_full=round(N2f, 3), ok=a_ok and b_ok and c_ok)


def main(out="results/check_law_horizon.txt"):
    rng = np.random.default_rng(20261127)
    res = {}
    for name, f in (("H1 lemma1", h1), ("H2 lemma2", h2), ("H3 band+proj", h3), ("H4 zone probs", h4),
                    ("H5 closed form", h5), ("H6 finite horizon", h6), ("H7 restricted lib", h7)):
        res[name] = f(rng)
    lines = [f"{k}: {v}" for k, v in res.items()]
    allok = all(v["ok"] for v in res.values())
    lines.append(f"ALL_OK {allok}")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return allok, {k: v["ok"] for k, v in res.items()}


if __name__ == "__main__":
    ok, per = main(*sys.argv[1:])
    print("ALL_OK", ok, per)
    sys.exit(0 if ok else 1)
