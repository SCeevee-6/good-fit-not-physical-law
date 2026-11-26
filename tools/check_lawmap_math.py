"""Numerical checks of the identities ROADMAP_v3 builds on.

1. min/max of a.xi over {xi : RSS(xi) <= tau} is a.xi_hat -/+ sqrt(eps a' G^-1 a), eps = tau - RSS(xi_hat).
2. Any sparse model with RSS <= tau lies in the full-library set, so its interval sits inside the full one.
3. Adding a row phi0 changes K(u) = a' G^-1 a by the Sherman-Morrison term.
4. A degree-d polynomial bounded by 1 on [-1, 1] stays below |T_d(x)| for |x| > 1 (Chebyshev).
5. With tau = (1 + delta) RSS(xi_hat) the interval width does not shrink with n. With a confidence
   threshold it shrinks like n^-1/2.
6. Two-sided growth of K(u) = phi(u)' G^-1 phi(u) outside a disk of data: T_d(g)^2/n <= K(u) <= N^2 T_d(g)^2.

Usage: python tools/check_lawmap_math.py [out_file]
"""
import sys
import numpy as np
from numpy.polynomial import Polynomial, chebyshev


def check_bound(rng, n=200, p=10, ndir=200_000):
    X = rng.normal(size=(n, p)); y = X @ rng.normal(size=p) + 0.3 * rng.normal(size=n)
    G = X.T @ X; xh = np.linalg.solve(G, X.T @ y)
    rss = lambda b: float(np.sum((X @ b - y) ** 2))
    tau = 1.1 * rss(xh); eps = tau - rss(xh)
    a = rng.normal(size=p); w = np.sqrt(eps * a @ np.linalg.solve(G, a))
    L = np.linalg.cholesky(np.linalg.inv(G))
    Z = rng.normal(size=(ndir, p)); Z /= np.linalg.norm(Z, axis=1)[:, None]
    sampled = float(((xh + np.sqrt(eps) * Z @ L.T) @ a).min())
    xs = xh - eps * np.linalg.solve(G, a) / w
    return dict(closed=float(a @ xh - w), sampled=sampled, argmin_on_boundary=rss(xs) / tau,
                ok=abs(a @ xs - (a @ xh - w)) < 1e-9 and sampled >= a @ xh - w - 1e-9
                and abs(rss(xs) / tau - 1) < 1e-9)


def check_sparse(rng, n=300, p=10, trials=500):
    X = rng.normal(size=(n, p)); bt = np.zeros(p); bt[[0, 2, 5]] = [1, -2, .5]
    y = X @ bt + 0.3 * rng.normal(size=n); G = X.T @ X; xh = np.linalg.solve(G, X.T @ y)
    tau = 1.1 * np.sum((X @ xh - y) ** 2); eps = tau - np.sum((X @ xh - y) ** 2)
    ok = 0
    for _ in range(trials):
        a = rng.normal(size=p); w = np.sqrt(eps * a @ np.linalg.solve(G, a))
        extra = rng.choice([1, 3, 4, 6, 7, 8, 9], size=rng.integers(0, 4), replace=False)
        S = [0, 2, 5] + list(extra); Xs = X[:, S]; bs = np.linalg.lstsq(Xs, y, rcond=None)[0]
        es = tau - np.sum((Xs @ bs - y) ** 2); As = a[S]
        ws = np.sqrt(es * As @ np.linalg.solve(Xs.T @ Xs, As))
        ok += (a @ xh - w <= As @ bs - ws + 1e-10) and (As @ bs + ws <= a @ xh + w + 1e-10)
    return dict(contained=int(ok), trials=trials, ok=ok == trials)


def check_rank_one(rng, p=10):
    X = rng.normal(size=(200, p)); G = X.T @ X; Gi = np.linalg.inv(G)
    a = rng.normal(size=p); f = rng.normal(size=p)
    sm = a @ Gi @ a - (a @ Gi @ f) ** 2 / (1 + f @ Gi @ f)
    direct = a @ np.linalg.solve(G + np.outer(f, f), a)
    return dict(abs_err=float(abs(sm - direct)), ok=abs(sm - direct) < 1e-10)


def check_chebyshev(rng, d=5, trials=2000, xs=(1.3, 2.0, 3.5)):
    grid = np.linspace(-1, 1, 4001); Td = chebyshev.Chebyshev.basis(d); worst = 0.0
    for _ in range(trials):
        P = Polynomial(rng.normal(size=d + 1)); P = P / np.abs(P(grid)).max()
        worst = max(worst, max(abs(P(x)) / abs(Td(x)) for x in xs))
    return dict(max_ratio=float(worst), ok=worst <= 1 + 1e-9)


def check_width_vs_n(rng, ns=(100, 1000, 10000), sigma=0.3, chi2=7.81):
    rows = []
    for n in ns:
        X = rng.normal(size=(n, 3)); y = X @ np.ones(3) + sigma * rng.normal(size=n); G = X.T @ X
        r = np.sum((X @ np.linalg.solve(G, X.T @ y) - y) ** 2); K = np.ones(3) @ np.linalg.solve(G, np.ones(3))
        rows.append((n, float(np.sqrt(0.1 * r * K)), float(np.sqrt(sigma ** 2 * chi2 * K))))
    ok = rows[-1][1] > 0.5 * rows[0][1] and rows[-1][2] < 0.2 * rows[0][2]
    return dict(rows=rows, ok=ok)


def check_horizon_bound(rng, d=4, n=400, ngrid=301):
    """6. Data in the unit disk K, full degree-d library in 2 variables, h^2/eps = K(u) = phi' G^-1 phi.
    Claim: T_d(g)^2 / n <= K(u) <= N^2 T_d(g)^2 for |u| = g >= 1, with N^2 = sup_K K."""
    expo = [(i, j) for i in range(d + 1) for j in range(d + 1 - i)]
    phi = lambda P: np.stack([P[:, 0] ** i * P[:, 1] ** j for i, j in expo], 1)
    r = np.sqrt(rng.uniform(size=n)); t = rng.uniform(0, 2 * np.pi, n)
    D = np.c_[r * np.cos(t), r * np.sin(t)]; Gi = np.linalg.inv(phi(D).T @ phi(D))
    Kf = lambda P: np.einsum("ij,jk,ik->i", phi(P), Gi, phi(P))
    g1 = np.linspace(-1, 1, ngrid); XX, YY = np.meshgrid(g1, g1); P = np.c_[XX.ravel(), YY.ravel()]
    P = P[np.hypot(P[:, 0], P[:, 1]) <= 1]; th = np.linspace(0, 2 * np.pi, 4000, endpoint=False)
    N2 = max(Kf(P).max(), Kf(np.c_[np.cos(th), np.sin(th)]).max())
    Td = chebyshev.Chebyshev.basis(d); lo_ok = hi_ok = True; worst_lo = worst_hi = 0.0
    for g in (1.0, 1.2, 1.5, 2.0, 3.0):
        U = g * np.c_[np.cos(th[::40]), np.sin(th[::40])]; k = Kf(U); T2 = Td(g) ** 2
        worst_lo = max(worst_lo, (T2 / n / k).max()); worst_hi = max(worst_hi, (k / (N2 * T2)).max())
    return dict(max_lower_ratio=float(worst_lo), max_upper_ratio=float(worst_hi), N2=float(N2),
                ok=worst_lo <= 1 + 1e-9 and worst_hi <= 1 + 1e-6)


def main(out="results/check_lawmap_math.txt"):
    rng = np.random.default_rng(20261126)
    res = {"1 bound": check_bound(rng), "2 sparse": check_sparse(rng), "3 rank-one": check_rank_one(rng),
           "4 chebyshev": check_chebyshev(rng), "5 width vs n": check_width_vs_n(rng),
           "6 horizon bound": check_horizon_bound(rng)}
    lines = [f"{k}: {v}" for k, v in res.items()]
    allok = all(v["ok"] for v in res.values())
    lines.append(f"ALL_OK {allok}")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return allok


if __name__ == "__main__":
    ok = main(*sys.argv[1:])
    print("ALL_OK", ok)
    sys.exit(0 if ok else 1)
