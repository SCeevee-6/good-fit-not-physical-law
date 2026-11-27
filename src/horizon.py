"""0F-3: two-sided bound on the half-width h outside the data (protocol 0F, ROADMAP_v3 item 7).

K_data = {(x - mu)' C^-1 (x - mu) <= R^2}, gauge g(x) = Mahalanobis(x) / R. For a law with weights c on a full
degree-d library:
  upper (any law):            h^2 <= eps N^2 T_d(g)^2,  N^2 = sup over K_data of K
  lower (single equation j):  h^2 >= eps lambda_min c_j^2 T_d(g)^2 / n
with lambda_min the smallest eigenvalue of the stacked residual covariance and n the stacked row count.
"""
import numpy as np
from numpy.polynomial import chebyshev

G_BAND = (1.25, 1.5, 2.0, 3.0, 4.0)
G_FINE = np.round(np.arange(1.0, 8.0001, 0.01), 4)
N_DIR = 200
N_SUP = 20000


def ellipsoid(Xtr):
    mu = Xtr.mean(0); C = np.cov(Xtr.T); W = np.linalg.cholesky(C)
    R = float(np.sqrt((np.linalg.solve(W, (Xtr - mu).T) ** 2).sum(0)).max())
    return mu, W, R


def sup_K(ds, c, mu, W, R, rng):
    """Lower estimate of sup over K_data of K = h^2 / eps: uniform sample plus Nelder-Mead from the best 20."""
    from scipy.optimize import minimize
    n = len(mu)
    V = rng.normal(size=(N_SUP, n)); V /= np.linalg.norm(V, axis=1)[:, None]
    V *= rng.uniform(size=(N_SUP, 1)) ** (1 / n)
    X = mu + R * V @ W.T
    _, h = ds.law(c, X, 0.0); k = h ** 2 / ds.eps
    best = float(np.nanmax(k))
    def neg(v):
        v = np.asarray(v); nv = np.linalg.norm(v)
        if nv > 1:
            v = v / nv
        _, hh = ds.law(c, (mu + R * W @ v)[None], 0.0)
        return -float(hh[0] ** 2 / ds.eps)
    for i in np.argsort(-k)[:20]:
        r = minimize(neg, V[i], method='Nelder-Mead', options=dict(maxiter=200, xatol=1e-6, fatol=1e-12))
        best = max(best, -float(r.fun))
    return best


def lambda_min(trajs, S, E, sigma, Xi_ref):
    """Smallest eigenvalue over trajectories of the stacked residual covariance used by StackedSet."""
    from lawmap import sg_mats, _jacobian_rows
    n = E.shape[1]; lm = np.inf
    for tr in trajs:
        X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S); J = _jacobian_rows(tr['Xs'], E, Xi_ref)
        mi, m = D.shape; B = np.zeros((n * mi, n * m))
        for j in range(n):
            for k in range(n):
                blk = -J[:, j, k][:, None] * H + (D if j == k else 0)
                B[j * mi:(j + 1) * mi, k * m:(k + 1) * m] = blk * sigma[k]
        lm = min(lm, float(np.linalg.eigvalsh(B @ B.T)[0]))
    return lm


def test_poly_denominator(ds, E, j, u, mu, W, R, rng):
    """||Z delta||^2 for delta = coefficients of p(x) = T_d(l(x)) in equation j, l(x) = u' W^-1 (x - mu) / R.
    Gives the exact lower bound K(x) >= c_j^2 T_d(g)^2 / ||Z delta||^2 along the ray mu + g R W u."""
    from lib import theta
    d = int(E.sum(1).max()); Td = chebyshev.Chebyshev.basis(d); n = len(mu)
    P = mu + rng.normal(size=(6 * E.shape[0], n)) * np.abs(W).sum(1) * R
    l = (np.linalg.solve(W, (P - mu).T).T @ u) / R
    coef = np.linalg.lstsq(theta(P, E), Td(l), rcond=None)[0]
    pk = int(ds.keep.sum()); dv = np.zeros(ds.n * pk); dv[j * pk:(j + 1) * pk] = coef[ds.keep]
    return float(((ds.Z @ dv) ** 2).sum())


def horizon_job(ds, trajs, S, E, sigma, Xi_ref, specs, rng):
    from systems import SYSTEMS  # noqa: F401
    Xtr = np.vstack([t['Xs'] for t in trajs]); mu, W, R = ellipsoid(Xtr)
    d = int(E.sum(1).max()); Td = chebyshev.Chebyshev.basis(d)
    nrows = ds.Z.shape[0]; lm = lambda_min(trajs, S, E, sigma, Xi_ref)
    U = rng.normal(size=(N_DIR, len(mu))); U /= np.linalg.norm(U, axis=1)[:, None]
    out = {}
    for sp in specs:
        c = np.array([float(x) for x in sp.c]); nz = np.nonzero(c)[0]
        if len(nz) != 1:
            continue
        j = int(nz[0]); N2 = sup_K(ds, sp.c, mu, W, R, rng)
        n_in = n_tot = 0; worst_lo = worst_hi = 0.0; widths = []; cens = 0
        n_in_x = 0; worst_lx = 0.0; widths_x = []; cens_x = 0
        for u in U:
            den = test_poly_denominator(ds, E, j, u, mu, W, R, rng)
            for g in G_BAND:
                x = mu + g * R * (W @ u)
                _, h = ds.law(sp.c, x[None], 0.0); h2 = float(h[0] ** 2)
                lo = ds.eps * lm * c[j] ** 2 * Td(g) ** 2 / nrows; hi = ds.eps * N2 * Td(g) ** 2
                lx = ds.eps * c[j] ** 2 * Td(g) ** 2 / den
                n_tot += 1; n_in += (lo <= h2 * (1 + 1e-9)) and (h2 <= hi * (1 + 1e-9))
                n_in_x += (lx <= h2 * (1 + 1e-6)) and (h2 <= hi * (1 + 1e-9))
                worst_lo = max(worst_lo, lo / h2); worst_hi = max(worst_hi, h2 / hi); worst_lx = max(worst_lx, lx / h2)
            Xr = mu + np.outer(G_FINE, R * (W @ u))
            m, h = ds.law(sp.c, Xr, sp.tol_strong); T2 = Td(G_FINE) ** 2
            up = np.sqrt(ds.eps * N2 * T2); lw = np.sqrt(ds.eps * lm * c[j] ** 2 * T2 / nrows)
            a = np.nonzero(np.abs(m) < up)[0]; b = np.nonzero(np.abs(m) < lw)[0]
            if len(a) and len(b):
                widths.append(float(np.log(G_FINE[b[0]] / G_FINE[a[0]])))
            else:
                cens += 1
            bx = np.nonzero(np.abs(m) < np.sqrt(ds.eps * c[j] ** 2 * T2 / den))[0]
            if len(a) and len(bx):
                widths_x.append(float(np.log(G_FINE[bx[0]] / G_FINE[a[0]])))
            else:
                cens_x += 1
        out[sp.name] = dict(n_points=n_tot, n_in_band=n_in, max_lower_ratio=worst_lo, max_upper_ratio=worst_hi,
                            N2=N2, lambda_min=lm, nrows=nrows, R=R, median_log_bracket=float(np.median(widths)) if widths else None,
                            n_bracket=len(widths), n_censored=cens,
                            n_in_band_exact=n_in_x, max_lower_exact_ratio=worst_lx,
                            median_log_bracket_exact=float(np.median(widths_x)) if widths_x else None,
                            n_bracket_exact=len(widths_x), n_censored_exact=cens_x)
    return out
