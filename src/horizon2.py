"""Post-hoc unit R2 (docs/protocol_0F3b.md): the 0F-3 horizon rays again, with the projected lower bound of
docs/law_horizon.md Theorem 4 added. The random stream is consumed exactly as in horizon.horizon_job, so rays, sup
estimates and test points are identical to 0F-3, and 0F-3's own numbers are reproduced as a check.

Projected bound along a ray with supporting functional l: the span of T_0(l(x)) .. T_d(l(x)) in equation j gives
    K(x) >= c_j^2 v' M^-1 v,  v_k = T_k(g),  M_kl = (Z delta_k) . (Z delta_l),
where delta_k are the library coefficients of T_k(l(x)). The 0F-3 'exact' bound is the single term k = d.
"""
import numpy as np
from numpy.polynomial import chebyshev
from horizon import G_BAND, G_FINE, N_DIR, ellipsoid, sup_K, lambda_min


def test_poly_gram(ds, E, j, u, mu, W, R, rng):
    """Gram matrix of Z delta_k, k = 0..d (same random points as horizon.test_poly_denominator)."""
    from lib import theta
    d = int(E.sum(1).max()); n = len(mu)
    P = mu + rng.normal(size=(6 * E.shape[0], n)) * np.abs(W).sum(1) * R
    l = (np.linalg.solve(W, (P - mu).T).T @ u) / R
    T = np.array([chebyshev.Chebyshev.basis(k)(l) for k in range(d + 1)]).T
    coef = np.linalg.lstsq(theta(P, E), T, rcond=None)[0]                     # (p, d+1)
    pk = int(ds.keep.sum()); Dv = np.zeros((ds.n * pk, d + 1)); Dv[j * pk:(j + 1) * pk] = coef[ds.keep]
    Zd = ds.Z @ Dv
    return Zd.T @ Zd


def horizon_job2(ds, trajs, S, E, sigma, Xi_ref, specs, rng):
    Xtr = np.vstack([t['Xs'] for t in trajs]); mu, W, R = ellipsoid(Xtr)
    d = int(E.sum(1).max()); Tk = [chebyshev.Chebyshev.basis(k) for k in range(d + 1)]
    nrows = ds.Z.shape[0]; lm = lambda_min(trajs, S, E, sigma, Xi_ref)
    U = rng.normal(size=(N_DIR, len(mu))); U /= np.linalg.norm(U, axis=1)[:, None]
    out = {}
    for sp in specs:
        c = np.array([float(x) for x in sp.c]); nz = np.nonzero(c)[0]
        if len(nz) != 1:
            continue
        j = int(nz[0]); N2 = sup_K(ds, sp.c, mu, W, R, rng)
        n_tot = 0; bad_proj = 0; worst_proj = 0.0
        cens_x = cens_p = 0; widths_x, widths_p, gap_x, gap_p = [], [], [], []
        for u in U:
            M = test_poly_gram(ds, E, j, u, mu, W, R, rng); Mi = np.linalg.pinv(M)
            den = float(M[d, d])
            for g in G_BAND:
                x = mu + g * R * (W @ u)
                _, h = ds.law(sp.c, x[None], 0.0); h2 = float(h[0] ** 2)
                v = np.array([T(g) for T in Tk]); lp = ds.eps * c[j] ** 2 * float(v @ Mi @ v)
                n_tot += 1; bad_proj += lp > h2 * (1 + 1e-6); worst_proj = max(worst_proj, lp / h2)
            Xr = mu + np.outer(G_FINE, R * (W @ u))
            m, h = ds.law(sp.c, Xr, sp.tol_strong)
            V = np.array([T(G_FINE) for T in Tk]).T
            lx = np.sqrt(ds.eps * c[j] ** 2 * V[:, d] ** 2 / den)
            lp = np.sqrt(ds.eps * c[j] ** 2 * np.einsum('gi,ij,gj->g', V, Mi, V))
            up = np.sqrt(ds.eps * N2 * V[:, d] ** 2)
            a = np.nonzero(np.abs(m) < up)[0]; t = np.nonzero(np.abs(m) < h)[0]
            bx = np.nonzero(np.abs(m) < lx)[0]; bp = np.nonzero(np.abs(m) < lp)[0]
            for b, cens_name, widths, gap in ((bx, 'x', widths_x, gap_x), (bp, 'p', widths_p, gap_p)):
                if len(a) and len(b):
                    widths.append(float(np.log(G_FINE[b[0]] / G_FINE[a[0]])))
                    if len(t):
                        gap.append(float(np.log(G_FINE[b[0]] / G_FINE[t[0]])))
                elif cens_name == 'x':
                    cens_x += 1
                else:
                    cens_p += 1
        out[sp.name] = dict(n_points=n_tot, proj_violations=int(bad_proj), max_proj_ratio=worst_proj, N2=N2,
                            n_censored_exact=cens_x, n_bracket_exact=len(widths_x),
                            median_log_bracket_exact=float(np.median(widths_x)) if widths_x else None,
                            median_log_gap_exact=float(np.median(gap_x)) if gap_x else None,
                            n_censored_proj=cens_p, n_bracket_proj=len(widths_p),
                            median_log_bracket_proj=float(np.median(widths_p)) if widths_p else None,
                            median_log_gap_proj=float(np.median(gap_p)) if gap_p else None)
    return out
