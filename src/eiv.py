"""Protocol 0G: confidence set from a bias-corrected (corrected-score) estimate for noisy regressors.

Smoothed states Xs = H X carry Gaussian noise u = H e with per-row variance v_tk = sigma_k^2 sum_l H_tl^2.
For the probabilists' Hermite polynomials with variance v (He_0 = 1, He_1 = y, He_{m+1} = y He_m - m v He_{m-1}),
E[He_m(x + u; v)] = x^m. Replacing every monomial of the Gram matrix and of the cross term by its Hermite version
gives estimating equations with expectation equal to their noise-free value:
    G~[a, b] = sum_t prod_k He_{a_k + b_k}(xs_tk; v_tk),   b~[a, j] = sum_t prod_k He_{a_k}(xs_tk; v_tk) dX_tj.
The cross term needs no correction because on interior rows of a symmetric Savitzky-Golay window sum_l H_tl D_tl = 0,
so (H e)_t and (D e)_t are independent Gaussians. The covariance is the delta-method sandwich
    V = (I (x) G~)^+ [Z' C Z] (I (x) G~)^+,  C = B B',
with B the linearised residual model of lawmap.StackedSet (Jacobian at the corrected estimate). Same interface as
StackedSet (Xi_hat, eps, rank, p, law, sample_boundary, quad_form), plus stat(Xi) for coverage.
"""
import numpy as np
from lib import theta


def hermite_table(y, v, m):
    """He_0..He_m of y (rows) with variances v (same shape). Returns (m + 1, len(y))."""
    out = np.empty((m + 1,) + y.shape); out[0] = 1.0
    if m >= 1:
        out[1] = y
    for k in range(1, m):
        out[k + 1] = y * out[k] - k * v * out[k - 1]
    return out


def corrected_moments(trajs, S, E, sigma):
    """G~ (p, p) and b~ (p, n) summed over trajectories, plus plain Theta' Theta for scaling."""
    from lawmap import sg_mats
    p, n = E.shape; d2 = 2 * int(E.sum(1).max())
    G = np.zeros((p, p)); b = np.zeros((p, n)); G0 = np.zeros((p, p))
    for tr in trajs:
        X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S)
        Xs, dX = np.asarray(tr['Xs'], float), np.asarray(tr['dX'], float)
        hv = (H ** 2).sum(1)                                   # per-row filter gain
        tabs = [hermite_table(Xs[:, k], sigma[k] ** 2 * hv, d2) for k in range(n)]
        Ga = np.ones((p, p, len(Xs)))
        for k in range(n):
            Ga *= tabs[k][E[:, k][:, None] + E[:, k][None, :]]
        G += Ga.sum(2)
        Tb = np.ones((p, len(Xs)))
        for k in range(n):
            Tb *= tabs[k][E[:, k]]
        b += Tb @ dX
        Th = theta(Xs, E); G0 += Th.T @ Th
    return G, b, G0


class CorrectedSet:
    """Reliability of a direction = generalized eigenvalue of the corrected Gram against the plain Gram (the share
    of that direction's signal left after removing the noise). Directions with reliability below kappa (default:
    non-positive) are treated as unidentified: infinite width wherever a law depends on them, no coverage claim.
    Development showed that larger kappa is not valid (the eigenbasis of low-reliability directions is unstable)."""

    def __init__(self, trajs, S, E, sigma, keep=None, alpha=0.05, kappa=1e-6, rcond=1e-12, null_tol=1e-8):
        from scipy.stats import chi2
        from lawmap import sg_mats, _jacobian_rows
        self.E = E; n = E.shape[1]; self.n = n; sig = np.asarray(sigma, float); self.sigma = sig; self.kappa = kappa
        self.keep = np.ones(E.shape[0], bool) if keep is None else np.asarray(keep, bool)
        Gf, bf, G0 = corrected_moments(trajs, S, E, sig)
        k = self.keep; G, bb, G0 = Gf[np.ix_(k, k)], bf[k], G0[np.ix_(k, k)]
        pk = int(k.sum()); self.p = pk
        col = np.sqrt(np.diag(G0)); col[col == 0] = 1.0; self.col = col
        Gs = G / col[:, None] / col[None, :]; G0s = G0 / col[:, None] / col[None, :]; self.G0s = G0s
        mu, W = np.linalg.eigh(G0s); k0 = mu > rcond * mu.max()
        Q = W[:, k0] / np.sqrt(mu[k0]); self.W0null = W[:, ~k0]
        lam, Rm = np.linalg.eigh(Q.T @ Gs @ Q); Uall = Q @ Rm      # U' G0s U = I, U' Gs U = diag(lam)
        kv = lam >= kappa
        self.reliab = lam; self.n_neg = int((lam < 0).sum()); self.n_drop = int((~kv).sum())
        self.min_eig_ratio = float(lam.min())
        self.U, self.lam, self.Ud = Uall[:, kv], lam[kv], Uall[:, ~kv]; self.r_g = int(kv.sum())
        cr = (self.U.T @ (bb / col[:, None])) / self.lam[:, None]  # identified coordinates (r, n)
        coef = self.U @ cr
        self.Xi_hat = np.zeros((E.shape[0], n)); self.Xi_hat[k] = coef / col[:, None]
        Sm = np.zeros((n * pk, n * pk))
        for tr in trajs:
            X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S)
            Th = theta(tr['Xs'], E)[:, k] / col; J = _jacobian_rows(tr['Xs'], E, self.Xi_hat)
            mi, m = D.shape; W_ = np.zeros((n * pk, n * m))
            for j in range(n):
                Bj = np.zeros((mi, n * m))
                for kk in range(n):
                    blk = -J[:, j, kk][:, None] * H
                    if kk == j:
                        blk = blk + D
                    Bj[:, kk * m:(kk + 1) * m] = blk * sig[kk]
                W_[j * pk:(j + 1) * pk] = Th.T @ Bj
            Sm += W_ @ W_.T
        P = np.kron(np.eye(n), self.U.T / self.lam[:, None])
        Vr = P @ Sm @ P.T; Vr = (Vr + Vr.T) / 2
        w, Qv = np.linalg.eigh(Vr); w = np.maximum(w, w.max() * 1e-14)
        self.Vr, self.Lr, self.Vr_inv = Vr, Qv * np.sqrt(w), (Qv / w) @ Qv.T
        self.rank = n * self.r_g; self.eps = float(chi2.ppf(1 - alpha, self.rank)); self.null_tol = null_tol
        self.dX = np.zeros((1, n))

    def _coords(self, a_scaled):
        """a (m, n*pk) scaled -> identified coefficients (m, n*r) and the relative weight on unidentified ones."""
        n, pk = self.n, self.p
        A = a_scaled.reshape(len(a_scaled), n, pk)
        cr = np.einsum('mjp,pr->mjr', A, self.U).reshape(len(a_scaled), -1)
        Nb = np.hstack([self.Ud, self.W0null])
        nn = np.linalg.norm(np.einsum('mjp,pr->mjr', A, Nb).reshape(len(a_scaled), -1), axis=1) if Nb.shape[1] else np.zeros(len(a_scaled))
        return cr, nn / np.maximum(np.linalg.norm(cr, axis=1) + nn, 1e-300)

    def vec(self, Xi):
        return (np.asarray(Xi)[self.keep] * self.col[:, None]).T.reshape(-1)      # scaled, equation-major

    def stat(self, Xi):
        d = (self.vec(self.Xi_hat) - self.vec(Xi)).reshape(self.n, self.p)
        cr = (d @ self.G0s @ self.U).reshape(-1)                         # c = U' G0s xi (U' G0s U = I)
        return float(cr @ self.Vr_inv @ cr)

    def law(self, c, X, tol):
        c = np.asarray([float(ci) for ci in c]); ph = theta(np.atleast_2d(X), self.E)
        g = (ph @ self.Xi_hat) @ c
        a = np.concatenate([cj * ph[:, self.keep] / self.col for cj in c], axis=1)
        cr, nn = self._coords(a)
        k = np.einsum('mi,ij,mj->m', cr, self.Vr, cr)
        k[nn > self.null_tol] = np.inf
        return g + tol, np.sqrt(self.eps * np.maximum(k, 0))

    def sample_boundary(self, k, rng):
        Z = rng.normal(size=(k, self.rank)); Z *= np.sqrt(self.eps) / np.linalg.norm(Z, axis=1)[:, None]
        dr = (Z @ self.Lr.T).reshape(k, self.n, self.r_g)
        d = np.einsum('kjr,pr->kjp', dr, self.U) / self.col          # (k, n, pk) original units
        out = np.repeat(self.Xi_hat[None], k, 0)
        out[:, self.keep, :] += d.transpose(0, 2, 1)
        return out

    def quad_form(self):
        """A with K = a' A a in the original coefficient basis (equation-major), identified part only."""
        T = np.kron(np.eye(self.n), self.U)                           # (n*pk, n*r)
        A = T @ self.Vr @ T.T
        cc = np.tile(self.col, self.n)
        return A / cc[:, None] / cc[None, :]
