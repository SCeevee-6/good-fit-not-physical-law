"""Law zones over the set of models that fit the training data about as well as the best one.

Set: E = {Xi : r(Xi) <= (1 + delta) r(Xi_hat)}, r = sum over equations of ||Theta Xi_j - dX_j||^2 / sd_j^2,
Xi_hat the least-squares fit on the (possibly reduced) library. Coefficients outside `keep` are fixed at 0.
Law at a state x: g(Xi, x) = sum_j c_j f_j(x) >= -tol. Over E:
    min, max g = g(Xi_hat, x) -/+ h(x),  h = sqrt(eps K(x)),  K = sum_j c_j^2 sd_j^2 phi(x)' (Theta'Theta)^-1 phi(x).
If phi(x) has a component in the null space of Theta (restricted to keep), K is infinite.
"""
import numpy as np
from lib import theta

DELTA_RMSE = 0.10                       # held-out RMSE within 10 % of best in every earlier gate
DELTA_RSS = (1 + DELTA_RMSE) ** 2 - 1   # same ratio on the residual sum of squares


class DataSet:
    def __init__(self, Xs, dX, deriv_sd, E, keep=None, rcond=1e-12, null_tol=1e-8):
        self.E = E
        self.keep = np.ones(E.shape[0], bool) if keep is None else np.asarray(keep, bool)
        self.sd = np.asarray(deriv_sd, float)
        Th = theta(Xs, E)[:, self.keep]
        self.col = np.linalg.norm(Th, axis=0); self.col[self.col == 0] = 1.0
        U, s, Vt = np.linalg.svd(Th / self.col, full_matrices=False)
        r = int((s > rcond * s[0]).sum())
        self.rank, self.p, self.cond = r, int(self.keep.sum()), float(s[0] / s[r - 1])
        self.s, self.V = s[:r], Vt[:r].T
        self.Vnull = np.linalg.svd(Vt[:r], full_matrices=True)[2][r:].T    # (p, p - r), works when rows < terms
        self.null_tol = null_tol
        coef = self.V @ ((U[:, :r].T @ dX) / self.s[:, None])          # scaled basis
        self.Xi_hat = np.zeros((E.shape[0], dX.shape[1]))
        self.Xi_hat[self.keep] = coef / self.col[:, None]
        self.Th, self.dX = Th, dX
        self.r_hat = self.r(self.Xi_hat)
        self.eps = DELTA_RSS * self.r_hat
        self.tau = self.r_hat + self.eps

    def r(self, Xi):
        res = self.Th @ np.asarray(Xi)[self.keep] - self.dX
        return float(((res / self.sd) ** 2).sum())

    def K(self, X):
        """phi' (Theta'Theta)^-1 phi at states X (m, n), in original units; inf where not identified."""
        ph = theta(np.atleast_2d(X), self.E)[:, self.keep] / self.col
        z = (ph @ self.V) / self.s
        k = (z ** 2).sum(1)
        if self.Vnull.shape[1]:
            nn = np.linalg.norm(ph @ self.Vnull, axis=1)
            k[nn > self.null_tol * np.linalg.norm(ph, axis=1)] = np.inf
        return k

    def law(self, c, X, tol):
        """m = g(Xi_hat, x) + tol and h at states X for law weights c."""
        c = np.asarray([float(ci) for ci in c])
        g = (theta(np.atleast_2d(X), self.E) @ self.Xi_hat) @ c
        w = float((c ** 2 * self.sd ** 2).sum())
        h = np.sqrt(self.eps * w * self.K(X))
        return g + tol, h


def zone(m, h):
    """1 certified, -1 refuted, 0 undecided."""
    z = np.zeros(m.shape, int)
    z[m >= h] = 1
    z[m + h < 0] = -1
    return z


def member_ratio(ds, Xi, c, X):
    """max over X of |g(Xi) - g(Xi_hat)| / h for one model; <= 1 for every model inside E."""
    c = np.asarray([float(ci) for ci in c])
    d = (theta(np.atleast_2d(X), ds.E) @ (np.asarray(Xi) - ds.Xi_hat)) @ c
    _, h = ds.law(c, X, 0.0)
    with np.errstate(divide='ignore', invalid='ignore'):
        q = np.where(np.isinf(h), 0.0, np.abs(d) / h)
    return float(np.nanmax(q))


def law_grid(spec, n1=2001, n2=201, n3=41, n4=15):
    """States on the law's domain: uniform grid in the query variables u, filtered by the domain."""
    from systems import K_FLOAT, in_domain
    k = len(spec.u_box); nper = {1: n1, 2: n2, 3: n3, 4: n4}[k]
    axes = [np.linspace(a, b, nper) for a, b in spec.u_box]
    U = [g.ravel() for g in np.meshgrid(*axes, indexing='ij')]
    ok = in_domain(spec, U)
    U = [u[ok] for u in U]
    pt = spec.point(U, K_FLOAT)
    m = U[0].shape[0]
    return np.stack([np.broadcast_to(np.asarray(p, float), (m,)) for p in pt], 1)
