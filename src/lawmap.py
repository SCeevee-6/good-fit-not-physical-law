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
    def __init__(self, Xs, dX, deriv_sd, E, keep=None, rcond=1e-12, null_tol=1e-8, Thw=None, dXw=None, eps=None):
        """Rashomon set by default. With Thw, dXw (whitened rows from `whiten`), deriv_sd = noise sd per state and
        eps = chi2 quantile, this is the GLS confidence set of protocol 0F."""
        self.E = E
        self.keep = np.ones(E.shape[0], bool) if keep is None else np.asarray(keep, bool)
        self.sd = np.asarray(deriv_sd, float)
        if Thw is not None:
            Th, dX = Thw[:, self.keep], dXw
        else:
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
        self.eps = DELTA_RSS * self.r_hat if eps is None else float(eps)
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


# ---------------------------------------------------------------- protocol 0F additions

def sg_mats(m, S):
    """Savitzky-Golay smoothing (H) and derivative (D) matrices for a series of length m, interior rows only."""
    from data import savgol_filter
    I = np.eye(m); w, p, dt = S.sg_window, S.sg_poly, S.dt
    H = savgol_filter(I, w, p, axis=0, mode='interp'); D = savgol_filter(I, w, p, deriv=1, delta=dt, axis=0, mode='interp')
    h = w // 2; sl = slice(h, m - h)
    return H[sl], D[sl], sl


def noise_sd(trajs, S):
    """Pooled estimate of the measurement noise sd per state from smoothing residuals, corrected by the exact
    expectation E||(S_int - H_int) e||^2 = sigma^2 tr(R'R) for iid noise."""
    num, den = 0.0, 0.0
    for tr in trajs:
        X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S)
        R = np.eye(len(X))[sl] - H
        num = num + ((R @ X) ** 2).sum(0); den += float((R ** 2).sum())
    return np.sqrt(num / den)


def whiten(trajs, S, E):
    """Stacked whitened library rows and derivatives: L^-1 Theta(Xs), L^-1 dX per trajectory, L L' = D D'."""
    Ts, Ys, n = [], [], 0
    for tr in trajs:
        X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S)
        C = D @ D.T; C += 1e-12 * np.trace(C) / len(C) * np.eye(len(C))
        L = np.linalg.cholesky(C)
        Ts.append(np.linalg.solve(L, theta(tr['Xs'], E))); Ys.append(np.linalg.solve(L, tr['dX'])); n += len(C)
    return np.vstack(Ts), np.vstack(Ys), n


def conf_eps(ds, alpha=0.05):
    from scipy.stats import chi2
    return float(chi2.ppf(1 - alpha, ds.rank * ds.dX.shape[1]))


def sample_boundary(ds, k, rng):
    """k models on the boundary of the set: Xi_hat + delta with sum_j delta_j' G_j delta_j = eps."""
    n = ds.dX.shape[1]; r = ds.rank
    Z = rng.normal(size=(k, r * n)); Z *= np.sqrt(ds.eps) / np.linalg.norm(Z, axis=1)[:, None]
    out = np.repeat(ds.Xi_hat[None], k, 0)
    for j in range(n):
        coef = (Z[:, j * r:(j + 1) * r] / ds.s) @ ds.V.T * ds.sd[j]       # scaled basis
        out[:, ds.keep, j] += coef / ds.col
    return out


def sample_check(ds, c, X, tol, k, rng):
    """Counts of sampled models that contradict the zone map: violating at certified points, obeying at refuted."""
    c = np.asarray([float(ci) for ci in c]); m, h = ds.law(c, X, tol); z = zone(m, h)
    B = ds.sample_boundary(k, rng) if hasattr(ds, 'sample_boundary') else sample_boundary(ds, k, rng)
    Th = theta(X, ds.E)
    g = np.einsum('mp,kpn,n->km', Th, B, c)
    bad_c = int(((g < -tol - 1e-9) & (z == 1)[None]).sum()); bad_r = int(((g >= -tol + 1e-9) & (z == -1)[None]).sum())
    return bad_c, bad_r


def _poly_coeffs(ds, c, tol):
    """Monomial coefficients (degree 2d basis) of M = g_hat + tol and P = M^2 - eps w K."""
    from lib import exponents
    n = ds.E.shape[1]; d = int(ds.E.sum(1).max()); E2 = exponents(n, 2 * d)
    idx = {tuple(e): i for i, e in enumerate(E2.tolist())}
    c = np.asarray([float(ci) for ci in c]); w = float((c ** 2 * ds.sd ** 2).sum()) if hasattr(ds, 'sd') else 1.0
    M = np.zeros(len(E2))
    for a, e in enumerate(ds.E.tolist()):
        M[idx[tuple(e)]] += ds.Xi_hat[a] @ c
    M[idx[(0,) * n]] += tol
    Ek = ds.E[ds.keep]; P = np.zeros(len(E2)); pk = len(Ek)
    if hasattr(ds, 'quad_form'):                       # stacked set: K = sum_{j,l} c_j c_l phi' A_jl phi
        A0 = ds.quad_form(); A = np.zeros((pk, pk)); w = 1.0
        for j in range(len(c)):
            for l in range(len(c)):
                A += c[j] * c[l] * A0[j * pk:(j + 1) * pk, l * pk:(l + 1) * pk]
    else:
        A = (ds.V / ds.s ** 2) @ ds.V.T; A = A / ds.col[:, None] / ds.col[None, :]
    nzM = [(i, M[i]) for i in np.nonzero(M)[0]]
    for i1, v1 in nzM:
        for i2, v2 in nzM:
            P[idx[tuple(np.array(E2[i1]) + np.array(E2[i2]))]] += v1 * v2
    for a in range(len(Ek)):
        for b in range(len(Ek)):
            P[idx[tuple(Ek[a] + Ek[b])]] -= ds.eps * w * A[a, b]
    return E2, M, P


def _certify_poly(spec, E2, coef, xmax, budget_s):
    """Exact check that sum coef_a x^a >= 0 on the spec's domain, after rounding to the 1e-6 grid with a proven slack."""
    import dataclasses
    from fractions import Fraction as Fr
    from verify import verify_spec_robust
    s = 1e3 / max(np.abs(coef).max(), 1e-300); q = np.round(coef * s, 6)
    mono_max = np.prod(np.asarray(xmax, float)[None, :] ** E2, axis=1)
    b = float(np.ceil((np.abs(q - coef * s) * mono_max).sum() * 1e6 + 1) / 1e6)
    q[np.where(E2.sum(1) == 0)[0][0]] -= b
    Xi = np.zeros((len(E2), E2.shape[1])); Xi[:, 0] = np.round(q, 6) + 0.0
    sp = dataclasses.replace(spec, c=(Fr(1),) + (Fr(0),) * (E2.shape[1] - 1), tol_strong=0.0, tol_weak=0.0)
    if len(spec.u_names) <= 3:
        r = verify_spec_robust(sp, Xi, E2, 'strong', timeout_ms=5000, budget_s=budget_s, prefer='bernstein' if len(spec.u_names) == 3 else 'z3')
    else:
        import dimtest
        r = dimtest.verify_nd(sp, Xi, E2, 'strong')
    return r.get('status', 'timeout'), b


def region_certificate(ds, spec, X, budget_s=60.0):
    """'verified' if every model in the set obeys the law on the whole domain (exact), 'counterexample' if the
    polynomial check finds a point, 'timeout' otherwise. X: grid states, used only for coordinate bounds."""
    E2, M, P = _poly_coeffs(ds, spec.c, spec.tol_strong)
    xmax = np.abs(X).max(0) * 1.01 + 1e-9
    st1, _ = _certify_poly(spec, E2, M, xmax, budget_s)
    if st1 != 'verified':
        return st1
    st2, _ = _certify_poly(spec, E2, P, xmax, budget_s)
    return st2


def plane_spec(spec):
    """SIR law restricted to S + I + R = 1: face laws on the segment (u = t), region laws on the triangle (u = S, I)."""
    import dataclasses
    from fractions import Fraction as Fr
    one = Fr(1)
    if len(spec.u_names) == 2:
        pt = lambda u, K, f=spec.point: f([u[0], 1 - u[0]], K)
        return dataclasses.replace(spec, u_names=('t',), point=pt, domain=[((one,), Fr(0)), ((-one,), -one)], u_box=((0, 1),))
    pt = lambda u, K, f=spec.point: f([u[0], u[1], 1 - u[0] - u[1]], K)
    return dataclasses.replace(spec, u_names=('S', 'I'), point=pt,
                               domain=[((one, Fr(0)), Fr(0)), ((Fr(0), one), Fr(0)), ((-one, -one), -one)], u_box=((0, 1), (0, 1)))


# ---------------------------------------------------------------- errors-in-variables confidence set (0F)

def _jacobian_rows(Xs, E, Xi):
    """J[t, j, k] = d f_j / d x_k at Xs[t] for f = theta(., E) @ Xi."""
    m, n = Xs.shape; J = np.zeros((m, Xi.shape[1], n))
    for k in range(n):
        Ek = E.copy(); ck = Ek[:, k].astype(float); Ek[:, k] = np.maximum(Ek[:, k] - 1, 0)
        J[:, :, k] = (theta(Xs, Ek) * ck) @ Xi
    return J


class StackedSet:
    """GLS confidence set for all equations jointly. Residual model, linearised in the state noise e (iid, sd sigma_k):
    dX_j - f_j(Xs) = sum_k B_jk e_k,  B_jk = delta_jk D - diag(J_jk) H,  J from Xi_ref.
    Same interface as DataSet for law(), r(), Xi_hat, eps, rank."""

    def __init__(self, trajs, S, E, sigma, Xi_ref, keep=None, alpha=0.05, rcond=1e-12, null_tol=1e-8):
        from scipy.stats import chi2
        from scipy.linalg import solve_triangular
        self.E = E; n = E.shape[1]; self.n = n; self.Xi_ref = np.asarray(Xi_ref, float)
        self.keep = np.ones(E.shape[0], bool) if keep is None else np.asarray(keep, bool)
        pk = int(self.keep.sum()); self.p = pk; sig = np.asarray(sigma, float); self.sigma = sig
        Zs, ys = [], []
        for tr in trajs:
            X = np.asarray(tr['X'], float); H, D, sl = sg_mats(len(X), S)
            Th = theta(tr['Xs'], E)[:, self.keep]; J = _jacobian_rows(tr['Xs'], E, Xi_ref)
            mi, m = D.shape
            B = np.zeros((n * mi, n * m))
            for j in range(n):
                for k in range(n):
                    blk = -J[:, j, k][:, None] * H
                    if j == k:
                        blk = blk + D
                    B[j * mi:(j + 1) * mi, k * m:(k + 1) * m] = blk * sig[k]
            C = B @ B.T; C += 1e-12 * np.trace(C) / len(C) * np.eye(len(C))
            L = np.linalg.cholesky(C)
            Zb = np.zeros((n * mi, n * pk))
            for j in range(n):
                Zb[j * mi:(j + 1) * mi, j * pk:(j + 1) * pk] = Th
            Zs.append(solve_triangular(L, Zb, lower=True)); ys.append(solve_triangular(L, tr['dX'].T.reshape(-1), lower=True))
        Z = np.vstack(Zs); y = np.concatenate(ys)
        self.col = np.linalg.norm(Z, axis=0); self.col[self.col == 0] = 1.0
        U, sv, Vt = np.linalg.svd(Z / self.col, full_matrices=False)
        r = int((sv > rcond * sv[0]).sum()); self.rank, self.cond = r, float(sv[0] / sv[r - 1])
        self.s, self.V = sv[:r], Vt[:r].T
        self.Vnull = np.linalg.svd(Vt[:r], full_matrices=True)[2][r:].T
        self.null_tol = null_tol
        b = (self.V @ ((U[:, :r].T @ y) / self.s)) / self.col
        self.Xi_hat = np.zeros((E.shape[0], n)); self.Xi_hat[self.keep] = b.reshape(n, pk).T
        self.Z, self.y = Z, y
        self.r_hat = self.r(self.Xi_hat)
        self.eps = float(chi2.ppf(1 - alpha, r)); self.tau = self.r_hat + self.eps
        self.dX = np.zeros((1, n))                     # shape only (interface)

    def vec(self, Xi):
        return np.asarray(Xi)[self.keep].T.reshape(-1)

    def r(self, Xi):
        return float(((self.Z @ self.vec(Xi) - self.y) ** 2).sum())

    def law(self, c, X, tol):
        c = np.asarray([float(ci) for ci in c]); ph = theta(np.atleast_2d(X), self.E)
        g = (ph @ self.Xi_hat) @ c
        a = np.concatenate([cj * ph[:, self.keep] for cj in c], axis=1) / self.col
        z = (a @ self.V) / self.s; k = (z ** 2).sum(1)
        if self.Vnull.shape[1]:
            nn = np.linalg.norm(a @ self.Vnull, axis=1); k[nn > self.null_tol * np.linalg.norm(a, axis=1)] = np.inf
        return g + tol, np.sqrt(self.eps * k)

    def sample_boundary(self, k, rng):
        Z = rng.normal(size=(k, self.rank)); Z *= np.sqrt(self.eps) / np.linalg.norm(Z, axis=1)[:, None]
        d = ((Z / self.s) @ self.V.T) / self.col
        out = np.repeat(self.Xi_hat[None], k, 0)
        out[:, self.keep, :] += d.reshape(k, self.n, self.p).transpose(0, 2, 1)
        return out

    def quad_form(self):
        """A[(j,a),(l,b)] with K = a' A a in the original coefficient basis (equation-major)."""
        A = (self.V / self.s ** 2) @ self.V.T
        return A / self.col[:, None] / self.col[None, :]


def stacked_set(trajs, S, E, sigma, keep=None, alpha=0.05):
    """One feasible-GLS step: Jacobian from the least-squares fit, then the confidence set."""
    keep_ = np.ones(E.shape[0], bool) if keep is None else keep
    Xs = np.vstack([t['Xs'] for t in trajs]); dX = np.vstack([t['dX'] for t in trajs])
    Th = theta(Xs, E)[:, keep_]; Xi0 = np.zeros((E.shape[0], dX.shape[1]))
    Xi0[keep_] = np.linalg.lstsq(Th, dX, rcond=None)[0]
    return StackedSet(trajs, S, E, sigma, Xi0, keep=keep_, alpha=alpha)

