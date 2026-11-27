"""0F-4: choosing new trajectories to settle a law (protocol 0F).

Arms: 'law' (greedy, minimise predicted undecided share of the target set), 'random', 'halton' (space-filling),
'dopt' (greedy log det). random / halton / dopt do not depend on the law, so one sequence per seed serves every law
of the system; 'law' runs once per law. Prediction for law and dopt: roll out the current fit from the candidate start,
whiten with the per-equation filter covariance sigma_j^2 D D' (the full errors-in-variables covariance is used for
every actual refit), and add the predicted rows to G.
"""
import time
import numpy as np

DESIGN = {'A': (4, ('A2',)), 'B': (3, ('B2S', 'B2R')), 'E4': (3, ('E2S', 'E2R'))}     # protocol 0F addendum 2 (H dropped)
SEEDS = tuple(range(200, 220))
N_CAND = 64
CAP = 30
TARGET = 1 / 20
MARGIN_TOL = 4.0          # target set: g_true(x) >= 4 tol_strong
GRID = dict(n1=501, n2=61, n3=25, n4=11)
ROLL_CAP_S = 3.0


def _rollout(Xi, E, t, x0, cap=ROLL_CAP_S):
    from scipy.integrate import solve_ivp
    from lib import make_rhs
    f0 = make_rhs(Xi, E); t0 = time.time()
    def f(tt, x):
        if time.time() - t0 > cap:
            raise TimeoutError
        return f0(tt, x)
    ev = lambda tt, x: 1e3 - np.max(np.abs(x)); ev.terminal = True
    try:
        sol = solve_ivp(f, (t[0], t[-1]), x0, t_eval=t, rtol=1e-6, atol=1e-9, method='RK45', events=ev, max_step=(t[-1] - t[0]) / 50)
    except Exception:
        return None
    if sol.y.shape[1] < 2 * 11 or not np.all(np.isfinite(sol.y)):
        return None
    return sol.y.T


def new_traj(S, ic, sigma, rng):
    from data import simulate_true, smooth
    t, Xc = simulate_true(S, ic)
    X = Xc + rng.normal(size=Xc.shape) * sigma
    Xs, dX, sl = smooth(X, S)
    return dict(X=X, Xs=Xs, dX=dX, sl=sl, ic=np.asarray(ic, float))


def target_points(S, names):
    from lawmap import law_grid
    out = {}
    f = S.f_true()
    for sp in S.specs:
        if sp.name not in names:
            continue
        X = law_grid(sp, **GRID)
        c = np.array([float(x) for x in sp.c])
        g = np.array([f(0.0, x) for x in X]) @ c
        out[sp.name] = (sp, X[g >= MARGIN_TOL * sp.tol_strong])
    return out


def cand_box(S):
    from lawmap import law_grid
    X = np.vstack([law_grid(sp, **GRID) for sp in S.specs])
    return np.maximum(X.min(0), 0.01), np.maximum(X.max(0), 0.02)


def undecided(ds, tp):
    from lawmap import zone
    out = {}
    for name, (sp, X) in tp.items():
        m, h = ds.law(sp.c, X, sp.tol_strong)
        out[name] = float((zone(m, h) == 0).mean()) if len(X) else 0.0
    return out


class Predictor:
    """Predicted effect of one more trajectory on the current set (fixed Xi_hat, eps)."""

    def __init__(self, ds, S, E, sigma, m_len):
        from lawmap import sg_mats
        self.ds, self.S, self.E, self.sigma = ds, S, E, np.asarray(sigma, float)
        H, D, sl = sg_mats(m_len, S); C = D @ D.T; C += 1e-12 * np.trace(C) / len(C) * np.eye(len(C))
        self.L = np.linalg.cholesky(C); self.sl = sl
        self.Ginv0 = None

    def rows(self, ic, t):
        from lib import theta
        from scipy.linalg import solve_triangular
        Xp = _rollout(self.ds.Xi_hat, self.E, t, np.asarray(ic, float))
        if Xp is None or len(Xp) != len(t):
            return None
        Th = theta(Xp[self.sl], self.E)[:, self.ds.keep]
        W = solve_triangular(self.L, Th, lower=True)
        n = self.ds.n; pk = W.shape[1]
        Z = np.zeros((n * W.shape[0], n * pk))
        for j in range(n):
            Z[j * W.shape[0]:(j + 1) * W.shape[0], j * pk:(j + 1) * pk] = W / self.sigma[j]
        return Z / self.ds.col                                      # scaled basis of the set

    def G(self):
        return (self.ds.V * self.ds.s ** 2) @ self.ds.V.T          # scaled-basis Gram on the identified subspace

    def undecided_after(self, Zn, tp):
        from lib import theta
        from lawmap import zone
        Gn = self.G() + Zn.T @ Zn; Gi = np.linalg.pinv(Gn)
        out = []
        for name, (sp, X) in tp.items():
            c = np.array([float(x) for x in sp.c]); ph = theta(X, self.E)
            g = (ph @ self.ds.Xi_hat) @ c + sp.tol_strong
            a = np.concatenate([cj * ph[:, self.ds.keep] for cj in c], axis=1) / self.ds.col
            h = np.sqrt(self.ds.eps * np.einsum('ij,jk,ik->i', a, Gi, a))
            out.append(float((zone(g, h) == 0).mean()) if len(X) else 0.0)
        return float(np.mean(out))

    def logdet_after(self, Zn):
        return float(np.linalg.slogdet(self.G() + Zn.T @ Zn + 1e-12 * np.eye(Zn.shape[1]))[1])


def run_sequence(S, E, base, sigma_noise, tp, laws, arm, cands, halton, rng, t_grid):
    """Add trajectories until every law in `laws` is below TARGET or CAP is reached. Returns per-law counts."""
    from lawmap import noise_sd, stacked_set
    trajs = list(base); need = {l: None for l in laws}; trace = []
    for k in range(CAP + 1):
        sig = noise_sd(trajs, S); ds = stacked_set(trajs, S, E, sig)
        u = undecided(ds, {l: tp[l] for l in laws}); trace.append(u)
        for l in laws:
            if need[l] is None and u[l] < TARGET:
                need[l] = k
        if all(v is not None for v in need.values()) or k == CAP:
            break
        if arm == 'random':
            ic = cands[rng.integers(len(cands))]
        elif arm == 'halton':
            ic = halton[k]
        else:
            P = Predictor(ds, S, E, sig, len(t_grid)); best, ic = None, None
            sub = {l: tp[l] for l in laws}
            for cnd in cands:
                Zn = P.rows(cnd, t_grid)
                if Zn is None:
                    continue
                v = P.undecided_after(Zn, sub) if arm == 'law' else -P.logdet_after(Zn)
                if best is None or v < best:
                    best, ic = v, cnd
            if ic is None:
                ic = cands[rng.integers(len(cands))]
        trajs.append(new_traj(S, ic, sigma_noise, rng))
    return {l: (CAP + 1 if need[l] is None else need[l]) for l in laws}, trace


def design_job(key, seed):
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from scipy.stats import qmc
    S = SYSTEMS[key]; d, laws = DESIGN[key]; E = S.E(d); t0 = time.time()
    cell = make_cell(S, NOISE[1], seed, 1); base = cell['train']; sigma_noise = np.asarray(cell['sigma'], float)
    tp = target_points(S, laws)
    lo, hi = cand_box(S)
    rng0 = np.random.default_rng([seed, 41]); cands = lo + (hi - lo) * rng0.uniform(size=(N_CAND, len(lo)))
    halton = qmc.scale(qmc.Halton(d=len(lo), scramble=True, seed=seed).random(CAP), lo, hi)
    t_grid = np.arange(0.0, S.T + S.dt / 2, S.dt)
    rec = dict(key=key, seed=seed, degree=d, laws=list(laws), n_target={l: int(len(tp[l][1])) for l in laws}, arms={})
    for arm in ('random', 'halton', 'dopt'):
        rng = np.random.default_rng([seed, 43, ['random', 'halton', 'dopt'].index(arm)])
        n, tr = run_sequence(S, E, base, sigma_noise, tp, laws, arm, cands, halton, rng, t_grid)
        rec['arms'][arm] = dict(need=n, trace=tr)
    for l in laws:
        rng = np.random.default_rng([seed, 47, laws.index(l)])
        n, tr = run_sequence(S, E, base, sigma_noise, tp, (l,), 'law', cands, halton, rng, t_grid)
        rec['arms'][f'law_{l}'] = dict(need=n, trace=tr)
    rec['seconds'] = round(time.time() - t0, 1)
    return rec
