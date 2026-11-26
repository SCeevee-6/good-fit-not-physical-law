"""Protocol 0C (docs/protocol_0C.md): robustness checks. No pandas.

0C-1 soft-penalty baseline on the 0B sample (noise index 1). 0C-2 noise sensitivity of the Gate A comparison (noise indices 0 and 2).
0C-3 independent sampling check of every model the verifier accepts. Arms S, C, Rnd are gate_a / gate_a2 code, unchanged.
"""
import os, sys, json, time, zlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np
import cvxpy as cp
import gate_a, gate_a2, gate_b, needles, oracle
from lib import theta, round_xi
from systems import SYSTEMS, spec_row, in_domain
from data import make_cell, stack
from ceg import MARGIN
from real_r import limit, StageTimeout

NOISE = gate_a.NOISE
SAMPLE_SEED = 20261127
N_PER_CELL = 4
SEEDS = tuple(range(20))
KEYS = ('A', 'Amis', 'B', 'H')
LAMS = (1e-2, 1.0, 1e2)                      # soft-penalty weights, frozen
GRID_N = 100                                 # constraint points per law, as the pilot's dense-grid arm
CHECK_N = 200000                             # independent sampling check, points per law
JOB_S, ARM_S = 900, 300                      # wall-clock limits per job and per arm run (Linux)
gate_a.STRUCT.setdefault('H', {'H1x', 'H1y'})
sysname = {'A': 'A', 'Amis': 'A', 'B': 'B', 'H': 'H'}


def cells():
    if os.environ.get('OC_QUICK'):             # pipeline test only: one cell per system
        return [(k, 0, 0) for k in KEYS]
    return [(k, ni, s) for ni in (0, 2) for k in KEYS for s in SEEDS]


def cell_id(k, ni, s):
    return f'{k}_n{ni}_s{s}'


def pool_job(key, ni, seed):
    t0 = time.time()
    S = SYSTEMS[sysname[key]]
    cell = make_cell(S, NOISE[ni], seed, ni)
    rows = gate_b.pool_rows(key, cell, seed)
    return dict(cell=cell_id(key, ni, seed), key=key, noise_idx=ni, seed=seed, n=len(rows), n_aes=sum(r['aes'] for r in rows),
                seconds=round(time.time() - t0, 1), rows=rows)


def draw_sample(pool_dir):
    rng = np.random.default_rng(SAMPLE_SEED)
    out, short = [], []
    for key, ni, seed in cells():
        P = json.load(open(os.path.join(pool_dir, cell_id(key, ni, seed) + '.json')))
        aes = [r for r in P['rows'] if r['aes']]
        if len(aes) < N_PER_CELL:
            short.append(P['cell'])
        for i in sorted(rng.choice(len(aes), size=min(N_PER_CELL, len(aes)), replace=False).tolist()):
            r = aes[i]
            out.append(dict(sample_id=f"{P['cell']}_c{r['cand']}", key=key, noise_idx=ni, seed=seed, cand=r['cand'], degree=r['degree'],
                            variant=r['variant'], threshold=r['threshold'], pool_iid=r['iid_rmse'], k_pool=r['k']))
    return out, short


# ------------------------------------------------------------------------------------------------ soft penalty
def dense_rows(S, E, M):
    """Rows (p, nf) and rhs of the pilot's dense grid, G(u) >= -MARGIN * tol, mapped to the free coefficients."""
    rows, rhs = [], []
    for sp in S.specs:
        dim = len(sp.u_names)
        m = {1: GRID_N, 2: int(np.sqrt(GRID_N)) + 1, 3: int(round(GRID_N ** (1 / 3))) + 1}[dim]
        mesh = np.meshgrid(*[np.linspace(lo, hi, m) for lo, hi in sp.u_box], indexing='ij')
        u = [g.ravel() for g in mesh]
        ok = in_domain(sp, u)
        u = [x[ok] for x in u]
        R = np.asarray(spec_row(sp, u, E)).reshape(len(u[0]), E.shape[0], -1)
        for r in R:
            rows.append(r @ M); rhs.append(-MARGIN * sp.tol_strong)
    return rows, rhs


def fit_soft(Th, dXf, allowed, th, rows, rhs, lb, eq, lam, max_iter=10):
    """gate_a2.fit_free2 with the constraint rows as a penalty lam * sum(max(0, rhs - r.z)^2) instead of hard constraints."""
    p, nf = Th.shape[1], dXf.shape[1]
    m = Th.shape[0]
    act = allowed.copy()
    scale = np.linalg.norm(Th, axis=0) / np.sqrt(m)
    scale[scale == 0] = 1.0
    Ths = Th / scale
    Xi = np.zeros((p, nf))
    for _ in range(max_iter):
        idx = [(j, f) for f in range(nf) for j in range(p) if act[j, f]]
        if not idx:
            return np.zeros((p, nf))
        z = cp.Variable(len(idx))
        obj, pos = 0, 0
        for f in range(nf):
            A = [j for j in range(p) if act[j, f]]
            if not A:
                continue
            Q, R = np.linalg.qr(Ths[:, A])
            obj = obj + cp.sum_squares(R @ z[pos:pos + len(A)] - Q.T @ dXf[:, f]) / m
            pos += len(A)
        Rm = np.array([[r[j, f] / scale[j] for (j, f) in idx] for r in rows])
        obj = obj + lam * cp.sum_squares(cp.pos(np.array(rhs) - Rm @ z)) + 1e-10 * cp.sum_squares(z)
        cons = []
        li = [k for k, (j, f) in enumerate(idx) if lb[j, f]]
        if li:
            cons.append(z[li] >= 0)
        if eq:
            Em = np.array([[r[j, f] / scale[j] for (j, f) in idx] for r in eq])
            Em = Em[np.abs(Em).sum(1) > 0]
            if len(Em):
                cons.append(Em @ z == 0)
        prob = cp.Problem(cp.Minimize(obj), cons)
        try:
            prob.solve(solver=cp.CLARABEL)
        except Exception:
            return None
        if prob.status not in ('optimal', 'optimal_inaccurate') or z.value is None:
            return None
        Xi = np.zeros((p, nf))
        for k, (j, f) in enumerate(idx):
            Xi[j, f] = z.value[k] / scale[j]
        new = act & (np.abs(Xi) >= th)
        for r in eq:
            j = int(np.where(r[:, 0] != 0)[0][0])
            keep = new[j, 0] or new[j, 1]
            new[j, 0] = keep and allowed[j, 0]; new[j, 1] = keep and allowed[j, 1]
        if (new == act).all():
            break
        act = new
    Xi = gate_a2.project_eq(Xi, eq)
    Xi[lb & (Xi < 0)] = 0.0
    return round_xi(Xi)


# ------------------------------------------------------------------------------------------------ arms
def sample_check(S, Xi, E, statuses, seed):
    """Independent check: uniform sample of each law's domain. Returns (n verified laws with a sampled point below -tol_strong, {law: sampled min / tol})."""
    rng = np.random.default_rng([SAMPLE_SEED, seed, 313])
    bad, mins = 0, {}
    for sp in S.specs:
        U = oracle._sample_domain(sp, CHECK_N, rng)
        mn = float(needles.g_vec(sp, Xi, E, U).min()) / sp.tol_strong
        mins[sp.name] = round(mn, 4)
        if statuses.get(sp.name) == 'verified' and mn < -1.0:
            bad += 1
    return bad, mins


def _setup(job):
    key, seed, ni = job['key'], int(job['seed']), int(job['noise_idx'])
    S = SYSTEMS[sysname[key]]
    cell = make_cell(S, NOISE[ni], seed, ni)
    E = gate_a.library(S, key, int(job['degree']))
    idx = dict(gate_a.make_variants(len(cell['train']), seed))[job['variant']]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    if key == 'B':
        allowed, M, fc, lb, eq = gate_a2.structure2(S, 'Bpos', E)
    else:
        allowed, M, fc = gate_a.structure(S, key, E)
        lb, eq = np.zeros(allowed.shape, bool), []
    return S, cell, E, Th, dX[:, fc], allowed, M, lb, eq


def _finish(job, S, cell, E, M, Xf, rec):
    if Xf is None:
        rec['status'] = 'infeasible'; return rec
    Xi = Xf @ M.T
    rec.update(gate_a.evaluate(S, job['key'], Xi, E, cell))
    rec['iid_ratio'] = rec['iid_rmse'] / float(job['pool_iid'])
    bad, mins = sample_check(S, Xi, E, json.loads(rec['statuses']), int(job['seed']))
    rec.update(unsound=bad, sampled_min=mins, status='ok')
    return rec


def arm_job(job):
    """Arms S, C, Rnd for one sampled candidate (0C-2)."""
    S, cell, E, Th, dXf, allowed, M, lb, eq = _setup(job)
    key = job['key']; seed, ni, cand = int(job['seed']), int(job['noise_idx']), int(job['cand'])
    rng = np.random.default_rng([seed, cand, ni, 5151])
    th = float(job['threshold'])
    out = []
    t0 = time.time()
    for arm in ('S', 'C', 'Rnd'):
        rec = dict(job, arm=arm)
        try:
            with limit(ARM_S):
                if key == 'B':
                    Xf, iters, nrows, stop = gate_a2.run_arm2(arm, S, Th, dXf, allowed, th, M, E, lb, eq, rng)
                else:
                    Xf, iters, nrows = gate_a.run_arm(arm, S, key, Th, dXf, allowed, th, M, E, rng); stop = None
                rec.update(iters=iters, n_rows=nrows, stop=stop)
                _finish(job, S, cell, E, M, Xf, rec)
        except StageTimeout:
            rec['status'] = 'timeout'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out


def soft_job(job):
    """Soft-penalty arms Soft1..Soft3 (lam = 1e-2, 1, 1e2) for one candidate of the 0B sample (0C-1)."""
    S, cell, E, Th, dXf, allowed, M, lb, eq = _setup(job)
    th = float(job['threshold'])
    rows, rhs = dense_rows(S, E, M)
    out = []
    t0 = time.time()
    for i, lam in enumerate(LAMS):
        rec = dict(job, arm=f'Soft{i + 1}', lam=lam, n_rows=len(rows))
        try:
            with limit(ARM_S):
                Xf = fit_soft(Th, dXf, allowed, th, rows, rhs, lb, eq, lam)
                _finish(job, S, cell, E, M, Xf, rec)
        except StageTimeout:
            rec['status'] = 'timeout'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out
