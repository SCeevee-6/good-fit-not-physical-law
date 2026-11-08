"""Protocol 0A (docs/protocol_0A.md): structure only (S), witness loop (C), random-counterexample loop (Rnd).
Everything here avoids pandas (the ceg-ode kernel cannot import it in this sandbox)."""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np
import cvxpy as cp
from lib import theta, round_xi
from systems import SYSTEMS, spec_G
from data import make_cell, stack
from sindy import deriv_rmse, stlsq
from verify import verify_spec_robust
from oracle import _sample_domain
from ceg import row_at, MARGIN

NOISE = (0.005, 0.02, 0.05)
DROP_AMIS = [(2, 0)]
M_B = np.array([[-1.0, -1.0], [1.0, 0.0], [0.0, 1.0]])   # Xi_full = Xi_free @ M.T ; free = (I-eq, R-eq); S-eq = -(I'+R')
STRUCT = {'A': {'A1x', 'A1y'}, 'Amis': {'A1x', 'A1y'}, 'B': {'B1I', 'B1R', 'B3p', 'B3m'}}
N_BOOT = 8


def make_variants(ntr, seed, n_boot=N_BOOT):
    rng = np.random.default_rng([seed, 777])
    return [('full', list(range(ntr)))] + [(f'boot{b}', list(rng.integers(0, ntr, ntr))) for b in range(n_boot)]


def library(S, key, degree):
    E = S.E(degree)
    if key == 'Amis':
        keep = np.array([tuple(e) not in DROP_AMIS for e in E.tolist()])
        E = E[keep]
    return E


def structure(S, key, E):
    """Allowed (library term, free equation) pairs and the map Xi_full = Xi_free @ M.T."""
    if key == 'B':
        allowed = np.column_stack([E[:, 1] >= 1, E[:, 2] >= 1])
        return allowed, M_B, [1, 2]
    allowed = np.column_stack([E[:, 0] >= 1, E[:, 1] >= 1])
    return allowed, np.eye(2), [0, 1]


def fit_free(Th, dXf, allowed, th, rows, rhs, max_iter=10):
    """Constrained STLSQ on the free equations: rows are linear constraints on Xi_free (p, nf); allowed masks the library."""
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
            qy = Q.T @ dXf[:, f]
            obj = obj + cp.sum_squares(R @ z[pos:pos + len(A)] - qy) / m
            pos += len(A)
        obj = obj + 1e-10 * cp.sum_squares(z)
        cons = []
        if rows:
            Rm = np.array([[r[j, f] / scale[j] for (j, f) in idx] for r in rows])
            cons.append(Rm @ z >= np.array(rhs))
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
        if (new == act).all():
            break
        act = new
    return round_xi(Xi)


def spec_family(S, key):
    struct = STRUCT[key]
    return ([s for s in S.specs if s.name in struct], [s for s in S.specs if s.name not in struct])


def run_arm(arm, S, key, Th, dXf, allowed, th, M, E, rng, budget=15, samples=2000):
    rows, rhs = [], []
    Xi_free = fit_free(Th, dXf, allowed, th, rows, rhs)
    iters = 0
    if arm == 'S' or Xi_free is None:
        return Xi_free, iters, 0
    for it in range(budget + 1):
        Xi_full = Xi_free @ M.T
        viol = []
        for sp in S.specs:
            if arm == 'C':
                r = verify_spec_robust(sp, Xi_full, E, 'strong', prefer='bernstein')
                if r['status'] == 'counterexample' and r.get('u') is not None:
                    viol.append((sp, r['u']))
            else:
                U = _sample_domain(sp, samples, rng)
                G = spec_G(sp, Xi_full, U, E)
                i = int(np.argmin(G))
                if G[i] < -sp.tol_strong:
                    viol.append((sp, [x[i] for x in U]))
        if not viol or it == budget:
            break
        for sp, u in viol:
            rr, b = row_at(sp, E, u)
            rows.append(rr @ M)
            rhs.append(b)
        iters += 1
        Xi_free = fit_free(Th, dXf, allowed, th, rows, rhs)
        if Xi_free is None:
            return None, iters, len(rows)
    return Xi_free, iters, len(rows)


def evaluate(S, key, Xi_full, E, cell):
    out = {}
    for sp in S.specs:
        r = verify_spec_robust(sp, Xi_full, E, 'strong', prefer='bernstein')
        out[sp.name] = r['status']
    sfam, nfam = spec_family(S, key)
    ok_s = [out[s.name] == 'verified' for s in sfam]
    ok_n = [out[s.name] == 'verified' for s in nfam]
    rm, _ = deriv_rmse(Xi_full, E, cell['iid'], cell['deriv_sd'])
    return dict(struct_adm=float(np.mean(ok_s)), nonstruct_adm=float(np.mean(ok_n)), iid_rmse=rm, k=int((Xi_full != 0).sum()),
                timeouts=sum(v == 'timeout' for v in out.values()), statuses=json.dumps(out))


def run_sample(job):
    job = dict(job, seed=int(job['seed']), degree=int(job['degree']), threshold=float(job['threshold']), cand=int(job['cand']), pool_iid=float(job['pool_iid']))
    rec0 = dict(job)
    key, seed, degree, variant, th = job['key'], job['seed'], job['degree'], job['variant'], job['threshold']
    sysname = 'B' if key == 'B' else 'A'
    S = SYSTEMS[sysname]
    cell = make_cell(S, NOISE[1], seed, 1)
    E = library(S, key, degree)
    idx = dict(make_variants(len(cell['train']), seed))[variant]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    allowed, M, free_cols = structure(S, key, E)
    dXf = dX[:, free_cols]
    rng = np.random.default_rng([seed, job['cand'], 5151])
    out = []
    t0 = time.time()
    for arm in ('S', 'C', 'Rnd'):
        Xi_free, iters, nrows = run_arm(arm, S, key, Th, dXf, allowed, th, M, E, rng)
        rec = dict(rec0, arm=arm, iters=iters, n_rows=nrows)
        if Xi_free is None:
            rec.update(status='infeasible')
        else:
            Xi_full = Xi_free @ M.T
            rec.update(evaluate(S, key, Xi_full, E, cell))
            rec['iid_ratio'] = rec['iid_rmse'] / job['pool_iid']
            rec['status'] = 'ok'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out


def sanity_reproduce(job):
    """Re-create the pool candidate with the pool's own STLSQ (no mask) and compare to the pilot CSV row."""
    from sindy import stlsq as _s
    key, seed, degree, variant, th = job['key'], job['seed'], job['degree'], job['variant'], job['threshold']
    S = SYSTEMS['B' if key == 'B' else 'A']
    cell = make_cell(S, NOISE[1], seed, 1)
    E = library(S, key, degree)
    idx = dict(make_variants(len(cell['train']), seed))[variant]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    Xi = _s(Th, dX, th)
    return float(np.abs(Xi).sum())
