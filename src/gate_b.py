"""Protocol 0B (docs/protocol_0B.md). Pools, sample draw, arms and needle jobs. No pandas.

Units: 0B-1 (system H), 0B-2 (needles), 0B-3 (Gate A over 20 seeds). The arms reuse src/gate_a.py and src/gate_a2.py unchanged.
"""
import os, sys, json, time, zlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
import numpy as np
from lib import theta
from systems import SYSTEMS, spec_G, in_domain
from data import make_cell, stack
from sindy import stlsq, deriv_rmse, BLOWUP
from scipy.integrate import solve_ivp
from lib import make_rhs

ROLL_WALL_S = 20.0      # per trajectory; a stiff model that has not escaped by then counts as a failed rollout


def _rollout(Xi, E, t, x0):
    """sindy.rollout with a wall-clock guard (same solver and settings)."""
    f0 = make_rhs(Xi, E)
    t_start = time.time()
    def f(tt, x):
        if time.time() - t_start > ROLL_WALL_S:
            raise TimeoutError
        return f0(tt, x)
    ev = lambda tt, x: BLOWUP - np.max(np.abs(x))
    ev.terminal = True
    try:
        sol = solve_ivp(f, (t[0], t[-1]), x0, t_eval=t, rtol=1e-6, atol=1e-9, method='RK45', events=ev, max_step=(t[-1] - t[0]) / 50)
    except Exception:
        return None
    if (not sol.success) or sol.y.shape[1] != len(t) or not np.all(np.isfinite(sol.y)):
        return None
    return sol.y.T


def rollout_rmse(Xi, E, trajs, state_sd, penalty=10.0):
    """sindy.rollout_rmse with the guarded rollout."""
    errs = []
    for tr in trajs:
        sl = tr['sl']
        Y = _rollout(Xi, E, tr['t'][sl], tr['Xs'][0])
        errs.append(penalty if Y is None else min(np.sqrt((((Y - tr['X'][sl]) / state_sd) ** 2).mean()), penalty))
    return float(np.mean(errs))
from verify import verify_spec_robust, exact_min_1d
import gate_a
from gate_a import NOISE, make_variants, library, structure, run_arm, evaluate

gate_a.STRUCT.setdefault('H', {'H1x', 'H1y'})      # H uses the A-type structure (x | x', y | y')
THRESHOLDS = np.logspace(-3, np.log10(0.3), 12)    # = pool.THRESHOLDS
SAMPLE_SEED = 20261122
N_PER_CELL = 4
SEEDS = tuple(range(20))
DEGREES = {'A': (2, 3, 4), 'Amis': (3, 4), 'B': (2, 3), 'H': (2, 3, 4)}


def sysname(key):
    return {'A': 'A', 'Amis': 'A', 'B': 'B', 'H': 'H'}[key]


def cells():
    """Fixed cell order. H at noise 1 and 2 (0B-1), A, Amis, B at noise 1 (0B-3). H noise 1 is shared with 0B-3."""
    out = [('H', ni, s) for ni in (1, 2) for s in SEEDS]
    out += [(k, 1, s) for k in ('A', 'Amis', 'B') for s in SEEDS]
    return out


def cell_id(key, ni, seed):
    return f'{key}_n{ni}_s{seed}'


# ------------------------------------------------------------------------------------------------ pool
def pool_rows(key, cell, seed):
    """The pilot pool loop (src/pool.py) without pandas and without the rollout/AIC columns."""
    S = SYSTEMS[sysname(key)]
    train = cell['train']
    variants = make_variants(len(train), seed)
    rows, seen = [], set()
    for d in DEGREES[key]:
        E = library(S, key, d)
        Efull = S.E(d)
        keep = np.array([tuple(e) in {tuple(x) for x in E.tolist()} for e in Efull.tolist()])
        for vname, idx in variants:
            Xs, dX = stack(train, idx)
            Th = theta(Xs, Efull)
            for th in THRESHOLDS:
                Xi = np.zeros((Efull.shape[0], dX.shape[1]))
                Xi[keep] = stlsq(Th[:, keep], dX, th) if not keep.all() else stlsq(Th, dX, th)
                kk = (d, Xi.tobytes())
                if kk in seen:
                    continue
                seen.add(kk)
                rm, _ = deriv_rmse(Xi, Efull, cell['iid'], cell['deriv_sd'])
                rows.append(dict(cand=len(rows), degree=d, variant=vname, threshold=float(th), k=int((Xi != 0).sum()),
                                 iid_rmse=rm, xi=Xi.tolist()))
    best = min(r['iid_rmse'] for r in rows)
    for r in rows:
        r['aes'] = bool(r['iid_rmse'] <= 1.10 * best)
        r['best_iid_rmse'] = best
    return rows


def law_volume(sp, Xi, E):
    """Share of the spec domain with G < -tol_strong. 1-D: exact, from the real roots of G + tol on the interval
    (G is a polynomial in u of degree <= deg E). 2-D/3-D: 10^6 uniform points, fixed seed."""
    tol = sp.tol_strong
    if len(sp.u_names) == 1:
        lo, hi = sp.u_box[0]
        deg = int(E.sum(1).max())
        z = lo + (hi - lo) * (1 - np.cos(np.pi * (np.arange(deg + 3) + 0.5) / (deg + 3))) / 2
        P = np.polyfit(z - lo, spec_G(sp, Xi, [z], E) + tol, deg)
        r = np.roots(P) if np.any(P) else np.array([])
        r = np.sort(np.concatenate([[0.0, hi - lo], [x.real for x in r if abs(x.imag) < 1e-12 and 0 < x.real < hi - lo]]))
        mid = (r[:-1] + r[1:]) / 2
        neg = np.polyval(P, mid) < 0
        return float(np.sum(np.diff(r)[neg]) / (hi - lo))
    from oracle import _sample_domain
    from needles import g_vec
    U = _sample_domain(sp, 10 ** 6, np.random.default_rng([SAMPLE_SEED, 61]))
    return float(np.mean(g_vec(sp, Xi, E, U) < -tol))


def pool_job(key, ni, seed):
    S = SYSTEMS[sysname(key)]
    t0 = time.time()
    cell = make_cell(S, NOISE[ni], seed, ni)
    rows = pool_rows(key, cell, seed)
    if key == 'H':                                   # 0B-1: verify every candidate, volume and depth of every violation
        for r in rows:
            Xi = np.array(r['xi']); E = S.E(r['degree'])
            st, tm, vol, depth = {}, {}, {}, {}
            for sp in S.specs:
                q = verify_spec_robust(sp, Xi, E, 'strong')
                st[sp.name], tm[sp.name] = q['status'], round(q['time'], 4)
                if q['status'] == 'counterexample':
                    vol[sp.name] = law_volume(sp, Xi, E)
                    depth[sp.name] = exact_min_1d(sp, Xi, E) / sp.tol_strong if len(sp.u_names) == 1 else q.get('G')
            r.update(status=st, z3_sec=tm, volume=vol, depth=depth, viol=any(v == 'counterexample' for v in st.values()),
                     timeout=any(v == 'timeout' for v in st.values()))
            if r['aes']:
                r['roll_ood'] = rollout_rmse(Xi, E, cell['ood'], cell['state_sd'])
    return dict(cell=cell_id(key, ni, seed), key=key, noise_idx=ni, seed=seed, n=len(rows), n_aes=sum(r['aes'] for r in rows),
                seconds=round(time.time() - t0, 1), rows=rows)


# ------------------------------------------------------------------------------------------------ sample
def draw_sample(pool_dir):
    """4 accuracy-equivalent candidates per cell, without replacement, one rng stream in the fixed cell order."""
    rng = np.random.default_rng(SAMPLE_SEED)
    out, short = [], []
    for key, ni, seed in cells():
        P = json.load(open(os.path.join(pool_dir, cell_id(key, ni, seed) + '.json')))
        aes = [r for r in P['rows'] if r['aes']]
        if len(aes) < N_PER_CELL:
            short.append(P['cell'])
        pick = sorted(rng.choice(len(aes), size=min(N_PER_CELL, len(aes)), replace=False).tolist())
        for i in pick:
            r = aes[i]
            out.append(dict(sample_id=f"{P['cell']}_c{r['cand']}", key=key, noise_idx=ni, seed=seed, cand=r['cand'],
                            degree=r['degree'], variant=r['variant'], threshold=r['threshold'], pool_iid=r['iid_rmse'], k_pool=r['k']))
    return out, short


# ------------------------------------------------------------------------------------------------ arms
def run_sample_H(job):
    """Arms S, C, Rnd on H. Same loop, budget, sampler and margin as 0A (gate_a.run_arm), plus out-of-range rollouts."""
    job = dict(job, seed=int(job['seed']), degree=int(job['degree']), threshold=float(job['threshold']), cand=int(job['cand']),
               pool_iid=float(job['pool_iid']), noise_idx=int(job['noise_idx']))
    S = SYSTEMS['H']
    seed, ni = job['seed'], job['noise_idx']
    cell = make_cell(S, NOISE[ni], seed, ni)
    E = library(S, 'H', job['degree'])
    idx = dict(make_variants(len(cell['train']), seed))[job['variant']]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    allowed, M, free_cols = structure(S, 'H', E)
    dXf = dX[:, free_cols]
    rng = np.random.default_rng([seed, job['cand'], ni, 5151])
    Xp = stlsq(Th, dX, job['threshold'])                   # the pool candidate itself (unstructured)
    ood0 = rollout_rmse(Xp, E, cell['ood'], cell['state_sd'])
    out = []
    t0 = time.time()
    for arm in ('S', 'C', 'Rnd'):
        Xi_free, iters, nrows = run_arm(arm, S, 'H', Th, dXf, allowed, job['threshold'], M, E, rng)
        rec = dict(job, arm=arm, iters=iters, n_rows=nrows, pool_ood=ood0)
        if Xi_free is None:
            rec.update(status='infeasible')
        else:
            Xi_full = Xi_free @ M.T
            rec.update(evaluate(S, 'H', Xi_full, E, cell))
            rec['iid_ratio'] = rec['iid_rmse'] / job['pool_iid']
            rec['roll_ood'] = rollout_rmse(Xi_full, E, cell['ood'], cell['state_sd'])
            rec['ood_ratio'] = rec['roll_ood'] / ood0
            rec['xi'] = json.dumps(Xi_full.tolist())
            rec['status'] = 'ok'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out


def arm_job(job, struct=''):
    if job['key'] == 'H':
        return run_sample_H(job)
    if job['key'] == 'B':
        import gate_a2
        return gate_a2.run_sample2(job, struct)
    return gate_a.run_sample(job)


# ------------------------------------------------------------------------------------------------ needles (0B-2)
def needle_bases(root):
    """True A (degree 4), true B (degree 3), and the first two arm-C models in sample order from 0A (A) and 0A2 (B-pos)
    that verify every spec of their system. A arm-C models are regenerated (0A did not store coefficients)."""
    import csv
    bases = []
    for k, d in (('A', 4), ('B', 3)):
        S = SYSTEMS[k]
        bases.append(dict(base=f'true_{k}', sys=k, degree=d, xi=S.true_Xi(d).tolist()))
    sample = list(csv.DictReader(open(os.path.join(root, 'results', '0A_sample.csv'))))
    got = 0
    for j in sample:
        if j['key'] != 'A' or got == 2:
            continue
        rec = json.load(open(os.path.join(root, 'results', '0A_jobs', j['sample_id'] + '.json')))
        c = [r for r in rec['records'] if r['arm'] == 'C' and r.get('status') == 'ok']
        if not c or any(v != 'verified' for v in json.loads(c[0]['statuses']).values()):
            continue
        Xi = regen_A_C(j)
        bases.append(dict(base=f"C_{j['sample_id']}", sys='A', degree=int(j['degree']), xi=Xi.tolist(), statuses_0A=c[0]['statuses']))
        got += 1
    got = 0
    for j in sample:
        if j['key'] != 'B' or got == 2:
            continue
        rec = json.load(open(os.path.join(root, 'results', '0A2_jobs', 'Bpos', j['sample_id'] + '.json')))
        c = [r for r in rec['records'] if r['arm'] == 'C' and r.get('status') == 'ok']
        if not c or any(v != 'verified' for v in json.loads(c[0]['statuses']).values()):
            continue
        bases.append(dict(base=f"Cpos_{j['sample_id']}", sys='B', degree=int(j['degree']), xi=json.loads(c[0]['xi'])))
        got += 1
    return bases


def regen_A_C(j):
    S = SYSTEMS['A']
    seed, d = int(j['seed']), int(j['degree'])
    cell = make_cell(S, NOISE[1], seed, 1)
    E = library(S, 'A', d)
    idx = dict(make_variants(len(cell['train']), seed))[j['variant']]
    Xs, dX = stack(cell['train'], idx)
    allowed, M, fc = structure(S, 'A', E)
    rng = np.random.default_rng([seed, int(j['cand']), 5151])
    Xf, _, _ = run_arm('C', S, 'A', theta(Xs, E), dX[:, fc], allowed, float(j['threshold']), M, E, rng)
    return Xf @ M.T


def needle_jobs(bases):
    import needles
    J = []
    for b in bases:
        S = SYSTEMS[b['sys']]
        for sp in S.specs:
            for ki, kap in enumerate(needles.KAPPAS):
                for vi, vt in enumerate(needles.VOLUMES):
                    J.append(dict(needle_id=f"{b['base']}__{sp.name}__k{ki}__v{vi}", base=b['base'], sys=b['sys'], degree=b['degree'],
                                  xi=b['xi'], spec=sp.name, kappa=kap, v_target=vt))
    return J


def needle_job(job):
    import needles
    S = SYSTEMS[job['sys']]
    sp = [q for q in S.specs if q.name == job['spec']][0]
    E = S.E(job['degree'])
    Xi = np.array(job['xi'])
    t0 = time.time()
    seed = zlib.crc32(job['needle_id'].encode())
    nd = needles.build_needle(sp, Xi, E, job['kappa'], job['v_target'], [SAMPLE_SEED, seed])
    out = {k: v for k, v in job.items() if k != 'xi'}
    if nd is None:
        out.update(status='unbuilt', seconds=round(time.time() - t0, 1))
        return out
    Xn = nd['Xi']
    r = verify_spec_robust(sp, Xn, E, 'strong', prefer='bernstein')
    hits = {str(n): needles.detect(sp, Xn, E, n, needles.N_SAMPLER_SEEDS, seed) for n in needles.BUDGETS}
    out.update(status='built', v=nd['v'], depth=nd['depth'], reached=nd['reached'], corner_k=nd['k'], power_d=nd['d'],
               a=nd['a'], t=nd['t'], verifier=r['status'], verifier_method=r.get('method'), verifier_sec=round(r['time'], 3),
               hits=hits, n_seeds=needles.N_SAMPLER_SEEDS, xi_needle=Xn.tolist(), seconds=round(time.time() - t0, 1))
    return out
