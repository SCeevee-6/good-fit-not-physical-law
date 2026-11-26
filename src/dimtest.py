"""Protocol 0D (docs/protocol_0D.md): does a dense constraint grid still enforce the laws at four states? Testbed E4 (SEIR).
Verifier for any number of states: z3 and a d-simplex Bernstein branch and bound (bern_nd.py), each witness re-checked exactly.
Arms S, C, Rnd, SoftA, SoftB on structure B-pos extended to four states. No pandas.
"""
import os, sys, json, time, zlib
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
from fractions import Fraction as Fr
import numpy as np
import data as data_mod, gate_a, gate_a2, gate_b, needles, oracle, robust, verify
from systems import System, Spec, SYSTEMS, K_FLOAT, spec_row, in_domain
from lib import theta
from data import make_cell, stack
from ceg import row_at
from sindy import deriv_rmse
from bern_nd import bernstein_nd
from real_r import limit, StageTimeout

NOISE_IDX = 1
SAMPLE_SEED = 20261128
SEEDS = tuple(range(20))
N_PER_CELL = 4
BETA, SIGMA, GAMMA = 0.6, 0.4, 0.2
ARMS = ('S', 'C', 'Rnd', 'SoftA', 'SoftB')
XCHECK_FRAC = 0.25                  # share of verifier calls on arm models that z3 also decides (agreement check)
M4 = np.array([[-1.0, -1.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])   # Xi_full = Xi_free @ M4.T, S' = -(E'+I'+R')
LAM_SOFT = 1.0
GRID_B_AXIS = 11
ARM_S = 900


def sys_E4():
    terms = [{(1, 0, 1, 0): -BETA}, {(1, 0, 1, 0): BETA, (0, 1, 0, 0): -SIGMA}, {(0, 1, 0, 0): SIGMA, (0, 0, 1, 0): -GAMMA}, {(0, 0, 1, 0): GAMMA}]
    def _ic(rng, lo, hi):
        N = rng.uniform(lo, hi)
        Ee, I0, R0 = rng.uniform(0.01, 0.1) * N, rng.uniform(0.01, 0.08) * N, rng.uniform(0.0, 0.1) * N
        return np.array([N - Ee - I0 - R0, Ee, I0, R0])
    tw, ts, cw, cs = 1e-3, 2e-2, 5e-3, 2e-2
    def dom(k):
        d = [(tuple(Fr(1) if j == i else Fr(0) for j in range(k)), Fr(0)) for i in range(k)]
        d.append((tuple(Fr(-1) for _ in range(k)), Fr(-3, 2)))
        return d
    box = lambda k: tuple((0, 1.5) for _ in range(k))
    def face(name, desc, c, zero, unames):
        def pt(u, K):
            out, it = [], iter(u)
            for q in range(4):
                out.append(0 * u[0] if q == zero else next(it))
            return out
        return Spec(name, desc, c, unames, pt, dom(3), tw, ts, box(3))
    ident = lambda u, K: [u[0], u[1], u[2], u[3]]
    specs = [
        face('E1S', "S'(0,E,I,R) >= 0", (Fr(1), Fr(0), Fr(0), Fr(0)), 0, ('E', 'I', 'R')),
        face('E1E', "E'(S,0,I,R) >= 0", (Fr(0), Fr(1), Fr(0), Fr(0)), 1, ('S', 'I', 'R')),
        face('E1I', "I'(S,E,0,R) >= 0", (Fr(0), Fr(0), Fr(1), Fr(0)), 2, ('S', 'E', 'R')),
        face('E1R', "R'(S,E,I,0) >= 0", (Fr(0), Fr(0), Fr(0), Fr(1)), 3, ('S', 'E', 'I')),
        Spec('E2S', "S' <= 0", (Fr(-1), Fr(0), Fr(0), Fr(0)), ('S', 'E', 'I', 'R'), ident, dom(4), tw, ts, box(4)),
        Spec('E2R', "R' >= 0", (Fr(0), Fr(0), Fr(0), Fr(1)), ('S', 'E', 'I', 'R'), ident, dom(4), tw, ts, box(4)),
        Spec('E3p', "S'+E'+I'+R' <= tol", (Fr(-1),) * 4, ('S', 'E', 'I', 'R'), ident, dom(4), cw, cs, box(4)),
        Spec('E3m', "S'+E'+I'+R' >= -tol", (Fr(1),) * 4, ('S', 'E', 'I', 'R'), ident, dom(4), cw, cs, box(4)),
    ]
    return System('E4_seir', ('S', 'E', 'I', 'R'), terms, (2, 3), 80.0, 0.25, 11, 3, lambda rng: _ic(rng, 0.7, 1.3), lambda rng: _ic(rng, 1.3, 1.5), specs)


SYSTEMS['E4'] = sys_E4()
data_mod.SYS_ID['E4_seir'] = 3
gate_b.DEGREES['E4'] = (2, 3)
_prev_sysname = gate_b.sysname
gate_b.sysname = lambda key: 'E4' if key == 'E4' else _prev_sysname(key)
STRUCT_LAWS = {'E1E', 'E1I', 'E1R', 'E3p', 'E3m'}
NONSTRUCT = ('E1S', 'E2S', 'E2R')


# ------------------------------------------------------------------------------------------------ verifier
def _point_map(spec):
    """For each state coordinate k: index m of the free variable it equals, or None if it is the constant 0."""
    d = len(spec.u_names)
    a = [np.array([0.11 * (j + 1)]) for j in range(d)]
    b = [np.array([0.07 * (j + 3)]) for j in range(d)]
    pa, pb = spec.point(a, K_FLOAT), spec.point(b, K_FLOAT)
    out = []
    for k in range(len(pa)):
        va, vb = float(np.asarray(pa[k]).ravel()[0]), float(np.asarray(pb[k]).ravel()[0])
        if va == 0.0 and vb == 0.0:
            out.append(None)
        else:
            m = [j for j in range(d) if np.isclose(va, a[j][0]) and np.isclose(vb, b[j][0])]
            assert len(m) == 1, spec.name
            out.append(m[0])
    return out


def law_poly(spec, Xi, E):
    d = len(spec.u_names)
    pm = _point_map(spec)
    g = np.asarray(Xi, float) @ np.array([float(c) for c in spec.c])
    poly = {}
    for j in range(E.shape[0]):
        if g[j] == 0.0:
            continue
        e = [0] * d; dead = False
        for k, m in enumerate(pm):
            if E[j, k] > 0:
                if m is None:
                    dead = True; break
                e[m] += int(E[j, k])
        if not dead:
            poly[tuple(e)] = poly.get(tuple(e), 0.0) + float(g[j])
    return poly


def verify_nd(spec, Xi, E, level='strong', order='bern', bern_s=60.0, z3_ms=3000):
    """order 'bern': Bernstein first, z3 on time-out. order 'z3': z3 first, Bernstein on time-out."""
    tol = spec.tol(level)
    d = len(spec.u_names)
    L = float(-spec.domain[-1][1])
    t0 = time.time()
    def check(u):
        uF = [Fr(float(x)) for x in u]
        return bool(verify._in_domain_exact(spec, uF) and verify._exact_G(spec, Xi, E, uF) < -Fr(tol).limit_denominator(10 ** 6))
    def run_bern():
        r = bernstein_nd(law_poly(spec, Xi, E), d, L, tol, check, time_budget=bern_s)
        return dict(spec=spec.name, status=r['status'], u=r['u'], time=time.time() - t0, method='bernstein_nd')
    def run_z3():
        r = verify.verify_spec(spec, Xi, E, level, z3_ms)
        return dict(spec=spec.name, status=r['status'], u=r.get('u'), time=time.time() - t0, method='z3')
    first, second = (run_bern, run_z3) if order == 'bern' else (run_z3, run_bern)
    r = first()
    return r if r['status'] != 'timeout' else second()


def evaluate_nd(S, Xi, E, cell, rng):
    """Statuses of all laws (Bernstein first); on a random XCHECK_FRAC of the laws z3 also runs (3 s) so agreement can be counted."""
    st, xc = {}, []
    for sp in S.specs:
        r = verify_nd(sp, Xi, E, 'strong', 'bern')
        st[sp.name] = r['status']
        if r['status'] in ('verified', 'counterexample') and rng.random() < XCHECK_FRAC:
            z = verify.verify_spec(sp, Xi, E, 'strong', 3000)['status']
            if z in ('verified', 'counterexample'):
                xc.append((sp.name, r['status'], z))
    ok_s = [st[s.name] == 'verified' for s in S.specs if s.name in STRUCT_LAWS]
    ok_n = [st[s.name] == 'verified' for s in S.specs if s.name not in STRUCT_LAWS]
    rm, _ = deriv_rmse(Xi, E, cell['iid'], cell['deriv_sd'])
    return dict(struct_adm=float(np.mean(ok_s)), nonstruct_adm=float(np.mean(ok_n)), iid_rmse=rm, k=int((Xi != 0).sum()),
                timeouts=sum(v == 'timeout' for v in st.values()), statuses=json.dumps(st), xcheck=xc,
                xcheck_disagree=sum(1 for _, b, z in xc if b != z))


# ------------------------------------------------------------------------------------------------ pool, sample, arms
CELLS = [('E4', NOISE_IDX, s) for s in SEEDS]


def cell_id(seed):
    return f'E4_n{NOISE_IDX}_s{seed}'


def pool_job(seed):
    t0 = time.time()
    S = SYSTEMS['E4']
    cell = make_cell(S, gate_a.NOISE[NOISE_IDX], seed, NOISE_IDX)
    rows = gate_b.pool_rows('E4', cell, seed)
    return dict(cell=cell_id(seed), seed=seed, n=len(rows), n_aes=sum(r['aes'] for r in rows), seconds=round(time.time() - t0, 1), rows=rows)


def draw_sample(pool_dir):
    rng = np.random.default_rng(SAMPLE_SEED)
    out, short = [], []
    for seed in SEEDS:
        if os.environ.get('OD_QUICK') and not os.path.exists(os.path.join(pool_dir, cell_id(seed) + '.json')):
            continue                          # pipeline test only (two pool cells)
        P = json.load(open(os.path.join(pool_dir, cell_id(seed) + '.json')))
        aes = [r for r in P['rows'] if r['aes']]
        if len(aes) < N_PER_CELL:
            short.append(P['cell'])
        for i in sorted(rng.choice(len(aes), size=min(N_PER_CELL, len(aes)), replace=False).tolist()):
            r = aes[i]
            out.append(dict(sample_id=f"{P['cell']}_c{r['cand']}", key='E4', noise_idx=NOISE_IDX, seed=seed, cand=r['cand'], degree=r['degree'],
                            variant=r['variant'], threshold=r['threshold'], pool_iid=r['iid_rmse'], k_pool=r['k']))
    return out, short


def grid_rows(S, E, per_axis):
    """Rows (p, 3) and rhs of a box grid filtered to the law domains, G >= -0.25 tol, mapped to the free coefficients.
    per_axis(dim) gives the points per axis for a law with dim free variables."""
    rows, rhs = [], []
    for sp in S.specs:
        dim = len(sp.u_names)
        m = per_axis(dim)
        mesh = np.meshgrid(*[np.linspace(lo, hi, m) for lo, hi in sp.u_box], indexing='ij')
        u = [g.ravel() for g in mesh]
        ok = in_domain(sp, u)
        u = [x[ok] for x in u]
        R = np.asarray(spec_row(sp, u, E)).reshape(len(u[0]), E.shape[0], -1)
        for r in R:
            rows.append(r @ M4); rhs.append(-0.25 * sp.tol_strong)
    return rows, rhs


GRID_A = lambda dim: {1: 100, 2: 11, 3: 6, 4: 4}[dim]        # the pilot's budget of about 100 points per law
GRID_B = lambda dim: GRID_B_AXIS                              # a fixed 11 points per axis in every law


def run_arm_nd(arm, S, Th, dXf, allowed, th, E, lb, rng, budget=15, samples=2000):
    """gate_a2.run_arm2 with the d-state verifier. arm C: verify_nd witnesses; arm Rnd: worst of 2000 random points."""
    rows, rhs = [], []
    Xi_free = gate_a2.fit_free2(Th, dXf, allowed, th, rows, rhs, lb, [])
    iters, stop = 0, None
    if arm == 'S' or Xi_free is None:
        return Xi_free, iters, 0, None
    for it in range(budget + 1):
        Xi_full = Xi_free @ M4.T
        viol = []
        for sp in S.specs:
            if arm == 'C':
                r = verify_nd(sp, Xi_full, E, 'strong', 'bern')
                if r['status'] == 'counterexample' and r.get('u') is not None:
                    viol.append((sp, r['u']))
            else:
                U = oracle._sample_domain(sp, samples, rng)
                G = needles.g_vec(sp, Xi_full, E, U)
                i = int(np.argmin(G))
                if G[i] < -sp.tol_strong:
                    viol.append((sp, [x[i] for x in U]))
        if not viol:
            stop = 'clean'; break
        if it == budget:
            stop = 'budget'; break
        for sp, u in viol:
            rr, b = row_at(sp, E, u)
            rows.append(rr @ M4); rhs.append(b)
        iters += 1
        Xi_free = gate_a2.fit_free2(Th, dXf, allowed, th, rows, rhs, lb, [])
        if Xi_free is None:
            return None, iters, len(rows), 'infeasible'
    return Xi_free, iters, len(rows), stop


def arm_job(job):
    seed = int(job['seed'])
    S = SYSTEMS['E4']
    cell = make_cell(S, gate_a.NOISE[NOISE_IDX], seed, NOISE_IDX)
    E = S.E(int(job['degree']))
    idx = dict(gate_a.make_variants(len(cell['train']), seed))[job['variant']]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    dXf = dX[:, [1, 2, 3]]
    allowed = np.ones((E.shape[0], 3), bool)
    lb = np.zeros((E.shape[0], 3), bool)
    for f in range(3):
        lb[:, f] = E[:, f + 1] == 0                                  # own-state-free monomials of the E, I, R equations have coefficient >= 0
    th = float(job['threshold'])
    rng = np.random.default_rng([seed, int(job['cand']), 5151])
    xrng = np.random.default_rng([SAMPLE_SEED, seed, int(job['cand']), 77])
    out, t0 = [], time.time()
    for arm in ARMS:
        rec = dict(job, arm=arm)
        try:
            with limit(ARM_S):
                if arm in ('S', 'C', 'Rnd'):
                    Xf, iters, nrows, stop = run_arm_nd(arm, S, Th, dXf, allowed, th, E, lb, rng)
                    rec.update(iters=iters, n_rows=nrows, stop=stop)
                else:
                    rows, rhs = grid_rows(S, E, GRID_A if arm == 'SoftA' else GRID_B)
                    Xf = robust.fit_soft(Th, dXf, allowed, th, rows, rhs, lb, [], LAM_SOFT)
                    rec.update(n_rows=len(rows), iters=0, stop=None)
                if Xf is None:
                    rec['status'] = 'infeasible'
                else:
                    Xi = Xf @ M4.T
                    rec.update(evaluate_nd(S, Xi, E, cell, xrng))
                    rec['iid_ratio'] = rec['iid_rmse'] / float(job['pool_iid'])
                    bad, mins = robust.sample_check(S, Xi, E, json.loads(rec['statuses']), seed)
                    rec.update(unsound=bad, sampled_min=mins, status='ok')
        except StageTimeout:
            rec['status'] = 'timeout'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out
