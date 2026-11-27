"""Protocol 0R (docs/protocol_0R.md): hare-lynx pelt series, same pipeline as 0B. No pandas.

Series -> normalise by the training mean -> Savitzky-Golay per segment -> pool (pilot STLSQ sweep) -> verify every candidate
-> arms S, C, Rnd on the selected candidates -> witness checks. One config = (SG window, shell level c), one JSON.
"""
import os, sys, csv, json, time, math, signal, contextlib
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[v] = '1'
from fractions import Fraction as Fr
import numpy as np
import systems
from systems import System, Spec, SYSTEMS, K_FLOAT
import data as data_mod
import gate_a, gate_b
from lib import theta
from sindy import stlsq, deriv_rmse, aic_score
from verify import verify_spec_robust, replay, exact_min_1d

SERIES = os.path.join(ROOT, 'data_raw', 'hare_lynx_pelts.csv')
TRAIN_END, VAL_END = 1895, 1915          # train 1845-1895, validation 1896-1915, extrapolation window 1916-1935
BLOCK = 5                                # block-bootstrap block length (points after smoothing)
N_ARM = 60                               # cap on accuracy-equivalent candidates sent through the arms
TW, TS = 1e-3, 2e-2                      # tolerances, as the pilot (normalised units)
VERIFY_S, ARM_S, WIT_S = 60, 180, 60      # per-candidate limits (addendum 1 to protocol 0R); no SIGALRM on Windows, so no limit there
PRIMARY = dict(window=7, c=3)
CONFIGS = [dict(window=w, c=c) for w in (7, 5, 9) for c in (3, 2, 4)]     # first one is the primary config
gate_b.DEGREES['R'] = (2, 3, 4)


def _sysname(key):
    return {'A': 'A', 'Amis': 'A', 'B': 'B', 'H': 'H', 'R': 'R'}[key]


gate_b.sysname = _sysname
gate_a.STRUCT['R'] = {'R1x', 'R1y'}


class StageTimeout(Exception):
    pass


@contextlib.contextmanager
def limit(seconds):
    """Wall-clock limit for a block of Python code (Linux). Native calls end when they return, z3 has its own 3 s limit."""
    if not hasattr(signal, 'SIGALRM'):
        yield; return
    def _h(sig, frm):
        raise StageTimeout()
    old = signal.signal(signal.SIGALRM, _h)
    signal.setitimer(signal.ITIMER_REAL, seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, old)


def selftest_limit(sec=2.0):
    """True if a busy Python loop is stopped by limit() within sec + 3 s. Always True where there is no SIGALRM."""
    if not hasattr(signal, 'SIGALRM'):
        return True
    t = time.time()
    try:
        with limit(sec):
            while True:
                sum(i * i for i in range(1000))
    except StageTimeout:
        return time.time() - t < sec + 3


_PROG = [None]


def progress(msg):
    if _PROG[0]:
        with open(_PROG[0], 'a') as f:
            f.write(f"{time.strftime('%H:%M:%S')} {msg}\n")


def load_series(path=SERIES):
    if path == 'synthetic_H':                      # off-grid test input only: annual samples of system H, 2 % noise
        from scipy.integrate import solve_ivp
        H = SYSTEMS['H']
        sol = solve_ivp(H.f_true(), (0, 90), [1.5, 0.8], t_eval=np.arange(91.0), rtol=1e-9, atol=1e-12)
        rng = np.random.default_rng(7)
        Z = sol.y.T * (1 + 0.02 * rng.normal(size=sol.y.T.shape))
        return np.arange(1845, 1936), Z[:, 0] * 53000, Z[:, 1] * 30000
    rows = list(csv.DictReader(open(path)))
    return (np.array([int(r['year']) for r in rows]), np.array([float(r['hare']) for r in rows]), np.array([float(r['lynx']) for r in rows]))


def _traj(years, X, S):
    Xs, dX, sl = data_mod.smooth(X, S)
    return dict(t=years.astype(float), X=X, Xs=Xs, dX=dX, sl=sl)


def build(cfg, series=SERIES):
    """Returns (S, cell, info). S is registered as SYSTEMS['R']."""
    yr, h, l = load_series(series)
    tr = yr <= TRAIN_END; va = (yr > TRAIN_END) & (yr <= VAL_END); te = yr > VAL_END
    sh, sl_ = h[tr].mean(), l[tr].mean()
    X = np.column_stack([h / sh, l / sl_])
    xm = float(np.ceil(X[tr, 0].max() * 20) / 20); ym = float(np.ceil(X[tr, 1].max() * 20) / 20)   # round up to 0.05
    xF, yF, C = Fr(xm).limit_denominator(100), Fr(ym).limit_denominator(100), Fr(cfg['c'])
    ymax, xmax = C * yF, C * xF
    specs = [
        Spec('R1x', 'hare cannot go negative: x\'(0,y) >= 0', (Fr(1), Fr(0)), ('y',), lambda u, K: [0 * u[0], u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -ymax)], TW, TS, ((0, float(ymax)),)),
        Spec('R1y', 'lynx cannot go negative: y\'(x,0) >= 0', (Fr(0), Fr(1)), ('x',), lambda u, K: [u[0], 0 * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -xmax)], TW, TS, ((0, float(xmax)),)),
        Spec('R2', f'trapping region: dV/dt <= 0 on V = x/{xF} + y/{yF} = {C}', (-1 / xF, -1 / yF), ('x',),
             lambda u, K: [u[0], K(C * yF) - K(yF / xF) * u[0]],
             [((Fr(1),), Fr(0)), ((Fr(-1),), -xmax)], TW, TS, ((0, float(xmax)),), shell_V=True),
    ]
    S = System('R_hare_lynx', ('x', 'y'), [{}, {}], (2, 3, 4), float(yr[-1] - yr[0]), 1.0, cfg['window'], 3,
               None, None, specs, n_train=0, n_iid=1, n_ood=1)
    SYSTEMS['R'] = S
    seg = lambda m: _traj(yr[m], X[m], S)
    T = seg(tr)
    n = T['Xs'].shape[0]; nb = max(1, n // BLOCK)
    cuts = [i * BLOCK for i in range(nb)] + [n]
    train = [dict(Xs=T['Xs'][cuts[i]:cuts[i + 1]], dX=T['dX'][cuts[i]:cuts[i + 1]]) for i in range(nb)]
    train[-1] = dict(Xs=T['Xs'][cuts[nb - 1]:], dX=T['dX'][cuts[nb - 1]:])
    cell = dict(train=train, iid=[seg(va)], ood=[seg(te)], state_sd=T['Xs'].std(0), deriv_sd=T['dX'].std(0), seed=0)
    cell['train_whole'] = T                     # unsplit training series with raw X (protocol 0F); not used by 0R
    vmax_obs = float((X[:, 0] / xm + X[:, 1] / ym).max())
    info = dict(cfg=cfg, scale_hare=float(sh), scale_lynx=float(sl_), xm=xm, ym=ym, n_train=int(tr.sum()), n_val=int(va.sum()), n_test=int(te.sum()),
                n_train_smoothed=int(n), n_blocks=nb, V_obs_max_all=vmax_obs,
                obs_box=[float(X[:, 0].min()), float(X[:, 0].max()), float(X[:, 1].min()), float(X[:, 1].max())])
    return S, cell, info


def _witness_state(sp, u):
    pt = sp.point([np.array([float(x)]) for x in u], K_FLOAT)
    return [float(np.asarray(p).ravel()[0]) for p in pt]


def witness_check(S, Xi, E, sp, q, info, T_roll=10.0):
    u = q['u']
    st = _witness_state(sp, u)
    box = info['obs_box']
    outside = not (box[0] <= st[0] <= box[1] and box[2] <= st[1] <= box[3])
    out = dict(spec=sp.name, state=st, G=None if q.get('G') is None else float(q['G']), outside_obs_box=bool(outside),
               replay_ok=bool(replay(sp, Xi, E, u)[0]))
    t = np.linspace(0.0, T_roll, 41)
    Y = gate_b._rollout(Xi, E, t, np.array(st))
    if Y is None:
        out.update(roll='failed_or_blowup')
    else:
        V = Y[:, 0] / info['xm'] + Y[:, 1] / info['ym']
        out.update(roll='ok', roll_min=float(Y.min()), roll_Vmax=float(V.max()))
    return out


def pool_stage(S, cell, info):
    t0 = time.time()
    rows = gate_b.pool_rows('R', cell, 0)
    if any(not np.isfinite(r['iid_rmse']) for r in rows):       # a NaN would poison min() and the equivalence rule
        for r in rows:
            r['iid_rmse'] = float(r['iid_rmse']) if np.isfinite(r['iid_rmse']) else float('inf')
        best = min(r['iid_rmse'] for r in rows)
        for r in rows:
            r['best_iid_rmse'] = best; r['aes'] = bool(r['iid_rmse'] <= 1.10 * best)
    progress(f'pool fitted: {len(rows)} candidates')
    E_of = {d: S.E(d) for d in (2, 3, 4)}
    for n_, r in enumerate(rows):
        Xi = np.array(r['xi']); E = E_of[r['degree']]
        st, vol, depth, tm = {}, {}, {}, {}
        r['limit_hits'] = []
        for sp in S.specs:
            try:
                with limit(VERIFY_S):
                    q = verify_spec_robust(sp, Xi, E, 'strong')
            except StageTimeout:
                st[sp.name], tm[sp.name] = 'timeout', float(VERIFY_S); r['limit_hits'].append(f'verify:{sp.name}')
                continue
            st[sp.name], tm[sp.name] = q['status'], round(q['time'], 4)
            if q['status'] != 'counterexample':
                continue
            try:
                with limit(VERIFY_S):
                    vol[sp.name] = gate_b.law_volume(sp, Xi, E)
                    depth[sp.name] = exact_min_1d(sp, Xi, E) / sp.tol_strong
            except (StageTimeout, Exception) as e:
                vol.setdefault(sp.name, None); depth.setdefault(sp.name, None); r['limit_hits'].append(f'volume:{sp.name}:{type(e).__name__}')
            try:
                with limit(WIT_S):
                    r.setdefault('wit', {})[sp.name] = witness_check(S, Xi, E, sp, q, info)
            except (StageTimeout, Exception) as e:
                r.setdefault('wit', {})[sp.name] = dict(spec=sp.name, error=type(e).__name__); r['limit_hits'].append(f'witness:{sp.name}:{type(e).__name__}')
        try:
            with limit(VERIFY_S):
                aic = aic_score(Xi, E, cell['train_flat'], cell['deriv_sd'])
        except (StageTimeout, Exception):
            aic = float('inf')
        r.update(status=st, z3_sec=tm, volume=vol, depth=depth, viol=any(v == 'counterexample' for v in st.values()),
                 timeout=any(v == 'timeout' for v in st.values()), aic=aic)
        if n_ % 20 == 0:
            progress(f'pool verified {n_ + 1}/{len(rows)}')
    progress(f'pool done, limit hits {sum(len(r["limit_hits"]) for r in rows)}')
    return rows, round(time.time() - t0, 1)


def select(rows):
    key = lambda r: (r['iid_rmse'] if np.isfinite(r['iid_rmse']) else float('inf'), r['k'], r['degree'], r['cand'])
    prim = min(rows, key=key)['cand']
    aic = min(rows, key=lambda r: (r['aic'] if np.isfinite(r['aic']) else float('inf'), r['k'], r['cand']))['cand']
    aes = sorted([r for r in rows if r['aes']], key=key)
    return prim, aic, [r['cand'] for r in aes]


def arms_stage(S, cell, info, rows, picks):
    byc = {r['cand']: r for r in rows}
    out = []
    for n_, cand in enumerate(picks):
        r = byc[cand]
        E = gate_a.library(S, 'R', r['degree'])
        idx = dict(gate_a.make_variants(len(cell['train']), 0))[r['variant']]
        Xs = np.vstack([cell['train'][i]['Xs'] for i in idx]); dX = np.vstack([cell['train'][i]['dX'] for i in idx])
        Th = theta(Xs, E)
        allowed, M, fc = gate_a.structure(S, 'R', E)
        rng = np.random.default_rng([cand, 5151])
        ood0 = r.get('roll_ood')
        if ood0 is None:
            ood0 = gate_b.rollout_rmse(np.array(r['xi']), E, cell['ood'], cell['state_sd'])
        for arm in ('S', 'C', 'Rnd'):
            rec = dict(cand=cand, arm=arm, pool_ood=ood0, pool_iid=r['iid_rmse'])
            try:
                with limit(ARM_S):
                    Xf, iters, nrows = gate_a.run_arm(arm, S, 'R', Th, dX[:, fc], allowed, r['threshold'], M, E, rng)
                    rec.update(iters=iters, n_rows=nrows)
                    if Xf is None:
                        rec['status'] = 'infeasible'
                    else:
                        Xi = Xf @ M.T
                        rec.update(gate_a.evaluate(S, 'R', Xi, E, cell))
                        rec['iid_ratio'] = rec['iid_rmse'] / r['iid_rmse']
                        rec['roll_ood'] = gate_b.rollout_rmse(Xi, E, cell['ood'], cell['state_sd'])
                        rec['ood_ratio'] = rec['roll_ood'] / ood0
                        rec['xi'] = Xi.tolist(); rec['status'] = 'ok'
            except StageTimeout:
                rec['status'] = 'timeout'
            out.append(rec)
        progress(f'arms {n_ + 1}/{len(picks)} (cand {cand}) ' + ','.join(a['status'] for a in out[-3:]))
    return out


def run_config(cfg, series=SERIES, name=None):
    t0 = time.time()
    if name:
        os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
        _PROG[0] = os.path.join(ROOT, 'logs', f'0R_progress_{name}.log')
    progress(f'start {cfg}')
    S, cell, info = build(cfg, series)
    cell['train_flat'] = [dict(Xs=np.vstack([t['Xs'] for t in cell['train']]), dX=np.vstack([t['dX'] for t in cell['train']]))]
    rows, tp = pool_stage(S, cell, info)
    prim, aic, aes = select(rows)
    E_of = {d: S.E(d) for d in (2, 3, 4)}
    aes_pick = aes[:N_ARM]
    for c in set([prim, aic] + aes_pick):          # extrapolation-window rollout of the unstructured candidate
        r = rows[c]
        if 'roll_ood' not in r:
            r['roll_ood'] = gate_b.rollout_rmse(np.array(r['xi']), E_of[r['degree']], cell['ood'], cell['state_sd'])
    picks = list(dict.fromkeys([prim, aic] + aes_pick))
    progress(f'selected: primary {prim}, aic {aic}, {len(aes)} AES, {len(picks)} to arms')
    arms = arms_stage(S, cell, info, rows, picks)
    return dict(cfg=cfg, info=info, n=len(rows), n_aes=len(aes), primary=prim, aic_pick=aic, aes=aes, arm_picks=picks,
                pool_seconds=tp, seconds=round(time.time() - t0, 1), rows=rows, arms=arms)
