"""Protocol 0E (docs/protocol_0E.md): COVID-19 reported-case SIR panel (JHU CSSE), same pipeline as 0B/0R. No pandas.

Per config (country, SG window): series -> states S, I, R (population fractions) -> split 60/20/20 by days -> SG per segment ->
pool (pilot STLSQ sweep, library degrees 2-3, structure-free) -> verify every candidate against the seven SIR laws of system B
-> arms S, C, Rnd under B-pos and B-fac on selected candidates. One config = one JSON.
"""
import os, sys, csv, json, time, math, dataclasses
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
import numpy as np
import real_r                                   # time limits, progress log, select(); registers system R, harmless here
from real_r import limit, StageTimeout, progress, select, VERIFY_S, ARM_S, WIT_S
import systems, data as data_mod, gate_a, gate_a2, gate_b, needles
from systems import SYSTEMS, K_FLOAT
from lib import theta
from sindy import aic_score
from verify import verify_spec_robust, replay
import oracle

JHU = os.path.join(ROOT, 'data_raw', 'jhu')
END = '8/4/21'
FRAC_TRAIN, FRAC_VAL = 0.6, 0.2
N_ARM = 20                                   # equivalent candidates sent through the arms (plus the primary and the AIC pick)
VOL_N = 200000                               # uniform points for the violation volume of a law
PRIMARY_W = 21
PANEL = json.load(open(os.path.join(JHU, 'panel.json')))['panel']
CONFIGS = [dict(slug=p['slug'], window=w) for w in (21, 15, 29) for p in PANEL]     # primary window first, all countries
gate_b.DEGREES['E'] = (2, 3)
gate_b.sysname = lambda key: 'B' if key == 'E' else real_r._sysname(key)      # pool_rows uses the module-level name map


def name(cfg):
    return f"{cfg['slug']}_w{cfg['window']}"


def _series(path):
    rd = list(csv.reader(open(path, encoding='utf8')))
    return rd[0], rd[1:]


def load_country(slug):
    p = [q for q in PANEL if q['slug'] == slug][0]
    out = {}
    for k in ('confirmed', 'recovered', 'deaths'):
        h, rows = _series(os.path.join(JHU, f'time_series_covid19_{k}_global.csv'))
        r = [x for x in rows if x[1] == p['country'] and x[0] == ''][0]
        out[k] = np.array([float(v) for v in r[4:]]); dates = h[4:]
    i0, i1 = dates.index(p['start']), dates.index(p['end'])
    sl = slice(i0, i1 + 1)
    C, R_, D = out['confirmed'][sl], out['recovered'][sl], out['deaths'][sl]
    N = p['population']
    X = np.column_stack([1 - C / N, (C - R_ - D) / N, (R_ + D) / N])
    return np.arange(len(C), dtype=float), X, dates[sl]


def synthetic_B():
    """Off-grid test input only: system B truth on a daily grid, one wave, 1 % relative noise."""
    from scipy.integrate import solve_ivp
    f = lambda t, u: np.array([-0.06 * u[0] * u[1], 0.06 * u[0] * u[1] - 0.025 * u[1], 0.025 * u[1]])
    sol = solve_ivp(f, (0, 499), [0.995, 0.005, 0.0], t_eval=np.arange(500.0), rtol=1e-10, atol=1e-13)
    rng = np.random.default_rng(11)
    return sol.t, np.clip(sol.y.T * (1 + 0.01 * rng.normal(size=sol.y.T.shape)), 0, None), None


def build(cfg, series=None):
    S0 = SYSTEMS['B']
    S = dataclasses.replace(S0, sg_window=cfg['window'], sg_poly=3, dt=1.0)       # same specs; SG settings for daily data
    t, X, dates = synthetic_B() if series == 'synthetic_B' else load_country(cfg['slug'])
    n = len(t); a = int(FRAC_TRAIN * n); b = a + int(FRAC_VAL * n)
    seg = lambda sl: dict(t=t[sl], X=X[sl], **dict(zip(('Xs', 'dX', 'sl'), data_mod.smooth(X[sl], S))))
    T, V, O = seg(slice(0, a)), seg(slice(a, b)), seg(slice(b, n))
    w = cfg['window']; m = T['Xs'].shape[0]; nb = max(1, m // w)
    cuts = [i * w for i in range(nb)] + [m]
    train = [dict(Xs=T['Xs'][cuts[i]:cuts[i + 1]], dX=T['dX'][cuts[i]:cuts[i + 1]]) for i in range(nb)]
    cell = dict(train=train, iid=[V], ood=[O], state_sd=T['Xs'].std(0), deriv_sd=T['dX'].std(0), seed=0,
                train_flat=[dict(Xs=T['Xs'], dX=T['dX'])], train_whole=T)   # train_whole: raw X too (protocol 0F); not used by 0E
    info = dict(cfg=cfg, n_days=n, n_train=a, n_val=b - a, n_ext=n - b, n_train_smoothed=int(m), n_blocks=nb,
                obs_box=[[float(X[:, j].min()), float(X[:, j].max())] for j in range(3)],
                start=None if dates is None else dates[0], end=None if dates is None else dates[-1],
                I_max=float(X[:, 1].max()), S_min=float(X[:, 0].min()))
    return S, cell, info


def vol3(sp, Xi, E):
    U = oracle._sample_domain(sp, VOL_N, np.random.default_rng([gate_b.SAMPLE_SEED, 61]))
    return float(np.mean(needles.g_vec(sp, Xi, E, U) < -sp.tol_strong))


def witness_check(sp, Xi, E, q, info):
    st = real_r._witness_state(sp, q['u'])
    box = info['obs_box']
    out = dict(spec=sp.name, state=st, G=None if q.get('G') is None else float(q['G']),
               outside_obs_box=bool(any(not (box[j][0] <= st[j] <= box[j][1]) for j in range(3))), replay_ok=bool(replay(sp, Xi, E, q['u'])[0]))
    Y = gate_b._rollout(Xi, E, np.linspace(0.0, 30.0, 31), np.array(st))
    if Y is None:
        out['roll'] = 'failed_or_blowup'
    else:
        out.update(roll='ok', roll_min=float(Y.min()), roll_sum_range=[float(Y.sum(1).min()), float(Y.sum(1).max())])
    return out


def pool_stage(S, cell, info):
    t0 = time.time()
    rows = gate_b.pool_rows('E', cell, 0)
    if any(not np.isfinite(r['iid_rmse']) for r in rows):
        for r in rows:
            r['iid_rmse'] = float(r['iid_rmse']) if np.isfinite(r['iid_rmse']) else float('inf')
        best = min(r['iid_rmse'] for r in rows)
        for r in rows:
            r['best_iid_rmse'] = best; r['aes'] = bool(r['iid_rmse'] <= 1.10 * best)
    progress(f'pool fitted: {len(rows)} candidates, {sum(r["aes"] for r in rows)} equivalent')
    E_of = {d: S.E(d) for d in (2, 3)}
    for n_, r in enumerate(rows):
        Xi = np.array(r['xi']); E = E_of[r['degree']]
        st, vol, depth, tm = {}, {}, {}, {}
        r['limit_hits'] = []
        for sp in S.specs:
            try:
                with limit(VERIFY_S):
                    q = verify_spec_robust(sp, Xi, E, 'strong', prefer='bernstein')
            except StageTimeout:
                st[sp.name], tm[sp.name] = 'timeout', float(VERIFY_S); r['limit_hits'].append(f'verify:{sp.name}'); continue
            st[sp.name], tm[sp.name] = q['status'], round(q['time'], 3)
            if q['status'] != 'counterexample':
                continue
            try:
                with limit(VERIFY_S):
                    vol[sp.name] = vol3(sp, Xi, E)
                    depth[sp.name] = None if q.get('G') is None else float(q['G']) / sp.tol_strong
            except (StageTimeout, Exception) as e:
                vol.setdefault(sp.name, None); depth.setdefault(sp.name, None); r['limit_hits'].append(f'volume:{sp.name}:{type(e).__name__}')
            if r['aes'] and n_ < 10 ** 9:
                try:
                    with limit(WIT_S):
                        r.setdefault('wit', {})[sp.name] = witness_check(sp, Xi, E, q, info)
                except (StageTimeout, Exception) as e:
                    r.setdefault('wit', {})[sp.name] = dict(spec=sp.name, error=type(e).__name__); r['limit_hits'].append(f'witness:{sp.name}:{type(e).__name__}')
        try:
            with limit(VERIFY_S):
                aic = aic_score(Xi, E, cell['train_flat'], cell['deriv_sd'])
        except (StageTimeout, Exception):
            aic = float('inf')
        r.update(status=st, z3_sec=tm, volume=vol, depth=depth, viol=any(v == 'counterexample' for v in st.values()),
                 timeout=any(v == 'timeout' for v in st.values()), aic=aic)
        if n_ % 10 == 0:
            progress(f'pool verified {n_ + 1}/{len(rows)}')
    progress(f'pool done, limit hits {sum(len(r["limit_hits"]) for r in rows)}')
    return rows, round(time.time() - t0, 1)


def arms_stage(S, cell, info, rows, picks):
    byc = {r['cand']: r for r in rows}
    out = []
    for n_, cand in enumerate(picks):
        r = byc[cand]
        E = gate_a.library(S, 'B', r['degree'])
        idx = dict(gate_a.make_variants(len(cell['train']), 0))[r['variant']]
        Xs = np.vstack([cell['train'][i]['Xs'] for i in idx]); dX = np.vstack([cell['train'][i]['dX'] for i in idx])
        Th = theta(Xs, E)
        ood0 = r.get('roll_ood')
        if ood0 is None:
            ood0 = gate_b.rollout_rmse(np.array(r['xi']), E, cell['ood'], cell['state_sd'])
        for struct in ('Bpos', 'Bfac'):
            allowed, M, fc, lb, eq = gate_a2.structure2(S, struct, E)
            rng = np.random.default_rng([cand, 5151])
            for arm in ('S', 'C', 'Rnd'):
                rec = dict(cand=cand, struct=struct, arm=arm, pool_ood=ood0, pool_iid=r['iid_rmse'])
                try:
                    with limit(ARM_S):
                        Xf, iters, nrows, stop = gate_a2.run_arm2(arm, S, Th, dX[:, fc], allowed, r['threshold'], M, E, lb, eq, rng)
                        rec.update(iters=iters, n_rows=nrows, stop=stop)
                        if Xf is None:
                            rec['status'] = 'infeasible'
                        else:
                            Xi = Xf @ M.T
                            rec.update(gate_a.evaluate(S, 'B', Xi, E, cell))
                            rec['iid_ratio'] = rec['iid_rmse'] / r['iid_rmse']
                            rec['roll_ood'] = gate_b.rollout_rmse(Xi, E, cell['ood'], cell['state_sd'])
                            rec['ood_ratio'] = rec['roll_ood'] / ood0
                            rec['xi'] = Xi.tolist(); rec['status'] = 'ok'
                except StageTimeout:
                    rec['status'] = 'timeout'
                out.append(rec)
        progress(f'arms {n_ + 1}/{len(picks)} (cand {cand}) ' + ','.join(a['status'] for a in out[-6:]))
    return out


def run_config(cfg, series=None, nm=None):
    t0 = time.time()
    if nm:
        os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
        real_r._PROG[0] = os.path.join(ROOT, 'logs', f'0E_progress_{nm}.log')
    progress(f'start {cfg}')
    S, cell, info = build(cfg, series)
    rows, tp = pool_stage(S, cell, info)
    prim, aic, aes = select(rows)
    E_of = {d: S.E(d) for d in (2, 3)}
    picks = list(dict.fromkeys([prim, aic] + aes[:N_ARM]))
    for c in picks:
        r = rows[c]
        if 'roll_ood' not in r:
            r['roll_ood'] = gate_b.rollout_rmse(np.array(r['xi']), E_of[r['degree']], cell['ood'], cell['state_sd'])
    progress(f'selected: primary {prim}, aic {aic}, {len(aes)} AES, {len(picks)} to arms')
    arms = arms_stage(S, cell, info, rows, picks)
    return dict(cfg=cfg, info=info, n=len(rows), n_aes=len(aes), primary=prim, aic_pick=aic, aes=aes, arm_picks=picks,
                pool_seconds=tp, seconds=round(time.time() - t0, 1), rows=rows, arms=arms)
