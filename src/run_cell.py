"""Full pipeline for one grid cell: pool -> verification (+cross-checks) -> oracles -> repair arms -> audit. Resumable (skips finished cells)."""
import os, sys, time, json, pickle
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[v] = '1'
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import pandas as pd

NOISE = (0.005, 0.02, 0.05)


def cell_id(key, noise_idx, seed):
    return f"{key}_n{noise_idx}_s{seed}"


def _log(msg, path):
    with open(path, 'a', encoding='utf-8') as f:
        f.write(time.strftime('%H:%M:%S ') + msg + '\n')


def run_cell(args):
    key, noise_idx, seed, outdir, logpath = args
    cid = cell_id(key, noise_idx, seed)
    outp = os.path.join(outdir, cid + '.pkl')
    if os.path.exists(outp):
        return cid, 'skipped'
    t_start = time.time()
    from systems import SYSTEMS
    import data, pool, oracle, ceg
    from lib import theta
    from verify import verify_model, verify_spec_robust, replay, exact_min_1d
    from bern import bernstein_verify
    sysname = 'A' if key in ('A', 'Amis') else 'B'
    S = SYSTEMS[sysname]
    drop = [(2, 0)] if key == 'Amis' else None
    degrees = (3, 4) if key == 'Amis' else None
    cell = data.make_cell(S, NOISE[noise_idx], seed, noise_idx)
    df = pool.build_pool(S, cell, seed, drop=drop, degrees=degrees).reset_index(drop=True)
    df.insert(0, 'cell', cid); df.insert(1, 'key', key); df.insert(2, 'noise_idx', noise_idx); df.insert(3, 'seed', seed); df.insert(4, 'cand', df.index)
    _log(f'{cid}: pool {len(df)} cands ({time.time()-t_start:.0f}s)', logpath)
    rng = np.random.default_rng([seed, noise_idx, 31337])

    # ---------------- verification (z3-first, 3 s limit, Bernstein fallback) + cross-checks ----------------
    spec_rows, cand_extra = [], []
    for i, row in df.iterrows():
        Xi = np.array(row['xi']).reshape(row['xi_shape']); E = S.E(int(row['degree']))
        t0 = time.time()
        res = verify_model(S, Xi, E)
        tv = time.time() - t0
        n_bern_chk = n_bern_dis = n_1d_chk = n_1d_dis = 0
        do_1d = (sysname == 'A') and (rng.random() < 0.25)
        for (spn, lv), r in res.items():
            sp_ = [s for s in S.specs if s.name == spn][0]
            rec = dict(cell=cid, cand=i, spec=spn, level=lv, status=r['status'], method=r.get('method', 'z3'), time=r.get('time'), G=r.get('G'),
                       witness_exact=r.get('witness_exact'), copied=(lv == 'weak' and res[(spn, 'strong')]['status'] == 'counterexample'))
            if lv == 'strong' and r['status'] == 'counterexample' and r.get('u') is not None:
                rec['replay_ok'] = bool(replay(sp_, Xi, E, r['u'])[0])
            if len(sp_.u_names) == 3 and r['status'] in ('verified', 'counterexample') and r.get('method') == 'z3' and not rec['copied']:
                br = bernstein_verify(sp_, Xi, E, lv, time_budget=20.0)
                if br['status'] != 'timeout':
                    n_bern_chk += 1; n_bern_dis += int(br['status'] != r['status'])
            if do_1d and len(sp_.u_names) == 1 and not rec['copied']:
                mn = exact_min_1d(sp_, Xi, E)
                n_1d_chk += 1; n_1d_dis += int((mn < -sp_.tol(lv)) != (r['status'] == 'counterexample'))
            spec_rows.append(rec)
        st = [res[(s.name, 'strong')]['status'] for s in S.specs]
        wk = [res[(s.name, 'weak')]['status'] for s in S.specs]
        cand_extra.append(dict(smt_viol_strong=any(x == 'counterexample' for x in st), smt_timeout_strong=any(x == 'timeout' for x in st) and not any(x == 'counterexample' for x in st),
                               smt_viol_weak=any(x == 'counterexample' for x in wk), smt_timeout_weak=any(x == 'timeout' for x in wk) and not any(x == 'counterexample' for x in wk),
                               viol_specs_strong=','.join(s.name for s, x in zip(S.specs, st) if x == 'counterexample'), verify_seconds=tv,
                               bern_checks=n_bern_chk, bern_disagree=n_bern_dis, exact1d_checks=n_1d_chk, exact1d_disagree=n_1d_dis))
    df = pd.concat([df, pd.DataFrame(cand_extra)], axis=1)
    _log(f'{cid}: verified ({time.time()-t_start:.0f}s) strong violators {int(df.smt_viol_strong.sum())}/{len(df)}', logpath)

    # ---------------- oracles: all AES + 25 % of non-AES ----------------
    run_or = df.aes.values | (rng.random(len(df)) < 0.25)
    orows = []
    for i, row in df[run_or].iterrows():
        Xi = np.array(row['xi']).reshape(row['xi_shape']); E = S.E(int(row['degree']))
        lo = oracle.local_oracle(S, Xi, E, seed=seed)
        tj = oracle.trajectory_oracle(S, Xi, E, seed=seed)
        r2 = np.random.default_rng([seed, i, 99])
        f200, s200 = oracle.random_falsify(S, Xi, E, 200, r2); f2000, s2000 = oracle.random_falsify(S, Xi, E, 2000, r2)
        orows.append(dict(cand=i, oracle_run=True, local_viol_strong=any(v['viol_strong'] for v in lo.values()), local_viol_weak=any(v['viol_weak'] for v in lo.values()),
                          local_min_json=json.dumps({k: v['min'] for k, v in lo.items()}), traj_unphysical=tj['unphysical'], traj_n_unphysical=tj['n_unphysical'],
                          traj_flags_json=json.dumps({k: tj[k] for k in ('blowup', 'negative', 'trap', 's_rise', 'r_fall', 'drift')}),
                          rf200=f200, rf200_score=s200, rf2000=f2000, rf2000_score=s2000))
    odf = pd.DataFrame(orows)
    df = df.merge(odf, on='cand', how='left'); df['oracle_run'] = df['oracle_run'].fillna(False).astype(bool)
    _log(f'{cid}: oracles on {int(run_or.sum())} ({time.time()-t_start:.0f}s)', logpath)

    # ---------------- repair arms on AES violators ----------------
    rep_rows = []
    viol = df[df.aes & df.smt_viol_strong]
    for i, row in viol.iterrows():
        rd = row.to_dict()
        try:
            arms = ceg.repair_candidate(S, cell, rd, seed, drop)
        except Exception as e:
            _log(f'{cid}: repair error cand {i}: {e}', logpath); continue
        for a in arms:
            a = dict(a); a.update(cell=cid, cand=i, orig_iid_rmse=row.iid_rmse, orig_k=row.k, best_iid_rmse=row.best_iid_rmse, orig_support_f1=row.support_f1)
            rep_rows.append(a)
    rdf = pd.DataFrame(rep_rows)
    _log(f'{cid}: repairs {len(rdf)} recs ({time.time()-t_start:.0f}s)', logpath)

    # ---------------- audit of CEG-repaired models (z3-first strong+weak, SMT-free local oracle) ----------------
    audit = []
    if len(rdf):
        sub = rdf[(rdf.arm == 'ceg_witness') & (rdf.status == 'repaired')].head(30)
        for _, a in sub.iterrows():
            Xi = np.array(a['xi']).reshape(df.loc[a['cand'], 'xi_shape']); E = S.E(int(df.loc[a['cand'], 'degree']))
            res = verify_model(S, Xi, E)
            lo = oracle.local_oracle(S, Xi, E, seed=seed)
            audit.append(dict(cell=cid, cand=a['cand'], audit_strong_all_verified=all(res[(s.name, 'strong')]['status'] == 'verified' for s in S.specs),
                              audit_weak_all_verified=all(res[(s.name, 'weak')]['status'] == 'verified' for s in S.specs),
                              audit_timeouts=sum(res[(s.name, 'strong')]['status'] == 'timeout' for s in S.specs),
                              audit_local_oracle_viol=any(v['viol_strong'] for v in lo.values())))
    adf = pd.DataFrame(audit)
    out = dict(cell=cid, cand=df, spec=pd.DataFrame(spec_rows), repair=rdf, audit=adf, seconds=time.time() - t_start)
    tmp = outp + '.tmp'
    with open(tmp, 'wb') as f: pickle.dump(out, f)
    os.replace(tmp, outp)
    _log(f'{cid}: DONE in {time.time()-t_start:.0f}s', logpath)
    return cid, 'done'
