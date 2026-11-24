"""Analysis for protocol 0R. Reads <jobs>/*.json, writes CSVs (no pandas). Gate thresholds are compared as exact fractions of counts.
usage (repo root): python src/analyze_0R.py [--jobs DIR] [--out DIR]
"""
import os, sys, csv, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
import numpy as np

ARGS = sys.argv[1:]
JOBS = ARGS[ARGS.index('--jobs') + 1] if '--jobs' in ARGS else os.path.join(ROOT, 'results', '0R_jobs')
OUT = ARGS[ARGS.index('--out') + 1] if '--out' in ARGS else os.path.join(ROOT, 'results')
LAWS = ('R1x', 'R1y', 'R2')


def wcsv(name, rows):
    rows = rows or [dict(empty=1)]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else float('nan')


def adm_all(r):
    return r.get('status') == 'ok' and all(v == 'verified' for v in json.loads(r['statuses']).values())


def summarize(J):
    rows = J['rows']; byc = {r['cand']: r for r in rows}
    aes = [byc[c] for c in J['aes']]
    nv = sum(r['viol'] for r in aes)
    arms = {(a['cand'], a['arm']): a for a in J['arms']}
    picks = J['arm_picks']
    aes_pick = [c for c in picks if byc[c]['aes']]
    out = dict(window=J['cfg']['window'], c=J['cfg']['c'], n_pool=J['n'], n_aes=J['n_aes'], n_aes_viol=nv,
               n_pool_viol=sum(r['viol'] for r in rows), n_timeout=sum(r['timeout'] for r in rows),
               primary=J['primary'], primary_viol=byc[J['primary']]['viol'], aic_pick=J['aic_pick'], aic_viol=byc[J['aic_pick']]['viol'],
               n_arm_aes=len(aes_pick), V_obs_max_all=J['info']['V_obs_max_all'], pool_seconds=J['pool_seconds'])
    for arm in ('S', 'C', 'Rnd'):
        A = [arms[(c, arm)] for c in aes_pick if (c, arm) in arms]
        out[f'f_{arm}'] = sum(adm_all(a) for a in A) / max(1, len(A))
        out[f'k_{arm}'] = sum(adm_all(a) for a in A)
        out[f'med_iid_{arm}'] = med([a.get('iid_ratio') for a in A if a['status'] == 'ok'])
        out[f'med_ood_{arm}'] = med([a.get('ood_ratio') for a in A if a['status'] == 'ok'])
        out[f'med_rows_{arm}'] = med([a.get('n_rows') for a in A if a['status'] == 'ok'])
    return out, aes_pick, arms, byc


def main():
    J = {}
    for f in glob.glob(os.path.join(JOBS, '*.json')):
        d = json.load(open(f))
        if d.get('error'):
            print('ERROR in', f, d['error'][-300:]); continue
        J[(d['cfg']['window'], d['cfg']['c'])] = d
    prim_key = (7, 3)
    summ, detail = [], {}
    for k in sorted(J):
        s, aes_pick, arms, byc = summarize(J[k]); summ.append(s); detail[k] = (aes_pick, arms, byc)
    wcsv('0R_summary.csv', summ)
    if prim_key not in J:
        print('primary config missing'); return
    P = J[prim_key]; aes_pick, arms, byc = detail[prim_key]
    S0 = [s for s in summ if (s['window'], s['c']) == prim_key][0]

    # candidate table for the primary config
    cand = []
    for c in P['arm_picks']:
        r = byc[c]
        row = dict(cand=c, role=('primary' if c == P['primary'] else '') + ('aic' if c == P['aic_pick'] else ''), aes=r['aes'], degree=r['degree'],
                   threshold=r['threshold'], k=r['k'], val_rmse=r['iid_rmse'], aic=r['aic'], test_rollout=r.get('roll_ood'), viol=r['viol'],
                   laws=','.join(l for l in LAWS if r['status'][l] == 'counterexample'),
                   **{f'vol_{l}': r['volume'].get(l, '') for l in LAWS}, **{f'depth_{l}': r['depth'].get(l, '') for l in LAWS})
        for arm in ('S', 'C', 'Rnd'):
            a = arms.get((c, arm))
            row[f'adm_{arm}'] = bool(a and adm_all(a))
            row[f'iid_ratio_{arm}'] = a.get('iid_ratio') if a else ''
            row[f'ood_ratio_{arm}'] = a.get('ood_ratio') if a else ''
            row[f'rows_{arm}'] = a.get('n_rows') if a else ''
        cand.append(row)
    wcsv('0R_candidates.csv', cand)

    # witnesses of the primary and AIC picks
    wit = []
    for c in dict.fromkeys([P['primary'], P['aic_pick']]):
        for l, w in (byc[c].get('wit') or {}).items():
            wit.append(dict(cand=c, role='primary' if c == P['primary'] else 'aic', **{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in w.items()}))
    wcsv('0R_witnesses.csv', wit)

    # gates, exact fractions
    gates = []
    pr = byc[P['primary']]
    gates.append(dict(rule='G-R1', stat='validation-selected model (min validation RMSE) breaks a law', value=pr['viol'],
                      laws=','.join(l for l in LAWS if pr['status'][l] == 'counterexample'), val_rmse=pr['iid_rmse'], k=pr['k'], degree=pr['degree'],
                      witness_outside_obs_box=','.join(f"{l}:{w['outside_obs_box']}" for l, w in (pr.get('wit') or {}).items()),
                      witness_replay_ok=','.join(f"{l}:{w['replay_ok']}" for l, w in (pr.get('wit') or {}).items()),
                      threshold='verifier counterexample at strong tolerance', passed=bool(pr['viol'])))
    n, k = S0['n_aes'], S0['n_aes_viol']
    gates.append(dict(rule='G-R2', stat='share of accuracy-equivalent candidates that break a law', count=k, n=n, value=k / max(1, n),
                      threshold='n >= 5 and share >= 1/5', passed=(n >= 5) and (5 * k >= n)))
    Cp = arms.get((P['primary'], 'C'))
    prim_ok = bool(Cp and adm_all(Cp) and Cp['iid_ratio'] <= 1.10 and Cp['ood_ratio'] <= 1.10)
    kc, na = S0['k_C'], S0['n_arm_aes']
    pooled_ok = (10 * kc >= 9 * na) and S0['med_iid_C'] <= 1.10 and S0['med_ood_C'] <= 1.10
    gates.append(dict(rule='G-R3', stat='C repairs the selected model and the pooled set at low cost',
                      primary_C_verified=bool(Cp and adm_all(Cp)), primary_iid_ratio=Cp.get('iid_ratio') if Cp else '', primary_ood_ratio=Cp.get('ood_ratio') if Cp else '',
                      pooled_count=kc, pooled_n=na, pooled_med_iid=S0['med_iid_C'], pooled_med_ood=S0['med_ood_C'],
                      threshold='primary C verified, ratios <= 1.10; pooled f(C) >= 9/10, medians <= 1.10', passed=bool(prim_ok and pooled_ok)))
    keep = gates[0]['passed'] and gates[1]['passed']
    gates.append(dict(rule='Verdict', stat='real-data claim retained (G-R1 and G-R2); cheap fix retained (also G-R3)',
                      value=('retained' if keep else ('inconclusive' if gates[0]['passed'] else 'dropped')), cheap_fix=bool(keep and gates[2]['passed']),
                      passed=bool(keep)))
    sens = [(s['window'], s['c'], s['primary_viol'], s['n_aes_viol'], s['n_aes']) for s in summ]
    gates.append(dict(rule='Sensitivity', stat='primary violates in every config / share of AES violating, per config',
                      value=all(x[2] for x in sens), detail=';'.join(f"w{w}c{c}:{int(v)}:{a}/{b}" for w, c, v, a, b in sens), passed=''))
    wcsv('0R_gate.csv', gates)
    for g in gates:
        print(g['rule'], g.get('passed'), g.get('value', ''))


if __name__ == '__main__':
    main()
