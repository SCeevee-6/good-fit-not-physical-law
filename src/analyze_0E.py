"""Analysis for protocol 0E. Reads <jobs>/*.json, writes CSVs (no pandas). Exact fractions for thresholds.
Country is the cluster. usage (repo root): python src/analyze_0E.py [--jobs DIR] [--out DIR]
"""
import os, sys, csv, json, glob
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
import numpy as np

ARGS = sys.argv[1:]
JOBS = ARGS[ARGS.index('--jobs') + 1] if '--jobs' in ARGS else os.path.join(ROOT, 'results', '0E_jobs')
OUT = ARGS[ARGS.index('--out') + 1] if '--out' in ARGS else os.path.join(ROOT, 'results')
NONSTRUCT = {'B1S', 'B2S', 'B2R'}
LAWS = ('B1S', 'B1I', 'B1R', 'B2S', 'B2R', 'B3p', 'B3m')
PRIMARY_W = 21
B, SEED = 2000, 20261126
K_MIN = 8                                # countries needed for a gate decision


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


def ok(r, fam=None):
    if r.get('status') != 'ok':
        return False
    st = json.loads(r['statuses'])
    return all(v == 'verified' for k, v in st.items() if fam is None or k in fam)


def boot(clusters, stat):
    rng = np.random.default_rng(SEED); m = len(clusters)
    reps = sorted(stat([clusters[i] for i in rng.integers(0, m, m)]) for _ in range(B))
    return reps[50], reps[1949]


def main():
    J, lost = {}, []
    for f in sorted(glob.glob(os.path.join(JOBS, '*.json'))):
        d = json.load(open(f)); cfg = d['cfg']
        if d.get('error') or d.get('timed_out'):
            lost.append(dict(slug=cfg['slug'], window=cfg['window'], state='error' if d.get('error') else 'timed_out', reason=(d.get('error') or d.get('reason') or '')[-200:])); continue
        J[(cfg['slug'], cfg['window'])] = d
    wcsv('0E_lost_configs.csv', lost)
    summ, cand_rows, wit_rows = [], [], []
    per = {}
    for (slug, w), d in sorted(J.items()):
        rows = d['rows']; byc = {r['cand']: r for r in rows}
        aes = [byc[c] for c in d['aes']]
        nv = sum(r['viol'] for r in aes)
        arms = {(a['cand'], a['struct'], a['arm']): a for a in d['arms']}
        pick_aes = [c for c in d['arm_picks'] if byc[c]['aes']]
        s = dict(country=slug, window=w, n_days=d['info']['n_days'], n_pool=d['n'], n_aes=d['n_aes'], n_aes_viol=nv, n_timeout=sum(r['timeout'] for r in rows),
                 primary_viol=byc[d['primary']]['viol'], aic_viol=byc[d['aic_pick']]['viol'], n_arm_aes=len(pick_aes), I_max=d['info']['I_max'], S_min=d['info']['S_min'],
                 pool_seconds=d['pool_seconds'])
        for st in ('Bpos', 'Bfac'):
            for arm in ('S', 'C', 'Rnd'):
                A = [arms[(c, st, arm)] for c in pick_aes if (c, st, arm) in arms]
                s[f'{st}_{arm}_all'] = sum(ok(a) for a in A); s[f'{st}_{arm}_nonstruct'] = sum(ok(a, NONSTRUCT) for a in A); s[f'{st}_{arm}_n'] = len(A)
                s[f'{st}_{arm}_med_iid'] = med([a.get('iid_ratio') for a in A if a['status'] == 'ok'])
                s[f'{st}_{arm}_med_ood'] = med([a.get('ood_ratio') for a in A if a['status'] == 'ok'])
        summ.append(s); per[(slug, w)] = (d, byc, arms, pick_aes)
        if w == PRIMARY_W:
            for c in d['arm_picks']:
                r = byc[c]
                row = dict(country=slug, cand=c, role=('primary' if c == d['primary'] else '') + ('aic' if c == d['aic_pick'] else ''), aes=r['aes'], degree=r['degree'], k=r['k'],
                           val_rmse=r['iid_rmse'], aic=r['aic'], test_rollout=r.get('roll_ood'), viol=r['viol'], laws=','.join(l for l in LAWS if r['status'][l] == 'counterexample'))
                for st in ('Bpos', 'Bfac'):
                    for arm in ('S', 'C', 'Rnd'):
                        a = arms.get((c, st, arm)); row[f'{st}_{arm}_adm'] = bool(a and ok(a)); row[f'{st}_{arm}_iid'] = a.get('iid_ratio') if a else ''; row[f'{st}_{arm}_ood'] = a.get('ood_ratio') if a else ''
                cand_rows.append(row)
            for c in (d['primary'], d['aic_pick']):
                for l, wv in (byc[c].get('wit') or {}).items():
                    wit_rows.append(dict(country=slug, cand=c, role='primary' if c == d['primary'] else 'aic', **{k: (json.dumps(v) if isinstance(v, list) else v) for k, v in wv.items()}))
    wcsv('0E_summary.csv', summ); wcsv('0E_candidates.csv', cand_rows); wcsv('0E_witnesses.csv', wit_rows)
    # per law, primary window, arm-set equivalent candidates pooled over countries: share satisfying each law in the pool and after each arm
    per_law = []
    for l in LAWS:
        row = dict(law=l, family='non-structural' if l in NONSTRUCT else 'structural')
        n = pv = 0; cnt = {}
        for (slug, w), (d, byc, arms, pick_aes) in per.items():
            if w != PRIMARY_W:
                continue
            for c in pick_aes:
                n += 1; pv += byc[c]['status'][l] == 'counterexample'
                for st in ('Bpos', 'Bfac'):
                    for arm in ('S', 'C', 'Rnd'):
                        a_ = arms.get((c, st, arm))
                        cnt[(st, arm)] = cnt.get((st, arm), 0) + bool(a_ and a_.get('status') == 'ok' and json.loads(a_['statuses'])[l] == 'verified')
        row.update(n=n, pool_violating=pv, pool_share=pv / max(1, n), **{f'{st}_{arm}_verified': v for (st, arm), v in cnt.items()})
        per_law.append(row)
    wcsv('0E_per_law.csv', per_law)

    P = [s for s in summ if s['window'] == PRIMARY_W]
    K = len(P); gates = []
    if K < K_MIN:
        gates.append(dict(rule='Verdict', value='undecided', stat=f'only {K} primary-window countries finished, need {K_MIN}', passed=''))
        wcsv('0E_gate.csv', gates); print('undecided: K =', K); return
    kp = sum(s['primary_viol'] for s in P)
    gates.append(dict(rule='G-E1', stat='countries whose validation-selected model breaks a law', count=kp, n=K, value=kp / K, threshold='>= 1/2 of countries', passed=2 * kp >= K))
    cl = [(s['n_aes_viol'], s['n_aes']) for s in P]
    k2, n2 = sum(a for a, _ in cl), sum(b for _, b in cl)
    lo, hi = boot(cl, lambda c: F(sum(a for a, _ in c), max(1, sum(b for _, b in c))))
    gates.append(dict(rule='G-E2', stat='share of equivalent candidates that break a law (pooled, country clusters)', count=k2, n=n2, value=k2 / max(1, n2), boot_lo=float(lo), boot_hi=float(hi),
                      threshold='>= 1/5 and boot_lo > 1/10', passed=(5 * k2 >= n2) and lo > F(1, 10)))
    kc = sum(s['Bpos_C_all'] for s in P); na = sum(s['Bpos_C_n'] for s in P)
    mi = med([x for s in P for x in [s['Bpos_C_med_iid']]]); mo = med([s['Bpos_C_med_ood'] for s in P])
    gates.append(dict(rule='G-E3', stat='B-pos C repairs the arm set at low cost (all seven laws)', count=kc, n=na, value=kc / max(1, na), med_iid_ratio=mi, med_ood_ratio=mo,
                      threshold='f(C) >= 9/10, medians <= 1.10 (median of country medians)', passed=(10 * kc >= 9 * na) and mi <= 1.10 and mo <= 1.10))
    keep = gates[0]['passed'] and gates[1]['passed']
    gates.append(dict(rule='Verdict', stat='second real series: claim retained (G-E1 and G-E2); cheap fix retained (also G-E3)', value='retained' if keep else ('inconclusive' if gates[0]['passed'] or gates[1]['passed'] else 'dropped'),
                      cheap_fix=bool(keep and gates[2]['passed']), passed=bool(keep)))
    # structure comparisons, non-structural family, country clusters
    for st in ('Bpos', 'Bfac'):
        cl2 = [(s[f'{st}_C_nonstruct'] - s[f'{st}_S_nonstruct'], s[f'{st}_C_nonstruct'] - s[f'{st}_Rnd_nonstruct'], s[f'{st}_C_n']) for s in P]
        n = sum(c[2] for c in cl2)
        for lab, j in (('C-S', 0), ('C-Rnd', 1)):
            lo, hi = boot(cl2, lambda c, j=j: F(sum(x[j] for x in c), max(1, sum(x[2] for x in c))))
            gates.append(dict(rule=f'{st} {lab}', stat='candidate-level difference in verified share, non-structural laws (report only)', count=sum(c[j] for c in cl2), n=n,
                              value=sum(c[j] for c in cl2) / max(1, n), boot_lo=float(lo), boot_hi=float(hi), passed=''))
    sens = [s for s in summ if s['window'] != PRIMARY_W]
    gates.append(dict(rule='Sensitivity', stat='primary breaks a law / equivalent share breaking, other windows (report only)',
                      value=';'.join(f"{s['country']}w{s['window']}:{int(s['primary_viol'])}:{s['n_aes_viol']}/{s['n_aes']}" for s in sens), passed=''))
    wcsv('0E_gate.csv', gates)
    for g in gates:
        print(g['rule'], g.get('passed'), g.get('value', ''))


if __name__ == '__main__':
    main()
