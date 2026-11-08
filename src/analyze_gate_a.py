"""Gate A analysis (protocol 0A + scoring note). Reads results/0A_jobs/*.json, writes
results/0A_runs.csv, results/0A_summary.csv, results/0A_gate.csv."""
import os, sys, csv, json, glob
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')
STRUCT = {'A': {'A1x', 'A1y'}, 'Amis': {'A1x', 'A1y'}, 'B': {'B1I', 'B1R', 'B3p', 'B3m'}}
ARMS = ('S', 'C', 'Rnd')


def load():
    recs, errs = [], []
    for p in sorted(glob.glob(os.path.join(R, '0A_jobs', '*.json'))):
        d = json.load(open(p))
        if d['error']:
            errs.append((d['sample_id'], d['error']))
        recs += d['records']
    return recs, errs


def score(r):
    if r.get('status') != 'ok':
        return dict(n_ok=0, n_tot=None, n_adm=0, s_adm=0, s_struct=None)
    st = json.loads(r['statuses']); s = STRUCT[r['key']]
    ns = [v == 'verified' for k, v in st.items() if k not in s]
    ss = [v == 'verified' for k, v in st.items() if k in s]
    return dict(n_ok=sum(ns), n_tot=len(ns), n_adm=int(all(ns)), s_adm=int(all(ss)), s_struct=ss)


def main():
    recs, errs = load()
    nspec = {'A': 1, 'Amis': 1, 'B': 3}
    rows = []
    for r in recs:
        sc = score(r)
        rows.append(dict(sample_id=r['sample_id'], key=r['key'], arm=r['arm'], status=r['status'], iters=r['iters'],
                         n_rows=r['n_rows'], nonstruct_all=sc['n_adm'], nonstruct_pass=sc['n_ok'], nonstruct_n=nspec[r['key']],
                         struct_all=sc['s_adm'], iid_ratio=r.get('iid_ratio', ''), k=r.get('k', ''),
                         timeouts=r.get('timeouts', ''), statuses=r.get('statuses', ''), seconds=r['seconds']))
    with open(os.path.join(R, '0A_runs.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    N = len({r['sample_id'] for r in rows})
    summ = []
    for grp in ('all', 'A', 'Amis', 'B'):
        for arm in ARMS:
            X = [r for r in rows if r['arm'] == arm and (grp == 'all' or r['key'] == grp)]
            ok = [r for r in X if r['status'] == 'ok']
            summ.append(dict(group=grp, arm=arm, n=len(X), n_ok=len(ok),
                             f_candidate=np.mean([r['nonstruct_all'] for r in X]),
                             f_pooled_spec=sum(r['nonstruct_pass'] for r in X) / sum(r['nonstruct_n'] for r in X),
                             struct_all=np.mean([r['struct_all'] for r in X]),
                             med_iid_ratio=np.median([r['iid_ratio'] for r in ok]) if ok else np.nan,
                             med_rows=np.median([r['n_rows'] for r in X]), med_iters=np.median([r['iters'] for r in X]),
                             timeouts=sum(int(r['timeouts'] or 0) for r in ok)))
    with open(os.path.join(R, '0A_summary.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
    a = {s['arm']: s for s in summ if s['group'] == 'all'}
    f = {k: a[k]['f_candidate'] for k in ARMS}; fp = {k: a[k]['f_pooled_spec'] for k in ARMS}
    ga1 = (f['C'] - f['S'] >= 0.20) and (a['C']['med_iid_ratio'] <= 1.10)
    ga2 = (f['C'] - f['Rnd'] >= 0.20)
    outcome = 1 if (ga1 and ga2) else (2 if not ga1 else 3)
    gate = [dict(criterion='G-A1', quantity='f(C)-f(S); med iid ratio C', value=f"{f['C']-f['S']:.4f}; {a['C']['med_iid_ratio']:.4f}",
                 threshold='>=0.20; <=1.10', passed=ga1, secondary_pooled=f"{fp['C']-fp['S']:.4f}"),
            dict(criterion='G-A2', quantity='f(C)-f(Rnd)', value=f"{f['C']-f['Rnd']:.4f}", threshold='>=0.20', passed=ga2,
                 secondary_pooled=f"{fp['C']-fp['Rnd']:.4f}"),
            dict(criterion='outcome', quantity='protocol 0A outcome', value=outcome, threshold='', passed='', secondary_pooled=''),
            dict(criterion='coverage', quantity='candidates with records; job errors', value=f"{N}; {len(errs)}", threshold='60; 0',
                 passed=(N == 60 and not errs), secondary_pooled='')]
    with open(os.path.join(R, '0A_gate.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(gate[0])); w.writeheader(); w.writerows(gate)
    return gate, summ, errs


if __name__ == '__main__':
    g, s, e = main()
    for x in g: print(x)
