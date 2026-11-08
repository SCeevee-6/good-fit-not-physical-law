"""Gate A re-decision (protocol 0A2). A and A-mis records from results/0A_jobs (0A, reused), B records from
results/0A2_jobs/<struct>. Writes results/0A2_runs.csv, 0A2_summary.csv, 0A2_per_law.csv, 0A2_gate.csv."""
import os, sys, csv, json, glob
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')
STRUCT = {'A': {'A1x', 'A1y'}, 'Amis': {'A1x', 'A1y'}, 'B': {'B1I', 'B1R', 'B3p', 'B3m'}}
NONSTRUCT = {'A': ['A2'], 'Amis': ['A2'], 'B': ['B1S', 'B2S', 'B2R']}
ARMS = ('S', 'C', 'Rnd')
BSTRUCTS = ('Bpos', 'Bfac')


def _load(pattern):
    recs, errs = [], []
    for p in sorted(glob.glob(pattern)):
        d = json.load(open(p))
        if d['error']:
            errs.append((d['sample_id'], d['error']))
        recs += d['records']
    return recs, errs


def row(r, struct):
    ok = r.get('status') == 'ok'
    st = json.loads(r['statuses']) if ok else {}
    ns = NONSTRUCT[r['key']]
    npass = sum(st.get(s) == 'verified' for s in ns)
    sfam = STRUCT[r['key']]
    return dict(sample_id=r['sample_id'], key=r['key'], struct=struct, arm=r['arm'], status=r['status'],
                iters=r['iters'], n_rows=r['n_rows'], stop=r.get('stop', ''),
                nonstruct_all=int(ok and npass == len(ns)), nonstruct_pass=npass, nonstruct_n=len(ns),
                struct_all=int(ok and all(st.get(s) == 'verified' for s in sfam)),
                iid_ratio=r.get('iid_ratio', '') if ok else '', k=r.get('k', ''),
                timeouts=r.get('timeouts', '') if ok else '', statuses=r.get('statuses', ''), seconds=r['seconds'])


def summarise(X, group, struct, arm):
    ok = [r for r in X if r['status'] == 'ok']
    return dict(group=group, struct=struct, arm=arm, n=len(X), n_ok=len(ok),
                f_candidate=float(np.mean([r['nonstruct_all'] for r in X])) if X else np.nan,
                f_pooled_spec=sum(r['nonstruct_pass'] for r in X) / max(1, sum(r['nonstruct_n'] for r in X)),
                struct_all=float(np.mean([r['struct_all'] for r in X])) if X else np.nan,
                med_iid_ratio=float(np.median([r['iid_ratio'] for r in ok])) if ok else np.nan,
                med_rows=float(np.median([r['n_rows'] for r in X])) if X else np.nan,
                med_iters=float(np.median([r['iters'] for r in X])) if X else np.nan,
                timeouts=sum(int(r['timeouts'] or 0) for r in ok), n_infeasible=sum(r['status'] == 'infeasible' for r in X),
                rnd_clean_but_violating=sum(1 for r in X if r['arm'] == 'Rnd' and r['stop'] == 'clean' and not r['nonstruct_all']))


def gate(rows):
    a = {arm: summarise([r for r in rows if r['arm'] == arm], 'all', '', arm) for arm in ARMS}
    f = {k: a[k]['f_candidate'] for k in ARMS}; fp = {k: a[k]['f_pooled_spec'] for k in ARMS}
    ga1 = (f['C'] - f['S'] >= 0.20) and (a['C']['med_iid_ratio'] <= 1.10)
    ga2 = (f['C'] - f['Rnd'] >= 0.20)
    return dict(dCS=f['C'] - f['S'], medC=a['C']['med_iid_ratio'], dCR=f['C'] - f['Rnd'], ga1=ga1, ga2=ga2,
                outcome=1 if (ga1 and ga2) else (2 if not ga1 else 3), pdCS=fp['C'] - fp['S'], pdCR=fp['C'] - fp['Rnd'],
                n=len({(r['key'], r['sample_id']) for r in rows}))


def main():
    rA, eA = _load(os.path.join(R, '0A_jobs', '*.json'))
    base = [row(r, '0A') for r in rA if r['key'] in ('A', 'Amis')]
    B0 = [row(r, 'B0A') for r in rA if r['key'] == 'B']
    Bs, errs = {}, list(eA)
    for s in BSTRUCTS:
        rec, e = _load(os.path.join(R, '0A2_jobs', s, '*.json'))
        Bs[s] = [row(r, s) for r in rec]; errs += e
    allrows = base + B0 + Bs['Bpos'] + Bs['Bfac']
    fields = list(allrows[0])
    with open(os.path.join(R, '0A2_runs.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(allrows)
    summ = []
    for grp, struct, X in [('A', '0A', [r for r in base if r['key'] == 'A']), ('Amis', '0A', [r for r in base if r['key'] == 'Amis']),
                           ('B', 'B0A', B0), ('B', 'Bpos', Bs['Bpos']), ('B', 'Bfac', Bs['Bfac'])]:
        for arm in ARMS:
            summ.append(summarise([r for r in X if r['arm'] == arm], grp, struct, arm))
    for s in BSTRUCTS:
        for arm in ARMS:
            summ.append(summarise([r for r in base + Bs[s] if r['arm'] == arm], 'all', 'A+Amis+' + s, arm))
    with open(os.path.join(R, '0A2_summary.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(summ[0])); w.writeheader(); w.writerows(summ)
    # per law: verified fraction and median held-out ratio of the models verified on that law
    per = []
    for grp, struct, X in [('A', '0A', [r for r in base if r['key'] == 'A']), ('Amis', '0A', [r for r in base if r['key'] == 'Amis']),
                           ('B', 'B0A', B0), ('B', 'Bpos', Bs['Bpos']), ('B', 'Bfac', Bs['Bfac'])]:
        for arm in ARMS:
            Y = [r for r in X if r['arm'] == arm]
            ok = [r for r in Y if r['status'] == 'ok']
            specs = sorted({k for r in ok for k in json.loads(r['statuses'])})
            for sp in specs:
                v = [json.loads(r['statuses'])[sp] == 'verified' for r in ok]
                per.append(dict(group=grp, struct=struct, arm=arm, spec=sp, family='structural' if sp in STRUCT[grp] else 'non-structural',
                                n=len(Y), verified=int(sum(v)), frac=sum(v) / len(Y) if Y else np.nan))
    with open(os.path.join(R, '0A2_per_law.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(per[0])); w.writeheader(); w.writerows(per)
    G = []
    for s in BSTRUCTS:
        g = gate(base + Bs[s])
        G.append(dict(decision='primary' if s == 'Bpos' else 'reported', b_structure=s, n_candidates=g['n'],
                      f_C_minus_f_S=round(g['dCS'], 4), med_iid_ratio_C=round(g['medC'], 4), G_A1=g['ga1'],
                      f_C_minus_f_Rnd=round(g['dCR'], 4), G_A2=g['ga2'], outcome=g['outcome'],
                      pooled_C_minus_S=round(g['pdCS'], 4), pooled_C_minus_Rnd=round(g['pdCR'], 4), job_errors=len(errs)))
    with open(os.path.join(R, '0A2_gate.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(G[0])); w.writeheader(); w.writerows(G)
    return G, summ, per, errs


if __name__ == '__main__':
    G, s, p, e = main()
    for x in G: print(x)
