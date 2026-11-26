"""Analysis for protocol 0D. CSVs only, exact fractions, seeds as bootstrap clusters (2000 replicates, rng 20261128).
usage (repo root): python src/analyze_0D.py
"""
import os, sys, csv, json, glob
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
import numpy as np
RES = os.path.join(ROOT, 'results'); JOBS = os.path.join(RES, '0D_jobs')
NONSTRUCT = ('E1S', 'E2S', 'E2R')
ARMS = ('S', 'C', 'Rnd', 'SoftA', 'SoftB')
B, SEED = 2000, 20261128


def wcsv(name, rows):
    rows = rows or [dict(empty=1)]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(RES, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def ok(r, fam):
    if r is None or r.get('status') != 'ok':
        return False
    st = json.loads(r['statuses'])
    return all(st[s] == 'verified' for s in (fam or st))


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else float('nan')


def boot(clusters, stat):
    rng = np.random.default_rng(SEED); m = len(clusters)
    reps = sorted(stat([clusters[i] for i in rng.integers(0, m, m)]) for _ in range(B))
    return float(reps[50]), float(reps[1949])


def main():
    recs = {}
    for f in glob.glob(os.path.join(JOBS, 'arms', '*.json')):
        d = json.load(open(f))
        if d.get('error'):
            print('ERROR', f, d['error'][-300:]); continue
        for r in d['records']:
            recs[(d['sample_id'], r['arm'])] = r
    ids = sorted({i for i, _ in recs})
    seed = {i: int(next(r for (j, a), r in recs.items() if j == i)['seed']) for i in ids}
    rows, per = [], {}
    for arm in ARMS:
        R = [recs.get((i, arm)) for i in ids]
        okr = [r for r in R if r and r.get('status') == 'ok']
        k = sum(ok(r, NONSTRUCT) for r in R); ka = sum(ok(r, None) for r in R)
        per[arm] = k
        rows.append(dict(arm=arm, n=len(ids), verified_nonstruct=k, share=k / max(1, len(ids)), verified_all=ka,
                         med_iid_ratio=med([r.get('iid_ratio') for r in okr]), med_rows=med([r.get('n_rows') for r in okr]),
                         n_not_ok=len(ids) - len(okr), n_timeouts_in_models=sum(r.get('timeouts', 0) for r in okr),
                         med_seconds=med([r.get('seconds') for r in okr])))
    wcsv('0D_arms.csv', rows)
    # per law, per arm
    laws = [r_ for r_ in json.loads(next(r for r in recs.values() if r.get('status') == 'ok')['statuses'])]
    pl = []
    for l in laws:
        row = dict(law=l, family='non-structural' if l in NONSTRUCT else 'structural')
        for arm in ARMS:
            row[arm] = sum(1 for i in ids if recs.get((i, arm)) and recs[(i, arm)].get('status') == 'ok' and json.loads(recs[(i, arm)]['statuses'])[l] == 'verified')
        pl.append(row)
    wcsv('0D_per_law.csv', pl)
    n = len(ids)
    fC = F(per['C'], n)
    mrC = med([recs[(i, 'C')].get('n_rows') for i in ids if recs.get((i, 'C')) and recs[(i, 'C')].get('status') == 'ok'])
    mrB = med([recs[(i, 'SoftB')].get('n_rows') for i in ids if recs.get((i, 'SoftB')) and recs[(i, 'SoftB')].get('status') == 'ok'])
    gates = []
    for arm in ('SoftA', 'SoftB'):
        cl = []
        for sd in range(20):
            G = [i for i in ids if seed[i] == sd]
            cl.append((sum(ok(recs.get((i, 'C')), NONSTRUCT) - ok(recs.get((i, arm)), NONSTRUCT) for i in G), len(G)))
        lo, hi = boot(cl, lambda c: F(sum(a for a, _ in c), max(1, sum(b for _, b in c))))
        gates.append(dict(rule=f'C - {arm}', stat='verified share on non-structural laws', value=float(fC - F(per[arm], n)), boot_lo=lo, boot_hi=hi, passed=''))
    a_ = fC - F(per['SoftA'], n) >= F(1, 10)
    b_ = fC - F(per['SoftB'], n) >= F(1, 10)
    gates.append(dict(rule='R-D1a', stat='C verifies >= 1/10 more than the pilot-budget grid (SoftA)', value=float(fC - F(per['SoftA'], n)), passed=bool(a_)))
    gates.append(dict(rule='R-D1b', stat='C verifies >= 1/10 more than the fine grid (SoftB)', value=float(fC - F(per['SoftB'], n)), passed=bool(b_)))
    gates.append(dict(rule='R-D1', stat='coverage advantage of the loop at four states (R-D1a or R-D1b)', value=bool(a_ or b_), passed=bool(a_ or b_)))
    # constraint rows against number of states, descriptive: 0B arms (C rows) and 0C soft arms (pilot-budget grid rows) for A and B, 0D for E4
    rb = []
    s0b = {r['sample_id']: r for r in csv.DictReader(open(os.path.join(RES, '0B_sample.csv')))}
    crow, srow = {'A': [], 'B': []}, {'A': [], 'B': []}
    for f in glob.glob(os.path.join(RES, '0B_jobs', 'arms', '*.json')):
        d = json.load(open(f)); k = s0b[d['sample_id']]['key']
        if k in crow and (k != 'B' or d['struct'] == 'Bpos'):
            crow[k] += [r.get('n_rows') for r in d['records'] if r['arm'] == 'C' and r.get('status') == 'ok']
    for f in glob.glob(os.path.join(JOBS, '..', '0C_jobs', 'soft', '*.json')):
        d = json.load(open(f)); k = s0b[d['sample_id']]['key']
        if k in srow:
            srow[k] += [r.get('n_rows') for r in d['records'] if r['arm'] == 'Soft2' and r.get('status') == 'ok']
    for k, ns in (('A', 2), ('B', 3)):
        rb.append(dict(system=k, states=ns, med_rows_C=med(crow[k]), med_rows_grid_pilot_budget=med(srow[k]), med_rows_grid_fine=''))
    rb.append(dict(system='E4', states=4, med_rows_C=mrC, med_rows_grid_pilot_budget=med([recs[(i, 'SoftA')].get('n_rows') for i in ids if recs.get((i, 'SoftA')) and recs[(i, 'SoftA')].get('status') == 'ok']), med_rows_grid_fine=mrB))
    wcsv('0D_rows_by_dimension.csv', rb)
    xc = dis = unsound = nm = 0
    for r in recs.values():
        if r.get('status') == 'ok':
            nm += 1; xc += len(r.get('xcheck', [])); dis += r.get('xcheck_disagree', 0); unsound += r.get('unsound', 0)
    gates.append(dict(rule='R-D2', stat='Bernstein vs z3 disagreements among laws both decided (must be 0)', count=dis, n=xc, value=dis, passed=dis == 0))
    gates.append(dict(rule='R-D3', stat='verified laws with a sampled point below -tol_strong (must be 0)', count=unsound, n=nm, value=unsound, passed=unsound == 0))
    # Gate A style, reported only
    for lab, a, b2 in (('C - S', 'C', 'S'), ('C - Rnd', 'C', 'Rnd')):
        cl = [(sum(ok(recs.get((i, a)), NONSTRUCT) - ok(recs.get((i, b2)), NONSTRUCT) for i in ids if seed[i] == sd), sum(1 for i in ids if seed[i] == sd)) for sd in range(20)]
        lo, hi = boot(cl, lambda c: F(sum(x for x, _ in c), max(1, sum(y for _, y in c))))
        gates.append(dict(rule=f'E4 {lab}', stat='verified share, non-structural laws (report only)', value=(per[a] - per[b2]) / n, boot_lo=lo, boot_hi=hi, passed=''))
    wcsv('0D_gate.csv', gates)
    for g in gates:
        print(g['rule'], g.get('passed'), g.get('value'))


if __name__ == '__main__':
    main()
