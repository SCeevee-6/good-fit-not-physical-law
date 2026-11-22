"""Analysis for protocol 0B. Reads results/0B_jobs and results/0B_sample.csv, writes CSVs only (no pandas).
Gate thresholds are compared as exact fractions of counts. Bootstrap: 2000 replicates, percentile interval
(sorted replicate 50 and 1949, 0-based), rng seed 20261122; replicate statistics are kept as exact fractions.
usage (repo root): python src/analyze_0B.py
"""
import os, sys, csv, json, glob, math
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import numpy as np

RES = os.path.join(ROOT, 'results'); JOBS = os.path.join(RES, '0B_jobs')
B = 2000
SEED = 20261122
NONSTRUCT = {'A': {'A2'}, 'Amis': {'A2'}, 'B': {'B1S', 'B2S', 'B2R'}, 'H': {'H2'}}


def wcsv(name, rows):
    if not rows:
        rows = [dict(empty=1)]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(RES, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (float('nan'),) * 2
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def boot(clusters, stat):
    """clusters: list of per-cluster payloads; stat(list of payloads) -> Fraction. Returns (lo, hi) Fractions."""
    rng = np.random.default_rng(SEED)
    m = len(clusters)
    reps = sorted(stat([clusters[i] for i in rng.integers(0, m, m)]) for _ in range(B))
    return reps[50], reps[1949]


def fx(x):
    return float(x) if isinstance(x, F) else x


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else float('nan')


def load_pool():
    P = {}
    for f in glob.glob(os.path.join(JOBS, 'pool', '*.json')):
        d = json.load(open(f)); P[d['cell']] = d
    return P


def load_arms():
    A = {}
    for f in glob.glob(os.path.join(JOBS, 'arms', '*.json')):
        d = json.load(open(f))
        for r in d['records']:
            A[(d['sample_id'], d['struct'] or '', r['arm'])] = r
    return A


def adm(r, fam):
    if r is None or r.get('status') != 'ok':
        return False
    st = json.loads(r['statuses'])
    return all(st[s] == 'verified' for s in st if fam is None or s in fam)


# ------------------------------------------------------------------------------------------------ 0B-1
def unit1(P, A, sample, gates):
    cells = [P[c] for c in sorted(P) if c.startswith('H_')]
    rows, per_cell = [], []
    for c in cells:
        aes = [r for r in c['rows'] if r['aes']]
        law = {s: sum(r['status'][s] == 'counterexample' for r in aes) for s in ('H1x', 'H1y', 'H2')}
        per_cell.append((sum(r['viol'] for r in aes), len(aes)))
        rows.append(dict(cell=c['cell'], noise_idx=c['noise_idx'], seed=c['seed'], n_pool=c['n'], n_aes=len(aes),
                         n_aes_viol=per_cell[-1][0], n_aes_timeout=sum(r['timeout'] for r in aes), n_pool_viol=sum(r['viol'] for r in c['rows']),
                         **{f'aes_viol_{k}': v for k, v in law.items()}, max_z3_sec=max(max(r['z3_sec'].values()) for r in c['rows']), seconds=c['seconds']))
    wcsv('0B1_pool_cells.csv', rows)
    k = sum(a for a, _ in per_cell); n = sum(b for _, b in per_cell)
    lo, hi = boot(per_cell, lambda cl: F(sum(a for a, _ in cl), max(1, sum(b for _, b in cl))))
    wl, wh = wilson(k, n)
    g1 = (5 * k >= n) and lo > F(1, 10)
    gates.append(dict(rule='G-B1', stat='fraction of accuracy-equivalent H candidates violating a law', count=k, n=n, value=k / n,
                      boot_lo=fx(lo), boot_hi=fx(hi), wilson_lo=wl, wilson_hi=wh, threshold='>= 1/5 and boot_lo > 1/10', passed=g1))

    # arms on the H sample (both noise levels)
    Hs = [s for s in sample if s['key'] == 'H']
    poolrow = {(c['cell'], r['cand']): r for c in cells for r in c['rows']}
    arm_rows, summ = [], []
    for s in Hs:
        cid = s['sample_id'].rsplit('_c', 1)[0]
        pr = poolrow[(cid, int(s['cand']))]
        arm_rows.append(dict(sample_id=s['sample_id'], noise_idx=s['noise_idx'], seed=s['seed'], arm='pool', status='ok',
                             adm_all=not pr['viol'] and not pr['timeout'], adm_H2=pr['status']['H2'] == 'verified',
                             **{f'st_{q}': pr['status'][q] for q in ('H1x', 'H1y', 'H2')}, iid_ratio=1.0, ood_ratio=1.0, iters=0, n_rows=0))
        for arm in ('S', 'C', 'Rnd'):
            r = A.get((s['sample_id'], '', arm))
            st = json.loads(r['statuses']) if r and r.get('status') == 'ok' else {}
            arm_rows.append(dict(sample_id=s['sample_id'], noise_idx=s['noise_idx'], seed=s['seed'], arm=arm, status=r['status'] if r else 'missing',
                                 adm_all=adm(r, None), adm_H2=adm(r, {'H2'}), **{f'st_{q}': st.get(q, '') for q in ('H1x', 'H1y', 'H2')},
                                 iid_ratio=r.get('iid_ratio') if r else None, ood_ratio=r.get('ood_ratio') if r else None,
                                 iters=r.get('iters') if r else None, n_rows=r.get('n_rows') if r else None))
    wcsv('0B1_arms.csv', arm_rows)
    for grp in ('pooled', '1', '2'):
        for arm in ('pool', 'S', 'C', 'Rnd'):
            R = [r for r in arm_rows if r['arm'] == arm and (grp == 'pooled' or str(r['noise_idx']) == grp)]
            ka = sum(r['adm_all'] for r in R); kh = sum(r['adm_H2'] for r in R); nn = len(R)
            summ.append(dict(noise=grp, arm=arm, n=nn, f_all=ka / max(nn, 1), f_all_lo=wilson(ka, nn)[0], f_all_hi=wilson(ka, nn)[1],
                             f_H2=kh / max(nn, 1), **{f'f_{q}': sum(r[f'st_{q}'] == 'verified' for r in R) / max(nn, 1) for q in ('H1x', 'H1y')},
                             med_iid_ratio=med([r['iid_ratio'] for r in R if r['status'] == 'ok']),
                             med_ood_ratio=med([r['ood_ratio'] for r in R if r['status'] == 'ok']),
                             med_rows=med([r['n_rows'] for r in R if r['status'] == 'ok'])))
    wcsv('0B1_summary.csv', summ)
    C = [r for r in arm_rows if r['arm'] == 'C']
    kc = sum(r['adm_all'] for r in C); mo = med([r['ood_ratio'] for r in C if r['status'] == 'ok'])
    gates.append(dict(rule='G-B2', stat='f(C) on H (all laws) and median OOD rollout ratio of C', count=kc, n=len(C), value=kc / max(len(C), 1),
                      med_ood_ratio_C=mo, threshold='f(C) >= 9/10 and median OOD ratio <= 1.10', passed=(10 * kc >= 9 * len(C)) and mo <= 1.10))
    # candidate volumes for G-B3
    vols = []
    for c in cells:
        for r in c['rows']:
            if r['viol']:
                vols.append(dict(cell=c['cell'], cand=r['cand'], aes=r['aes'], k=r['k'], degree=r['degree'],
                                 laws=','.join(sorted(r['volume'])), vol_max=max(r['volume'].values()),
                                 **{f'vol_{q}': r['volume'].get(q, '') for q in ('H1x', 'H1y', 'H2')},
                                 **{f'depth_{q}': r['depth'].get(q, '') for q in ('H1x', 'H1y', 'H2')}))
    wcsv('0B1_violation_volumes.csv', vols)
    return vols


# ------------------------------------------------------------------------------------------------ 0B-2
def unit2(vols, gates):
    N = [json.load(open(f)) for f in glob.glob(os.path.join(JOBS, 'needles', '*.json'))]
    rows = []
    for d in sorted(N, key=lambda x: x['needle_id']):
        r = {k: d.get(k) for k in ('needle_id', 'base', 'sys', 'degree', 'spec', 'kappa', 'v_target', 'status', 'v', 'depth', 'reached',
                                    'corner_k', 'power_d', 'verifier', 'verifier_method', 'verifier_sec', 'n_seeds', 'seconds', 'error')}
        for n, h in (d.get('hits') or {}).items():
            r[f'hits_{n}'] = h
            r[f'ref_detect_{n}'] = 1 - (1 - d['v']) ** int(n)
        rows.append(r)
    wcsv('0B2_needles.csv', rows)
    use = [r for r in rows if r['status'] == 'built' and r['v'] and r['v'] > 0]
    fits = []
    v50 = {}
    for n in ('2000', '20000', '200000'):
        etas = np.exp(np.linspace(np.log(0.01), np.log(100), 20001))
        ll = np.zeros_like(etas)
        for r in use:
            S_, h, v = r['n_seeds'], r[f'hits_{n}'], r['v']
            a = int(n) * etas * math.log1p(-v)                 # log miss probability
            ll += (S_ - h) * a + (h * np.log(-np.expm1(a)) if h else 0.0)
        eta = float(etas[int(np.argmax(ll))])
        v50[n] = -math.expm1(math.log(0.5) / (int(n) * eta))
        miss = sum(r['n_seeds'] - r[f'hits_{n}'] for r in use); tot = sum(r['n_seeds'] for r in use)
        ref = sum(r['n_seeds'] * (1 - r[f'ref_detect_{n}']) for r in use)
        fits.append(dict(budget=int(n), n_needles=len(use), eta_ml=eta, v50=v50[n], v50_reference=-math.expm1(math.log(0.5) / int(n)),
                         observed_misses=miss, reference_misses=ref, trials=tot))
    bins = [(10.0 ** (-e - 0.5), 10.0 ** (-e + 0.5)) for e in range(2, 7)]
    for lo, hi in bins:
        R = [r for r in use if lo <= r['v'] < hi]
        for n in ('2000', '20000', '200000'):
            fits.append(dict(budget=int(n), bin_lo=lo, bin_hi=hi, n_needles=len(R),
                             obs_miss_rate=(sum(r['n_seeds'] - r[f'hits_{n}'] for r in R) / max(1, sum(r['n_seeds'] for r in R))) if R else '',
                             ref_miss_rate=(sum(r['n_seeds'] * (1 - r[f'ref_detect_{n}']) for r in R) / max(1, sum(r['n_seeds'] for r in R))) if R else ''))
    wcsv('0B2_detector_fit.csv', fits)
    aesv = [v for v in vols if v['aes']]
    below = sum(v['vol_max'] < v50['20000'] for v in aesv)
    gates.append(dict(rule='G-B3', stat='share of accuracy-equivalent H violators with volume < v50(20k)', count=below, n=len(aesv),
                      value=below / max(1, len(aesv)), v50_20k=v50['20000'], threshold='>= 1/20', passed=20 * below >= len(aesv),
                      verifier_found=sum(r['verifier'] == 'counterexample' for r in use), needles_built=len(use)))


# ------------------------------------------------------------------------------------------------ 0B-3
def unit3(A, sample, gates):
    keys = [('A', 'A', ''), ('Amis', 'Amis', ''), ('Bpos', 'B', 'Bpos'), ('Bfac', 'B', 'Bfac'), ('H', 'H', '')]
    per = {}                                        # (label, seed) -> list of (S, C, Rnd) admissibility
    rows = []
    for lab, key, struct in keys:
        for s in sample:
            if s['key'] != key or (key == 'H' and str(s['noise_idx']) != '1'):
                continue
            fam = NONSTRUCT[key]
            a = {arm: adm(A.get((s['sample_id'], struct, arm)), fam) for arm in ('S', 'C', 'Rnd')}
            r = {arm: A.get((s['sample_id'], struct, arm)) for arm in ('S', 'C', 'Rnd')}
            per.setdefault((lab, int(s['seed'])), []).append(a)
            rows.append(dict(label=lab, sample_id=s['sample_id'], seed=s['seed'], **{f'adm_{k}': v for k, v in a.items()},
                             **{f'iid_ratio_{k}': (r[k] or {}).get('iid_ratio') for k in r}, **{f'status_{k}': (r[k] or {}).get('status', 'missing') for k in r}))
    wcsv('0B3_candidates.csv', rows)
    out = []
    groups = {lab: [lab] for lab, _, _ in keys}
    groups['A+Amis+Bpos'] = ['A', 'Amis', 'Bpos']
    for g, labs in groups.items():
        clusters = [sum((per.get((l, sd), []) for l in labs), []) for sd in range(20)]
        clusters = [c for c in clusters if c]
        allc = sum(clusters, [])
        n = len(allc)
        f = {arm: sum(c[arm] for c in allc) for arm in ('S', 'C', 'Rnd')}
        def diff(x, y):
            return lambda cl: F(sum(c[x] - c[y] for cc in cl for c in cc), max(1, sum(len(cc) for cc in cl)))
        cs = boot(clusters, diff('C', 'S')); cr = boot(clusters, diff('C', 'Rnd'))
        med_c = med([r['iid_ratio_C'] for r in rows if r['label'] in labs and r['status_C'] == 'ok'])
        rec = dict(group=g, n=n, **{f'f_{a}': f[a] / max(1, n) for a in f}, **{f'f_{a}_wlo': wilson(f[a], n)[0] for a in f},
                   **{f'f_{a}_whi': wilson(f[a], n)[1] for a in f}, d_CS=(f['C'] - f['S']) / max(1, n), d_CS_lo=fx(cs[0]), d_CS_hi=fx(cs[1]),
                   d_CRnd=(f['C'] - f['Rnd']) / max(1, n), d_CRnd_lo=fx(cr[0]), d_CRnd_hi=fx(cr[1]), med_iid_ratio_C=med_c)
        out.append(rec)
        if g == 'A+Amis+Bpos':
            gates.append(dict(rule='Replication', stat='pooled A+Amis+Bpos: C-Rnd interval reaches below 1/5, C-S interval above 0',
                              n=n, d_CS_lo=fx(cs[0]), d_CS_hi=fx(cs[1]), d_CRnd_lo=fx(cr[0]), d_CRnd_hi=fx(cr[1]),
                              threshold='d_CRnd_lo < 1/5 and d_CS_lo > 0', passed=(cr[0] < F(1, 5)) and (cs[0] > 0)))
    wcsv('0B3_summary.csv', out)


def main():
    P = load_pool(); A = load_arms()
    sample = list(csv.DictReader(open(os.path.join(RES, '0B_sample.csv'))))
    gates = []
    vols = unit1(P, A, sample, gates)
    unit2(vols, gates)
    unit3(A, sample, gates)
    wcsv('0B_gate.csv', gates)
    for g in gates:
        print(g['rule'], 'PASS' if g['passed'] else 'FAIL')


if __name__ == '__main__':
    main()
