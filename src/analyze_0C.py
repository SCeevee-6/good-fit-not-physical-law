"""Analysis for protocol 0C. CSVs only, exact fractions for thresholds, seeds as bootstrap clusters (2000 replicates, rng 20261127).
usage (repo root): python src/analyze_0C.py
"""
import os, sys, csv, json, glob
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
import numpy as np
RES = os.path.join(ROOT, 'results'); JOBS = os.path.join(RES, '0C_jobs')
NONSTRUCT = {'A': {'A2'}, 'Amis': {'A2'}, 'B': {'B1S', 'B2S', 'B2R'}, 'H': {'H2'}}
B, SEED = 2000, 20261127


def wcsv(name, rows):
    rows = rows or [dict(empty=1)]
    keys = []
    for r in rows:
        keys += [k for k in r if k not in keys]
    with open(os.path.join(RES, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def adm(r, key):
    if r is None or r.get('status') != 'ok':
        return False
    st = json.loads(r['statuses'])
    return all(st[s] == 'verified' for s in NONSTRUCT[key])


def adm_all(r):
    return r is not None and r.get('status') == 'ok' and all(v == 'verified' for v in json.loads(r['statuses']).values())


def med(xs):
    xs = [x for x in xs if x is not None and np.isfinite(x)]
    return float(np.median(xs)) if xs else float('nan')


def boot(clusters, stat):
    rng = np.random.default_rng(SEED); m = len(clusters)
    reps = sorted(stat([clusters[i] for i in rng.integers(0, m, m)]) for _ in range(B))
    return float(reps[50]), float(reps[1949])


def load(dirn):
    out = {}
    for f in glob.glob(os.path.join(JOBS, dirn, '*.json')):
        d = json.load(open(f))
        if d.get('error'):
            print('ERROR', f, d['error'][-200:]); continue
        for r in d['records']:
            out[(d['sample_id'], r['arm'])] = r
    return out


def load_0b():
    sample = {r['sample_id']: r for r in csv.DictReader(open(os.path.join(RES, '0B_sample.csv')))}
    out = {}
    for f in glob.glob(os.path.join(RES, '0B_jobs', 'arms', '*.json')):
        d = json.load(open(f)); s = sample[d['sample_id']]
        if s['key'] == 'B' and d['struct'] != 'Bpos':
            continue
        for r in d['records']:
            out[(d['sample_id'], r['arm'])] = r
    return sample, out


def main():
    unsound, rows_unsound = 0, []
    # ---- 0C-1: soft penalty on the 0B sample, noise 1
    s0b, a0b = load_0b(); soft = load('soft')
    ids = [i for i, s in s0b.items() if s['key'] in ('A', 'Amis', 'B', 'H') and s['noise_idx'] == '1']
    arms = ['S', 'C', 'Rnd', 'Soft1', 'Soft2', 'Soft3']
    get = lambda i, a: (a0b if a in ('S', 'C', 'Rnd') else soft).get((i, a))
    t1 = []
    for grp in ('A', 'Amis', 'B', 'H', 'pooled'):
        G = [i for i in ids if grp == 'pooled' or s0b[i]['key'] == grp]
        for arm in arms:
            R = [(i, get(i, arm)) for i in G]
            k = sum(adm(r, s0b[i]['key']) for i, r in R); ka = sum(adm_all(r) for i, r in R)
            ok = [r for i, r in R if r and r.get('status') == 'ok']
            t1.append(dict(group=grp, arm=arm, n=len(G), verified_nonstruct=k, share=k / max(1, len(G)), verified_all=ka,
                           med_iid_ratio=med([r.get('iid_ratio') for r in ok]), n_not_ok=len(G) - len(ok)))
    wcsv('0C1_soft.csv', t1)
    P = {r['arm']: r for r in t1 if r['group'] == 'pooled'}
    fC = F(P['C']['verified_nonstruct'], P['C']['n'])
    need = []
    for arm in ('Soft1', 'Soft2', 'Soft3'):
        fs = F(P[arm]['verified_nonstruct'], P[arm]['n'])
        need.append(dict(arm=arm, f_soft=float(fs), f_C=float(fC), gap=float(fC - fs), med_iid=P[arm]['med_iid_ratio'],
                         loop_needed_vs_this=bool((fC - fs >= F(1, 10)) or P[arm]['med_iid_ratio'] > 1.10)))
    # per-seed clusters for the C - Soft gap
    gates = []
    for r in need:
        arm = r['arm']
        cl = []
        for sd in range(20):
            G = [i for i in ids if int(s0b[i]['seed']) == sd]
            cl.append((sum(adm(get(i, 'C'), s0b[i]['key']) - adm(get(i, arm), s0b[i]['key']) for i in G), len(G)))
        lo, hi = boot(cl, lambda c: F(sum(a for a, _ in c), max(1, sum(b for _, b in c))))
        gates.append(dict(rule=f'R-C1 {arm}', stat='C minus soft penalty, verified share on non-structural laws, pooled', value=r['gap'], boot_lo=lo, boot_hi=hi,
                          f_soft=r['f_soft'], med_iid_soft=r['med_iid'], loop_needed=r['loop_needed_vs_this'], passed=''))
    gates.append(dict(rule='R-C1', stat='counterexample loop needed over soft penalty: for every lambda, gap >= 1/10 or median validation ratio > 1.10',
                      value=all(r['loop_needed_vs_this'] for r in need), passed=all(r['loop_needed_vs_this'] for r in need)))

    # ---- 0C-2: noise sensitivity
    samp = {r['sample_id']: r for r in csv.DictReader(open(os.path.join(RES, '0C_sample.csv')))}
    ar = load('arms')
    t2, ref = [], []
    def level(lab, ids_, getter, keyof, seedof):
        for grp, keys in (('A+Amis+B', ('A', 'Amis', 'B')), ('A', ('A',)), ('Amis', ('Amis',)), ('B', ('B',)), ('H', ('H',))):
            G = [i for i in ids_ if keyof(i) in keys]
            if not G:
                continue
            f = {a: sum(adm(getter(i, a), keyof(i)) for i in G) for a in ('S', 'C', 'Rnd')}
            cl = []
            for sd in range(20):
                Gs = [i for i in G if seedof(i) == sd]
                cl.append((sum(adm(getter(i, 'C'), keyof(i)) - adm(getter(i, 'S'), keyof(i)) for i in Gs), sum(adm(getter(i, 'C'), keyof(i)) - adm(getter(i, 'Rnd'), keyof(i)) for i in Gs), len(Gs)))
            cs = boot(cl, lambda c: F(sum(x[0] for x in c), max(1, sum(x[2] for x in c)))); cr = boot(cl, lambda c: F(sum(x[1] for x in c), max(1, sum(x[2] for x in c))))
            ok = [getter(i, 'C') for i in G if getter(i, 'C') and getter(i, 'C').get('status') == 'ok']
            row = dict(noise=lab, group=grp, n=len(G), f_S=f['S'] / len(G), f_C=f['C'] / len(G), f_Rnd=f['Rnd'] / len(G), d_CS=(f['C'] - f['S']) / len(G), d_CS_lo=cs[0], d_CS_hi=cs[1],
                       d_CRnd=(f['C'] - f['Rnd']) / len(G), d_CRnd_lo=cr[0], d_CRnd_hi=cr[1], med_iid_C=med([r.get('iid_ratio') for r in ok]))
            t2.append(row)
            if grp == 'A+Amis+B':
                gates.append(dict(rule=f'R-C2 noise {lab}', stat='pooled A+Amis+B-pos: C-S lower > 0 and C-Rnd lower < 1/5', value=row['d_CRnd'], d_CS_lo=cs[0], d_CS_hi=cs[1], d_CRnd_lo=cr[0], d_CRnd_hi=cr[1],
                                  passed=(cs[0] > 0) and (cr[0] < 0.2)))
    for ni, lab in ((0, '0.5%'), (2, '5%')):
        ids_ = [i for i, s in samp.items() if int(s['noise_idx']) == ni]
        level(lab, ids_, lambda i, a: ar.get((i, a)), lambda i: samp[i]['key'], lambda i: int(samp[i]['seed']))
    level('2% (0B)', [i for i in ids], lambda i, a: a0b.get((i, a)), lambda i: s0b[i]['key'], lambda i: int(s0b[i]['seed']))
    wcsv('0C2_noise.csv', t2)

    # ---- 0C-3: independent sampling check
    for src_, recs in (('arms', ar), ('soft', soft)):
        for (sid, arm), r in recs.items():
            if r.get('status') == 'ok' and r.get('unsound', 0):
                unsound += r['unsound']; rows_unsound.append(dict(src=src_, sample_id=sid, arm=arm, unsound=r['unsound'], sampled_min=json.dumps(r.get('sampled_min'))))
    n_models = sum(1 for recs in (ar, soft) for r in recs.values() if r.get('status') == 'ok')
    wcsv('0C3_unsound.csv', rows_unsound)
    gates.append(dict(rule='R-C3', stat='verified laws where a 200,000-point uniform sample finds G < -tol_strong (must be 0)', count=unsound, n=n_models, value=unsound, passed=unsound == 0))
    wcsv('0C_gate.csv', gates)
    for g in gates:
        print(g['rule'], g.get('passed'), round(g['value'], 3) if isinstance(g.get('value'), float) else g.get('value'))


if __name__ == '__main__':
    main()
