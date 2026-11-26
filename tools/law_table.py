"""Law taxonomy table: one row per (system, law). Tabulates results that already exist; runs nothing new.

Sources: results/pilot_candidates.csv (A, A-mis, B pools, noise index 1), results/0B_jobs/pool (H pools),
results/0B_jobs/arms and results/0B_sample.csv (arms S, C, Rnd on A, A-mis, B-pos, B-fac, H; seeds 0-19),
results/0R_jobs/w7_c3.json (R pool and arms, primary config).
Shares are exact counts; classes use exact fractions. Descriptive only, not a gate.
usage (repo root): python tools/law_table.py            -> results/law_taxonomy.csv
"""
import os, sys, csv, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')
csv.field_size_limit(10 ** 9)
LAWS = {'A': ['A1x', 'A1y', 'A2'], 'Amis': ['A1x', 'A1y', 'A2'], 'Bpos': ['B1S', 'B1I', 'B1R', 'B2S', 'B2R', 'B3p', 'B3m'],
        'Bfac': ['B1S', 'B1I', 'B1R', 'B2S', 'B2R', 'B3p', 'B3m'], 'H': ['H1x', 'H1y', 'H2'], 'R': ['R1x', 'R1y', 'R2']}
STRUCT_LAWS = {'A': {'A1x', 'A1y'}, 'Amis': {'A1x', 'A1y'}, 'Bpos': {'B1I', 'B1R', 'B3p', 'B3m'}, 'Bfac': {'B1I', 'B1R', 'B3p', 'B3m'},
               'H': {'H1x', 'H1y'}, 'R': {'R1x', 'R1y'}}
POOL_KEY = {'A': 'A', 'Amis': 'Amis', 'Bpos': 'B', 'Bfac': 'B', 'H': 'H', 'R': 'R'}


def pool_counts():
    """(system key) -> (n accuracy-equivalent candidates, {law: n violating at strong tolerance})."""
    out = {}
    cnt = {}
    for r in csv.DictReader(open(os.path.join(R, 'pilot_candidates.csv'))):
        if r['noise_idx'] != '1' or r['aes'] != 'True':
            continue
        k = r['key']
        n, c = cnt.setdefault(k, [0, {}])
        cnt[k][0] += 1
        for l in filter(None, r['viol_specs_strong'].split(',')):
            c[l] = c.get(l, 0) + 1
    out.update({k: tuple(v) for k, v in cnt.items()})
    n, c = 0, {}
    vols = {}
    for f in glob.glob(os.path.join(R, '0B_jobs', 'pool', 'H_*.json')):
        for r in json.load(open(f))['rows']:
            if r['aes']:
                n += 1
                for l, s in r['status'].items():
                    c[l] = c.get(l, 0) + (s == 'counterexample')
                for l, v in r['volume'].items():
                    if v is not None:
                        vols.setdefault(l, []).append(v)
    out['H'] = (n, c); out['H_vol'] = vols
    J = json.load(open(os.path.join(R, '0R_jobs', 'w7_c3.json')))
    n, c, vols = 0, {}, {}
    for r in J['rows']:
        if r['aes']:
            n += 1
            for l, s in r['status'].items():
                c[l] = c.get(l, 0) + (s == 'counterexample')
            for l, v in r['volume'].items():
                if v is not None:
                    vols.setdefault(l, []).append(v)
    out['R'] = (n, c); out['R_vol'] = vols
    return out, J


def arm_counts(J):
    """(label, arm) -> (n, {law: n verified}). Candidates whose arm run is not ok count as not verified."""
    sample = {r['sample_id']: r for r in csv.DictReader(open(os.path.join(R, '0B_sample.csv')))}
    acc = {}
    for f in glob.glob(os.path.join(R, '0B_jobs', 'arms', '*.json')):
        d = json.load(open(f))
        s = sample[d['sample_id']]
        key, struct = s['key'], d['struct']
        lab = {'A': 'A', 'Amis': 'Amis', 'H': 'H', 'B': struct}[key]
        if key == 'H' and s['noise_idx'] != '1':
            continue
        for rec in d['records']:
            n, c = acc.setdefault((lab, rec['arm']), [0, {}])
            acc[(lab, rec['arm'])][0] += 1
            if rec.get('status') == 'ok':
                for l, v in json.loads(rec['statuses']).items():
                    c[l] = c.get(l, 0) + (v == 'verified')
    by = {r['cand']: r for r in J['rows']}
    for a in J['arms']:
        if a['cand'] not in J['aes'] or a['cand'] not in J['arm_picks']:
            continue
        n, c = acc.setdefault(('R', a['arm']), [0, {}])
        acc[('R', a['arm'])][0] += 1
        if a.get('status') == 'ok':
            for l, v in json.loads(a['statuses']).items():
                c[l] = c.get(l, 0) + (v == 'verified')
    return acc


def klass(pool_v, pool_n, s, c, rnd, n):
    if 20 * pool_v <= pool_n:
        return 'rarely violated, no fix needed'
    if s == n:
        return 'enforced by structure'
    if 10 * c >= 9 * n:
        return 'needs the counterexample loop'
    return 'not fixed by C'


def main():
    pools, J = pool_counts()
    acc = arm_counts(J)
    rows = []
    for lab, laws in LAWS.items():
        pn, pc = pools[POOL_KEY[lab]]
        for l in laws:
            row = dict(system=lab, law=l, family='structural' if l in STRUCT_LAWS[lab] else 'non-structural', n_pool_aes=pn,
                       pool_violating=pc.get(l, 0), pool_share=pc.get(l, 0) / max(1, pn))
            for arm in ('S', 'C', 'Rnd'):
                n, c = acc[(lab, arm)]
                row[f'n_{arm}'] = n; row[f'{arm}_verified'] = c.get(l, 0); row[f'{arm}_share'] = c.get(l, 0) / max(1, n)
            row['C_minus_Rnd'] = row['C_verified'] - row['Rnd_verified']
            row['pool_ok_share'] = 1 - row['pool_share']                    # share of equivalent pool candidates that satisfy the law
            row['S_minus_pool_ok'] = row['S_share'] - row['pool_ok_share']  # negative: structure alone made the law fail more often
            row['class'] = klass(row['pool_violating'], pn, row['S_verified'], row['C_verified'], row['Rnd_verified'], row['n_C'])
            vols = pools.get(f'{lab}_vol', {}).get(l)
            row['median_violation_volume'] = sorted(vols)[len(vols) // 2] if vols else ''
            rows.append(row)
    with open(os.path.join(R, 'law_taxonomy.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(len(rows), 'rows')


if __name__ == '__main__':
    main()
