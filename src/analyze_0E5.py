"""Protocol 0E5 tables. python src/analyze_0E5.py -> results/0E5_calib.csv, 0E5_zones.csv, 0E5_design.csv, 0E5_gate.csv."""
import os, json, glob, csv
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); RES = os.path.join(ROOT, 'results'); J = os.path.join(RES, '0E5_jobs')
GATE_SEEDS = set(range(100, 200)); CHECK_SEEDS = set(range(100, 120))


def load(phase):
    rs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(J, phase, '*.json')))]
    return [r for r in rs if 'error' not in r], [r for r in rs if 'error' in r]


def write(name, rows):
    if rows:
        with open(os.path.join(RES, name), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def ratio_stats(P, seed=20261126):
    r = np.array(P, float); med = float(np.median(r))
    bs = np.median(r[np.random.default_rng(seed).integers(0, len(r), (10000, len(r)))], axis=1)
    return med, float(np.quantile(bs, 0.025))


def main():
    cal, cbad = load('calib'); des, dbad = load('design'); g = []
    crow, zrow = [], []
    for r in cal:
        for d, R in r['degrees'].items():
            crow.append(dict(seed=r['seed'], degree=int(d), **{f'{t}_{k}': R[t][k] for t in ('est', 'known') if t in R for k in ('stat', 'eps', 'rank', 'covered', 'n_neg', 'min_reliability')}))
            for law, z in R.get('zones_conf', {}).items():
                chk = R.get('sample_check', {}).get(law, {})
                zrow.append(dict(seed=r['seed'], degree=int(d), law=law, certified=z['certified'], refuted=z['refuted'], undecided=z['undecided'],
                                 bad=(chk['bad_certified'] + chk['bad_refuted']) if chk else '', certificate=R.get('certificates', {}).get(law, {}).get('status', '')))
    write('0E5_calib.csv', crow); write('0E5_zones.csv', zrow)
    for d in (2, 3):
        rr = [x for x in crow if x['degree'] == d and x['seed'] in GATE_SEEDS]; c = sum(x['est_covered'] for x in rr); ck = sum(x['known_covered'] for x in rr)
        nn = sum(x['est_n_neg'] > 0 for x in rr)
        g.append(dict(rule=f'{"G-E5a" if d == 2 else "R-E5a"} E5 degree {d} coverage, 0G set (estimated sigma; known sigma)', value=f'{c}/{len(rr)}; {ck}/{len(rr)}',
                      threshold='>= 90/100' if d == 2 else 'report', result=('pass' if Fr(c, 100) >= Fr(9, 10) else 'fail') if d == 2 else '',
                      note=f'seeds with a non-positive reliability: {nn}'))
    chk = [x for x in zrow if x['seed'] in CHECK_SEEDS and x['bad'] != '']; bad = sum(int(x['bad']) for x in chk)
    g.append(dict(rule='G-E5b sampled models contradicting the zone map', value=str(bad), threshold='0', result='pass' if chk and bad == 0 else 'fail', note=f'{len(chk)} checks'))
    cert = [x for x in zrow if x['certificate']]
    g.append(dict(rule='region certificates (5-simplex verifier)', value=', '.join(f"{s}: {sum(x['certificate'] == s for x in cert)}" for s in sorted({x['certificate'] for x in cert})) or 'none attempted',
                  threshold='report', result='', note=''))
    drows, pairs = [], []
    for dj in des:
        for l in dj['laws']:
            nl = dj['arms'][f'law_{l}']['need'][l]; nr = dj['arms']['random']['need'][l]
            drows.append(dict(key='E5', seed=dj['seed'], law=l, n_target=dj['n_target'][l], law_design=nl, random=nr, halton=dj['arms']['halton']['need'][l],
                              dopt=dj['arms']['dopt']['need'][l], start_undecided=dj['arms']['random']['trace'][0][l]))
            if nl > 0:
                pairs.append((nr / nl, nl, dj['arms']['dopt']['need'][l]))
    write('0E5_design.csv', drows)
    e4 = [x for x in csv.DictReader(open(os.path.join(RES, '0F_design.csv'))) if x['key'] == 'E4' and int(x['law_design']) > 0]
    e4r = [int(x['random']) / int(x['law_design']) for x in e4]
    if pairs:
        med, lo = ratio_stats([p[0] for p in pairs]); ml = float(np.median([p[1] for p in pairs])); md = float(np.median([p[2] for p in pairs]))
        m4, lo4 = ratio_stats(e4r)
        g.append(dict(rule='G-E5c design: median random/law (bootstrap low), predicted >= 1.5 with low > 1', value=f'{med:.3f} ({lo:.3f})', threshold='>= 1.5, > 1',
                      result='pass' if med >= 1.5 and lo > 1 else 'fail', note=f'{len(pairs)} pairs; law median {ml} vs D-optimal {md}'))
        g.append(dict(rule='R-E5c E4 (0F-4) for comparison: median random/law (bootstrap low)', value=f'{m4:.3f} ({lo4:.3f})', threshold='report', result='',
                      note=f'{len(e4r)} pairs; E5 minus E4 median ratio {med - m4:+.3f}'))
        for l in ('E2S', 'E2R'):
            rr = [x for x in drows if x['law'] == l]
            if rr:
                g.append(dict(rule=f'E5 {l}: median trajectories law / dopt / random / halton (censored at 31)', threshold='report', result='',
                              value=' / '.join(str(float(np.median([x[a] for x in rr]))) for a in ('law_design', 'dopt', 'random', 'halton')),
                              note=f"start undecided median {float(np.median([x['start_undecided'] for x in rr])):.3f}; censored law/random {sum(x['law_design'] == 31 for x in rr)}/{sum(x['random'] == 31 for x in rr)}"))
    g.append(dict(rule='jobs with errors or time limits', value=f'{len(cbad)} calib, {len(dbad)} design', threshold='report', result='', note=''))
    write('0E5_gate.csv', g)
    for x in g:
        print(x['rule'], '|', x['value'], '|', x['result'], '|', x['note'])


if __name__ == '__main__':
    main()
