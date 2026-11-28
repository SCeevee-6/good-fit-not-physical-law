"""Protocol 0H tables. python src/analyze_0H.py -> results/0H_calib.csv, results/0H_gate.csv (thresholds as exact fractions)."""
import os, json, glob, csv
from fractions import Fraction as Fr
from statistics import median
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); RES = os.path.join(ROOT, 'results')
ARMS = ('base', 'oracle', 'known', 'w15', 'w31')


def main():
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RES, '0H_jobs', '*.json')))]
    bad = [r for r in recs if 'error' in r]; recs = [r for r in recs if 'error' not in r]
    rows = []
    for r in recs:
        for d, R in r['degrees'].items():
            row = dict(noise_idx=r['noise_idx'], seed=r['seed'], degree=int(d), nu=R['base']['bias_noncentrality'], rank=R['base']['rank'])
            for a in ARMS:
                row[f'{a}_stat'] = R[a]['stat']; row[f'{a}_eps'] = R[a]['eps']; row[f'{a}_covered'] = R[a]['covered']
            rows.append(row)
    with open(os.path.join(RES, '0H_calib.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    sel = lambda ni, d: [x for x in rows if x['noise_idx'] == ni and x['degree'] == d]
    cov = lambda ni, d, a: sum(x[f'{a}_covered'] for x in sel(ni, d))
    g = []
    o = [cov(0, d, 'oracle') for d in (2, 3, 4)]
    g.append(dict(rule='G-H1 oracle coverage at 0.5 %, degrees 2/3/4', value='/'.join(map(str, o)), threshold='>= 90 each',
                  result='pass' if all(Fr(c, 100) >= Fr(9, 10) for c in o) else 'fail'))
    b = [cov(0, d, 'base') for d in (2, 3, 4)]
    g.append(dict(rule='G-H2 base coverage at 0.5 %, degrees 2/3/4', value='/'.join(map(str, b)), threshold='< 90 at one or more',
                  result='pass' if any(Fr(c, 100) < Fr(9, 10) for c in b) else 'fail'))
    m0 = [median(x['nu'] for x in sel(0, d)) for d in (2, 3, 4)]; m1 = [median(x['nu'] for x in sel(1, d)) for d in (2, 3, 4)]
    g.append(dict(rule='G-H3 median nu at 0.5 % / 2 %, degrees 2/3/4', value='/'.join(f'{a:.2f}:{c:.2f}' for a, c in zip(m0, m1)),
                  threshold='ratio >= 4 each', result='pass' if all(c > 0 and a >= 4 * c for a, c in zip(m0, m1)) else 'fail'))
    for ni in (0, 1):
        for a in ARMS:
            g.append(dict(rule=f'coverage {a}, noise index {ni}, degrees 2/3/4', value='/'.join(str(cov(ni, d, a)) for d in (2, 3, 4)),
                          threshold='report', result=''))
    g.append(dict(rule='jobs with errors', value=str(len(bad)), threshold='report', result=''))
    with open(os.path.join(RES, '0H_gate.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(g[0])); w.writeheader(); w.writerows(g)
    for x in g:
        print(x['rule'], '|', x['value'], '|', x['result'])


if __name__ == '__main__':
    main()
