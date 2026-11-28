"""Post-hoc unit R2 tables. python src/analyze_0F3b.py -> results/0F3b_horizon.csv, results/0F3b_gate.csv, figures/fig_0F3b.png."""
import os, json, glob, csv
from statistics import median
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); RES = os.path.join(ROOT, 'results')


def main():
    recs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(RES, '0F3b_jobs', '*.json')))]
    bad = [r for r in recs if 'error' in r]; recs = [r for r in recs if 'error' not in r]
    old = {(x['key'], int(x['seed']), int(x['degree']), x['law']): x for x in csv.DictReader(open(os.path.join(RES, '0F_horizon.csv')))}
    rows = []
    for r in recs:
        for d, laws in r['degrees'].items():
            for law, v in laws.items():
                o = old.get((r['key'], r['seed'], int(d), law))
                rows.append(dict(key=r['key'], seed=r['seed'], degree=int(d), law=law, **v,
                                 old_n_censored_exact=int(o['n_censored_exact']) if o else ''))
    with open(os.path.join(RES, '0F3b_horizon.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    g = []
    viol = sum(x['proj_violations'] for x in rows); pts = sum(x['n_points'] for x in rows)
    g.append(dict(rule='projected bound below the exact h^2 at every band point', value=f'{viol} violations in {pts}', threshold='0', result='pass' if viol == 0 else 'fail'))
    rep = [x for x in rows if x['old_n_censored_exact'] != '']; mism = sum(x['n_censored_exact'] != x['old_n_censored_exact'] for x in rep)
    g.append(dict(rule='0F-3 censoring counts reproduced', value=f'{len(rep) - mism} of {len(rep)} rows equal', threshold='all', result='pass' if rep and mism == 0 else 'fail'))
    for key in sorted({x['key'] for x in rows}):
        for d in sorted({x['degree'] for x in rows if x['key'] == key}):
            sel = [x for x in rows if x['key'] == key and x['degree'] == d]
            rays = sum(x['n_censored_exact'] + x['n_bracket_exact'] for x in sel)
            cx = sum(x['n_censored_exact'] for x in sel) / rays; cp = sum(x['n_censored_proj'] for x in sel) / rays
            gx = [x['median_log_gap_exact'] for x in sel if x['median_log_gap_exact'] is not None]
            gp = [x['median_log_gap_proj'] for x in sel if x['median_log_gap_proj'] is not None]
            g.append(dict(rule=f'{key} degree {d}: censored share exact -> projected (predicted <= half)', value=f'{cx:.3f} -> {cp:.3f}', threshold='cp <= cx/2',
                          result='held' if cp <= cx / 2 else 'did not hold'))
            g.append(dict(rule=f'{key} degree {d}: median log(g_lower / g_true), exact vs projected', value=f"{median(gx):.3f} vs {median(gp):.3f}" if gx and gp else 'n/a',
                          threshold='report', result=''))
    g.append(dict(rule='jobs with errors', value=str(len(bad)), threshold='report', result=''))
    with open(os.path.join(RES, '0F3b_gate.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(g[0])); w.writeheader(); w.writerows(g)
    for x in g:
        print(x['rule'], '|', x['value'], '|', x['result'])


if __name__ == '__main__':
    main()
