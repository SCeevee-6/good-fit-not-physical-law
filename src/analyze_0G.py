"""Protocol 0G tables. python src/analyze_0G.py
Writes results/0G_calib.csv, 0G_zones_calib.csv, 0G_real.csv, 0G_gate.csv. Gate thresholds are compared as exact
fractions. 0F numbers for the paired comparison are read from results/0F_calib.csv, 0F_zones_calib.csv, 0F_real.csv.
"""
import os, json, csv, glob
from fractions import Fraction as Fr
from statistics import median
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = os.path.join(ROOT, 'results'); J = os.path.join(RES, '0G_jobs')
GATE_SEEDS = set(range(100, 200)); CHECK_SEEDS = set(range(100, 120))
if os.environ.get('OG_QUICK') == '1':                 # pipeline test only
    GATE_SEEDS = CHECK_SEEDS = {10, 11}
NEED = Fr(9, 10)


def load(phase):
    ok, bad = [], []
    for f in sorted(glob.glob(os.path.join(J, phase, '*.json'))):
        r = json.load(open(f))
        (bad if 'error' in r else ok).append(r)
    return ok, bad


def write(name, rows):
    if not rows:
        return
    with open(os.path.join(RES, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    cal, cal_bad = load('calib'); real, real_bad = load('real')
    crow, zrow = [], []
    for r in cal:
        for d, R in r['degrees'].items():
            row = dict(key=r['key'], noise_idx=r['noise_idx'], seed=r['seed'], degree=int(d))
            for tag in ('est', 'known'):
                for k in ('eps', 'rank', 'stat', 'covered', 'n_neg', 'n_drop', 'min_reliability', 'max_abs_err'):
                    row[f'{tag}_{k}'] = R.get(tag, {}).get(k, '')
            crow.append(row)
            for law, z in R.get('zones_conf', {}).items():
                chk = R.get('sample_check', {}).get(law, {})
                zrow.append(dict(key=r['key'], seed=r['seed'], degree=int(d), law=law, certified=z['certified'], refuted=z['refuted'],
                                 undecided=z['undecided'], h_inf=z['h_inf'], all_certified=z['all_certified'],
                                 bad_certified=chk.get('bad_certified', ''), bad_refuted=chk.get('bad_refuted', ''),
                                 certificate=R.get('certificates', {}).get(law, {}).get('status', '')))
    write('0G_calib.csv', crow); write('0G_zones_calib.csv', zrow)
    rrow = []
    for r in real:
        for d, R in r['degrees'].items():
            for law, z in R['zones_conf'].items():
                rrow.append(dict(name=r['name'], degree=int(d), law=law, rank=R['rank'], n_neg=R['n_neg'], certified=z['certified'],
                                 refuted=z['refuted'], undecided=z['undecided'], h_inf=z['h_inf'],
                                 certificate=R['certificates'].get(law, '')))
    write('0G_real.csv', rrow)

    def cov(key, d, ni=1, tag='est', src=crow, prefix=''):
        rows = [x for x in src if x['key'] == key and int(x['degree']) == d and int(x['noise_idx']) == ni and int(x['seed']) in GATE_SEEDS]
        c = sum(str(x[f'{tag}_covered']) == 'True' for x in rows)
        ratio = median(float(x[f'{tag}_stat']) / float(x[f'{tag}_rank']) for x in rows) if rows else None
        return c, len(rows), ratio
    f_cal = list(csv.DictReader(open(os.path.join(RES, '0F_calib.csv'))))
    g = []
    c, n, rt = cov('E4', 2)
    g.append(dict(rule='G-G1 E4 degree 2 coverage', value=f'{c}/{n}', threshold='>= 9/10', result='pass' if n and Fr(c, 100) >= NEED else 'fail',
                  note=f'median stat/rank {rt:.3f}' if rt else ''))
    for key, ds in (('A', (2, 3, 4)), ('B', (2, 3))):
        for d in ds:
            c, n, rt = cov(key, d)
            c0, n0, _ = cov(key, d, src=f_cal)
            g.append(dict(rule=f'G-G2 {key} degree {d} coverage', value=f'{c}/{n}', threshold='>= 9/10',
                          result='pass' if n and Fr(c, 100) >= NEED else 'fail', note=f'0F {c0}/{n0}; median stat/rank {rt:.3f}' if rt else ''))
    chk = [x for x in zrow if int(x['seed']) in CHECK_SEEDS and x['bad_certified'] != '']
    bad = sum(int(x['bad_certified']) + int(x['bad_refuted']) for x in chk)
    g.append(dict(rule='G-G3 sampled models contradicting the zone map', value=str(bad), threshold='0',
                  result='pass' if chk and bad == 0 else 'fail', note=f'{len(chk)} system x degree x law checks'))
    c, n, rt = cov('E4', 3); c0, n0, _ = cov('E4', 3, src=f_cal)
    g.append(dict(rule='R-G1 E4 degree 3 coverage (predicted < 9/10)', value=f'{c}/{n}', threshold='report',
                  result='held' if n and Fr(c, 100) < NEED else 'did not hold', note=f'0F {c0}/{n0}; median stat/rank {rt:.3f}' if rt else ''))
    f_z = list(csv.DictReader(open(os.path.join(RES, '0F_zones_calib.csv'))))
    for key in ('E4', 'A', 'B'):
        for d in sorted({int(x['degree']) for x in zrow if x['key'] == key}):
            gz = [float(x['certified']) for x in zrow if x['key'] == key and int(x['degree']) == d and int(x['seed']) in GATE_SEEDS]
            fz = [float(x['certified']) for x in f_z if x['key'] == key and int(x['degree']) == d and x['set'] == 'conf'
                  and int(x['noise_idx']) == 1 and int(x['seed']) in GATE_SEEDS]
            g.append(dict(rule=f'R-G2 {key} degree {d} mean certified share', value=f'{sum(gz) / len(gz):.3f}' if gz else '',
                          threshold='report', result='', note=f'0F {sum(fz) / len(fz):.3f}' if fz else ''))
    for ni in (0, 2):
        for key in ('E4', 'A', 'B'):
            for d in sorted({int(x['degree']) for x in crow if x['key'] == key}):
                c, n, rt = cov(key, d, ni=ni); c0, n0, _ = cov(key, d, ni=ni, src=f_cal)
                g.append(dict(rule=f'R-G3 {key} degree {d} coverage, noise index {ni}', value=f'{c}/{n}', threshold='report', result='',
                              note=f'0F {c0}/{n0}'))
    for key in ('E4', 'A', 'B'):
        for d in sorted({int(x['degree']) for x in crow if x['key'] == key}):
            c, n, _ = cov(key, d, tag='known')
            g.append(dict(rule=f'R-G3 {key} degree {d} coverage, known sigma', value=f'{c}/{n}', threshold='report', result='', note=''))
    f_r = list(csv.DictReader(open(os.path.join(RES, '0F_real.csv'))))
    for d in (2, 3):
        def med_ref(rows, setcol=None):
            names = sorted({x['name'] for x in rows if x['name'].startswith('covid_') and int(x['degree']) == d})
            per = []
            for nm in names:
                v = [float(x['refuted']) for x in rows if x['name'] == nm and int(x['degree']) == d and (setcol is None or x['set'] == setcol)]
                if v:
                    per.append(sum(v) / len(v))
            return (median(per), len(per)) if per else (None, 0)
        m1, k1 = med_ref(rrow); m0, k0 = med_ref(f_r, 'conf')
        g.append(dict(rule=f'R-G4 COVID degree {d} median refuted share (predicted >= 0.10 at degree 2)', value=f'{m1:.3f}' if m1 is not None else '',
                      threshold='report', result=('held' if m1 is not None and m1 >= 0.10 else 'did not hold') if d == 2 else '',
                      note=f'{k1} countries; 0F {m0:.3f}' if m0 is not None else ''))
    g.append(dict(rule='jobs with errors or time limits', value=f'{len(cal_bad)} calib, {len(real_bad)} real', threshold='report', result='', note=''))
    write('0G_gate.csv', g)
    for x in g:
        print(x['rule'], '|', x['value'], '|', x['result'], '|', x['note'])


if __name__ == '__main__':
    main()
