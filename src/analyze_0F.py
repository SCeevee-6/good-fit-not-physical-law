"""Tables and gates for 0F-1 and 0F-2. Writes results/0F_gate.csv, 0F_calib.csv, 0F_zones_calib.csv, 0F_real.csv,
0F_witnesses.csv and prints the gate lines. Gate comparisons use exact fractions of counts."""
import os, sys, csv, json, glob
from fractions import Fraction as Fr
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(ROOT, 'results')
J = os.path.join(RES, '0F_jobs')


def write(name, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(RES, name), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)


def main():
    C = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(J, 'calib', '*.json')))]
    errs = [c for c in C if 'error' in c]; C = [c for c in C if 'error' not in c]
    cal, zc, gate = [], [], []
    cov = defaultdict(list); covk = defaultdict(list); bad = 0; n_chk = 0; certs = defaultdict(int)
    for c in C:
        for d, R in c['degrees'].items():
            row = dict(key=c['key'], noise_idx=c['noise_idx'], seed=c['seed'], degree=int(d), **{f'est_{k}': v for k, v in R['est'].items()})
            if 'known' in R:
                row.update({f'known_{k}': v for k, v in R['known'].items()}); covk[(c['key'], int(d), c['noise_idx'])].append(R['known']['covered'])
            cal.append(row); cov[(c['key'], int(d), c['noise_idx'])].append(R['est']['covered'])
            for setname in ('zones_conf', 'zones_rashomon'):
                for law, z in R[setname].items():
                    zc.append(dict(key=c['key'], noise_idx=c['noise_idx'], seed=c['seed'], degree=int(d), set=setname[6:], law=law, **z))
            for law, v in R.get('sample_check', {}).items():
                bad += v['bad_certified'] + v['bad_refuted']; n_chk += 1
            for law, v in R.get('certificates', {}).items():
                certs[v['status']] += 1
    gate.append(dict(rule='G-F0 sampled-model contradictions', value=bad, n=n_chk, threshold='0', result='pass' if bad == 0 and n_chk > 0 else 'FAIL'))
    gate.append(dict(rule='G-F0 region certificates (reported)', value=json.dumps(dict(certs)), n=sum(certs.values()), threshold='-', result='reported'))
    for (k, d, ni), v in sorted(cov.items()):
        if ni != 1:
            continue
        f = Fr(sum(v), len(v)); ok = len(v) >= 90 and f >= Fr(9, 10)
        gate.append(dict(rule=f'G-F1 coverage {k} degree {d}', value=f'{sum(v)}/{len(v)}', n=len(v), threshold='>= 9/10, n >= 90',
                         result='pass' if ok else ('FAIL' if len(v) >= 90 else 'too few seeds')))
    for (k, d, ni), v in sorted(cov.items()):
        if ni != 1:
            gate.append(dict(rule=f'coverage (reported) {k} degree {d} noise {ni}', value=f'{sum(v)}/{len(v)}', n=len(v), threshold='-', result='reported'))
    for (k, d, ni), v in sorted(covk.items()):
        gate.append(dict(rule=f'coverage known sigma (reported) {k} degree {d}', value=f'{sum(v)}/{len(v)}', n=len(v), threshold='-', result='reported'))
    gate.append(dict(rule='calib jobs with errors or time limit', value=len(errs), n=len(errs) + len(C), threshold='-', result='reported'))
    # 0F-2
    Rr = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(J, 'real', '*.json')))]
    real, wit = [], []
    low = {}
    for r in Rr:
        if 'error' in r:
            gate.append(dict(rule=f'real job error {r.get("job")}', value=r['error'][:80], n=1, threshold='-', result='error')); continue
        dmin = min(int(d) for d in r['degrees'])
        for d, R in r['degrees'].items():
            for setname in ('zones_conf', 'zones_rashomon'):
                for law, z in R[setname].items():
                    nr = R['near'].get(law, {})
                    real.append(dict(name=r['name'], degree=int(d), set=setname[6:], law=law, rank=R['rank'], p=R['p'], **z,
                                     n_near=nr.get('n_near'), certified_near=nr.get('certified_near') if setname == 'zones_conf' else None,
                                     certificate=R['certificates'].get(law) if setname == 'zones_conf' else None))
            if int(d) == dmin:
                low[r['name']] = float(np.mean([z['undecided'] for z in R['zones_conf'].values()]))
        for w in r['witnesses']:
            wit.append(dict(name=r['name'], **w))
    if low:
        hl = low.get('hare_lynx'); cv = [v for k, v in low.items() if k.startswith('covid_')]
        gate.append(dict(rule='G-F2 hare-lynx mean undecided share, lowest degree (prediction >= 1/2)', value=None if hl is None else round(hl, 4), n=1, threshold='>= 0.5', result='reported'))
        if cv:
            gate.append(dict(rule='G-F2 COVID median over countries (prediction >= 1/2)', value=round(float(np.median(cv)), 4), n=len(cv), threshold='>= 0.5', result='reported'))
    write('0F_calib.csv', cal); write('0F_zones_calib.csv', zc); write('0F_real.csv', real); write('0F_witnesses.csv', wit); write('0F_gate.csv', gate)
    for g in gate:
        print(g['rule'], '|', g['value'], '|', g['result'])


if __name__ == '__main__':
    main()
