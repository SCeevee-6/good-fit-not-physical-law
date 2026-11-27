"""F0: law zone maps for stored fits (Rashomon set on the training data, no gates).

Writes results/F0_zones.csv (one row per dataset x degree x law), results/F0_witness_zones.csv,
results/F0_maps.npz (grids for the figure). Usage: python src/run_F0.py [n_seeds]
"""
import os, sys, csv, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from systems import SYSTEMS
from data import make_cell, stack
from gate_a import NOISE, library
import dimtest                                   # registers SYSTEMS['E4']
import real_r, real_e
from lawmap import DataSet, zone, member_ratio, law_grid, DELTA_RSS

RES = os.path.join(ROOT, 'results')
NEAR = 0.25                                      # "near data": within 0.25 state sd of a training point


def near_mask(X, Xtr, sd, chunk=2000):
    out = np.zeros(len(X), bool); A = Xtr / sd
    for i in range(0, len(X), chunk):
        B = X[i:i + chunk] / sd
        d2 = (B ** 2).sum(1)[:, None] + (A ** 2).sum(1)[None, :] - 2 * B @ A.T
        out[i:i + chunk] = d2.min(1) <= NEAR ** 2
    return out


def datasets(n_seeds):
    for key in ('A', 'Amis', 'B', 'H', 'E4'):
        S = SYSTEMS['A' if key == 'Amis' else key]
        for seed in range(n_seeds):
            cell = make_cell(S, NOISE[1], seed, 1)
            sub = '0D_jobs' if key == 'E4' else '0B_jobs'
            pool = json.load(open(os.path.join(RES, sub, 'pool', f'{key}_n1_s{seed}.json')))
            yield dict(group='synthetic', name=f'{key}_s{seed}', key=key, S=S, cell=cell, rows=pool['rows'], primary=None)
    S, cell, info = real_r.build(real_r.PRIMARY)
    j = json.load(open(os.path.join(RES, '0R_jobs', 'w7_c3.json')))
    yield dict(group='real', name='hare_lynx', key='R', S=S, cell=cell, rows=j['rows'], primary=j['primary'])
    for p in real_e.PANEL:
        S, cell, info = real_e.build(dict(slug=p['slug'], window=real_e.PRIMARY_W))
        j = json.load(open(os.path.join(RES, '0E_jobs', f"{p['slug']}_w{real_e.PRIMARY_W}.json")))
        yield dict(group='real', name=f"covid_{p['slug']}", key='B', S=S, cell=cell, rows=j['rows'], primary=j['primary'])


def run_dataset(D, maps):
    S, cell = D['S'], D['cell']
    Xtr, dXtr = stack(cell['train'])
    out, sets = [], {}
    for d in sorted({r['degree'] for r in D['rows']}):
        Efull = S.E(d)
        Ekeep = {tuple(e) for e in library(S, D['key'], d).tolist()}
        keep = np.array([tuple(e) in Ekeep for e in Efull.tolist()])
        ds = DataSet(Xtr, dXtr, cell['deriv_sd'], Efull, keep); sets[d] = ds
        mem = [r for r in D['rows'] if r['degree'] == d and r['aes']]
        inE = [r for r in mem if ds.r(np.array(r['xi'])) <= ds.tau]
        for sp in S.specs:
            X = law_grid(sp)
            m, h = ds.law(sp.c, X, sp.tol_strong)
            z = zone(m, h); nr = near_mask(X, Xtr, cell['state_sd'])
            ratio = max([member_ratio(ds, np.array(r['xi']), sp.c, X) for r in inE], default=float('nan'))
            region_cert = bool((z == 1).all())
            contra = sum(1 for r in inE if region_cert and r.get('status', {}).get(sp.name) == 'counterexample')
            hfin = h[np.isfinite(h)]
            out.append(dict(group=D['group'], dataset=D['name'], key=D['key'], degree=d, law=sp.name, n_grid=len(X),
                            certified=float((z == 1).mean()), refuted=float((z == -1).mean()), undecided=float((z == 0).mean()),
                            n_near=int(nr.sum()), certified_near=float((z[nr] == 1).mean()) if nr.any() else '',
                            undecided_near=float((z[nr] == 0).mean()) if nr.any() else '',
                            h_inf_share=float(np.isinf(h).mean()), h_median=float(np.median(hfin)) if hfin.size else '',
                            rank=ds.rank, p=ds.p, cond=ds.cond, eps=ds.eps, n_train=len(Xtr),
                            n_aes=len(mem), n_aes_in_E=len(inE), max_member_ratio=ratio,
                            region_certified=region_cert, contradictions=contra))
            if D['name'] in ('A_s0', 'hare_lynx', 'H_s0'):
                maps[f"{D['name']}|{d}|{sp.name}"] = dict(X=X, m=m, h=h, z=z)
        if D['name'] in ('A_s0', 'hare_lynx', 'H_s0'):
            maps[f"{D['name']}|{d}|train"] = dict(X=Xtr)
            ext = 1.05 * max(b for sp in S.specs for _, b in sp.u_box)
            g = np.linspace(0, ext, 241)
            GX, GY = np.meshgrid(g, g); P = np.c_[GX.ravel(), GY.ravel()]
            maps[f"{D['name']}|{d}|plane"] = dict(X=P, hx=np.sqrt(ds.eps * cell['deriv_sd'][0] ** 2 * ds.K(P)))
    return out, sets


def witness_rows(D, sets, wit):
    out = []
    for w in wit:
        r = D['rows'][int(w['cand'])]; ds = sets[r['degree']]
        sp = [s for s in D['S'].specs if s.name == w['spec']][0]
        x = np.array(json.loads(w['state']), float)[None, :]
        m, h = ds.law(sp.c, x, sp.tol_strong)
        out.append(dict(dataset=D['name'], cand=w['cand'], role=w['role'], law=w['spec'], degree=r['degree'],
                        state=w['state'], cand_in_E=ds.r(np.array(r['xi'])) <= ds.tau,
                        m=float(m[0]), h=float(h[0]), zone=int(zone(m, h)[0])))
    return out


def write(path, rows):
    with open(path, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)


def main(n_seeds=5):
    t0 = time.time(); allrows, wrows, maps = [], [], {}
    wR = list(csv.DictReader(open(os.path.join(RES, '0R_witnesses.csv'))))
    wE = list(csv.DictReader(open(os.path.join(RES, '0E_witnesses.csv'))))
    for D in datasets(n_seeds):
        rows, sets = run_dataset(D, maps); allrows += rows
        if D['name'] == 'hare_lynx':
            wrows += witness_rows(D, sets, wR)
        elif D['name'].startswith('covid_'):
            wrows += witness_rows(D, sets, [w for w in wE if f"covid_{w['country']}" == D['name']])
        print(f"{D['name']}: {len(rows)} rows, {time.time() - t0:.0f}s", flush=True)
    write(os.path.join(RES, 'F0_zones.csv'), allrows)
    write(os.path.join(RES, 'F0_witness_zones.csv'), wrows)
    np.savez_compressed(os.path.join(RES, 'F0_maps.npz'),
                        **{f'{k}|{f}': v for k, dct in maps.items() for f, v in dct.items()})
    print('delta_rss', DELTA_RSS, 'rows', len(allrows), 'witnesses', len(wrows), f'{time.time() - t0:.0f}s')


if __name__ == '__main__':
    main(*[int(a) for a in sys.argv[1:]])
