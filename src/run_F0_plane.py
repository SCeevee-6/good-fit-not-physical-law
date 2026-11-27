"""F0, COVID on the plane S + I + R = 1.

The JHU fractions satisfy S + I + R = 1 exactly, so the data cannot identify any polynomial that is a multiple
of (1 - S - I - R). Off that plane every law is undecided for that reason alone (F0_zones.csv, h_inf_share = 1).
Here the law domains are cut down to the plane: face laws on the segment of the face that lies on it, region
laws on the triangle. Writes results/F0_covid_plane.csv. Usage: python src/run_F0_plane.py
"""
import os, sys, csv, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from systems import K_FLOAT, in_domain
from data import stack
from gate_a import library
import real_e
from lawmap import DataSet, zone
from run_F0 import near_mask, write

RES = os.path.join(ROOT, 'results')


def plane_points(spec, n1=2001, n2=201):
    if len(spec.u_names) == 2:                       # face law: u = (a, b) with a + b = 1 on the plane
        t = np.linspace(0, 1, n1); U = [t, 1 - t]
    else:                                            # region law: u = (S, I, R), R = 1 - S - I
        g = np.linspace(0, 1, n2); S, I = [a.ravel() for a in np.meshgrid(g, g, indexing='ij')]
        ok = S + I <= 1 + 1e-12; U = [S[ok], I[ok], np.clip(1 - S[ok] - I[ok], 0, None)]
    ok = in_domain(spec, U); U = [u[ok] for u in U]
    pt = spec.point(U, K_FLOAT); m = U[0].shape[0]
    return np.stack([np.broadcast_to(np.asarray(p, float), (m,)) for p in pt], 1)


def main():
    rows = []
    for p in real_e.PANEL:
        S, cell, info = real_e.build(dict(slug=p['slug'], window=real_e.PRIMARY_W))
        j = json.load(open(os.path.join(RES, '0E_jobs', f"{p['slug']}_w{real_e.PRIMARY_W}.json")))
        Xtr, dXtr = stack(cell['train'])
        for d in sorted({r['degree'] for r in j['rows']}):
            Efull = S.E(d); keep = np.ones(Efull.shape[0], bool)
            ds = DataSet(Xtr, dXtr, cell['deriv_sd'], Efull, keep)
            for sp in S.specs:
                X = plane_points(sp)
                assert np.abs(X.sum(1) - 1).max() < 1e-9
                m, h = ds.law(sp.c, X, sp.tol_strong); z = zone(m, h)
                nr = near_mask(X, Xtr, cell['state_sd'])
                rows.append(dict(country=p['slug'], degree=d, law=sp.name, n_grid=len(X),
                                 certified=float((z == 1).mean()), refuted=float((z == -1).mean()),
                                 undecided=float((z == 0).mean()), h_inf_share=float(np.isinf(h).mean()),
                                 n_near=int(nr.sum()), certified_near=float((z[nr] == 1).mean()) if nr.any() else '',
                                 rank=ds.rank, p=ds.p, I_max=info['I_max'], S_min=info['S_min']))
        print(p['slug'], flush=True)
    write(os.path.join(RES, 'F0_covid_plane.csv'), rows)
    print('rows', len(rows))


if __name__ == '__main__':
    main()
