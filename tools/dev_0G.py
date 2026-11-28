"""Development runs for protocol 0G (seeds outside every gate sample). Not a gate.

python tools/dev_0G.py <systems comma> <seed_lo> <seed_hi> <noise_idx> [out_csv] [kappa]
Per (system, seed, degree): stat/q of the 0F stacked set and of the corrected set (estimated and true sigma),
coverage, eigen diagnostics of the corrected Gram matrix, and two sanity numbers (Xs = H X, sum H*D = 0).
"""
import os, sys, csv, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'src'))


def main(systems, lo, hi, ni, out=None, kappa='0.5'):
    kappa = float(kappa)
    import dimtest  # noqa: F401
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from lawmap import stacked_set, noise_sd, sg_mats
    from eiv import CorrectedSet
    rows = []
    for key in systems.split(','):
        S = SYSTEMS[key]
        for seed in range(int(lo), int(hi) + 1):
            t0 = time.time()
            cell = make_cell(S, NOISE[int(ni)], seed, int(ni)); trajs = cell['train']
            sh = noise_sd(trajs, S); st = np.asarray(cell['sigma'], float)
            X = np.asarray(trajs[0]['X'], float); H, D, sl = sg_mats(len(X), S)
            xs_err = float(np.abs(H @ X - trajs[0]['Xs']).max()); hd = float(np.abs((H * D).sum(1)).max())
            for d in S.degrees:
                E = S.E(d); xt = S.true_Xi(d)
                ds0 = stacked_set(trajs, S, E, sh)
                r = dict(sys=key, seed=seed, d=d, xs_err=xs_err, hd=hd,
                         q0=ds0.rank, ratio_0F=(ds0.r(xt) - ds0.r_hat) / ds0.eps)
                for tag, sg in (('est', sh), ('true', st)):
                    cs = CorrectedSet(trajs, S, E, sg, kappa=kappa)
                    r[f'ratio_G_{tag}'] = cs.stat(xt) / cs.eps
                    if tag == 'est':
                        r.update(qG=cs.rank, n_neg=cs.n_neg, n_drop=cs.n_drop, min_eig=cs.min_eig_ratio,
                                 err_0F=float(np.abs(ds0.Xi_hat - xt).max()), err_G=float(np.abs(cs.Xi_hat - xt).max()))
                r['sec'] = round(time.time() - t0, 1)
                rows.append(r)
                print(key, seed, d, f"0F {r['ratio_0F']:.2f}  G {r['ratio_G_est']:.2f} / true-sigma {r['ratio_G_true']:.2f}  neg {r['n_neg']}", flush=True)
    if out:
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


if __name__ == '__main__':
    main(*sys.argv[1:])
