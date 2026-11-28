"""Development diagnostic for protocol 0G: coverage over noise redraws on one fixed set of clean trajectories.
Not a gate. python tools/dev_0G_mc.py <system> <seed> <degree> <reps> <kappa> [noise_idx]
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(HERE), 'src'))


def main(key, seed, d, reps, kappa, ni=1):
    import dimtest  # noqa: F401
    from systems import SYSTEMS
    from data import make_cell, smooth
    from gate_a import NOISE
    from lawmap import noise_sd
    from eiv import CorrectedSet
    S = SYSTEMS[key]; seed, d, reps, ni, kappa = int(seed), int(d), int(reps), int(ni), float(kappa)
    cell = make_cell(S, NOISE[ni], seed, ni); trajs = cell['train']; sig = np.asarray(cell['sigma'], float)
    E = S.E(d); xt = S.true_Xi(d); rng = np.random.default_rng([seed, 777]); ratio, drop = [], []
    for r in range(reps):
        new = []
        for tr in trajs:
            X = tr['X_clean'] + rng.normal(size=tr['X_clean'].shape) * sig
            Xs, dX, sl = smooth(X, S); new.append(dict(X=X, Xs=Xs, dX=dX))
        cs = CorrectedSet(new, S, E, noise_sd(new, S), kappa=kappa)
        ratio.append(cs.stat(xt) / cs.eps); drop.append(cs.n_drop)
    ratio = np.array(ratio)
    print(f"{key} seed {seed} d {d} kappa {kappa}: coverage {np.mean(ratio <= 1):.3f} over {reps}, median stat/eps {np.median(ratio):.2f}, "
          f"90th pct {np.quantile(ratio, 0.9):.2f}, dropped directions per draw {np.mean(drop):.2f} (min {min(drop)}, max {max(drop)})")


if __name__ == '__main__':
    main(*sys.argv[1:])
