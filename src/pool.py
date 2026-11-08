"""Candidate pool: STLSQ sweeps over degree x threshold x bootstrap, with the conventional validation battery."""
import time
import numpy as np
import pandas as pd
from lib import theta, support
from data import stack
from sindy import stlsq, deriv_rmse, rollout_rmse, aic_score

THRESHOLDS = np.logspace(-3, np.log10(0.3), 12)
N_BOOT = 8


def make_variants(ntr, seed, n_boot=N_BOOT):
    rng = np.random.default_rng([seed, 777])
    return [('full', list(range(ntr)))] + [(f'boot{b}', list(rng.integers(0, ntr, ntr))) for b in range(n_boot)]


def build_pool(S, cell, seed, n_boot=N_BOOT, thresholds=THRESHOLDS, drop=None, degrees=None):
    train = cell['train']
    ntr = len(train)
    variants = make_variants(ntr, seed, n_boot)
    rows, seen = [], set()
    drop = [tuple(x) for x in (drop or [])]
    for d in (degrees or S.degrees):
        E = S.E(d)
        Ttrue = S.true_Xi(d)
        keep = np.array([tuple(e) not in drop for e in E.tolist()])
        for vname, idx in variants:
            Xs, dX = stack(train, idx)
            Th = theta(Xs, E)
            for th in thresholds:
                if keep.all():
                    Xi = stlsq(Th, dX, th)
                else:
                    Xi = np.zeros((E.shape[0], dX.shape[1]))
                    Xi[keep] = stlsq(Th[:, keep], dX, th)
                key = (d, Xi.tobytes())
                if key in seen:
                    continue
                seen.add(key)
                iid_rmse, iid_r2 = deriv_rmse(Xi, E, cell['iid'], cell['deriv_sd'])
                tr_rmse, _ = deriv_rmse(Xi, E, train, cell['deriv_sd'])
                roll_iid = rollout_rmse(Xi, E, cell['iid'], cell['state_sd'])
                roll_ood = rollout_rmse(Xi, E, cell['ood'], cell['state_sd'])
                aic = aic_score(Xi, E, train, cell['deriv_sd'])
                sup, sup_t = Xi != 0, Ttrue != 0
                tp = int((sup & sup_t).sum()); fp = int((sup & ~sup_t).sum()); fn = int((~sup & sup_t).sum())
                f1 = 2 * tp / max(2 * tp + fp + fn, 1)
                rows.append(dict(degree=d, variant=vname, threshold=float(th), k=int(sup.sum()),
                                 iid_rmse=iid_rmse, iid_r2=iid_r2, train_rmse=tr_rmse, roll_iid=roll_iid, roll_ood=roll_ood, aic=aic,
                                 support_f1=f1, exact_support=bool((sup == sup_t).all()),
                                 coef_err=float(np.linalg.norm(Xi - Ttrue) / np.linalg.norm(Ttrue)),
                                 arm='misspec' if drop else 'wellspec', xi=Xi.ravel().tolist(), xi_shape=Xi.shape))
    df = pd.DataFrame(rows)
    best = df['iid_rmse'].min()
    df['best_iid_rmse'] = best
    df['aes'] = df['iid_rmse'] <= 1.10 * best
    return df
