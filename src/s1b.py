"""Protocol S1b (docs/protocol_S1b.md): zones of the energy law over the data-consistent set of linear closures.
python src/s1b.py  (reads results/S1_samples.npz) -> results/S1b_zones.csv, results/S1b_gate.csv, figures/fig_S1b.png

Closure regression b_t = C a_t + e_t (rows are snapshots). The law at resolved state a: the closure dissipates,
a' C a <= 0. It is linear in vec(C) with weights a (x) a (row-major). Covariance of vec(C_hat): Newey-West sandwich
with Bartlett weights over the snapshot sequence, V = (I (x) G^-1) S (I (x) G^-1), G = A' A, scores s_t = e_t (x) a_t.
Set: (c - c_hat)' V^+ (c - c_hat) <= chi2_{rank, 0.95}. Zones at a: certified (all closures dissipate) if m + h <= 0,
refuted (all inject) if m - h > 0, undecided otherwise; m = a' C_hat a, h = sqrt(eps (a (x) a)' V (a (x) a)).
Region certificate for "every closure dissipates at every a", exact in rationals: lambda_max(sym C_hat) + sqrt(eps lambda_max(V)) <= 0.
"""
import os, sys, csv
import numpy as np
from fractions import Fraction
from scipy.stats import chi2
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
KS = (1, 2, 3); LAG = 50; LAG_SENS = (20, 100); N_DIR = 20000; N_SAMPLE = 200; SEED = 11


def nw_cov(A, B, lag):
    T, m = A.shape; G = A.T @ A; Gi = np.linalg.inv(G)
    Ct = Gi @ (A.T @ B); E = B - A @ Ct                         # Ct = C'
    S_t = np.einsum('ti,tj->tij', E, A).reshape(T, m * m)        # e_t (x) a_t, row i of C first
    S = S_t.T @ S_t
    for l in range(1, lag + 1):
        w = 1 - l / (lag + 1); X = S_t[l:].T @ S_t[:-l]; S += w * (X + X.T)
    P = np.kron(np.eye(m), Gi)
    return Ct.T, P @ S @ P.T


def zones(C, V, eps, X, batch=1000):
    m = np.einsum('ni,ij,nj->n', X, C, X); q = np.empty(len(X))
    for s in range(0, len(X), batch):
        W = np.einsum('ni,nj->nij', X[s:s + batch], X[s:s + batch]).reshape(-1, X.shape[1] ** 2)
        q[s:s + batch] = np.einsum('nk,nk->n', W @ V, W)
    h = np.sqrt(eps * np.maximum(q, 0))
    return np.where(m + h <= 0, 1, np.where(m - h > 0, -1, 0)), m, h


def main():
    d = np.load(os.path.join(ROOT, 'results', 'S1_samples.npz')); rng = np.random.default_rng(SEED)
    rows, gate = [], []
    for K in KS:
        A, B = d[f'A{K}'], d[f'B{K}']; mdim = A.shape[1]
        dirs = rng.normal(size=(N_DIR, mdim)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        data_dirs = A / np.linalg.norm(A, axis=1)[:, None]
        for lag in (LAG,) + LAG_SENS:
            C, V = nw_cov(A, B, lag); V = (V + V.T) / 2
            w, Q = np.linalg.eigh(V); keep = w > w.max() * 1e-12; rank = int(keep.sum()); eps = float(chi2.ppf(0.95, rank))
            for name, X in (('sphere', dirs), ('data', data_dirs)):
                z, m, h = zones(C, V, eps, X)
                rows.append(dict(K=K, lag=lag, where=name, n=len(X), certified=float((z == 1).mean()), refuted=float((z == -1).mean()),
                                 undecided=float((z == 0).mean()), rank=rank, eps=eps))
            if lag == LAG:
                L = Q[:, keep] * np.sqrt(w[keep]); Zs = rng.normal(size=(N_SAMPLE, rank)); Zs *= np.sqrt(eps) / np.linalg.norm(Zs, axis=1)[:, None]
                bad = 0
                for X in (dirs[:2000], data_dirs):
                    z, _, _ = zones(C, V, eps, X)
                    Cs = C[None] + (Zs @ L.T).reshape(N_SAMPLE, mdim, mdim)
                    g = np.einsum('ni,sij,nj->sn', X, Cs, X)
                    bad += int(((g > 1e-12) & (z == 1)[None]).sum() + ((g <= -1e-12) & (z == -1)[None]).sum())
                lam_c = float(np.linalg.eigvalsh((C + C.T) / 2).max()); bound = lam_c + float(np.sqrt(eps * w.max()))
                exact = Fraction(lam_c) + Fraction(float(np.sqrt(eps * w.max()))) <= 0
                gate.append(dict(rule=f'K={K}: sampled closures contradicting the zone map (200 boundary closures)', value=str(bad), threshold='0', result='pass' if bad == 0 else 'fail'))
                gate.append(dict(rule=f'K={K}: region certificate, every closure dissipates everywhere', value=f'lambda_max(sym C_hat) {lam_c:.3f}; bound {bound:.3f}', threshold='bound <= 0',
                                 result='verified' if exact else ('impossible (centre already indefinite)' if lam_c > 0 else 'not verified')))
    with open(os.path.join(ROOT, 'results', 'S1b_zones.csv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
    for r in rows:
        if r['lag'] == LAG:
            pred = (r['refuted'] >= 0.5) if r['where'] == 'data' else (r['undecided'] >= 0.5)
            gate.append(dict(rule=f"K={r['K']}, {r['where']}: zone shares (predicted: {'refuted >= 1/2' if r['where'] == 'data' else 'undecided >= 1/2'})",
                             value=f"certified {r['certified']:.3f}, refuted {r['refuted']:.3f}, undecided {r['undecided']:.3f}", threshold='report',
                             result='held' if pred else 'did not hold'))
    with open(os.path.join(ROOT, 'results', 'S1b_gate.csv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(gate[0])); wr.writeheader(); wr.writerows(gate)
    for x in gate:
        print(x['rule'], '|', x['value'], '|', x['result'])
    figure(rows)
    d.close()


def figure(rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.3), constrained_layout=True, sharey=True)
    cols = {'certified': '#2c6fb3', 'refuted': '#e08214', 'undecided': '0.75'}
    for a, where, title in ((ax[0], 'data', 'At observed resolved states'), (ax[1], 'sphere', 'At random directions on the energy shell')):
        sel = [r for r in rows if r['where'] == where and r['lag'] == LAG]; x = np.arange(len(sel)); bottom = np.zeros(len(sel))
        for zname in ('certified', 'refuted', 'undecided'):
            v = np.array([r[zname] for r in sel]); a.bar(x, v, 0.6, bottom=bottom, color=cols[zname], label=zname if where == 'data' else None); bottom += v
        a.set_xticks(x); a.set_xticklabels([f"K = {r['K']}" for r in sel]); a.set_title(title, loc='left')
    ax[0].set_ylabel('share of states'); ax[0].legend(frameon=False, loc='lower left', fontsize=6)
    for a, l in zip(ax, 'ab'):
        a.text(-0.18, 1.08, l, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom')
    fig.savefig(os.path.join(ROOT, 'figures', 'fig_S1b.png'), dpi=300)


if __name__ == '__main__':
    main()
