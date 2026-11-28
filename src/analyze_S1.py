"""Protocol S1 analysis: closure fit per truncation, exact energy-stability test, gates, figure.
python src/analyze_S1.py  (reads results/S1_run.json and S1_samples.npz)
Writes results/S1_closure.csv, results/S1_gate.csv, figures/fig_S1.png.
"""
import os, sys, json, csv
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import ns_galerkin as G       # noqa: E402
from ns_dns import Grid       # noqa: E402

KS = (1, 2, 3)


def main():
    rec = json.load(open(os.path.join(ROOT, 'results', 'S1_run.json')))
    d = np.load(os.path.join(ROOT, 'results', 'S1_samples.npz')); g = Grid(rec['N'], rec['nu'])
    rows, gate = [], []
    for K in KS:
        A, B = d[f'A{K}'], d[f'B{K}']
        C, r2 = G.fit_closure(A, B); Q = (C + C.T) / 2
        ev = np.linalg.eigvalsh(Q); ok_exact, step, piv = G.nsd_exact(Q); ok_eig = bool(ev.max() <= 0)
        idx = G.half_modes(g, K); k2 = np.concatenate([g.K2[idx], g.K2[idx]])
        true_rate = np.einsum('ti,ti->t', A, B); model_rate = np.einsum('ti,ij,tj->t', A, C, A)
        visc = -2 * rec['nu'] * np.einsum('ti,i,ti->t', A, k2, A) / 2       # d/dt |a|^2/2 from viscosity: -nu sum k^2 a_i^2
        rows.append(dict(K=K, dof=A.shape[1], samples=A.shape[0], fit_r2=r2, lambda_max=float(ev.max()), lambda_min=float(ev.min()),
                         n_positive=int((ev > 0).sum()), nsd_exact=ok_exact, nsd_eig=ok_eig, exact_step=step,
                         closure_rate_true_mean=float(true_rate.mean()), closure_rate_model_mean=float(model_rate.mean()),
                         viscous_rate_mean=float(visc.mean()), share_samples_closure_injects=float((true_rate > 0).mean())))
        np.save(os.path.join(ROOT, 'results', f'S1_Q{K}.npy'), Q)
    with open(os.path.join(ROOT, 'results', 'S1_closure.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    b1 = max(rec['budget_residual_whole'], rec['budget_residual_max_window'])
    gate.append(dict(rule='S-G1 energy budget relative residual (whole run, worst 1-unit window)', value=f"{rec['budget_residual_whole']:.2e}, {rec['budget_residual_max_window']:.2e}",
                     threshold='< 1e-3', result='pass' if b1 < 1e-3 else 'fail'))
    b2 = max(rec['advection_neutrality_max_rel'].values())
    gate.append(dict(rule='S-G2 truncated advection energy rate, worst relative over 1000 states per K', value=f'{b2:.2e}', threshold='< 1e-10',
                     result='pass' if b2 < 1e-10 else 'fail'))
    agree = all(r['nsd_exact'] == r['nsd_eig'] for r in rows)
    gate.append(dict(rule='S-G3 exact LDL and eigenvalue test agree on every K', value='; '.join(f"K={r['K']}: exact {r['nsd_exact']}, eig {r['nsd_eig']}" for r in rows),
                     threshold='agree', result='pass' if agree else 'fail'))
    for r in rows:
        gate.append(dict(rule=f"outcome K={r['K']}: closure energy-stable (Q <= 0)", value=f"{r['nsd_exact']}; lambda_max {r['lambda_max']:.3e}; {r['n_positive']} of {r['dof']} eigenvalues > 0; fit R2 {r['fit_r2']:.3f}",
                         threshold='report', result=''))
    with open(os.path.join(ROOT, 'results', 'S1_gate.csv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(gate[0])); w.writeheader(); w.writerows(gate)
    for x in gate:
        print(x['rule'], '|', x['value'], '|', x['result'])
    figure(rec, d, rows)


def figure(rec, d, rows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3), constrained_layout=True)
    t = np.arange(len(d['E'])) * rec['dt'] * 100
    ax[0].plot(t, d['E'], color='0.2', lw=1.0); ax[0].axvline(rec['t_spin'], color='0.7', lw=0.6, ls=':')
    ax[0].set_title('DNS energy', loc='left')
    ax[0].set_xlabel('time'); ax[0].set_ylabel('kinetic energy')
    cols = {1: '#9ecae1', 2: '#4292c6', 3: '#08306b'}
    for r in rows:
        Q = np.load(os.path.join(ROOT, 'results', f"S1_Q{r['K']}.npy")); ev = np.sort(np.linalg.eigvalsh(Q))[::-1]
        ax[1].plot(np.arange(1, len(ev) + 1), ev, 'o-', ms=2.5, lw=0.8, color=cols[r['K']], label=f"K = {r['K']} ({r['dof']} modes)")
    ax[1].axhline(0, color='0.5', lw=0.6); ax[1].set_xlabel('eigenvalue rank'); ax[1].set_ylabel('eigenvalue of Q')
    ax[1].set_title('Closure energy form Q = sym(C)', loc='left'); ax[1].legend(frameon=False)
    x = np.arange(len(rows)); wd = 0.36
    ax[2].bar(x - wd / 2, [r['closure_rate_true_mean'] for r in rows], wd, color='#4292c6', label='closure (true flux)')
    ax[2].bar(x + wd / 2, [r['viscous_rate_mean'] for r in rows], wd, color='0.6', label='resolved viscosity')
    ax[2].axhline(0, color='0.4', lw=0.6); ax[2].set_xticks(x); ax[2].set_xticklabels([f"K = {r['K']}" for r in rows])
    ax[2].set_ylabel('mean energy rate'); ax[2].set_title('Resolved-mode energy budget', loc='left'); ax[2].legend(frameon=False, fontsize=6)
    for a, l in zip(ax, 'abc'):
        a.text(-0.2, 1.08, l, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom')
    os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
    fig.savefig(os.path.join(ROOT, 'figures', 'fig_S1.png'), dpi=300)


if __name__ == '__main__':
    main()
