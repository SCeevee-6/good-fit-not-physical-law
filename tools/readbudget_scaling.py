"""Measurements needed to decide a law sign with probability p, as a function of the margin, under three models
(docs/readbudget_link.md). Exact formulas, no simulation.

 fixed library (Theorem 5 of docs/law_horizon.md): n >= sigma^2 Kbar(u) (sqrt q + z_p)^2 / mu0^2         ~ mu0^-2
 slope-bounded class, readout noise (exploratory/local_theorem.md): n >= (4 L sigma z_p / m^2)^2          ~ m^-4
 slope-bounded class, accumulating noise (exploratory/mixed_theorem.md): K >= 3 L eta^2 z_p^2 / m^3        ~ m^-3
Fitted log-log slopes are printed and written to results/readbudget_scaling.txt; figure figures/fig_readbudget.png.
Usage: python tools/readbudget_scaling.py
"""
import os
import numpy as np
from numpy.polynomial import chebyshev
from scipy.stats import chi2, norm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
P = 0.9; ZP = norm.ppf(P)


def n_fixed(mu0, d=3, g=1.5, sigma=0.05):
    """1-D library of degree d, data at Chebyshev nodes (Kbar exact), law point at gauge g."""
    Kbar = 1 + 2 * sum(chebyshev.chebval(g, [0] * k + [1]) ** 2 for k in range(1, d + 1))
    q = chi2.ppf(0.95, d + 1)
    return sigma ** 2 * Kbar * (np.sqrt(q) + ZP) ** 2 / mu0 ** 2


def n_readout(m, L=1.0, sigma=0.0125):
    return (4 * L * sigma * ZP / m ** 2) ** 2


def k_accum(m, L=1.0, eta=0.025):
    return 3 * L * eta ** 2 * ZP ** 2 / m ** 3


def main():
    m = np.array([0.4, 0.2, 0.1, 0.05])
    rows = {'fixed library': n_fixed(m), 'slope-bounded, readout noise': n_readout(m), 'slope-bounded, accumulating noise': k_accum(m)}
    out = [f'p = {P}; margins {m.tolist()}']
    for k, v in rows.items():
        s = np.polyfit(np.log(m), np.log(v), 1)[0]
        out.append(f'{k}: needed {np.round(v, 2).tolist()}, log-log slope {s:.3f}, cost of halving the margin x{v[1] / v[0]:.2f}')
    open(os.path.join(ROOT, 'results', 'readbudget_scaling.txt'), 'w').write('\n'.join(out) + '\n')
    print('\n'.join(out))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'legend.fontsize': 7,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(3.4, 2.6), constrained_layout=True)
    mm = np.geomspace(0.04, 0.5, 100)
    styles = {'fixed library': ('#1f5fa8', '-', n_fixed), 'slope-bounded, readout noise': ('0.35', '--', n_readout),
              'slope-bounded, accumulating noise': ('0.6', ':', k_accum)}
    for k, (c, ls, f) in styles.items():
        y = f(mm) / f(0.4)
        ax.loglog(mm, y, color=c, ls=ls, lw=1.4 if k == 'fixed library' else 1.0)
        lab = {'fixed library': 'fixed library (-2)', 'slope-bounded, readout noise': 'slope-bounded, readout (-4)',
               'slope-bounded, accumulating noise': 'slope-bounded, accumulating (-3)'}[k]
        if k == 'fixed library':
            ax.text(0.06, f(0.06) / f(0.4) / 6, lab, color=c, fontsize=7, va='top')
        else:
            ax.text(0.0415, y[0] * 1.25, lab, color=c, fontsize=7, va='bottom')
    ax.set_xlabel('law margin'); ax.set_ylabel('measurements needed\n(relative to margin 0.4)')
    ax.set_title('Halving the margin: 4x in a fixed library', loc='left', fontsize=8)
    ax.set_xlim(0.035, 0.5); ax.set_xticks([0.05, 0.1, 0.2, 0.4]); ax.set_xticklabels(['0.05', '0.1', '0.2', '0.4']); ax.minorticks_off()
    ax.set_ylim(0.3, 3e4)
    fig.savefig(os.path.join(ROOT, 'figures', 'fig_readbudget.png'), dpi=300)


if __name__ == '__main__':
    main()
