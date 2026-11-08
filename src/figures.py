# Figure code for the pilot report. Executed with exec() inside a kernel where the figure-style helpers
# (apply_figure_style, panel_letter) are defined; reads the analysis tables produced by analyze.py.
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
import analyze as an

COL = {'A': '#4C72B0', 'B': '#55A868', 'A-mis': '#C44E52'}
LAB = {'A': 'Predator–prey', 'B': 'SIR', 'A-mis': 'Predator–prey,\nmissing x² term'}


def _bbox_check(fig, name):
    r = fig.canvas.get_renderer()
    ticks = set(t for a in fig.axes for t in a.get_xticklabels(which='both') + a.get_yticklabels(which='both'))
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible() and t not in ticks]
    ov = [(a.get_text()[:20], b.get_text()[:20]) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    print(name, 'text overlaps:', ov)


def fig_gap(cand, g1, path):
    apply_figure_style()
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(wspace=0.32))
    ax = axs[0]
    rhos = [0.005, 0.02, 0.05]
    for sysn in ('A', 'B', 'A-mis'):
        for subset, mk, filled in (('AES', 'o', True), ('AES-parsimonious', 's', False)):
            s = g1[(g1.system == sysn) & (g1.subset == subset)].sort_values('rho')
            if s.empty: continue
            x = np.arange(3)[[rhos.index(r) for r in s.rho]] + {'A': -0.08, 'B': 0.0, 'A-mis': 0.08}[sysn]
            ax.plot(x, s.frac, color=COL[sysn], lw=1 if filled else 0.7, ls='-' if filled else ':', marker=mk, ms=4.5, mfc=COL[sysn] if filled else 'white', mec=COL[sysn], zorder=3 if filled else 2)
            if filled:
                ax.vlines(x, s.seed_min, s.seed_max, color=COL[sysn], lw=0.8, alpha=0.6)
    ax.axhline(0.2, color='k', lw=0.6, ls='--'); ax.text(2.35, 0.215, 'G1 threshold 20 %', fontsize=6, ha='right')
    ax.set_xticks(range(3)); ax.set_xticklabels(['0.5 %', '2 %', '5 %']); ax.set_xlabel('Measurement noise'); ax.set_ylabel('Fraction violating a law')
    ax.set_ylim(-0.03, 1.05); ax.set_xlim(-0.4, 2.4)
    ax.set_title('Many equally accurate models break a law', loc='left')
    for sysn, yy in (('A', 0.30), ('B', 0.78), ('A-mis', 0.97)):
        pass
    h = [plt.Line2D([], [], color=COL[s], marker='o', ms=4, lw=1, label=LAB[s].replace('\n', ' ')) for s in ('A', 'B', 'A-mis')]
    h += [plt.Line2D([], [], color='0.4', marker='o', ms=4, lw=1, label='all accuracy-equivalent'), plt.Line2D([], [], color='0.4', marker='s', ms=4, mfc='white', lw=0.7, ls=':', label='only k ≤ k_true+4 terms')]
    fig.legend(handles=h, loc='lower center', ncol=5, bbox_to_anchor=(0.5, -0.08), fontsize=5.5, handlelength=1.6, columnspacing=1.2)
    panel_letter(ax, 'a')
    # b: violation probability by number of terms (pooled, well-specified + mis arms separately coloured)
    ax = axs[1]
    a = cand[cand.aes & ~cand.undecided].copy(); a['dk'] = a.k - a.ktrue
    bins = [-6, 0, 2, 4, 8, 16, 100]; names = ['≤0', '1–2', '3–4', '5–8', '9–16', '>16']
    a['kb'] = pd.cut(a.dk, bins=bins, labels=names, right=True)
    for sysn in ('A', 'B', 'A-mis'):
        g = a[a.system == sysn].groupby('kb', observed=True).agg(p=('inadmissible', 'mean'), n=('inadmissible', 'size')).reindex(names)
        x = np.arange(len(names)) + {'A': -0.1, 'B': 0.0, 'A-mis': 0.1}[sysn]
        ok = g.n.fillna(0) >= 5
        ax.plot(x[ok.values], g.p[ok], color=COL[sysn], marker='o', ms=4, lw=1)
    ax.set_xticks(range(len(names))); ax.set_xticklabels(names); ax.set_xlabel('Non-zero terms beyond the true model'); ax.set_ylabel('Fraction violating a law')
    ax.set_ylim(-0.03, 1.05); ax.set_title('Violations concentrate in larger models', loc='left'); ax.margins(x=0.06)
    panel_letter(ax, 'b')
    fig.savefig(path, dpi=300, bbox_inches='tight'); _bbox_check(fig, 'fig_gap'); plt.close(fig)


def _roc(y, s):
    y = np.asarray(y, bool); s = np.asarray(s, float)
    thr = np.unique(s)[::-1]
    P, N = y.sum(), (~y).sum()
    tpr = [0.0]; fpr = [0.0]
    for t in thr:
        tpr.append((s[y] >= t).sum() / P); fpr.append((s[~y] >= t).sum() / N)
    return np.array(fpr), np.array(tpr)


def fig_detection(cand, sel, path):
    apply_figure_style()
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), gridspec_kw=dict(wspace=0.38, width_ratios=[1, 1, 1.05]))
    a = cand[cand.aes & cand.oracle_run & ~cand.undecided & (cand.system != 'A-mis')].copy()
    for sc in an.SCORES: a['pr_' + sc] = a.groupby('cell')[sc].rank(pct=True)
    styles = [('iid_rmse', 'held-out error', '#8c8c8c', '-'), ('aic', 'AIC', '#DD8452', '-'), ('roll_ood', 'out-of-range rollout*', '#937860', '--'), ('rf200_score', 'random testing (200 pts)', '#8172B3', '-')]
    for ax, sysn, letter in zip(axs[:2], ('A', 'B'), 'ab'):
        g = a[a.system == sysn]
        for sc, lab, c, ls in styles:
            f, t = _roc(g.inadmissible, g['pr_' + sc]); ax.plot(f, t, color=c, ls=ls, lw=1.1, label=lab)
        vt = g[g.inadmissible].smt_viol_strong.mean(); vf = (g.smt_viol_strong & ~g.inadmissible).sum() / max((~g.inadmissible).sum(), 1)
        ax.plot([vf], [vt], color='#C44E52', marker='*', ms=7, ls='', label='SMT verifier')
        ax.plot([0, 1], [0, 1], color='0.8', lw=0.6)
        ax.set_xlabel('False-positive rate'); ax.set_ylabel('True-positive rate'); ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.04)
        ax.set_title(LAB[sysn].replace('\n', ' ') + f' (n = {len(g)})', loc='left'); panel_letter(ax, letter)
    axs[0].legend(loc='lower right', fontsize=5.5, handlelength=1.5, borderaxespad=0.3)
    ax = axs[2]
    order = ['A', 'B', 'A-mis']; w = 0.36
    for j, (col, lab, hatch) in enumerate((('validation_pick_inadmissible', 'best held-out error', None), ('aic_pick_inadmissible', 'best AIC', None))):
        v = [sel[sel.system == s][col].mean() for s in order]
        ax.bar(np.arange(3) + (j - 0.5) * w, v, w * 0.97, color=[COL[s] for s in order], alpha=1.0 if j == 0 else 0.55, edgecolor='white', lw=0.5)
        for xi, vi in zip(np.arange(3) + (j - 0.5) * w, v): ax.text(xi, vi + 0.02, f'{vi:.2f}', ha='center', fontsize=5)
    ax.set_xticks(range(3)); ax.set_xticklabels([LAB[s] for s in order], fontsize=5.5)
    ax.set_ylabel('Fraction of cells where pick\nviolates a law'); ax.set_ylim(0, 1.15)
    ax.set_title('The model validation picks', loc='left'); panel_letter(ax, 'c')
    ax.text(0.02, 0.97, 'solid: held-out error, light: AIC', transform=ax.transAxes, fontsize=5.5, va='top')
    fig.savefig(path, dpi=300, bbox_inches='tight'); _bbox_check(fig, 'fig_detection'); plt.close(fig)


def fig_repair(rep, path):
    apply_figure_style()
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.8), gridspec_kw=dict(wspace=0.42))
    arms = [('ceg_witness', 'CEG'), ('ceg_worst', 'CEG +\nworst pt'), ('dense_constraints', 'dense\ngrid'), ('random_constraints', 'random\npoints')]
    ax = axs[0]
    for j, sysn in enumerate(('A', 'B', 'A-mis')):
        v = [rep[(rep.system == sysn) & (rep.arm == a)].repaired.mean() for a, _ in arms] + [0.0]
        x = np.arange(5) + (j - 1) * 0.26
        ax.bar(x, v, 0.26, color=COL[sysn], label=LAB[sysn].replace('\n', ' '), edgecolor='white', lw=0.4)
    ax.set_xticks(range(5)); ax.set_xticklabels([l for _, l in arms] + ['reject\nonly'], fontsize=5.5)
    ax.set_ylabel('Violating models repaired'); ax.set_ylim(0, 1.45); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.axhline(0.8, color='k', lw=0.6, ls='--'); ax.text(4.45, 0.82, 'G4: 80 %', fontsize=5.5, ha='right')
    ax.legend(loc='upper left', fontsize=5.2, handlelength=1.0, bbox_to_anchor=(0.0, 1.0), ncol=1, borderaxespad=0.1); ax.set_title('CEG repairs nearly all violators', loc='left')
    panel_letter(ax, 'a')
    ax = axs[1]
    arms_b = arms[:3]
    for j, (a, l) in enumerate(arms_b):
        for k, sysn in enumerate(('A', 'B', 'A-mis')):
            d = rep[(rep.system == sysn) & (rep.arm == a)].n_rows.dropna()
            if len(d) == 0: continue
            x = j + (k - 1) * 0.26
            ax.vlines(x, d.quantile(0.25), d.quantile(0.75), color=COL[sysn], lw=3, alpha=0.5); ax.plot(x, d.median(), 'o', color=COL[sysn], ms=3.5)
    ax.set_yscale('log'); ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in arms_b], fontsize=5.5)
    ax.set_ylabel('Constraint rows (median, IQR)'); ax.set_title('CEG needs few constraints', loc='left'); panel_letter(ax, 'b')
    ax = axs[2]
    for k, sysn in enumerate(('A', 'B', 'A-mis')):
        d = rep[(rep.system == sysn) & (rep.arm == 'ceg_witness') & rep.repaired].rmse_ratio.dropna()
        if len(d) == 0: continue
        jit = np.random.default_rng(k).uniform(-0.15, 0.15, len(d)); dd = d.clip(upper=1.6)
        ax.plot(np.full(len(d), k) + jit, dd, '.', color=COL[sysn], ms=2.5, alpha=0.35)
        ax.text(k, 1.67, f'{(d > 1.1).mean():.0%} above\nthe limit', ha='center', va='bottom', fontsize=5, color=COL[sysn])
        ax.plot([k - 0.22, k + 0.22], [d.median()] * 2, color='k', lw=1.2)
    ax.axhline(1.0, color='0.6', lw=0.6); ax.axhline(1.10, color='k', lw=0.6, ls='--'); ax.text(-0.45, 1.12, 'G4 limit: +10 %', fontsize=5.5, ha='left')
    ax.set_xticks(range(3)); ax.set_xticklabels([LAB[s] for s in ('A', 'B', 'A-mis')], fontsize=5.2)
    ax.set_ylim(0.9, 1.86); ax.set_ylabel('Held-out error after ÷ before repair'); ax.set_title('Median accuracy cost is nil', loc='left'); ax.set_xlim(-0.5, 2.5); panel_letter(ax, 'c')
    fig.savefig(path, dpi=300, bbox_inches='tight'); _bbox_check(fig, 'fig_repair'); plt.close(fig)
