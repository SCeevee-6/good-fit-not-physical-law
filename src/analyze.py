"""Pilot analysis: G1-G4, sensitivity subsets, strata. Reads results/pilot_*.csv (produced by combine.py)."""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import rankdata, binomtest

HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
RHO = {0: 0.005, 1: 0.02, 2: 0.05}
KTRUE = {'A': 5, 'A-mis': 5, 'B': 4}
SCORES = ['iid_rmse', 'roll_iid', 'roll_ood', 'aic', 'rf200_score', 'rf2000_score']
CONVENTIONAL = ['iid_rmse', 'roll_iid', 'roll_ood', 'aic']
PRACTICAL = ['iid_rmse', 'roll_iid', 'aic']
TOL_STRONG = 0.02


def auroc(y, s):
    y = np.asarray(y, bool); s = np.asarray(s, float)
    if y.all() or (~y).all() or len(y) < 4:
        return np.nan
    r = rankdata(s); n1 = y.sum(); n0 = (~y).sum()
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def tpr_at_fpr(y, s, fpr=0.05):
    y = np.asarray(y, bool); s = np.asarray(s, float)
    if y.all() or (~y).all():
        return np.nan
    thr = np.quantile(s[~y], 1 - fpr)
    return float((s[y] > thr).mean())


def boot_ci(vals, n=2000, seed=0):
    v = np.asarray([x for x in vals if np.isfinite(x)])
    if len(v) < 2:
        return (np.nan, np.nan)
    rng = np.random.default_rng(seed)
    m = [rng.choice(v, len(v)).mean() for _ in range(n)]
    return float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))


def load():
    cand = pd.read_csv(os.path.join(RES, 'pilot_candidates.csv'))
    spec = pd.read_csv(os.path.join(RES, 'pilot_spec_queries.csv'))
    rep = pd.read_csv(os.path.join(RES, 'pilot_repairs.csv'))
    aud = pd.read_csv(os.path.join(RES, 'pilot_repair_audit.csv'))
    cand['system'] = cand.key.map({'A': 'A', 'Amis': 'A-mis', 'B': 'B'}); cand['rho'] = cand.noise_idx.map(RHO)
    for c in ('aes', 'smt_viol_strong', 'smt_timeout_strong', 'smt_viol_weak', 'oracle_run', 'local_viol_strong', 'local_viol_weak', 'traj_unphysical', 'rf200', 'rf2000', 'exact_support'):
        cand[c] = cand[c].fillna(False).astype(bool)
    # replay validity per candidate (all strong counterexamples replay)
    cx = spec[(spec.level == 'strong') & (spec.status == 'counterexample')].copy()
    cx['replay_ok'] = cx.replay_ok.astype(bool)
    # amended replay (post hoc, documented): a strict-replay failure that is a confirmed finite-time blow-up counts as confirmed behaviour
    dp = os.path.join(RES, 'replay_failure_diagnosis.csv')
    if os.path.exists(dp):
        dg = pd.read_csv(dp); dg = dg[dg.blowup_before_0p05.astype(bool) & (dg.change_at_1e5 < 0)][['cell', 'cand', 'spec']].assign(amended=True)
        cx = cx.merge(dg, on=['cell', 'cand', 'spec'], how='left'); cx['replay_ok_strict'] = cx.replay_ok; cx['replay_ok'] = cx.replay_ok | cx.amended.fillna(False).astype(bool)
    rp = cx.groupby(['cell', 'cand']).replay_ok.all().rename('replay_all').reset_index()
    cand = cand.merge(rp, on=['cell', 'cand'], how='left'); cand['replay_all'] = cand.replay_all.fillna(True).astype(bool)
    # label: inadmissible = (SMT strong cex & replay-validated) OR local-oracle strong violation; undecided if timeout & not inadmissible
    inad = (cand.smt_viol_strong & cand.replay_all) | (cand.oracle_run & cand.local_viol_strong)
    cand['inadmissible'] = inad
    cand['undecided'] = (~inad) & cand.smt_timeout_strong
    cand['admissible'] = (~inad) & (~cand.undecided)
    # severity from local oracle minima
    def sev(js):
        if not isinstance(js, str): return np.nan
        d = json.loads(js); return float(-min(d.values()) / TOL_STRONG)
    cand['severity'] = cand.local_min_json.apply(sev)
    # per-cell derived subsets
    g = cand.groupby('cell')
    cand['best_roll_iid'] = g.roll_iid.transform('min')
    cand['aes_strict'] = (cand.iid_rmse <= 1.05 * cand.best_iid_rmse) & (cand.roll_iid <= 1.25 * cand.best_roll_iid)
    cand['ktrue'] = cand.system.map(KTRUE)
    cand['aes_parsimonious'] = cand.aes & (cand.k <= cand.ktrue + 4)
    return cand, spec, rep, aud


def g1_table(cand):
    rows = []
    for (sysn, rho), g in cand.groupby(['system', 'rho']):
        for subset, mask in (('AES', g.aes), ('AES-strict', g.aes_strict), ('AES-parsimonious', g.aes_parsimonious)):
            s = g[mask & ~g.undecided]
            per_seed = s.groupby('seed').inadmissible.mean()
            rows.append(dict(system=sysn, rho=rho, subset=subset, n_cells=g.cell.nunique(), n=len(s), n_inadmissible=int(s.inadmissible.sum()), frac=float(s.inadmissible.mean()) if len(s) else np.nan,
                             seed_min=float(per_seed.min()) if len(per_seed) else np.nan, seed_max=float(per_seed.max()) if len(per_seed) else np.nan))
        # weak tolerance (SMT-only)
        s = g[g.aes]; rows.append(dict(system=sysn, rho=rho, subset='AES (weak tol, SMT only)', n_cells=g.cell.nunique(), n=len(s), n_inadmissible=int(s.smt_viol_weak.sum()), frac=float(s.smt_viol_weak.mean()), seed_min=np.nan, seed_max=np.nan))
    return pd.DataFrame(rows)


def selection_table(cand):
    rows = []
    for (sysn, rho), g in cand.groupby(['system', 'rho']):
        res = dict(system=sysn, rho=rho, n_cells=g.cell.nunique(), val_pick_inad=[], aic_pick_inad=[], any_admissible_in_aes=[], guided_ratio=[])
        for cell, c in g.groupby('cell'):
            c = c[~c.undecided]
            v = c.loc[c.iid_rmse.idxmin()]; a = c.loc[c.aic.idxmin()]
            res['val_pick_inad'].append(bool(v.inadmissible)); res['aic_pick_inad'].append(bool(a.inadmissible))
            adm = c[c.aes & c.admissible]
            res['any_admissible_in_aes'].append(len(adm) > 0)
            res['guided_ratio'].append(float(adm.iid_rmse.min() / c.iid_rmse.min()) if len(adm) else np.nan)
        rows.append(dict(system=sysn, rho=rho, n_cells=res['n_cells'], validation_pick_inadmissible=np.mean(res['val_pick_inad']), aic_pick_inadmissible=np.mean(res['aic_pick_inad']),
                         admissible_model_exists_in_AES=np.mean(res['any_admissible_in_aes']), guided_pick_rmse_ratio_median=np.nanmedian(res['guided_ratio']) if np.isfinite(res['guided_ratio']).any() else np.nan))
    return pd.DataFrame(rows)


def detection_tables(cand):
    """Per-cell AUROC / TPR@5%FPR within the AES, restricted to cells where oracle ran for the whole AES."""
    rows = []
    for cell, c in cand.groupby('cell'):
        a = c[c.aes & c.oracle_run & ~c.undecided]
        if len(a) < 8: continue
        y = a.inadmissible.values
        for sc in SCORES:
            rows.append(dict(cell=cell, system=a.system.iloc[0], rho=a.rho.iloc[0], seed=a.seed.iloc[0], score=sc, n=len(a), n_pos=int(y.sum()), auroc=auroc(y, a[sc]), tpr5=tpr_at_fpr(y, a[sc])))
    d = pd.DataFrame(rows)
    out = []
    for (sysn, sc), g in d.groupby(['system', 'score']):
        lo, hi = boot_ci(g.auroc.values); lo5, hi5 = boot_ci(g.tpr5.values)
        out.append(dict(system=sysn, score=sc, n_cells_defined=int(g.auroc.notna().sum()), n_cells=len(g), auroc_mean=g.auroc.mean(), auroc_lo=lo, auroc_hi=hi, tpr5_mean=g.tpr5.mean(), tpr5_lo=lo5, tpr5_hi=hi5))
    return d, pd.DataFrame(out)


def strata_table(cand):
    """Pooled AUROC of each score within complexity strata, using within-cell percentile ranks of the score."""
    rows = []
    a = cand[cand.aes & cand.oracle_run & ~cand.undecided & (cand.system != 'A-mis')].copy()
    for sc in SCORES:
        a['pr_' + sc] = a.groupby('cell')[sc].rank(pct=True)
    def stratum(r):
        d = r.k - r.ktrue
        return 'k <= k_true+4' if d <= 4 else ('k_true+5 .. k_true+12' if d <= 12 else 'k > k_true+12')
    a['stratum'] = a.apply(stratum, axis=1)
    for (sysn, st), g in a.groupby(['system', 'stratum']):
        for sc in SCORES:
            rows.append(dict(system=sysn, stratum=st, n=len(g), n_pos=int(g.inadmissible.sum()), score=sc, auroc=auroc(g.inadmissible, g['pr_' + sc])))
    return pd.DataFrame(rows)


def severity_table(cand):
    a = cand[cand.aes & cand.oracle_run & cand.inadmissible & cand.severity.notna()].copy()
    bins = [1.0, 1.5, 3, 10, np.inf]; labels = ['1-1.5x tol', '1.5-3x', '3-10x', '>10x']
    a['sev_bin'] = pd.cut(a.severity, bins=bins, labels=labels, right=False)
    rows = []
    for (sysn, b), g in a.groupby(['system', 'sev_bin'], observed=True):
        rows.append(dict(system=sysn, severity=str(b), n=len(g), smt_detect=g.smt_viol_strong.mean(), random200_detect=g.rf200.mean(), random2000_detect=g.rf2000.mean(), traj_unphysical=g.traj_unphysical.mean()))
    return pd.DataFrame(rows)


def slack_table(cand):
    a = cand[cand.oracle_run & ~cand.undecided]
    rows = []
    for (sysn,), g in a.groupby(['system']):
        v = g[g.smt_viol_strong]; ok = g[~g.smt_viol_strong]; okw = g[~g.smt_viol_weak]
        rows.append(dict(system=sysn, n_viol=len(v), traj_unphysical_given_viol=v.traj_unphysical.mean(), n_verified_strong=len(ok), traj_unphysical_given_strong_verified=ok.traj_unphysical.mean(),
                         n_verified_weak=len(okw), traj_unphysical_given_weak_verified=okw.traj_unphysical.mean() if len(okw) else np.nan))
    return pd.DataFrame(rows)


def g3_table(cand, spec):
    st = spec[spec.level == 'strong']
    cx = st[st.status == 'counterexample']
    out = dict(n_candidates=len(cand), n_strong_queries=len(st), strong_counterexamples=len(cx), replay_validated_frac=float(cx.replay_ok.astype(bool).mean()), exact_witness_frac=float(cx.witness_exact.astype(bool).mean()),
               decided_frac=float((st.status != 'timeout').mean()), z3_alone_decided_frac=float((st.method == 'z3').mean()), fallback_bernstein=int((st.method == 'bernstein_after_z3_timeout').sum()), fallback_split=int((st.method == 'z3_split').sum()),
               bernstein_cross_checks=int(cand.bern_checks.sum()), bernstein_disagreements=int(cand.bern_disagree.sum()), exact1d_cross_checks=int(cand.exact1d_checks.sum()), exact1d_disagreements=int(cand.exact1d_disagree.sum()))
    r = cand[cand.oracle_run]
    out['smt_missed_vs_local_oracle'] = int((r.local_viol_strong & ~r.smt_viol_strong).sum()); out['oracle_missed_vs_smt'] = int((r.smt_viol_strong & ~r.local_viol_strong).sum())
    out['verifier_tpr_vs_label'] = float(r[r.inadmissible].smt_viol_strong.mean())
    return pd.DataFrame([out]).T.rename(columns={0: 'value'})


def g4_tables(cand, rep, aud):
    rep = rep.copy()
    rep = rep.merge(cand[['cell', 'cand', 'system', 'rho', 'seed', 'iid_rmse']].rename(columns={'iid_rmse': 'orig_iid_rmse2'}), on=['cell', 'cand'], how='left')
    rep['repaired'] = rep.status == 'repaired'
    rep['rmse_ratio'] = rep.iid_rmse / rep.orig_iid_rmse
    rep['rmse_vs_best'] = rep.iid_rmse / rep.best_iid_rmse
    arm = rep.groupby(['system', 'arm']).agg(n=('repaired', 'size'), repaired=('repaired', 'mean'), infeasible=('status', lambda s: (s == 'infeasible').mean()), iters_median=('iters', 'median'), rows_median=('n_rows', 'median'),
                                             seconds_median=('seconds', 'median')).reset_index()
    for sysn_arm, g in rep[rep.repaired].groupby(['system', 'arm']):
        m = (arm.system == sysn_arm[0]) & (arm.arm == sysn_arm[1])
        arm.loc[m, 'rmse_ratio_median'] = g.rmse_ratio.median(); arm.loc[m, 'rmse_vs_best_median'] = g.rmse_vs_best.median(); arm.loc[m, 'support_f1_median'] = g.support_f1.median(); arm.loc[m, 'k_median'] = g.k.median()
    # paired CEG vs random
    w = rep[rep.arm == 'ceg_witness'].set_index(['cell', 'cand']).repaired; r = rep[rep.arm == 'random_constraints'].set_index(['cell', 'cand']).repaired
    j = pd.concat([w.rename('ceg'), r.rename('rand')], axis=1).dropna()
    b = int((j.ceg & ~j.rand).sum()); c = int((~j.ceg & j.rand).sum())
    p = binomtest(b, b + c, 0.5).pvalue if b + c > 0 else np.nan
    paired = pd.DataFrame([dict(n_pairs=len(j), ceg_only=b, random_only=c, ceg_rate=float(j.ceg.mean()), random_rate=float(j.rand.mean()), diff_pp=100 * float(j.ceg.mean() - j.rand.mean()), mcnemar_p=p)])
    audit = pd.DataFrame([dict(n=len(aud), strong_z3_verified=float(aud.audit_strong_all_verified.mean()), weak_z3_verified=float(aud.audit_weak_all_verified.mean()), local_oracle_violation=float(aud.audit_local_oracle_viol.mean()), z3_timeouts=int(aud.audit_timeouts.sum()))])
    return arm, paired, audit, rep
