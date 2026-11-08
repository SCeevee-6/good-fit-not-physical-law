
import pickle, sys, pandas as pd, numpy as np
pd.set_option('display.width', 200)
for path in sys.argv[1:]:
    d = pickle.load(open(path,'rb')); c, r, a, sp = d['cand'], d['repair'], d['audit'], d['spec']
    print('=====', d['cell'], '| %.0f s' % d['seconds'])
    print('cands', len(c), 'AES', int(c.aes.sum()), 'strong violators', int(c.smt_viol_strong.sum()), '| in AES', int((c.aes & c.smt_viol_strong).sum()), '| strong timeouts', int(c.smt_timeout_strong.sum()))
    if len(r):
        g = r.groupby('arm').agg(n=('status','size'), repaired=('status', lambda s: (s=='repaired').mean()), infeasible=('status', lambda s: (s=='infeasible').mean()),
                                 failed=('status', lambda s: s.isin(['failed','failed_max_iter']).mean()), iters=('iters','median'), rows=('n_rows','median'), sec=('seconds','median')).round(3)
        print(g.to_string())
        ok = r[(r.arm=='ceg_witness') & (r.status=='repaired')]
        if len(ok): print('CEG repaired n=%d: iid_rmse/original median %.3f | /best-in-cell median %.3f | k median %d | support F1 median %.2f (orig %.2f)' % (len(ok), (ok.iid_rmse/ok.orig_iid_rmse).median(), (ok.iid_rmse/ok.best_iid_rmse).median(), ok.k.median(), ok.support_f1.median(), ok.orig_support_f1.median()))
    if len(a): print('audit n=%d:' % len(a), a[['audit_strong_all_verified','audit_weak_all_verified','audit_local_oracle_viol']].mean().round(3).to_dict(), 'audit timeouts', int(a.audit_timeouts.sum()))
    print('cross-checks: bernstein %d (disagree %d) | exact-1D %d (disagree %d)' % (c.bern_checks.sum(), c.bern_disagree.sum(), c.exact1d_checks.sum(), c.exact1d_disagree.sum()))
    cx = sp[(sp.level=='strong') & (sp.status=='counterexample')]
    print('strong cex n=%d | replay ok %.4f | exact witness %.4f' % (len(cx), cx.replay_ok.mean(), cx.witness_exact.mean()))
    print('methods:', sp[sp.level=='strong'].method.value_counts().to_dict())
