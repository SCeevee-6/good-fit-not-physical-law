
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, r"C:\Users\shree\Downloads\CEG_ODE_Verification\src")
from systems import SYSTEMS
from verify import verify_spec_robust
from lib import make_rhs
from scipy.integrate import solve_ivp
R = r"C:\Users\shree\Downloads\CEG_ODE_Verification\results"
cand = pd.read_csv(R + r"\pilot_candidates.csv"); spec = pd.read_csv(R + r"\pilot_spec_queries.csv")
cx = spec[(spec.level == 'strong') & (spec.status == 'counterexample')]
bad = cx[~cx.replay_ok.astype(bool)]
rows = []
for _, b in bad.iterrows():
    c = cand[(cand.cell == b.cell) & (cand.cand == b.cand)].iloc[0]
    S = SYSTEMS['A']; E = S.E(int(c.degree)); Xi = np.array(json.loads(c.xi)).reshape(json.loads(c.xi_shape))
    sp = [s for s in S.specs if s.name == b.spec][0]
    r = verify_spec_robust(sp, Xi, E, 'strong')
    pt = [float(p) for p in sp.point([np.array(x) for x in r['u']], lambda fr: float(fr))]
    f = make_rhs(Xi, E); out = {}
    for h in (0.05, 0.01, 0.001):
        try:
            sol = solve_ivp(f, (0, h), np.array(pt), rtol=1e-12, atol=1e-14, method='LSODA')
            out[h] = (bool(sol.success), float(sol.t[-1]))
        except Exception as e:
            out[h] = (False, float('nan'))
    blow = False
    try:
        ev = lambda t, x: 1e6 - np.max(np.abs(x)); ev.terminal = True
        sol2 = solve_ivp(f, (0, 0.05), np.array(pt), rtol=1e-10, atol=1e-12, method='RK45', events=ev)
        blow = len(sol2.t_events[0]) > 0 or sol2.status != 0
    except Exception:
        blow = True
    # behavioural violation at tiny horizon: c.x change
    cvec = np.array([float(x) for x in sp.c]); sol3 = solve_ivp(f, (0, 1e-5), np.array(pt), rtol=1e-12, atol=1e-14, method='LSODA')
    d_small = float(cvec @ sol3.y[:, -1] - cvec @ np.array(pt))
    rows.append(dict(cell=b.cell, cand=b.cand, spec=b.spec, G=r['G'], pt0=pt[0], pt1=pt[1], reached_0p05=bool(out[0.05][1] >= 0.05 - 1e-9), success_0p05=out[0.05][0], blowup_before_0p05=blow, t_end_0p05=out[0.05][1], change_at_1e5=d_small))
df = pd.DataFrame(rows); df.to_csv(R + r'\replay_failure_diagnosis.csv', index=False); pd.set_option('display.width', 200)
print(df.describe(include='all').T[['count','mean','min','max']].round(4).to_string())
print('integration reached the 0.05 horizon:', int(df.reached_0p05.sum()), 'of', len(df), '| blow-up (|x|>1e6) before 0.05:', int(df.blowup_before_0p05.sum()), '| violation visible at t=1e-5 (change<0):', int((df.change_at_1e5 < 0).sum()))
print(df[['spec','G','pt0','pt1','t_end_0p05','blowup_before_0p05','change_at_1e5']].head(8).round(4).to_string())
