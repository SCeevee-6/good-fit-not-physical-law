"""Check that each 0A structural mask keeps every term of the true model."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from systems import SYSTEMS
import gate_a
for key in ("A","Amis","B"):
    S = SYSTEMS["B" if key=="B" else "A"]
    for d in S.degrees:
        E = gate_a.library(S, key, d); Xi = S.true_Xi(d)
        if key=="Amis":
            E0 = S.E(d); keep=[tuple(e) not in gate_a.DROP_AMIS for e in E0.tolist()]; Xi = Xi[np.array(keep)]
        allowed, M, fc = gate_a.structure(S, key, E)
        out = [(tuple(E[j]), f) for f_i, f in enumerate(fc) for j in range(len(E)) if Xi[j, f] != 0 and not allowed[j, f_i]]
        print(key, d, "true terms excluded by mask:", out)

# 0A2 structures: the true B model must be feasible, not only unmasked.
import gate_a2
from verify import verify_spec_robust
S = SYSTEMS["B"]
for st in gate_a2.STRUCTS:
    for d in S.degrees:
        E = S.E(d); Xi = S.true_Xi(d)
        allowed, M, fc, lb, eq = gate_a2.structure2(S, st, E)
        Xf = Xi[:, fc]
        excl = [(tuple(E[j]), fc[f]) for f in range(2) for j in range(len(E)) if Xf[j, f] != 0 and not allowed[j, f]]
        lbv = [(tuple(E[j]), fc[f], Xf[j, f]) for f in range(2) for j in range(len(E)) if lb[j, f] and Xf[j, f] < 0]
        eqr = max([abs(float((r * Xf).sum())) for r in eq], default=0.0)
        rec = float(np.abs(Xf @ M.T - Xi).max())
        ver = {sp.name: verify_spec_robust(sp, Xi, E, "strong")["status"] for sp in S.specs}
        print(st, d, "excluded:", excl, "| lb violated:", lbv, "| max eq residual:", eqr,
              "| max |Xi_free M^T - Xi_true|:", rec, "| n_eq:", len(eq), "| true model:", ver)
