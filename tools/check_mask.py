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
