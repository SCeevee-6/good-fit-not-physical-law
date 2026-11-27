"""Integrate signed auxiliary power before and after a crossing."""
import csv
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

TERMS = ['xv3','sine','cubic','even']

def integral(beta, damping, t):
    if t == 0:
        return np.zeros(4)
    def rhs(_, state):
        x,v = state[:2]
        return [v,-x-beta*x**3-damping*v,x*v**3,np.sin(3*x)*v*v,x**3*v*v,x*x*v*v]
    sol = solve_ivp(rhs,(0,t),[0,np.sqrt(1.5),0,0,0,0],rtol=1e-12,atol=1e-14,method='DOP853')
    assert sol.success
    return sol.y[2:,-1]

def main():
    rows = []
    for beta in (0,2):
        scale = integral(beta,0,.5)
        for damping in (0,.2):
            for offset in (0,.025,.25):
                left,right = offset-.25,offset+.25
                ratios = (integral(beta,damping,right)-integral(beta,damping,left))/scale
                if damping==offset==0:
                    assert max(abs(ratios[:3]))<1e-8
                    assert ratios[3]>.1
                rows.extend([beta,damping,offset,term,ratio] for term,ratio in zip(TERMS,ratios))
    with Path(__file__).with_name('reflection.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['stiffness','damping','center_offset','term','relative_meter_shift'])
        writer.writerows(rows)

if __name__=='__main__':
    main()
