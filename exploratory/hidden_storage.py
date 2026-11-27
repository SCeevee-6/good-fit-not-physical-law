"""Reproduce the two passive realizations and the delay-bank resolution test."""
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent


def components(tau, error):
    state = [0., 0.]
    port_error = balance_error = 0.
    for start,end,voltage in [(0,1,2.),(1,3,1.)]:
        def rhs(t,y):
            i,vc = y
            return [(voltage-i)/tau,(voltage-vc)/(tau*(1+error))]
        sol = solve_ivp(rhs,(start,end),state,rtol=1e-11,atol=1e-13,
                        t_eval=np.linspace(start,end,501))
        assert sol.success
        i,vc = sol.y
        di,dvc = rhs(0,sol.y)
        current = i+voltage-vc
        rate = tau*i*di+tau*(1+error)*vc*dvc
        balance = voltage*current-i*i-(voltage-vc)**2
        port_error = max(port_error,max(abs(current-voltage))/2)
        balance_error = max(balance_error,max(abs(rate-balance)))
        state = sol.y[:,-1]
    return rate[0],port_error,balance_error


def power(t):
    return 1.2+.65*np.sin(3*t)+.35*np.cos(7*t)


def main():
    with (ROOT/'hidden_components.csv').open('w',newline='') as f:
        out = csv.writer(f)
        out.writerow(['tau','C_fraction_error','rate_at_step_plus','peak_relative_port_error','max_balance_residual'])
        for tau in (.2,3.):
            for e in (0.,.01):
                row = components(tau,e)
                assert row[2]<1e-10
                out.writerow([tau,e,*row])
    rng = np.random.default_rng(91273)
    delays = 2*rng.beta(2,4,80)
    weights = rng.dirichlet(np.ones(len(delays)))
    times = np.linspace(3,10,701)
    p = power(times)
    actual = p-power(times[:,None]-delays)@weights
    with (ROOT/'hidden_delay.csv').open('w',newline='') as f:
        out = csv.writer(f)
        out.writerow(['bins','mean_width','certified_sign_fraction','containment_fraction'])
        for n in (1,4,8,16,32):
            edges = np.linspace(0,2,n+1)
            # The measured calibration step is the aggregate terminal heat CDF.
            cdf = (delays<=edges[:,None])@weights
            masses = np.diff(cdf)
            mids = .5*(edges[:-1]+edges[1:])
            # Rigorous smooth-wave extrema enclosures, using a subgrid plus
            # Lipschitz padding. Includes every admissible within-bin delay.
            offsets = np.linspace(-1/n,1/n,33)
            values = power(times[:,None,None]-mids[None,:,None]-offsets)
            pad = (3*.65+7*.35)/(32*n)
            lo = p-(values.max(2)+pad)@masses
            hi = p-(values.min(2)-pad)@masses
            covered = (actual>=lo)&(actual<=hi)
            assert covered.all()
            out.writerow([n,np.mean(hi-lo),np.mean((lo>0)|(hi<0)),np.mean(covered)])

if __name__=='__main__':
    main()
