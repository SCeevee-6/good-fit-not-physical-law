"""Run with python exploratory/crossover.py. No calibration on evaluation truths."""
import csv
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ndtr

from study import GRID_SIZE, NEG, PARAMS, PHASES, POS, SIGMA, TIMES, meter


def readings(params, s, phase, coupling, extra, term="xv"):
    a, b, c = params.T
    n = len(a)
    def rhs(t, state):
        x, v, _ = state.reshape(3, n)
        perturbation = {"xv": lambda: x*v, "x4": lambda: x**4,
                        "radial": lambda: (x*x+v*v-1.5)**2}[term]()
        gain = a+b*(x*x+v*v-1)+c*x*x+extra*perturbation
        return np.concatenate([v,-x-.2*v+coupling*v*gain,v*v*gain])
    initial = np.repeat(np.sqrt(s)*np.array([np.sin(phase),np.cos(phase),0]),n)
    sol = solve_ivp(rhs,[0,2],initial,t_eval=TIMES,rtol=1e-9,atol=1e-11)
    assert sol.success
    return sol.y[2*n:]


def design(arm):
    pos, neg = POS[:GRID_SIZE], NEG[:GRID_SIZE]
    best = None
    for s in (1.,1.25,1.5):
        for phase in PHASES:
            trace = meter(PARAMS[:GRID_SIZE],s,phase,.3)
            if arm == 'repeat':
                features = trace[:,1:2]
                cov = np.array([[SIGMA**2/4]])
            else:
                features = trace
                cov = SIGMA**2*(.4*np.eye(4)+.6)
            w = np.linalg.solve(cov,features[pos].mean(0)-features[neg].mean(0))
            w /= np.sqrt(w@cov@w)
            projection = features@w
            low, high = projection[pos].min(), projection[neg].max()
            candidate = ((low-high)/2,s,phase,w,(low+high)/2)
            if best is None or candidate[0]>best[0]:
                best = candidate
    return best


def main():
    designs = {arm:design(arm) for arm in ('trace','repeat')}
    switch = ((designs['repeat'][0]/designs['trace'][0])**2-1)/3
    print('Predicted correlation:',switch)
    rng = np.random.default_rng(94371)
    params = rng.uniform([-.22,-.1,-.35],[-.18,.8,.35],(1500,3))
    margin = params[:,0]+.5*params[:,1]
    mask = abs(margin)>=.08
    sign = np.sign(margin[mask])
    rows = []
    for coupling in (.1,.5):
        for extra in (0.,.25):
            for arm,(gap,s,p,w,threshold) in designs.items():
                obs = readings(params[mask],s,p,coupling,extra)
                if not extra:
                    np.testing.assert_allclose(obs[:3],meter(params[mask][:3],s,p,coupling),atol=1e-7)
                features = obs[:,1:2] if arm=='repeat' else obs
                for rho in np.linspace(0,1,11):
                    scale = np.sqrt(1+3*rho) if arm=='repeat' else 1
                    accuracy = ndtr(sign*(features@w-threshold)/scale)
                    rows.append([coupling,extra,arm,s,p,rho,ndtr(gap/scale),
                                 min(accuracy[sign==k].mean() for k in (-1,1)),
                                 accuracy.min(),sum(sign>0),sum(sign<0)])
    with Path(__file__).with_name('crossover.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow('coupling omitted_gain arm s phase reset_correlation predicted_floor worst_sign_mean worst_sample positive_n negative_n'.split())
        writer.writerows(rows)


if __name__=='__main__':
    main()
