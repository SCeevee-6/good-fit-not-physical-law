"""Run from the repository root: python exploratory/study.py"""
import csv
from pathlib import Path

import numpy as np
from scipy.special import ndtr
from scipy.stats import chi2

OUT = Path(__file__).resolve().parent
TIMES = (.25, .5, 1., 2.)
SIGMA = .02
PHASES = np.arange(8) * np.pi / 8
PARAMS = np.array([(a,b,c) for a in (-.22,-.20,-.18)
                   for b in np.linspace(-.1,.8,10)
                   for c in np.linspace(-.35,.35,5)])
GRID_SIZE = len(PARAMS)
OFF_GRID = np.random.default_rng(5091).uniform([-.22,-.1,-.35],[-.18,.8,.35],(120,3))
PARAMS = np.vstack([PARAMS, OFF_GRID])
MARGIN = PARAMS[:,0] + .5 * PARAMS[:,1]
POS = MARGIN >= .08 - 1e-12
NEG = MARGIN <= -.08 + 1e-12


def meter(params, radius_squared, phase, coupling, dt=.01):
    """Vectorized RK4. Meter integrates unscaled channel power."""
    a,b,c = params.T
    def rhs(z):
        x,v,_ = z
        gain = a + b*(x*x + v*v - 1) + c*x*x
        return np.array([v,-x-.2*v+coupling*v*gain,v*v*gain])
    z = np.tile(np.sqrt(radius_squared)*np.array([np.sin(phase),np.cos(phase),0.]),
                (len(params),1)).T
    samples = []
    for step in range(1,round(2/dt)+1):
        k1 = rhs(z)
        k2 = rhs(z+dt*k1/2)
        k3 = rhs(z+dt*k2/2)
        k4 = rhs(z+dt*k3)
        z = z+dt*(k1+2*k2+2*k3+k4)/6
        if step in [round(t/dt) for t in TIMES]:
            samples.append(z[2].copy())
    return np.stack(samples,axis=1)


def detector(means, positive, negative):
    """Project along whitened class centroids; use extreme projected means."""
    direction = means[positive].mean(0)-means[negative].mean(0)
    direction /= max(np.linalg.norm(direction),1e-15)
    scores = means @ direction
    low, high = scores[positive].min(), scores[negative].max()
    return (low-high)/2, scores-(low+high)/2


def run():
    rows = []
    for coupling in (0.,.3,.6):
        trajectories = np.array([[meter(PARAMS,s,p,coupling) for p in PHASES]
                                 for s in (1.,1.25,1.5)])
        refined = meter(PARAMS,1.5,PHASES[0],coupling,.005)
        assert np.max(abs(refined-trajectories[2,0])) < 1e-7
        train = trajectories[0,:, :,0].T
        for rho in (0.,.6):
            cov = SIGMA**2*((1-rho)*np.eye(4)+rho*np.ones((4,4)))
            white = np.linalg.inv(np.linalg.cholesky(cov)).T
            schedules = []
            for s,group in enumerate(trajectories):
                for phase,trace in enumerate(group):
                    schedules.append(('trace',s,phase,trace @ white))
                    schedules.append(('repeat',s,phase,trace[:,1:2]*2/SIGMA))
            rng = np.random.default_rng(20261127)
            # Select a single trace on the full development grid before any trial data.
            fixed_index = max((i for i,item in enumerate(schedules) if item[0]=='trace'),
                              key=lambda i: detector(schedules[i][3],POS & (np.arange(len(PARAMS))<GRID_SIZE),
                                                     NEG & (np.arange(len(PARAMS))<GRID_SIZE))[0])
            for trial in range(240):
                population = (np.arange(len(PARAMS)) < GRID_SIZE) == (trial < 120)
                truth = rng.choice(np.flatnonzero(population & (POS if trial%2 else NEG)))
                obs = train[truth]+rng.normal(0,SIGMA,8)
                admitted = (((train-obs)/SIGMA)**2).sum(1)<=chi2.ppf(.95,8)
                covered = int(admitted[truth]) if trial < 120 else -1
                admitted[GRID_SIZE:] = False
                pos,neg = admitted & POS, admitted & NEG
                common = [coupling,rho,trial,truth,covered,pos.sum(),neg.sum()]
                for policy in ('fixed_repeat','fixed_trace','development_trace','repeat','trace','choose'):
                    indices = [i for i,(kind,s,p,_) in enumerate(schedules)
                               if (policy=='development_trace' and i==fixed_index) or policy=='choose' or policy==kind or
                               (policy=='fixed_'+kind and s==2 and p==0)]
                    if pos.any() and neg.any():
                        options = [(detector(schedules[i][3],pos,neg),i) for i in indices]
                        (gap,scores),i = max(options,key=lambda o:o[0][0])
                        acc = ndtr(scores[truth]*(1 if POS[truth] else -1)) if gap>0 else .5
                        predicted = ndtr(max(gap,0))
                    else:
                        i = indices[0]
                        predicted = float('nan')
                        acc = float(pos.any()==POS[truth]) if pos.any() or neg.any() else .5
                    kind,s,p,_ = schedules[i]
                    rows.append(common+[policy,kind,s,p,predicted,acc])
    with (OUT/'trials.csv').open('w',newline='') as f:
        writer = csv.writer(f)
        writer.writerow('coupling rho trial truth covered npos nneg policy arm energy_index phase_index conditional_floor true_accuracy'.split())
        writer.writerows(rows)


if __name__ == '__main__':
    run()
