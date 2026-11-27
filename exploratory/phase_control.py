"""Post hoc control: four half-unit resets split between two start phases."""
import csv
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.special import ndtr
from crossover import readings
from study import GRID_SIZE, NEG, PARAMS, PHASES, POS, SIGMA, meter


def bank(params, k, extra, term="xv"):
    return np.column_stack([readings(params,1.5,p,k,extra,term)[:,1] for p in PHASES])


def main():
    rng = np.random.default_rng(1738)
    truth = rng.uniform([-.22,-.1,-.35],[-.18,.8,.35],(2500,3))
    margin = truth[:,0]+.5*truth[:,1]
    truth, sign = truth[abs(margin)>=.08], np.sign(margin[abs(margin)>=.08])
    rows = []
    for robust in (False,True):
        variants = (-.25,0,.25) if robust else (0.,)
        dev = np.vstack([bank(PARAMS[:GRID_SIZE],.3,d) for d in variants])
        pos,neg = [np.tile(z[:GRID_SIZE],len(variants)) for z in (POS,NEG)]
        for mode in ('single','pair'):
            best = None
            options = [(i,) for i in range(8)] if mode=='single' else combinations(range(8),2)
            for idx in options:
                x = dev[:,idx]
                w = x[pos].mean(0)-x[neg].mean(0)
                w /= np.linalg.norm(w)*SIGMA*np.sqrt(len(idx)/4)
                lo,hi = (x@w)[pos].min(),(x@w)[neg].max()
                if best is None or lo-hi>best[0]:
                    best = lo-hi,idx,w,(hi+lo)/2
            gap,idx,w,mid = best
            for k in (.1,.5):
                for term,d in [(term,d) for term in ("xv","x4","radial") for d in (-.35,.35)]+[("xv",0.)]:
                    score = (bank(truth,k,d,term)[:,idx]@w-mid)*sign
                    acc = ndtr(score)
                    rows.append([robust,mode,idx,k,term,d,ndtr(gap/2),
                                 min(acc[sign==i].mean() for i in (-1,1)),acc.min()])
    with Path(__file__).with_name('phase_control.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow('augmented_design schedule phase_indices coupling term extra_gain design_floor worst_sign_mean worst_sample'.split())
        writer.writerows(rows)


if __name__=='__main__':
    main()
