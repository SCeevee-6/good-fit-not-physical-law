"""Continuous-time lower and upper bounds for minimax local sign recovery."""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import lsq_linear
from scipy.special import ndtr

T, M, L = 2., .35, 1.


def area(u, v, dt):
    # Exact absolute integral of a linear spline, including its zero crossings.
    same = u*v >= 0
    trap = dt*(abs(u)+abs(v))/2
    cross = dt*(u*u+v*v)/(2*np.maximum(abs(u)+abs(v),1e-300))
    return np.where(same,trap,cross).sum()


def bounds(times, sigma, cells=64):
    times = np.array(times)
    durations = times+T
    # For p(t)=p(0)+integral_t^0 d(s) ds, |d|<=L.
    breaks = np.unique(np.r_[-T,times,0.])
    nodes = np.unique(np.concatenate([np.linspace(a,b,cells+1)
                          for a,b in zip(breaks[:-1],breaks[1:])]))
    kernel = np.minimum(nodes[:,None]+T,durations)
    A = (kernel[:-1]+kernel[1:]).T*np.diff(nodes)/2
    mat = np.column_stack((durations,A))
    fit = lsq_linear(mat,-M*durations,bounds=(np.r_[0.,np.full(A.shape[1],-L)],
                           np.r_[np.inf,np.full(A.shape[1],L)]),method='bvls',tol=1e-12)
    assert fit.success
    mu = mat@fit.x+M*durations
    upper = np.linalg.norm(mu)/sigma
    if upper < 1e-9:
        return .5,ndtr(upper),np.zeros(len(times))
    w = mu/np.linalg.norm(mu)
    curve = np.minimum(breaks[:,None]+T,durations)@w
    gap = M*(durations@w)-L*area(curve[:-1],curve[1:],np.diff(breaks))
    lower = ndtr(max(0,gap/sigma)) if w@durations>=0 else .5
    assert lower <= ndtr(upper)+1e-10
    return lower,ndtr(upper),w


def main():
    designs = [(f'pair_{h:.2f}',[-h,0]) for h in np.arange(.05,1.01,.05)]
    best = max(designs,key=lambda d:bounds(d[1],.0125)[0])
    designs += [('endpoint',[0]),('early_only',best[1][:1]),
                ('uniform_four',[-1.5,-1.,-.5,0]),
                ('nested_four',[-1.5,-1.,*best[1]])]
    with Path(__file__).with_name('local_power.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['design','sigma','times','accuracy_lower','accuracy_upper','weights'])
        for name,times in designs:
            for noise in (.0125,.05):
                lo,hi,w = bounds(times,noise)
                writer.writerow([name,noise,str([float(t) for t in times]),lo,hi,w.tolist()])
                if name==best[0] or not name.startswith('pair'):
                    print(name,noise,round(lo,6),round(hi,6))
    h = -best[1][0]
    print('plain difference and short gate',*[ndtr((M*h-L*h*h/2)/(.0125*k)) for k in (np.sqrt(2),1.)])

if __name__=='__main__':
    main()
