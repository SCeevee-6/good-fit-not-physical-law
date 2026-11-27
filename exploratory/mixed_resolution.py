"""Run from the repository root with Python, numpy and scipy."""
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
from scipy.special import ndtr

M,L,T,SIGMA = .35,1.,2.,.0125
H = M/L


def cut(n, eta):
    a = 4*L*SIGMA**2/(n*eta**2)
    return M*M/(M+a+np.sqrt(a*(2*M+a)))


def universal(n, eta, c):
    drift_energy = (c**3/3+c*c*(M-c))/L
    residual_area = (M-c)**2/(2*L)
    return np.sqrt(drift_energy/eta**2+n*(residual_area/(2*SIGMA))**2)


def linear_score(times, weights, eta):
    mean = (M*times+L*times**2/2)@weights
    covariance = eta**2*np.minimum.outer(T+times,T+times)+SIGMA**2*np.eye(len(times))
    return mean/np.sqrt(weights@covariance@weights)


def spaced(n,eta):
    times = np.linspace(-T,0,n)
    dt = np.diff(times)
    D = np.diff(np.eye(n),axis=0)
    cov = np.diag(eta**2*dt)+SIGMA**2*(D@D.T)
    chol = np.linalg.cholesky(cov)
    increments = dt*(M+L*(times[:-1]+times[1:])/2)
    # Best nonnegative increment weights for this grid.
    w,_ = nnls(chol.T,np.linalg.solve(chol,increments))
    return max(0,linear_score(times,D.T@w,eta)) if w.any() else 0.


def main():
    rows = []
    for eta in (.025,.075,.15):
        for n in (2,4,8,16,64):
            c = cut(n,eta)
            k = n//2
            starts = -H+c/L*(np.arange(k)+.5)/k
            times = np.r_[starts,np.zeros(k)]
            weights = np.r_[-np.ones(k),np.ones(k)]/k
            z = linear_score(times,weights,eta)
            bound = universal(n,eta,c)
            assert z<=bound+1e-10
            base = linear_score(np.r_[np.full(k,-H),np.zeros(k)],weights,eta)
            rows.append([eta,n,c,ndtr(z),ndtr(bound),ndtr(base),ndtr(spaced(n,eta)),
                         ndtr(np.sqrt(M**3/(3*L))/eta)])
            if n==8 and eta==.075:
                rng = np.random.default_rng(82349)
                # Simulate the entire read vector on the positive extremal ramp.
                cov = eta**2*np.minimum.outer(T+times,T+times)+SIGMA**2*np.eye(n)
                readings = rng.multivariate_normal(M*times+L*times**2/2,cov,100000)
                print('Monte Carlo',np.mean(readings@weights>0), 'analytic',ndtr(z))
    with Path(__file__).with_name('mixed_resolution.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow('eta n clip schedule upper two_times uniform_nonnegative drift_ceiling'.split())
        writer.writerows(rows)

if __name__=='__main__':
    main()
