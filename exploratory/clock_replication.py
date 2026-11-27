"""Run with no arguments to download pinned data, or supply Cs, cable, GPS paths."""
import csv
import gzip
import hashlib
from pathlib import Path
import sys
import urllib.request

import numpy as np
from scipy.optimize import nnls
from scipy.special import ndtri
from scipy.stats import binom

from clock_transfer import REV, ROOT, calibration

FILES = [
    ('cesium','Cs5071A/5071A_phase.txt.gz','3092941352deef03c5c966fd9f3fcb5d24f0d465'),
    ('cable','Keysight53230A_ti_noise_floor/tic_phase.txt','da2f77c48a5c58ba03162c8793115c2cb201fc2d'),
    ('gps','gps/gps_1pps_phase_data.txt.gz','d1c29a45f00807cfd12327a0ca213dfc73c0bd31'),
]


def read(name, sha, local):
    url = f'https://raw.githubusercontent.com/aewallin/allantools/{REV}/tests/{name}'
    data = Path(local).read_bytes() if local else urllib.request.urlopen(url,timeout=60).read()
    if hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()!=sha:
        raise ValueError(name)
    if name.endswith('.gz'):
        data = gzip.decompress(data)
    series = np.loadtxt(data.splitlines())*1e9
    return [s[:len(s)//128*128].reshape(-1,128) for s in np.array_split(series,3)]


def tolerance(scores):
    n = len(scores)
    # P(at least k samples in the lowest population decile) >= .975.
    k = int(binom.ppf(.025,n,.1))
    if k<1:
        raise ValueError('insufficient calibration blocks')
    sorted_scores = np.sort(scores)
    lo,hi = sorted_scores[k-1],sorted_scores[-k]
    return (hi+lo)/2, np.nextafter((hi-lo)/2,np.inf)


def main(paths):
    if paths and len(paths)!=3:
        raise ValueError('supply all three paths or none')
    rows = []
    for i,(record,name,sha) in enumerate(FILES):
        train,cal,test = read(name,sha,paths[i] if paths else None)
        trend,eta2,sigma2 = calibration(train)
        for h in (4,16,64):
            X,C,Y = [np.diff(a[:,-h-1:],axis=1)-trend for a in (train,cal,test)]
            ramp = (np.arange(h)+.5)/h
            D = np.diff(np.eye(h+1),axis=0)
            cov = eta2*np.eye(h)+sigma2*D@D.T
            w0 = np.linalg.solve(cov,ramp)
            chol = np.linalg.cholesky(.2*cov+.8*np.cov(X,rowvar=False))
            w1 = nnls(chol.T,np.linalg.solve(chol,ramp))[0]
            for method,w in [('gaussian',w0),('tails',w0),('learned',w1)]:
                w = w/(w@ramp)
                assert w.min()>=-1e-10
                if method=='gaussian':
                    shift,center = ndtri(.9)*np.sqrt(w@cov@w),0.
                else:
                    center,shift = tolerance(C@w)
                errors = Y@w-center
                pos,neg = np.mean(errors>-shift),np.mean(errors<shift)
                quarters = [min(np.mean(q>-shift),np.mean(q<shift)) for q in np.array_split(errors,4)]
                rows.append([record,h,method,len(train),len(cal),len(test),shift,center,pos,neg,min(pos,neg),*quarters])
    with (ROOT/'clock_replication.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow('record h method train_n cal_n test_n m_ns_per_s threshold positive negative worst q1 q2 q3 q4'.split())
        writer.writerows(rows)

if __name__=='__main__':
    main(sys.argv[1:])
