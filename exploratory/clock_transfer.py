"""Usage: python exploratory/clock_transfer.py [5071A_phase.txt.gz]."""
import csv
import gzip
import hashlib
import sys
import urllib.request
from pathlib import Path
import numpy as np
from scipy.optimize import nnls
from scipy.special import ndtr, ndtri

REV = '32ad8fce6be9f51dda4ab407931f71503c50d101'
URL = f'https://raw.githubusercontent.com/aewallin/allantools/{REV}/tests/Cs5071A/5071A_phase.txt.gz'
BLOB = '3092941352deef03c5c966fd9f3fcb5d24f0d465'
ROOT = Path(__file__).resolve().parent


def load_blocks(path=None):
    raw = Path(path).read_bytes() if path else urllib.request.urlopen(URL, timeout=60).read()
    git_header = f'blob {len(raw)}\0'.encode()
    if hashlib.sha1(git_header+raw).hexdigest()!=BLOB:
        raise ValueError('unexpected upstream data')
    y = np.loadtxt(gzip.decompress(raw).splitlines())*1e9  # nanoseconds
    return [a[:len(a)//128*128].reshape(-1,128) for a in np.array_split(y,3)]


def calibration(train):
    trend = np.median((train[:,-1]-train[:,0])/127)
    lags = np.array([1,2,4,8,16,32])
    variances = []
    for lag in lags:
        d = train[:,lag:]-train[:,:-lag]
        variances.append((np.median(abs(d-np.median(d)))/ndtri(.75))**2)
    A = np.column_stack((lags,2*np.ones(len(lags))))/np.array(variances)[:,None]
    eta2,sigma2 = nnls(A,np.ones(len(lags)))[0]
    return trend,eta2,sigma2


def accuracy(noise, shift):
    return np.array([np.mean(noise+shift>0),np.mean(noise-shift<0)])


def main(path=None):
    train,cal,test = load_blocks(path)
    trend,eta2,sigma2 = calibration(train)
    print('blocks',len(train),len(cal),len(test),'trend, eta², sigma² (ns units)',trend,eta2,sigma2)
    rows = []
    for h in (4,16,64):
        diff = np.diff(np.eye(h+1),axis=0)
        cov = eta2*np.eye(h)+sigma2*diff@diff.T
        ramp = (np.arange(h)+.5)/h
        matched = np.linalg.solve(cov,ramp)
        m = ndtri(.9)/np.sqrt(matched@ramp)
        for name,w in [('endpoint',np.ones(h)),('triangle',ramp),('mixed',matched)]:
            shift = m*(w@ramp)
            forecast = ndtr(shift/np.sqrt(w@cov@w))
            scores = [(np.diff(x[:,-h-1:],axis=1)-trend)@w for x in (cal,test)]
            ref,held = [accuracy(s,shift) for s in scores]
            quarters = [min(accuracy(s,shift)) for s in np.array_split(scores[1],4)]
            rows.append([h,name,m,forecast,min(ref),*held,held.mean(),min(held),*quarters])
    with (ROOT/'clock_transfer.csv').open('w',newline='') as out:
        writer = csv.writer(out)
        writer.writerow('h method m_ns_per_s model comparator_worst positive negative balanced worst first second third fourth'.split())
        writer.writerows(rows)

if __name__=='__main__':
    main(*sys.argv[1:])
