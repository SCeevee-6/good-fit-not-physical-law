"""Protocol 0H (docs/protocol_0H.md): is A's under-coverage at 0.5 % noise caused by smoothing bias?
python src/run_0H.py [workers]  -> results/0H_jobs/A_n<ni>_s<seed>.json
Env OH_QUICK=1: development seeds 30-31 (pipeline test only; outputs are discarded).
Arms per degree: base (0G corrected set), oracle (filter bias delta = D x_clean - f(H x_clean) removed from dX),
w15 and w31 (Savitzky-Golay window 15 or 31 instead of 21), known (true sigma). Every arm reports stat, eps and
coverage. The base arm also reports the bias noncentrality nu = c' V^-1 c, where c is the shift of the identified
coordinates caused by delta alone.
"""
import os, sys, json, time, signal, dataclasses
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
QUICK = os.environ.get('OH_QUICK') == '1'
RES = os.path.join(ROOT, 'results', '0H_jobs')
SEEDS = (30, 31) if QUICK else tuple(range(200, 300))
NOISE_IDX = (0, 1)
WINDOWS = (15, 31)
JOB_LIMIT_S = 1800


def resmooth(trajs, S):
    from data import smooth
    out = []
    for tr in trajs:
        Xs, dX, sl = smooth(tr['X'], S)
        out.append(dict(tr, Xs=Xs, dX=dX, sl=sl))
    return out


def filter_bias(trajs, S, E, xt):
    """delta per trajectory (interior rows, n states): D x_clean - Theta(H x_clean) xi_true."""
    from lawmap import sg_mats
    from lib import theta
    out = []
    for tr in trajs:
        Xc = np.asarray(tr['X_clean'], float); H, D, sl = sg_mats(len(Xc), S)
        out.append(D @ Xc - theta(H @ Xc, E) @ xt)
    return out


def noncentrality(cs, trajs, E, delta):
    from lib import theta
    n, pk = cs.n, cs.p; v = np.zeros(n * pk)
    for tr, dl in zip(trajs, delta):
        Th = theta(tr['Xs'], E)[:, cs.keep] / cs.col
        for j in range(n):
            v[j * pk:(j + 1) * pk] += Th.T @ dl[:, j]
    c = np.kron(np.eye(n), cs.U.T / cs.lam[:, None]) @ v
    return float(c @ cs.Vr_inv @ c)


def job(args):
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from lawmap import noise_sd
    from eiv import CorrectedSet
    ni, seed = args
    S = SYSTEMS['A']; t0 = time.time()
    cell = make_cell(S, NOISE[ni], seed, ni); trajs = cell['train']
    sig_true = np.asarray(cell['sigma'], float); sig_hat = noise_sd(trajs, S)
    rec = dict(key='A', noise_idx=ni, seed=seed, sigma_true=sig_true.tolist(), sigma_hat=sig_hat.tolist(), degrees={})
    for d in S.degrees:
        E = S.E(d); xt = S.true_Xi(d); R = {}
        def put(tag, cs):
            st = cs.stat(xt)
            R[tag] = dict(stat=float(st), eps=cs.eps, rank=cs.rank, covered=bool(st <= cs.eps), n_neg=cs.n_neg)
            return cs
        base = put('base', CorrectedSet(trajs, S, E, sig_hat))
        delta = filter_bias(trajs, S, E, xt)
        R['base']['bias_noncentrality'] = noncentrality(base, trajs, E, delta)
        orc = [dict(tr, dX=tr['dX'] - dl) for tr, dl in zip(trajs, delta)]
        put('oracle', CorrectedSet(orc, S, E, sig_hat))
        put('known', CorrectedSet(trajs, S, E, sig_true))
        for w in WINDOWS:
            Sw = dataclasses.replace(S, sg_window=w); tw = resmooth(trajs, Sw)
            put(f'w{w}', CorrectedSet(tw, Sw, E, noise_sd(tw, Sw)))
        rec['degrees'][str(d)] = R
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def _alarm(sig, frm):
    raise TimeoutError('time limit')


def _worker(a):
    path = os.path.join(RES, f'A_n{a[0]}_s{a[1]}.json')
    if os.path.exists(path):
        return path, 'skip'
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, _alarm); signal.alarm(JOB_LIMIT_S)
    try:
        rec, st = job(a), 'ok'
    except Exception as e:                                  # recorded, not hidden
        rec, st = dict(job=list(a), error=repr(e)[:500]), 'error'
    finally:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
    os.makedirs(RES, exist_ok=True)
    json.dump(rec, open(path, 'w'))
    return path, st


def jobs():
    return [(ni, s) for ni in NOISE_IDX for s in SEEDS]


def main(workers=4):
    import multiprocessing as mp
    js = jobs(); t0 = time.time(); bad = 0
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    log = open(os.path.join(ROOT, 'logs', '0H.log'), 'a')
    with mp.get_context('spawn').Pool(int(workers)) as pool:
        for i, (p, st) in enumerate(pool.imap_unordered(_worker, js)):
            bad += st not in ('ok', 'skip'); log.write(f'{time.time() - t0:.0f}s {i + 1}/{len(js)} {os.path.basename(p)} {st}\n'); log.flush()
    print('jobs', len(js), 'not ok', bad, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main(*(sys.argv[1:2] or [4]))
