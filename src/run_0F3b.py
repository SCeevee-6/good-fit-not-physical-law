"""Post-hoc unit R2 (docs/protocol_0F3b.md): 0F-3 horizon rays with the Theorem 4 projected bound.
python src/run_0F3b.py [workers] -> results/0F3b_jobs/<sys>_s<seed>.json. Env OF3B_QUICK=1: seeds 0-1 (pipeline test).
"""
import os, sys, json, time, signal
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import run_0F as F                                  # noqa: E402
QUICK = os.environ.get('OF3B_QUICK') == '1'
RES = os.path.join(ROOT, 'results', '0F3b_jobs')
SEEDS = (0, 1) if QUICK else F.HORIZON_SEEDS


def job(a):
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from lawmap import noise_sd, stacked_set
    from horizon2 import horizon_job2
    F._setup()
    key, seed = a; S = SYSTEMS[key]; t0 = time.time()
    cell = make_cell(S, NOISE[1], seed, 1); trajs = cell['train']; sig = noise_sd(trajs, S)
    rng = np.random.default_rng([seed, 31])
    rec = dict(key=key, seed=seed, degrees={})
    for d in S.degrees:
        E = S.E(d); ds = stacked_set(trajs, S, E, sig, alpha=F.ALPHA)
        rec['degrees'][str(d)] = horizon_job2(ds, trajs, S, E, sig, ds.Xi_ref, S.specs, rng)
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def _worker(a):
    path = os.path.join(RES, f'{a[0]}_s{a[1]}.json')
    if os.path.exists(path):
        return path, 'skip'
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, F._alarm); signal.alarm(F.JOB_LIMIT_S)
    try:
        rec, st = job(a), 'ok'
    except Exception as e:                                  # recorded, not hidden
        rec, st = dict(job=list(a), error=repr(e)[:500]), 'error'
    finally:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
    os.makedirs(RES, exist_ok=True); json.dump(rec, open(path, 'w'), default=F._np_json)
    return path, st


def jobs():
    return [(k, s) for k in F.HORIZON_SYS for s in SEEDS]


def main(workers=4):
    import multiprocessing as mp
    js = jobs(); t0 = time.time(); bad = 0
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True); log = open(os.path.join(ROOT, 'logs', '0F3b.log'), 'a')
    with mp.get_context('spawn').Pool(int(workers)) as pool:
        for i, (p, st) in enumerate(pool.imap_unordered(_worker, js)):
            bad += st not in ('ok', 'skip'); log.write(f'{time.time() - t0:.0f}s {i + 1}/{len(js)} {os.path.basename(p)} {st}\n'); log.flush()
    print('jobs', len(js), 'not ok', bad, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main(*(sys.argv[1:2] or [4]))
