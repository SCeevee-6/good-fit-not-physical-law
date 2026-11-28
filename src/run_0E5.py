"""Protocol 0E5 (docs/protocol_0E5.md): five-state testbed E5 (src/e5.py).
python src/run_0E5.py <calib|design> [workers] -> results/0E5_jobs/<phase>/<name>.json
Env OE5_QUICK=1: development seeds 10-11 for calib, seed 0 with a cap of 2 for design (pipeline test only).
calib: the 0G corrected set (src/eiv.py) at degrees 2 and 3, noise index 1, seeds 100-199, estimated and true sigma;
zones, sample checks (seeds 100-119) and region certificates as in 0G. design: the 0F-4 comparison (src/design.py)
on E5 at degree 3 for E2S and E2R, seeds 200-219 (the 0F-4 seeds), cap 30 added trajectories.
"""
import os, sys, json, time, signal
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
QUICK = os.environ.get('OE5_QUICK') == '1'
import run_0G as G0                                       # noqa: E402  (calib_job, reused with E5)
import run_0F as F                                        # noqa: E402
RES = os.path.join(ROOT, 'results', '0E5_jobs')
CALIB_SEEDS = (10, 11) if QUICK else tuple(range(100, 200))
G0.CHECK_SEEDS = (10, 11) if QUICK else tuple(range(100, 120))
LIMITS = {'calib': 1800, 'design': 6 * 3600}


def calib(job):
    import e5  # noqa: F401
    return G0.calib_job(job)


def design(job):
    import e5  # noqa: F401
    import design as D
    F._setup()
    if QUICK:
        D.CAP = 2
    return D.design_job(*job)


def jobs(phase):
    if phase == 'calib':
        return [('E5', 1, s) for s in CALIB_SEEDS]
    import design as D
    return [('E5', s) for s in ((0,) if QUICK else D.SEEDS)]


def _worker(args):
    phase, job = args
    name = f'{job[0]}_n{job[1]}_s{job[2]}' if phase == 'calib' else f'{job[0]}_s{job[1]}'
    path = os.path.join(RES, phase, name + '.json')
    if os.path.exists(path):
        return path, 'skip'
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, F._alarm); signal.alarm(LIMITS[phase])
    try:
        rec = (calib if phase == 'calib' else design)(job); st = 'ok'
    except F.Limit:
        rec, st = dict(job=list(job), error='time limit'), 'limit'
    except Exception as e:                                  # recorded, not hidden
        rec, st = dict(job=list(job), error=repr(e)[:500]), 'error'
    finally:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
    os.makedirs(os.path.dirname(path), exist_ok=True); json.dump(rec, open(path, 'w'), default=F._np_json)
    return path, st


def main(phase, workers=4):
    import multiprocessing as mp
    js = [(phase, j) for j in jobs(phase)]; t0 = time.time(); bad = 0
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True); log = open(os.path.join(ROOT, 'logs', f'0E5_{phase}.log'), 'a')
    with mp.get_context('spawn').Pool(int(workers)) as pool:
        for i, (p, st) in enumerate(pool.imap_unordered(_worker, js)):
            bad += st not in ('ok', 'skip'); log.write(f'{time.time() - t0:.0f}s {i + 1}/{len(js)} {os.path.basename(p)} {st}\n'); log.flush()
    print(phase, 'jobs', len(js), 'not ok', bad, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], *(sys.argv[2:3] or [4]))
