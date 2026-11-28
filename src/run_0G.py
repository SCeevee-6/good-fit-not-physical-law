"""Protocol 0G (docs/protocol_0G.md): bias-corrected confidence set (src/eiv.py) for noisy regressors.
python src/run_0G.py calib [workers]  -> results/0G_jobs/calib/<sys>_n<ni>_s<seed>.json
python src/run_0G.py real  [workers]  -> results/0G_jobs/real/<name>.json   (report only)
Env OG_QUICK=1: development seeds 10-11, one country (pipeline test only; outputs are discarded).
Resumable: a job whose JSON exists is skipped. Helpers (grids, zone shares, real series) are imported from run_0F.
"""
import os, sys, json, time, signal
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
QUICK = os.environ.get('OG_QUICK') == '1'
if QUICK:
    os.environ['OF_QUICK'] = '1'
import run_0F as F                                   # noqa: E402
RES = os.path.join(ROOT, 'results', '0G_jobs')
CALIB_SYS = ('E4', 'A', 'B')
CALIB_SEEDS = (10, 11) if QUICK else tuple(range(100, 200))
CHECK_SEEDS = (10, 11) if QUICK else tuple(range(100, 120))
N_SAMPLE = 200
ALPHA = 0.05
JOB_LIMIT_S = 1800


def calib_job(job):
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from lawmap import noise_sd, law_grid, sample_check, region_certificate
    from eiv import CorrectedSet
    F._setup()
    key, ni, seed = job
    S = SYSTEMS[key]; t0 = time.time()
    cell = make_cell(S, NOISE[ni], seed, ni); trajs = cell['train']
    sig_hat = noise_sd(trajs, S); sig_true = np.asarray(cell['sigma'], float)
    grids = {sp.name: law_grid(sp, **F.GRID) for sp in S.specs} if ni == 1 else None
    rng = np.random.default_rng([seed, ni, 17])
    rec = dict(key=key, noise_idx=ni, seed=seed, sigma_true=sig_true.tolist(), sigma_hat=sig_hat.tolist(), degrees={})
    for d in S.degrees:
        E = S.E(d); xt = S.true_Xi(d); R = {}
        for tag, sig in (('est', sig_hat), ('known', sig_true)) if ni == 1 else (('est', sig_hat),):
            ds = CorrectedSet(trajs, S, E, sig, alpha=ALPHA)
            stat = ds.stat(xt)
            R[tag] = dict(eps=ds.eps, rank=ds.rank, stat=float(stat), covered=bool(stat <= ds.eps), p=ds.p,
                          n_neg=ds.n_neg, n_drop=ds.n_drop, min_reliability=ds.min_eig_ratio,
                          max_abs_err=float(np.abs(ds.Xi_hat - xt).max()))
            if tag == 'est' and ni == 1:
                R['zones_conf'] = F.zone_block(ds, S.specs, grids)
                if seed in CHECK_SEEDS:
                    chk, cert = {}, {}
                    for sp in S.specs:
                        bc, br = sample_check(ds, sp.c, grids[sp.name], sp.tol_strong, N_SAMPLE, rng)
                        chk[sp.name] = dict(bad_certified=bc, bad_refuted=br)
                        if R['zones_conf'][sp.name]['all_certified']:
                            t1 = time.time()
                            cert[sp.name] = dict(status=region_certificate(ds, sp, grids[sp.name], budget_s=30.0), sec=round(time.time() - t1, 2))
                    R['sample_check'] = chk; R['certificates'] = cert
        rec['degrees'][str(d)] = R
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def real_job(name):
    from lawmap import noise_sd, law_grid, plane_spec, region_certificate, zone
    from run_F0 import near_mask
    from data import stack
    from eiv import CorrectedSet
    F._setup()
    D = [x for x in F.real_datasets() if x[0] == name][0]
    _, S, cell, pool_path, plane, slug = D
    t0 = time.time(); T = cell['train_whole']; trajs = [T]
    specs = [plane_spec(sp) for sp in S.specs] if plane else list(S.specs)
    grids = {sp.name: law_grid(sp, n1=2001, n2=201, n3=41, n4=15) for sp in specs}
    sig_hat = noise_sd(trajs, S); Xtr, _ = stack(cell['train'])
    pool = json.load(open(pool_path))
    rec = dict(name=name, plane=plane, sigma_hat=sig_hat.tolist(), degrees={})
    for d in sorted({r['degree'] for r in pool['rows']}):
        E = S.E(d)
        ds = CorrectedSet(trajs, S, E, sig_hat, alpha=ALPHA)
        R = dict(eps=ds.eps, rank=ds.rank, p=ds.p, n_neg=ds.n_neg, n_drop=ds.n_drop, min_reliability=ds.min_eig_ratio,
                 zones_conf=F.zone_block(ds, specs, grids), near={}, certificates={})
        for sp in specs:
            X = grids[sp.name]; nr = near_mask(X, Xtr, cell['state_sd'])
            m, h = ds.law(sp.c, X, sp.tol_strong); z = zone(m, h)
            R['near'][sp.name] = dict(n_near=int(nr.sum()), certified_near=float((z[nr] == 1).mean()) if nr.any() else None,
                                      undecided_near=float((z[nr] == 0).mean()) if nr.any() else None)
            if R['zones_conf'][sp.name]['all_certified']:
                R['certificates'][sp.name] = region_certificate(ds, sp, X, budget_s=60.0)
        rec['degrees'][str(d)] = R
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def jobs(phase):
    if phase == 'calib':
        return [(k, ni, s) for ni in (1, 0, 2) for k in CALIB_SYS for s in CALIB_SEEDS]
    return [x[0] for x in F.real_datasets()]


def _worker(args):
    phase, job = args
    name = f'{job[0]}_n{job[1]}_s{job[2]}' if phase == 'calib' else job
    path = os.path.join(RES, phase, name + '.json')
    if os.path.exists(path):
        return path, 'skip'
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, F._alarm); signal.alarm(JOB_LIMIT_S)
    try:
        rec = (calib_job if phase == 'calib' else real_job)(job); st = 'ok'
    except F.Limit:
        rec, st = dict(job=list(job) if isinstance(job, tuple) else job, error='time limit'), 'limit'
    except Exception as e:                                  # recorded, not hidden
        rec, st = dict(job=list(job) if isinstance(job, tuple) else job, error=repr(e)[:500]), 'error'
    finally:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(rec, open(path, 'w'), default=F._np_json)
    return path, st


def main(phase, workers=4):
    import multiprocessing as mp
    js = [(phase, j) for j in jobs(phase)]
    t0 = time.time(); n_ok = n_bad = 0
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    log = open(os.path.join(ROOT, 'logs', f'0G_{phase}.log'), 'a')
    with mp.get_context('spawn').Pool(int(workers)) as pool:
        for i, (path, st) in enumerate(pool.imap_unordered(_worker, js)):
            n_ok += st in ('ok', 'skip'); n_bad += st not in ('ok', 'skip')
            log.write(f'{time.time() - t0:.0f}s {i + 1}/{len(js)} {os.path.basename(path)} {st}\n'); log.flush()
    print(phase, 'jobs', len(js), 'ok', n_ok, 'not ok', n_bad, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], *(sys.argv[2:3] or [4]))
