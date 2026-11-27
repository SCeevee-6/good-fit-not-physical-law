"""Protocol 0F, units 0F-1 (calibration, implementation checks) and 0F-2 (real series).

python src/run_0F.py calib [workers]   -> results/0F_jobs/calib/<sys>_n<ni>_s<seed>.json
python src/run_0F.py real  [workers]   -> results/0F_jobs/real/<name>.json
python src/run_0F.py horizon [workers] -> results/0F_jobs/horizon/<sys>_s<seed>.json  (0F-3)
Env OF_QUICK=1: development seeds 0-1 and hare-lynx plus one country (pipeline test only; outputs are discarded).
Resumable: a job whose JSON exists is skipped.
"""
import os, sys, json, time, csv, signal
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

RES = os.path.join(ROOT, 'results', '0F_jobs')
CALIB_SYS = ('A', 'B', 'E4')
CALIB_SEEDS = tuple(range(100, 200))
CHECK_SEEDS = tuple(range(100, 120))      # G-F0 sampling and region certificates, noise index 1 only
N_SAMPLE = 200
ALPHA = 0.05
GRID = dict(n1=501, n2=61, n3=25, n4=11)
JOB_LIMIT_S = 1800
QUICK = os.environ.get('OF_QUICK') == '1'
if QUICK:                                   # pipeline test on development seeds, never the gate sample
    CALIB_SEEDS = (0, 1); CHECK_SEEDS = (0, 1)


class Limit(Exception):
    pass


def _alarm(sig, frm):
    raise Limit()


def _setup():
    import dimtest  # noqa: F401  registers SYSTEMS['E4']


def zone_block(ds, specs, grids, tol_attr='tol_strong'):
    from lawmap import zone
    out = {}
    for sp in specs:
        X = grids[sp.name]
        m, h = ds.law(sp.c, X, getattr(sp, tol_attr)); z = zone(m, h)
        out[sp.name] = dict(n=len(X), certified=float((z == 1).mean()), refuted=float((z == -1).mean()),
                            undecided=float((z == 0).mean()), h_inf=float(np.isinf(h).mean()),
                            all_certified=bool((z == 1).all()))
    return out


def conf_set(trajs, S, E, sig):
    """Errors-in-variables GLS confidence set (lawmap.StackedSet), one feasible-GLS step from least squares."""
    from lawmap import stacked_set
    ds = stacked_set(trajs, S, E, sig, alpha=ALPHA)
    return ds, int(ds.Z.shape[0])


def calib_job(job):
    from systems import SYSTEMS
    from data import make_cell, stack
    from gate_a import NOISE
    from lawmap import DataSet, noise_sd, law_grid, sample_check, region_certificate
    _setup()
    key, ni, seed = job
    S = SYSTEMS[key]; t0 = time.time()
    cell = make_cell(S, NOISE[ni], seed, ni); trajs = cell['train']
    sig_hat = noise_sd(trajs, S); sig_true = np.asarray(cell['sigma'], float)
    Xtr, dXtr = stack(trajs)
    grids = {sp.name: law_grid(sp, **GRID) for sp in S.specs}
    rng = np.random.default_rng([seed, ni, 7])
    rec = dict(key=key, noise_idx=ni, seed=seed, sigma_true=sig_true.tolist(), sigma_hat=sig_hat.tolist(), degrees={})
    for d in S.degrees:
        E = S.E(d); xt = S.true_Xi(d); R = {}
        for tag, sig in (('est', sig_hat), ('known', sig_true)) if ni == 1 else (('est', sig_hat),):
            ds, nrows = conf_set(trajs, S, E, sig)
            stat = ds.r(xt) - ds.r_hat
            R[tag] = dict(eps=ds.eps, q=ds.rank, stat=float(stat), covered=bool(stat <= ds.eps),
                          rank=ds.rank, p=ds.p, cond=ds.cond, nrows=nrows)
            if tag == 'est':
                R['zones_conf'] = zone_block(ds, S.specs, grids)
                if ni == 1 and seed in CHECK_SEEDS:
                    chk, cert = {}, {}
                    for sp in S.specs:
                        bc, br = sample_check(ds, sp.c, grids[sp.name], sp.tol_strong, N_SAMPLE, rng)
                        chk[sp.name] = dict(bad_certified=bc, bad_refuted=br)
                        if R['zones_conf'][sp.name]['all_certified']:
                            t1 = time.time()
                            cert[sp.name] = dict(status=region_certificate(ds, sp, grids[sp.name], budget_s=30.0), sec=round(time.time() - t1, 2))
                    R['sample_check'] = chk; R['certificates'] = cert
        dsr = DataSet(Xtr, dXtr, cell['deriv_sd'], E)
        R['zones_rashomon'] = zone_block(dsr, S.specs, grids)
        rec['degrees'][str(d)] = R
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def real_datasets():
    import real_r, real_e
    S, cell, info = real_r.build(real_r.PRIMARY)
    yield 'hare_lynx', S, cell, os.path.join(ROOT, 'results', '0R_jobs', 'w7_c3.json'), False, None
    for p in (real_e.PANEL[:1] if QUICK else real_e.PANEL):
        S, cell, info = real_e.build(dict(slug=p['slug'], window=real_e.PRIMARY_W))
        yield f"covid_{p['slug']}", S, cell, os.path.join(ROOT, 'results', '0E_jobs', f"{p['slug']}_w{real_e.PRIMARY_W}.json"), True, p['slug']


def real_job(name):
    from data import stack
    from lawmap import DataSet, noise_sd, law_grid, plane_spec, region_certificate, zone
    from run_F0 import near_mask
    _setup()
    D = [x for x in real_datasets() if x[0] == name][0]
    _, S, cell, pool_path, plane, slug = D
    t0 = time.time()
    T = cell['train_whole']; trajs = [T]
    specs = [plane_spec(sp) for sp in S.specs] if plane else list(S.specs)
    grids = {sp.name: law_grid(sp, n1=2001, n2=201, n3=41, n4=15) for sp in specs}
    sig_hat = noise_sd(trajs, S); Xtr, dXtr = stack(cell['train'])
    pool = json.load(open(pool_path))
    rec = dict(name=name, plane=plane, sigma_hat=sig_hat.tolist(), degrees={}, witnesses=[])
    sets = {}
    for d in sorted({r['degree'] for r in pool['rows']}):
        E = S.E(d)
        ds, nrows = conf_set(trajs, S, E, sig_hat); sets[d] = ds
        dsr = DataSet(Xtr, dXtr, cell['deriv_sd'], E)
        R = dict(eps=ds.eps, rank=ds.rank, p=ds.p, cond=ds.cond, nrows=nrows,
                 zones_conf=zone_block(ds, specs, grids), zones_rashomon=zone_block(dsr, specs, grids), near={}, certificates={})
        for sp in specs:
            X = grids[sp.name]; nr = near_mask(X, Xtr, cell['state_sd'])
            m, h = ds.law(sp.c, X, sp.tol_strong); z = zone(m, h)
            R['near'][sp.name] = dict(n_near=int(nr.sum()), certified_near=float((z[nr] == 1).mean()) if nr.any() else None,
                                      undecided_near=float((z[nr] == 0).mean()) if nr.any() else None)
            if R['zones_conf'][sp.name]['all_certified']:
                R['certificates'][sp.name] = region_certificate(ds, sp, X, budget_s=60.0)
        rec['degrees'][str(d)] = R
    wfile = 'F0_witness_zones.csv'
    for w in csv.DictReader(open(os.path.join(ROOT, 'results', wfile))):
        if w['dataset'] != name:
            continue
        sp = [s_ for s_ in S.specs if s_.name == w['law']][0]
        x = np.array(json.loads(w['state']), float)[None, :]; ds = sets[int(w['degree'])]
        m, h = ds.law(sp.c, x, sp.tol_strong)
        rec['witnesses'].append(dict(cand=w['cand'], law=w['law'], degree=int(w['degree']), state=w['state'],
                                     on_plane=bool(abs(x.sum() - 1) < 1e-9) if plane else None, zone=int(zone(m, h)[0]),
                                     m=float(m[0]), h=float(h[0])))
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


HORIZON_SYS = ('A', 'B', 'H', 'E4')
HORIZON_SEEDS = tuple(range(100, 120))


def horizon_run(job):
    from systems import SYSTEMS
    from data import make_cell
    from gate_a import NOISE
    from lawmap import noise_sd, stacked_set
    from horizon import horizon_job
    _setup()
    key, seed = job; S = SYSTEMS[key]; t0 = time.time()
    cell = make_cell(S, NOISE[1], seed, 1); trajs = cell['train']; sig = noise_sd(trajs, S)
    rng = np.random.default_rng([seed, 31])
    rec = dict(key=key, seed=seed, degrees={})
    for d in S.degrees:
        E = S.E(d); ds = stacked_set(trajs, S, E, sig, alpha=ALPHA)
        rec['degrees'][str(d)] = horizon_job(ds, trajs, S, E, sig, ds.Xi_ref, S.specs, rng)
    rec['seconds'] = round(time.time() - t0, 1)
    return rec


def _worker(args):
    phase, job = args
    name = {'calib': lambda j: f'{j[0]}_n{j[1]}_s{j[2]}', 'real': lambda j: j, 'horizon': lambda j: f'{j[0]}_s{j[1]}'}[phase](job)
    path = os.path.join(RES, phase, name + '.json')
    if os.path.exists(path):
        return path, 'skip'
    if hasattr(signal, 'SIGALRM'):
        signal.signal(signal.SIGALRM, _alarm); signal.alarm(JOB_LIMIT_S)
    try:
        rec = {'calib': calib_job, 'real': real_job, 'horizon': horizon_run}[phase](job)
        st = 'ok'
    except Limit:
        rec, st = dict(job=list(job) if isinstance(job, tuple) else job, error='time limit'), 'limit'
    except Exception as e:                                  # recorded, not hidden
        rec, st = dict(job=list(job) if isinstance(job, tuple) else job, error=repr(e)[:500]), 'error'
    finally:
        if hasattr(signal, 'SIGALRM'):
            signal.alarm(0)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(rec, open(path, 'w'))
    return path, st


def jobs(phase):
    if phase == 'calib':
        seeds = CALIB_SEEDS
        # primary noise first, so anything lost to a time limit is a sensitivity job
        return [(k, ni, s) for ni in (1, 0, 2) for k in CALIB_SYS for s in seeds]
    if phase == 'horizon':
        return [(k, s) for k in HORIZON_SYS for s in ((0, 1) if QUICK else HORIZON_SEEDS)]
    return [x[0] for x in real_datasets()]


def main(phase, workers=4):
    import multiprocessing as mp
    js = [(phase, j) for j in jobs(phase)]
    t0 = time.time(); n_ok = n_bad = 0
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    log = open(os.path.join(ROOT, 'logs', f'0F_{phase}.log'), 'a')
    with mp.get_context('spawn').Pool(int(workers)) as pool:
        for i, (path, st) in enumerate(pool.imap_unordered(_worker, js)):
            n_ok += st in ('ok', 'skip'); n_bad += st not in ('ok', 'skip')
            log.write(f'{time.time() - t0:.0f}s {i + 1}/{len(js)} {os.path.basename(path)} {st}\n'); log.flush()
    print(phase, 'jobs', len(js), 'ok', n_ok, 'not ok', n_bad, f'{time.time() - t0:.0f}s', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], *(sys.argv[2:3] or [4]))
