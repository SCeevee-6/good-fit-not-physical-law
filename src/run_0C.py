"""Driver for protocol 0C. One JSON per job, resumable.
phases: pool    -> results/0C_jobs/pool/<cell>.json     (160 cells, noise 0 and 2)
        sample  -> results/0C_sample.csv                (drawn once)
        arms    -> results/0C_jobs/arms/<sample_id>.json   (S, C, Rnd at noise 0 and 2)
        soft    -> results/0C_jobs/soft/<sample_id>.json   (Soft1-3 on the 0B sample, noise 1)
usage (repo root): python src/run_0C.py <phase> [n_workers]      log: logs/0C_<phase>.log
"""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, 'results', '0C_jobs')
SAMPLE = os.path.join(ROOT, 'results', '0C_sample.csv')
SAMPLE_0B = os.path.join(ROOT, 'results', '0B_sample.csv')


def _write(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, default=float)
    os.replace(tmp, path)


def w_pool(c):
    import robust
    nm = robust.cell_id(*c)
    try:
        _write(os.path.join(OUT, 'pool', nm + '.json'), robust.pool_job(*c)); return nm, None
    except Exception as e:
        return nm, repr(e)


def w_job(a):
    phase, job = a
    import robust
    p = os.path.join(OUT, phase, job['sample_id'] + '.json')
    try:
        recs, err = (robust.arm_job(job) if phase == 'arms' else robust.soft_job(job)), None
    except Exception as e:
        import traceback
        recs, err = [], traceback.format_exc()[-800:]
    _write(p, dict(sample_id=job['sample_id'], error=err, records=recs))
    return job['sample_id'], err


def todo(phase):
    import robust
    if phase == 'pool':
        J = robust.cells()
        return [c for c in J if not os.path.exists(os.path.join(OUT, 'pool', robust.cell_id(*c) + '.json'))], len(J), w_pool
    if phase == 'arms':
        J = list(csv.DictReader(open(SAMPLE)))
        order = {'B': 0, 'H': 1, 'Amis': 2, 'A': 3}
        J.sort(key=lambda j: order[j['key']])
    else:
        J = [r for r in csv.DictReader(open(SAMPLE_0B)) if r['key'] in robust.KEYS and r['noise_idx'] == '1']
        order = {'B': 0, 'H': 1, 'Amis': 2, 'A': 3}
        J.sort(key=lambda j: order[j['key']])
        if os.environ.get('OC_QUICK'):
            J = [next(j for j in J if j['key'] == k) for k in robust.KEYS]
    return [(phase, j) for j in J if not os.path.exists(os.path.join(OUT, phase, j['sample_id'] + '.json'))], len(J), w_job


def sample():
    import robust
    if os.path.exists(SAMPLE):
        print('sample exists, not redrawn'); return
    rows, short = robust.draw_sample(os.path.join(OUT, 'pool'))
    with open(SAMPLE, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(ROOT, 'results', '0C_sample_short_cells.txt'), 'w') as f:
        f.write('\n'.join(short) + ('\n' if short else 'none\n'))
    print(f'sample {len(rows)} rows, cells with < 4 equivalent candidates: {len(short)}')


def main(phase, n=5):
    if phase == 'sample':
        return sample()
    import multiprocessing as mp
    os.makedirs(os.path.join(OUT, phase), exist_ok=True); os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    J, N, fn = todo(phase)
    log = os.path.join(ROOT, 'logs', f'0C_{phase}.log')
    with open(log, 'a') as lg:
        lg.write(f"{time.strftime('%H:%M:%S')} start {phase}: {len(J)} to run of {N}\n"); lg.flush()
        done, nerr = N - len(J), 0
        with mp.get_context('spawn').Pool(n) as pool:
            for jid, err in pool.imap_unordered(fn, J):
                done += 1; nerr += err is not None
                lg.write(f"{time.strftime('%H:%M:%S')} {jid} {'ERR ' + err if err else 'done'} ({done}/{N})\n"); lg.flush()
        lg.write(f"{time.strftime('%H:%M:%S')} finished {phase}, errors {nerr}\n")
    print(f'{phase}: {done}/{N} done, errors {nerr}')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5)
