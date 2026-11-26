"""Driver for protocol 0D. phases: pool -> results/0D_jobs/pool/<cell>.json (20 cells), sample -> results/0D_sample.csv,
arms -> results/0D_jobs/arms/<sample_id>.json (S, C, Rnd, SoftA, SoftB). Resumable.
usage (repo root): python src/run_0D.py <phase> [n_workers]      log: logs/0D_<phase>.log
"""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, 'results', '0D_jobs')
SAMPLE = os.path.join(ROOT, 'results', '0D_sample.csv')


def _write(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, default=float)
    os.replace(tmp, path)


def w_pool(seed):
    import dimtest
    nm = dimtest.cell_id(seed)
    try:
        _write(os.path.join(OUT, 'pool', nm + '.json'), dimtest.pool_job(seed)); return nm, None
    except Exception:
        import traceback
        return nm, traceback.format_exc()[-600:]


def w_arm(job):
    import dimtest
    p = os.path.join(OUT, 'arms', job['sample_id'] + '.json')
    try:
        recs, err = dimtest.arm_job(job), None
    except Exception:
        import traceback
        recs, err = [], traceback.format_exc()[-800:]
    _write(p, dict(sample_id=job['sample_id'], error=err, records=recs))
    return job['sample_id'], err


def todo(phase):
    import dimtest
    if phase == 'pool':
        J = list(dimtest.SEEDS)
        return [s for s in J if not os.path.exists(os.path.join(OUT, 'pool', dimtest.cell_id(s) + '.json'))], len(J), w_pool
    J = list(csv.DictReader(open(SAMPLE)))
    if os.environ.get('OD_QUICK'):
        J = J[:2]
    return [j for j in J if not os.path.exists(os.path.join(OUT, 'arms', j['sample_id'] + '.json'))], len(J), w_arm


def sample():
    import dimtest
    if os.path.exists(SAMPLE):
        print('sample exists, not redrawn'); return
    rows, short = dimtest.draw_sample(os.path.join(OUT, 'pool'))
    with open(SAMPLE, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(ROOT, 'results', '0D_sample_short_cells.txt'), 'w') as f:
        f.write('\n'.join(short) + ('\n' if short else 'none\n'))
    print(f'sample {len(rows)} rows, cells with < 4 equivalent candidates: {len(short)}')


def main(phase, n=5):
    if phase == 'sample':
        return sample()
    import multiprocessing as mp
    os.makedirs(os.path.join(OUT, phase), exist_ok=True); os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    J, N, fn = todo(phase)
    if phase == 'pool' and os.environ.get('OD_QUICK'):
        J = J[:2]; N = 2
    with open(os.path.join(ROOT, 'logs', f'0D_{phase}.log'), 'a') as lg:
        lg.write(f"{time.strftime('%H:%M:%S')} start {phase}: {len(J)} to run of {N}\n"); lg.flush()
        done, nerr = N - len(J), 0
        with mp.get_context('spawn').Pool(min(n, max(1, len(J)))) as pool:
            for jid, err in pool.imap_unordered(fn, J):
                done += 1; nerr += err is not None
                lg.write(f"{time.strftime('%H:%M:%S')} {jid} {'ERR ' + err if err else 'done'} ({done}/{N})\n"); lg.flush()
        lg.write(f"{time.strftime('%H:%M:%S')} finished {phase}, errors {nerr}\n")
    print(f'{phase}: {done}/{N} done, errors {nerr}')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 5)
