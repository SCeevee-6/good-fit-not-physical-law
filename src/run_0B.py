"""Driver for protocol 0B. One JSON per job, resumable (finished JSONs are skipped).

phases: pool    -> results/0B_jobs/pool/<cell>.json          (100 cells)
        sample  -> results/0B_sample.csv                      (drawn once, never redrawn)
        arms    -> results/0B_jobs/arms/<sample_id>[__<struct>].json
        needles -> results/0B_jobs/needles/<needle_id>.json
usage (repo root): python src/run_0B.py <phase> [n_workers]      log: logs/0B_<phase>.log
"""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, 'results', '0B_jobs')
SAMPLE = os.path.join(ROOT, 'results', '0B_sample.csv')


def _write(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, default=float)
    os.replace(tmp, path)


def w_pool(c):
    import gate_b
    p = os.path.join(OUT, 'pool', gate_b.cell_id(*c) + '.json')
    try:
        _write(p, gate_b.pool_job(*c)); return gate_b.cell_id(*c), None
    except Exception as e:
        return gate_b.cell_id(*c), repr(e)


def w_arm(job):
    import gate_b
    struct = job.pop('struct', '')
    p = os.path.join(OUT, 'arms', job['sample_id'] + (f'__{struct}' if struct else '') + '.json')
    try:
        recs, err = gate_b.arm_job(job, struct), None
    except Exception as e:
        recs, err = [], repr(e)
    _write(p, dict(sample_id=job['sample_id'], struct=struct, error=err, records=recs))
    return os.path.basename(p), err


def w_needle(job):
    import gate_b
    p = os.path.join(OUT, 'needles', job['needle_id'] + '.json')
    try:
        rec, err = gate_b.needle_job(job), None
    except Exception as e:
        rec, err = dict(needle_id=job['needle_id'], status='error'), repr(e)
    rec['error'] = err
    _write(p, rec)
    return job['needle_id'], err


def todo(phase):
    import gate_b
    if phase == 'pool':
        J = gate_b.cells()
        return [c for c in J if not os.path.exists(os.path.join(OUT, 'pool', gate_b.cell_id(*c) + '.json'))], len(J), w_pool
    if phase == 'arms':
        rows = list(csv.DictReader(open(SAMPLE)))
        J = []
        for r in rows:
            for s in (('Bpos', 'Bfac') if r['key'] == 'B' else ('',)):
                J.append(dict(r, struct=s))
        order = {'B': 0, 'H': 1, 'Amis': 2, 'A': 3}
        J.sort(key=lambda j: order[j['key']])
        f = lambda j: os.path.join(OUT, 'arms', j['sample_id'] + (f"__{j['struct']}" if j['struct'] else '') + '.json')
        return [j for j in J if not os.path.exists(f(j))], len(J), w_arm
    if phase == 'needles':
        J = gate_b.needle_jobs(gate_b.needle_bases(ROOT))
        J.sort(key=lambda j: -len([s for s in [j['spec']] if s.startswith('B2') or s.startswith('B3')]))
        return [j for j in J if not os.path.exists(os.path.join(OUT, 'needles', j['needle_id'] + '.json'))], len(J), w_needle
    raise ValueError(phase)


def sample():
    import gate_b
    if os.path.exists(SAMPLE):
        print('sample exists, not redrawn'); return
    rows, short = gate_b.draw_sample(os.path.join(OUT, 'pool'))
    with open(SAMPLE, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    with open(os.path.join(ROOT, 'results', '0B_sample_short_cells.txt'), 'w') as f:
        f.write('\n'.join(short) + ('\n' if short else 'none\n'))
    print(f'sample {len(rows)} rows, cells with < 4 equivalent candidates: {len(short)}')


def main(phase, n=5):
    if phase == 'sample':
        return sample()
    import multiprocessing as mp
    os.makedirs(os.path.join(OUT, phase), exist_ok=True)
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    J, N, fn = todo(phase)
    log = os.path.join(ROOT, 'logs', f'0B_{phase}.log')
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
