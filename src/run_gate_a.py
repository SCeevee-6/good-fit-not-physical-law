"""Driver for protocols 0A and 0A2. One JSON per job (resumable).
0A : results/0A_jobs/<sample_id>.json, log logs/gate_0A_run.log (frozen behaviour, default).
0A2: results/0A2_jobs/<struct>/<sample_id>.json, B only, struct in (Bpos, Bfac); log logs/gate_0A2_run.log."""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, 'results', '0A_jobs'); LOG = os.path.join(ROOT, 'logs', 'gate_0A_run.log')
OUT2 = os.path.join(ROOT, 'results', '0A2_jobs'); LOG2 = os.path.join(ROOT, 'logs', 'gate_0A2_run.log')


def _path(job):
    if job.get('struct'):
        return os.path.join(OUT2, job['struct'], job['sample_id'] + '.json')
    return os.path.join(OUT, job['sample_id'] + '.json')


def work(job):
    p = _path(job)
    try:
        if job.get('struct'):
            import gate_a2
            recs = gate_a2.run_sample2({k: v for k, v in job.items() if k != 'struct'}, job['struct'])
        else:
            import gate_a
            recs = gate_a.run_sample(job)
        err = None
    except Exception as e:
        recs, err = [], repr(e)
    tmp = p + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(dict(sample_id=job['sample_id'], struct=job.get('struct', ''), error=err, records=recs), f, default=float)
    os.replace(tmp, p)
    return (job.get('struct', '') + '/' + job['sample_id']).lstrip('/'), err


def jobs(protocol='0A'):
    with open(os.path.join(ROOT, 'results', '0A_sample.csv')) as f:
        J = list(csv.DictReader(f))
    if protocol == '0A2':                        # B only, one job per (candidate, structure)
        J = [dict(j, struct=s) for s in ('Bpos', 'Bfac') for j in J if j['key'] == 'B']
    else:
        order = {'B': 0, 'Amis': 1, 'A': 2}      # slowest first
        J.sort(key=lambda j: order[j['key']])
    return [j for j in J if not os.path.exists(_path(j))], len(J)


def main(n=6, protocol='0A'):
    """protocol '0A' (default, frozen behaviour) or '0A2' (docs/protocol_0A2.md)."""
    import multiprocessing as mp
    out, log = (OUT2, LOG2) if protocol == '0A2' else (OUT, LOG)
    for s in (('Bpos', 'Bfac') if protocol == '0A2' else ('',)):
        os.makedirs(os.path.join(out, s), exist_ok=True)
    os.makedirs(os.path.dirname(log), exist_ok=True)
    todo, N = jobs(protocol)
    with open(log, 'a') as lg:
        lg.write(f"{time.strftime('%H:%M:%S')} start {protocol}: {len(todo)} to run of {N}\n"); lg.flush()
        done = N - len(todo)
        with mp.get_context('spawn').Pool(n) as pool:
            for sid, err in pool.imap_unordered(work, todo):
                done += 1
                lg.write(f"{time.strftime('%H:%M:%S')} {sid} {'ERR '+err if err else 'done'} ({done}/{N})\n"); lg.flush()
        lg.write(f"{time.strftime('%H:%M:%S')} finished\n")
    return done
