"""Driver for protocol 0A. One JSON per sample in results/0A_jobs/ (resumable); log in logs/gate_0A_run.log."""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(ROOT, 'results', '0A_jobs'); LOG = os.path.join(ROOT, 'logs', 'gate_0A_run.log')


def work(job):
    import gate_a
    p = os.path.join(OUT, job['sample_id'] + '.json')
    try:
        recs = gate_a.run_sample(job)
        err = None
    except Exception as e:
        recs, err = [], repr(e)
    tmp = p + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(dict(sample_id=job['sample_id'], error=err, records=recs), f, default=float)
    os.replace(tmp, p)
    return job['sample_id'], err


def jobs():
    with open(os.path.join(ROOT, 'results', '0A_sample.csv')) as f:
        J = list(csv.DictReader(f))
    order = {'B': 0, 'Amis': 1, 'A': 2}          # slowest first
    J.sort(key=lambda j: order[j['key']])
    return [j for j in J if not os.path.exists(os.path.join(OUT, j['sample_id'] + '.json'))], len(J)


def main(n=6):
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    todo, N = jobs()
    with open(LOG, 'a') as lg:
        lg.write(f"{time.strftime('%H:%M:%S')} start: {len(todo)} to run of {N}\n"); lg.flush()
        done = N - len(todo)
        with mp.get_context('spawn').Pool(n) as pool:
            for sid, err in pool.imap_unordered(work, todo):
                done += 1
                lg.write(f"{time.strftime('%H:%M:%S')} {sid} {'ERR '+err if err else 'done'} ({done}/{N})\n"); lg.flush()
        lg.write(f"{time.strftime('%H:%M:%S')} finished\n")
    return done
