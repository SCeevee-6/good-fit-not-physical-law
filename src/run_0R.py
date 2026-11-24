"""Driver for protocol 0R. One JSON per config in <out>/ (resumable). Config 0 is the primary, the rest are sensitivity runs.
usage (repo root): python src/run_0R.py [n_workers] [--series PATH|synthetic_H] [--out DIR]     log: logs/0R_run.log
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def name(cfg):
    return f"w{cfg['window']}_c{cfg['c']}"


def work(a):
    cfg, series, out = a
    import real_r
    p = os.path.join(out, name(cfg) + '.json')
    try:
        rec, err = real_r.run_config(cfg, series), None
    except Exception as e:
        import traceback
        rec, err = dict(cfg=cfg), traceback.format_exc()[-1500:]
    rec['error'] = err
    tmp = p + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(rec, f, default=float)
    os.replace(tmp, p)
    return name(cfg), err


def main():
    import multiprocessing as mp, real_r
    args = sys.argv[1:]
    series = args[args.index('--series') + 1] if '--series' in args else real_r.SERIES
    out = args[args.index('--out') + 1] if '--out' in args else os.path.join(ROOT, 'results', '0R_jobs')
    n = int(args[0]) if args and args[0].isdigit() else 3
    os.makedirs(out, exist_ok=True); os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    only = int(args[args.index('--only') + 1]) if '--only' in args else len(real_r.CONFIGS)
    todo = [(c, series, out) for c in real_r.CONFIGS[:only] if not os.path.exists(os.path.join(out, name(c) + '.json'))]
    with open(os.path.join(ROOT, 'logs', '0R_run.log'), 'a') as lg:
        lg.write(f"{time.strftime('%H:%M:%S')} start: {len(todo)} to run of {len(real_r.CONFIGS)}, series {series}\n"); lg.flush()
        nerr = 0
        with mp.get_context('spawn').Pool(min(n, max(1, len(todo)))) as pool:
            for nm, err in pool.imap_unordered(work, todo):
                nerr += err is not None
                lg.write(f"{time.strftime('%H:%M:%S')} {nm} {'ERR ' + err if err else 'done'}\n"); lg.flush()
        lg.write(f"{time.strftime('%H:%M:%S')} finished, errors {nerr}\n")
    print(f'0R: {len(todo)} run, errors {nerr}')


if __name__ == '__main__':
    main()
