"""Driver for protocol 0E (with addendum 1). Each config runs in its own subprocess, killed after CFG_LIMIT_S.
One JSON per config in <out>/ (resumable). A killed config is written as {timed_out: true} and reported, not retried.
usage (repo root): python src/run_0E.py [n_workers] [--series synthetic_B] [--out DIR] [--only N] [--cfg-limit SEC] [--total-limit SEC]
log: logs/0E_run.log, per-config progress: logs/0E_progress_<name>.log
"""
import os, sys, json, time, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def name(cfg):
    return f"{cfg['slug']}_w{cfg['window']}"


def _arg(args, key, default):
    return args[args.index(key) + 1] if key in args else default


def _write(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, default=float)
    os.replace(tmp, path)


def one(args):
    import real_e
    nm = _arg(args, '--one', None)
    cfg = [c for c in real_e.CONFIGS if name(c) == nm][0]
    out = _arg(args, '--out', os.path.join(ROOT, 'results', '0E_jobs'))
    series = _arg(args, '--series', None)
    try:
        rec, err = real_e.run_config(cfg, series, nm), None
    except Exception:
        import traceback
        rec, err = dict(cfg=cfg), traceback.format_exc()[-1500:]
    rec['error'] = err
    _write(os.path.join(out, nm + '.json'), rec)


def main():
    import real_e
    args = sys.argv[1:]
    if '--one' in args:
        return one(args)
    series = _arg(args, '--series', None)
    out = _arg(args, '--out', os.path.join(ROOT, 'results', '0E_jobs'))
    n = int(args[0]) if args and args[0].isdigit() else 3
    only = int(_arg(args, '--only', len(real_e.CONFIGS)))
    cfg_limit = float(_arg(args, '--cfg-limit', 6000)); total_limit = float(_arg(args, '--total-limit', 23400))
    os.makedirs(out, exist_ok=True); os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    todo = [c for c in real_e.CONFIGS[:only] if not os.path.exists(os.path.join(out, name(c) + '.json'))]
    t_start = time.time(); running = {}; nerr = nkill = 0
    lg = open(os.path.join(ROOT, 'logs', '0E_run.log'), 'a')
    def log(m):
        lg.write(f"{time.strftime('%H:%M:%S')} {m}\n"); lg.flush()
    log(f'start: {len(todo)} to run of {len(real_e.CONFIGS)}, series {series}, cfg limit {cfg_limit:.0f} s, total {total_limit:.0f} s')
    while todo or running:
        out_of_time = time.time() - t_start > total_limit
        while todo and len(running) < n and not out_of_time:
            c = todo.pop(0)
            fh = open(os.path.join(ROOT, 'logs', f'0E_stdout_{name(c)}.log'), 'w')
            pr = subprocess.Popen([sys.executable, os.path.abspath(__file__), '--one', name(c), '--series', series if series else '', '--out', out],
                                  stdout=fh, stderr=subprocess.STDOUT, cwd=ROOT)
            running[name(c)] = (pr, time.time(), fh, c); log(f'launched {name(c)}')
        if out_of_time:
            for c in todo:
                _write(os.path.join(out, name(c) + '.json'), dict(cfg=c, error=None, timed_out=True, reason='total time limit before launch'))
                log(f'{name(c)} not launched (total limit)')
            todo = []
        for nm, (pr, t0, fh, c) in list(running.items()):
            rc = pr.poll()
            over = time.time() - t0 > cfg_limit or (out_of_time and time.time() - t_start > total_limit + 600)
            if rc is None and over:
                pr.kill(); pr.wait(); rc = 'killed'; nkill += 1
                tail = ''
                pg = os.path.join(ROOT, 'logs', f'0E_progress_{nm}.log')
                if os.path.exists(pg):
                    tail = ''.join(open(pg).readlines()[-5:])
                _write(os.path.join(out, nm + '.json'), dict(cfg=c, error=None, timed_out=True, reason=f'config limit {cfg_limit:.0f} s', progress_tail=tail))
            if rc is not None:
                fh.close(); del running[nm]
                ok = os.path.exists(os.path.join(out, nm + '.json'))
                if rc not in (0, 'killed') or not ok:
                    nerr += 1
                    if not ok:
                        _write(os.path.join(out, nm + '.json'), dict(cfg=c, error=f'worker exit {rc}'))
                log(f'{nm} {"killed at limit" if rc == "killed" else ("done" if rc == 0 else f"exit {rc}")} ({time.time() - t0:.0f} s)')
        time.sleep(2)
    log(f'finished, errors {nerr}, killed {nkill}'); lg.close()
    print(f'0R: errors {nerr}, killed {nkill}')


if __name__ == '__main__':
    main()
