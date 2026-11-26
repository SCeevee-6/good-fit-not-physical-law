"""Push tools/kaggle_0B.py as a private Kaggle script kernel, or fetch its output.
usage: python tools/kaggle_launch.py push <mode> <commit> <user>    (mode: main | needles | R | E | C | D)
       python tools/kaggle_launch.py status <mode> <user>
       python tools/kaggle_launch.py fetch <mode> <user> <outdir>
Needs the kaggle CLI and its token in the environment.
"""
import os, sys, json, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def slug(mode):
    return {'R': 'ceg-0r', 'E': 'ceg-0e', 'C': 'ceg-0c', 'D': 'ceg-0d'}.get(mode, f'ceg-0b-{mode}')


def kg(*a):
    r = subprocess.run([sys.executable, '-m', 'kaggle', *a], capture_output=True, text=True)
    print((r.stdout + r.stderr).strip()[-1500:])
    return r.returncode


def main():
    cmd, mode = sys.argv[1], sys.argv[2]
    if cmd == 'push':
        commit, user = sys.argv[3], sys.argv[4]
        d = tempfile.mkdtemp()
        if mode in ('C', 'D'):
            src = open(os.path.join(HERE, f'kaggle_0{mode}.py'), encoding='utf-8').read().replace("COMMIT = 'REPLACE'", f"COMMIT = '{commit}'")
            if len(sys.argv) > 5 and sys.argv[5] == 'test':
                src = src.replace('QUICK = False', 'QUICK = True')
        elif mode in ('R', 'E'):
            src = open(os.path.join(HERE, f'kaggle_0{mode}.py'), encoding='utf-8').read().replace("COMMIT = 'REPLACE'", f"COMMIT = '{commit}'")
            if len(sys.argv) > 5 and sys.argv[5] == 'test':
                ts = 'synthetic_H' if mode == 'R' else 'synthetic_B'
                src = src.replace("SERIES_ARGS = []", f"SERIES_ARGS = ['--series', '{ts}', '--only', '1']")
        else:
            src = open(os.path.join(HERE, 'kaggle_0B.py'), encoding='utf-8').read()
            src = src.replace("MODE = 'main'", f"MODE = '{mode}'").replace("COMMIT = '9f1cf045'", f"COMMIT = '{commit}'")
        open(os.path.join(d, 'run.py'), 'w', encoding='utf-8', newline='\n').write(src)
        meta = dict(id=f'{user}/{slug(mode)}', title=slug(mode), code_file='run.py', language='python', kernel_type='script',
                    is_private=True, enable_gpu=False, enable_tpu=False, enable_internet=True,
                    dataset_sources=[], competition_sources=[], kernel_sources=[])
        json.dump(meta, open(os.path.join(d, 'kernel-metadata.json'), 'w'))
        sys.exit(kg('kernels', 'push', '-p', d))
    user = sys.argv[3]
    if cmd == 'status':
        sys.exit(kg('kernels', 'status', f'{user}/{slug(mode)}'))
    if cmd == 'fetch':
        out = sys.argv[4]; os.makedirs(out, exist_ok=True)
        sys.exit(kg('kernels', 'output', f'{user}/{slug(mode)}', '-p', out))


if __name__ == '__main__':
    main()
