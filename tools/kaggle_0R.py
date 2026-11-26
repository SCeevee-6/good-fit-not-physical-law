"""Run protocol 0R on a Kaggle CPU machine. The launcher (tools/kaggle_launch.py) sets COMMIT.
Code comes from the frozen commit on GitHub; packages are pinned to the versions of the local ceg-ode env.
Output: /kaggle/working/0R.zip (results/0R_jobs, logs/0R_run.log).
"""
import os, sys, subprocess, time, zipfile, io, urllib.request
COMMIT = 'REPLACE'
REPO = 'SCeevee-6/good-fit-not-physical-law'
PINS = ['numpy==2.4.6', 'scipy==1.17.1', 'z3-solver==5.1.0.0', 'cvxpy==1.9.3', 'clarabel==0.11.1', 'sympy==1.14.0']
T0 = time.time()
SERIES_ARGS = []        # the launcher can set ['--series', 'synthetic_H', '--only', '1'] for a test run


def sh(cmd, **kw):
    print('$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


lib = '/tmp/pylib'
sh([sys.executable, '-m', 'pip', 'install', '-q', '--target', lib, *PINS])
os.environ['PYTHONPATH'] = lib
with urllib.request.urlopen(urllib.request.Request(f'https://api.github.com/repos/{REPO}/zipball/{COMMIT}', headers={'User-Agent': 'kaggle-run'}), timeout=300) as r:
    zipfile.ZipFile(io.BytesIO(r.read())).extractall('/tmp/src')
root = os.path.join('/tmp/src', os.listdir('/tmp/src')[0])
n = os.cpu_count() or 4
print('cpus', n, 'root', root, flush=True)
sh([sys.executable, '-c', 'import sys; sys.path.insert(0, "src"); import real_r; print("alarm self-test", real_r.selftest_limit())'], cwd=root)
try:
    sh([sys.executable, 'src/run_0R.py', str(n)] + SERIES_ARGS, cwd=root)
except subprocess.CalledProcessError as e:
    print('run_0R failed', e, flush=True)
finally:
    print('run_0R ended', f'{time.time() - T0:.0f} s', flush=True)
out = '/kaggle/working/0R.zip'
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for dp, _, fs in os.walk(os.path.join(root, 'results', '0R_jobs')):
        for f in fs:
            p = os.path.join(dp, f); zf.write(p, os.path.relpath(p, root))
    for f in os.listdir(os.path.join(root, 'logs')):
        if f.startswith('0R_'):
            zf.write(os.path.join(root, 'logs', f), 'logs/' + f)
print('wrote', out, f'total {time.time() - T0:.0f} s', flush=True)
