"""Run protocol S1 (Navier-Stokes Galerkin closure) on a Kaggle CPU machine. The launcher (tools/kaggle_launch.py) sets COMMIT (and QUICK for a pipeline test).
Phases: calib -> real. Code comes from the frozen commit on GitHub, packages pinned as in 0B.
Output: /kaggle/working/S1.zip (results/S1_run.json, results/S1_samples.npz).
"""
import os, sys, subprocess, time, zipfile, io, urllib.request
COMMIT = 'REPLACE'
QUICK = False
PHASES = ('all',)
REPO = 'SCeevee-6/good-fit-not-physical-law'
PINS = ['numpy==2.4.6', 'scipy==1.17.1', 'z3-solver==5.1.0.0', 'cvxpy==1.9.3', 'clarabel==0.11.1', 'sympy==1.14.0']
T0 = time.time()


def sh(cmd, **kw):
    print('$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


lib = '/tmp/pylib'
sh([sys.executable, '-m', 'pip', 'install', '-q', '--target', lib, *PINS])
os.environ['PYTHONPATH'] = lib
if QUICK:
    os.environ['S1_QUICK'] = '1'
with urllib.request.urlopen(urllib.request.Request(f'https://api.github.com/repos/{REPO}/zipball/{COMMIT}', headers={'User-Agent': 'kaggle-run'}), timeout=300) as r:
    zipfile.ZipFile(io.BytesIO(r.read())).extractall('/tmp/src')
root = os.path.join('/tmp/src', os.listdir('/tmp/src')[0])
n = os.cpu_count() or 4
print('cpus', n, 'root', root, flush=True)
try:
    for ph in PHASES:
        t = time.time()
        sh([sys.executable, 'src/run_S1.py'], cwd=root)
        print(f'phase {ph} {time.time() - t:.0f} s', flush=True)
except subprocess.CalledProcessError as e:
    print('run_S1 failed', e, flush=True)
finally:
    out = '/kaggle/working/S1.zip'
    with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
        for f in ('S1_run.json', 'S1_samples.npz'):
            p = os.path.join(root, 'results', f)
            if os.path.exists(p):
                zf.write(p, 'results/' + f)
    print('wrote', out, f'total {time.time() - T0:.0f} s', flush=True)
