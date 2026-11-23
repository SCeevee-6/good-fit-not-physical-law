"""Run protocol 0B phases on a Kaggle CPU machine. The launcher (tools/kaggle_launch.py) sets MODE and COMMIT.
main:    pool -> sample -> arms        needles: needles
Code comes from the frozen commit on GitHub; packages are pinned to the versions of the local ceg-ode env.
Output: /kaggle/working/0B_<MODE>.zip (results/0B_jobs/..., results/0B_sample*.{csv,txt}, logs/0B_*.log).
"""
import os, sys, subprocess, time, zipfile, io, urllib.request, shutil
MODE = 'main'
COMMIT = '9f1cf045'
REPO = 'SCeevee-6/good-fit-not-physical-law'
PINS = ['numpy==2.4.6', 'scipy==1.17.1', 'z3-solver==5.1.0.0', 'cvxpy==1.9.3', 'clarabel==0.11.1', 'sympy==1.14.0']
T0 = time.time()


def sh(cmd, **kw):
    print('$', ' '.join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


lib = '/tmp/pylib'                       # venv has no ensurepip on the Kaggle image; pinned packages go first on the path
sh([sys.executable, '-m', 'pip', 'install', '-q', '--target', lib, *PINS])
os.environ['PYTHONPATH'] = lib
py = sys.executable
sh([py, '-c', 'import numpy, scipy, z3, cvxpy, clarabel, sympy, sys; print(sys.version, numpy.__version__, scipy.__version__, z3.get_version_string(), cvxpy.__version__)'])
url = f'https://api.github.com/repos/{REPO}/zipball/{COMMIT}'
with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'kaggle-run'}), timeout=300) as r:
    z = zipfile.ZipFile(io.BytesIO(r.read()))
z.extractall('/tmp/src')
root = os.path.join('/tmp/src', os.listdir('/tmp/src')[0])
n = os.cpu_count() or 4
print('cpus', n, 'root', root, flush=True)
phases = ['pool', 'sample', 'arms'] if MODE == 'main' else ['needles']
for ph in phases:
    t = time.time()
    sh([py, 'src/run_0B.py', ph, str(n)], cwd=root)
    print(f'phase {ph} {time.time() - t:.0f} s', flush=True)
out = f'/kaggle/working/0B_{MODE}.zip'
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED) as zf:
    for base in ('results/0B_jobs', 'logs'):
        for dp, _, fs in os.walk(os.path.join(root, base)):
            for f in fs:
                if base == 'logs' and not f.startswith('0B_'):
                    continue
                p = os.path.join(dp, f); zf.write(p, os.path.relpath(p, root))
    for f in ('results/0B_sample.csv', 'results/0B_sample_short_cells.txt'):
        if os.path.exists(os.path.join(root, f)):
            zf.write(os.path.join(root, f), f)
print('wrote', out, f'total {time.time() - T0:.0f} s', flush=True)
