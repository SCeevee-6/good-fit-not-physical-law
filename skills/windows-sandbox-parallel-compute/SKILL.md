---
name: windows-sandbox-parallel-compute
description: "What works and what is blocked for parallel/long-running Python on the native Windows sandbox (application control blocks joblib, sklearn, Parquet; subprocess pipes fail), with working patterns for pools, background cells, polling and killing a stuck kernel."
---

# windows-sandbox-parallel-compute

## blocked here
- `import joblib` and `import sklearn` (asyncio `_overlapped` DLL denied). Write AUROC with `scipy.stats.rankdata`.
- `pyarrow` Parquet (`_parquet` DLL denied). Use CSV.
- `python.exe` is not on PATH in the PowerShell tool and fails with DLL errors if called directly. Use the Python tool.
- Children started with `subprocess` from a kernel cannot create named pipes, so `multiprocessing` inside them fails. Plain scripts (no pools) are fine.

## works
- `multiprocessing.get_context('spawn').Pool(n)` created inside the kernel cell, workers in an importable module on `sys.path` (not defined in the cell). Set `OMP/OPENBLAS/MKL_NUM_THREADS=1` at the top of the worker module.
- Foreground cells are capped at 600 s. Use `background=true`, make jobs resumable (write one pickle per job, skip if present), log to a file.
- Poll the log from another environment (a PowerShell cell in the `r` environment) with `Start-Sleep -Seconds 540` then `Get-Content log -Tail 3`. Subprocesses are not tracked by `wait_for_notification`.
- Pickles with Arrow-backed strings do not load in a kernel that imported pandas before pyarrow was installed. Merge in a fresh interpreter (subprocess) and write CSV.

## stuck kernel
A cell stuck in native code (z3, cvxpy) ignores `exec_interrupt`. Find the python process with large CPU time (`Get-Process python`) and `Stop-Process -Id`; the kernel restarts on the next call and in-memory state is gone. Keep everything important in files.

## budget (pilot)
45 grid cells on 7 workers took 1h51m. SIR cells at 5% noise took up to 55 min each (3-D queries, large models). Put the slowest jobs first in the queue.
