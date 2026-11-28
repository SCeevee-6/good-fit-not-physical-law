# Protocol S1 result: energy stability of linear closures in a Galerkin truncation

Run on Kaggle 2026-11-28 from commit 46f63f2d (protocol `docs/protocol_S1.md`, frozen 2026-11-08, Addendum 1 written before the run), 286 s for one DNS. Outputs: `results/S1_run.json`, `results/S1_samples.npz`; `src/analyze_S1.py` writes `results/S1_gate.csv`, `results/S1_closure.csv`, `results/S1_Q{1,2,3}.npy` and `figures/fig_S1.png`.

## Gates

| rule | value | threshold | result |
|---|---|---|---|
| S-G1 energy budget residual (whole run; worst 1-unit window) | 2.4e-10; 1.0e-7 | < 1e-3 | pass |
| S-G2 truncated advection energy rate (worst relative, 1000 states per K) | 5.8e-16 | < 1e-10 | pass |
| S-G3 exact rational test and eigenvalues agree | agree at K = 1, 2, 3 | agree | pass |

The DNS is trustworthy as a data source, and the truncated advection does no work on the resolved energy, so any energy change in the truncated model comes from viscosity or the closure.

## Outcome

| K | resolved modes | fit R^2 | Q energy-stable? | lambda_max(Q) | eigenvalues > 0 | mean closure energy rate | mean resolved viscous rate |
|---|---|---|---|---|---|---|---|
| 1 | 8 | 0.46 | no | 0.359 | 5 of 8 | +0.00153 | -0.00099 |
| 2 | 24 | 0.96 | no | 1.24 | 14 of 24 | +0.00447 | -0.00424 |
| 3 | 48 | 0.999 | no | 5.02 | 29 of 48 | +0.01095 | -0.01082 |

No fitted linear closure is energy-stable at any truncation, as Addendum 1 predicted. The reason is physical, not a fitting failure. The forcing sits at |k| = 4, outside every truncation, so the resolved modes have no direct energy source. The true subgrid flux puts energy into them: its rate is positive in 82.6 % of snapshots at K = 1 and in all of them at K = 2 and 3. At K = 3 it balances resolved viscosity to about 1 %. A closure that is negative semidefinite in energy would remove the only source and starve the resolved flow. The least-squares closure reproduces the mean energy rate of the true flux exactly (a property of the fit), and it is indefinite because of that.

So for this flow, "the closure is dissipative" is the wrong law to demand. The law that holds is the budget, closure input = resolved dissipation (+ the change in resolved energy). The exact test answers the question it was built for, and the answer is no at every K.

## Limits

- The K = 1 modes are not statistically steady over the sampling window. Their energy grows from 0.010 to 0.064, at a rate (5.3e-4) that matches closure input minus viscous loss (5.4e-4). This is the large-scale growth typical of forced 2D flow. At K = 2 and 3 the drift is smaller (2.3e-4 and 1.3e-4).
- The linear closure explains only 46 % of the flux variance at K = 1, against 96 % and 99.9 % at K = 2 and 3.
- One run, one seed, one viscosity, linear closures only, and 2D. No statement about 3D or about singularities.
