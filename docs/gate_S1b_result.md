# Protocol S1b result: where the data decide the energy law for a closure

Run locally 2026-11-28 from the S1 snapshots, after S1b was frozen at 7b87805d. `src/s1b.py` writes `results/S1b_zones.csv`, `results/S1b_gate.csv` and `figures/fig_S1b.png`.

## Rules

| rule | K = 1 | K = 2 | K = 3 |
|---|---|---|---|
| sampled closures contradicting the zone map (must be 0) | 0 | 0 | 0 |
| region certificate "every closure dissipates everywhere" | impossible: centre indefinite (lambda_max 0.36) | impossible (1.24) | impossible (5.02) |
| observed states: refuted share (predicted >= 1/2) | **0.252, did not hold** | 0.754, held | 1.000, held |
| random shell directions: undecided share (predicted >= 1/2) | 1.000, held | 0.996, held | 0.998, held |

Certified shares are 0 at observed states and at most 0.00015 on the shell. Lag sensitivity (Newey-West lag 20 and 100 against 50) moves the refuted share at observed states by at most 0.05 (K = 2: 0.808 / 0.745 against 0.754) and leaves the picture unchanged.

## What it says

- At K = 3 the data refute the dissipation law at every observed resolved state: every closure in the 95 % set injects energy there. At K = 2 this holds at three quarters of the observed states.
- Away from the observed states (random directions on the energy shell) the data decide almost nothing, as in every other system in the project. The set is wide in most of the 2304 coefficient directions at K = 3.
- At K = 1 only a quarter of the observed states are refuted. With 8 resolved modes the linear closure explains 46 % of the flux variance (S1), so the residual is large and the set is wide. The prediction failed there for that reason.
- No region certificate is possible, because the best-fit closure already breaks the law. Here the law is refuted on the data rather than certified anywhere, and the exact checker reports that directly.

Together with S1: for this truncation the correct statement is not "the closure is dissipative" but "the data show the subgrid scales feeding the resolved ones wherever the flow actually goes". The data are silent elsewhere.

## Limits

- The probability wording is nominal. There is no true linear closure to check coverage against, and the snapshots are serially correlated (Newey-West, lag chosen in advance).
- One run, one seed, 2D, linear closures.
