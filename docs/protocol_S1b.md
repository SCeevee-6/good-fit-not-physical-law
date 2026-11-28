# Protocol S1b: zones of the energy law over the data-consistent set of linear closures

Written 2026-11-28 after the S1 result (`docs/S1_result.md`) and before this analysis runs. Frozen at the commit that adds this file. It reuses the S1 snapshots (`results/S1_samples.npz`, 1000 snapshots at K = 1, 2, 3) and runs locally in seconds.

## Question

S1 showed that the least-squares linear closure is not energy-stable at any K, because the subgrid scales are the only energy source for the resolved modes. In the language of the main project, the next question is where the data decide this. At which resolved states does every closure consistent with the data inject energy (the dissipation law refuted), where does every one dissipate (certified), and where can the data not tell?

## Set and zones

Code: `src/s1b.py`.

- Regression of the subgrid flux on the resolved state: b_t = C a_t + e_t.
- Covariance of vec(C_hat): Newey-West sandwich with Bartlett weights and lag 50 snapshots (5 time units), because snapshots 0.1 apart are serially correlated. Lags 20 and 100 are reported as sensitivity.
- Set: chi2 ellipsoid at 95 % on the identified directions. The law at a is a' C a <= 0 (the closure dissipates there). It is linear in vec(C) with weights a (x) a, so zones follow ROADMAP_v3 items 2-3.
- Zones are evaluated at 20,000 random unit directions on the energy shell and at the 1000 observed resolved states (normalised).
- Region certificate for "every closure in the set dissipates at every state": lambda_max(sym C_hat) + sqrt(eps lambda_max(V)) <= 0, checked in rationals. It is sufficient, since |a' (C - C_hat) a| <= sqrt(eps (a (x) a)' V (a (x) a)) <= sqrt(eps lambda_max(V)) for |a| = 1.

## Rules

| rule | kind | pass rule |
|---|---|---|
| sampled closures contradicting the zone map | check | 0, from 200 closures on the set boundary at 2000 sphere directions and all data directions, per K |
| region certificate | report | expected impossible: S1 found lambda_max(sym C_hat) > 0 at every K, so the centre already violates the law |
| zone shares at observed states | prediction | refuted >= 1/2 at every K (S1: the true flux injects in 83-100 % of snapshots) |
| zone shares at random sphere directions | prediction | undecided >= 1/2 at every K (most directions are far from the data) |

## Limits

- No calibration check is possible. There is no true linear closure, since the real subgrid flux is not linear in a, so the probability wording of the set is nominal.
- One run, one seed, lag chosen in advance without tuning.
