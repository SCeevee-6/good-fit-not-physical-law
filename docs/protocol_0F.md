# Protocol 0F: law zones over a confidence set

Written 2026-11-26, frozen at the commit that adds this file, before any run that counts. It covers 0F-1 (implementation checks and calibration), 0F-2 (real series), 0F-3 (horizon) and 0F-4 (measurement design). The rules for all four are fixed here. Code for 0F-3 and 0F-4 is written after this freeze and before their runs, under these rules; any change to a rule is a dated addendum made before the run it affects. Background: `docs/ROADMAP_v3.md` (items 1-7, Addendum 1) and `docs/F0_result.md`.

## What was seen before the freeze

All of this is listed so the reader can judge what is blind.

1. F0 (Rashomon set, all stored fits, both real series): `docs/F0_result.md`.
2. Smoke test of an earlier version of the set (per-equation GLS, no regressor noise) on calibration seed 100, systems A and B. The truth was outside the set: stat/q 3.2-4.4. That version is not used.
3. Development seeds 0-9 and 7, which are outside every gate sample, used to find the cause. Smoothing bias is negligible (r(true)/q <= 0.024 on noise-free data). With noise-free regressors stat/q is 1.08-1.37. With the errors-in-variables set below, A and B give median stat/q 0.77-1.12 with 8-9 of 10 covered. E4 does not calibrate (stat/q 3.4 at degree 2, 436 at degree 3). With the true Jacobian it is still 2.8 and 4.2, so the cause is bias in the estimate, not the covariance. I expect G-F1 to fail on E4.
4. Smoke test of the 0F-2 code on hare-lynx and Turkey. I saw the zone shares (hare-lynx degree 2: undecided 0.62, refuted 0.18; Turkey degree 2: undecided 0.46, refuted 0.10, B1R region-certified). The local pipeline test of the full runner (development seeds 0-1, outputs deleted) also ran the real phase on hare-lynx and Argentina, the first country in the panel (Argentina, lowest degree: mean undecided 0.11). The G-F2 prediction below was written after this and is not blind for hare-lynx, Turkey or Argentina.

## Set

Model per trajectory, linearised in the iid Gaussian state noise e_k (sd sigma_k):

dX_j - f_j(Xs) = sum_k B_jk e_k,  B_jk = delta_jk D - diag(J_jk) H

- D and H are the Savitzky-Golay derivative and smoothing matrices, interior rows (`lawmap.sg_mats`).
- J is the Jacobian of the least-squares fit at Xs. There is one feasible-GLS step.
- All equations are stacked. The covariance is B diag(sigma^2) B'. The fit is GLS after Cholesky whitening (`lawmap.StackedSet`).
- Confidence set: r(xi) - r(xi_hat) <= chi2_{q, 0.95}, with q the rank of the whitened design. It is simultaneous over every law and every point.
- sigma_k is estimated from the smoothing residuals with the exact hat-matrix correction (`lawmap.noise_sd`). "Known" sigma (synthetic only) is reported next to it.
- Rashomon set as in F0, reported for comparison.

Zones, as in ROADMAP_v3 item 3:
- certified: m >= h;
- refuted: m + h < 0;
- undecided: otherwise.

m is the law value of the GLS fit plus the strong tolerance, and h = sqrt(eps K).

Region certificate (`lawmap.region_certificate`): M = g + tol >= 0 and M^2 - eps K >= 0 on the law domain. Both are checked exactly with z3, or Bernstein for 3-variable domains (`bern_nd` for 4), after rounding to the 1e-6 grid with a proven slack subtracted. It is attempted only where the grid says certified everywhere. The certificate is exact for the computed set; the set itself comes from floating-point linear algebra.

Law domains: as in `systems.py`. For COVID they are restricted to S + I + R = 1 (`lawmap.plane_spec`), per Addendum 1. Grids for 0F-1 use 501 / 61 / 25 / 11 points per axis for 1-4 query variables. For 0F-2 they use 2001 / 201 / 41 / 15.

## 0F-1: implementation checks and calibration

Jobs:
- systems A, B, E4 (truth in the library);
- noise index 1 (2 %) primary, 0 and 2 reported;
- seeds 100-199, every library degree of the system;
- `python src/run_0F.py calib`.

G-F0, implementation check. At noise index 1 and seeds 100-119, for every system, degree and law, 200 models are sampled on the set's boundary.
- The count of sampled models that break the law at a certified grid point, plus the count that obey it at a refuted point, must be 0 in total. Pass: total = 0. Fail: bug, stop.
- Region certificates are reported as attempted, verified, counterexample and timeout. A counterexample means the grid missed a point. It is reported and is not a failure.

G-F1, calibration. Coverage is the share of seeds with r(xi_true) - r(xi_hat) <= eps, using estimated sigma at noise index 1, over 100 seeds.
- Decided separately for each (system, degree). Pass: coverage >= 0.90.
- Probability wording ("with 95 % confidence") is used only for (system, degree) cells that pass.
- If every cell of a system fails, that system's zones are reported as a set-membership statement without a probability.
- Known-sigma coverage and the other noise levels are reported, not gated.

## 0F-2: real series

- Hare-lynx at the primary 0R config, and the 10 COVID countries at the primary 0E window, on the plane.
- Training data only. The set is built from the unsplit training series.
- Every library degree in the stored pool.
- `python src/run_0F.py real`.

G-F2 is descriptive, reported either way. Prediction, written after item 4 above: at the lowest degree, the mean undecided share over laws is >= 1/2 for hare-lynx and for the COVID panel median. Also reported:
- refuted shares;
- near-data shares;
- region certificates;
- the zone of every stored 0R and 0E witness.

Because real noise is not iid and the library is surely wrong, zones on real series are statements about the set and not about the true system.

## 0F-3: horizon

Systems A, B, H and E4 at noise index 1, seeds 100-119, every degree. Laws are single-equation laws (one nonzero weight).

K_data is the ellipsoid {(x - mu)' C^-1 (x - mu) <= R^2}:
- mu and C are the training-state mean and covariance;
- R is the largest Mahalanobis radius of a training point;
- the gauge is g(x) = Mahalanobis(x) / R.

For 200 random directions from mu and g in {1.25, 1.5, 2, 3, 4}, check the band. Here lambda_min is the smallest eigenvalue of the stacked covariance, and n is the number of stacked rows. N^2 is the largest K over 20,000 points sampled in K_data, plus a local refinement from the best 20.

eps lambda_min / n * c_j^2 T_d(g)^2 <= h^2 <= eps N^2 T_d(g)^2

Rule G-F3:
- the band holds at >= 99 % of the points (100 % is expected; the margin is for N^2 being estimated from below);
- report the median width log(g_hi / g_lo) of the horizon bracket;
- if a closed-form radius is proved before the run, it is added by addendum with the rule "median abs log(measured / predicted) <= log 1.25".

## 0F-4: measurement design

- Systems and laws: A (degree 2, A2), B (degree 2, B2S, B2R), H (degree 2, H2), E4 (degree 2, E2S, E2R).
- Seeds 200-219. Noise index 1. The set is the 0F confidence set with estimated sigma.
- Target set: law-domain grid points where the true law value g_true(x) >= 4 tol_strong, i.e. a true margin of at least 5 tol. The truth is used only to define the target set and to simulate new data.
- Start from the stored training trajectories.
- Each step adds one trajectory of the system's training length and sampling, with the same noise level. The arm chooses its initial condition from 64 candidates, drawn uniformly in the bounding box of the law-domain states with every coordinate >= 0.01 (rng seed fixed per seed).
- After each step: simulate the truth, add noise, smooth, refit the set, re-estimate sigma.

Arms:
1. law-design (greedy): predict the new trajectory by rolling out the current fit, update the design, and choose the candidate that minimises the predicted undecided share of the target set at fixed xi_hat and eps;
2. random: a uniformly random candidate;
3. space-filling: candidates in Halton order (scrambled, seed fixed);
4. D-optimal (greedy): maximise log det of the updated design.

Outcome per seed and law:
- the number of added trajectories until the undecided share of the target set is < 1/20;
- the cap is 30, and an arm that has not reached the target at the cap counts as 31;
- a (seed, law) already below 1/20 at the start is dropped and counted.

Rule G-F4, on all (seed, law) pairs pooled. Pass if all three hold:
- the median of random / law-design is >= 2;
- the paired bootstrap 95 % interval of that median (10,000 resamples of pairs, rng 20261126) has its low end > 1;
- the median count of law-design is <= that of D-optimal.

If any fails, the design claim is dropped (ROADMAP_v3 kill rule).

## Not changed by results

Seeds, grids, alpha, tolerances, the 1/20 target, the cap, and the arms. The time limit per 0F-1 job is 30 minutes. A job lost to it is recorded and not rerun, and gates need at least 90 of 100 seeds per cell.

## Addendum 1 (2026-11-27, before any counted 0F-3 run)

On development seed 7 (system A), the lambda_min lower bound in 0F-3 is valid but useless. lambda_min of the Savitzky-Golay residual covariance is tiny, so the bound is below 1e-3 of h^2, and the horizon bracket almost never closes by g = 8.

Added: a second, exact lower bound. It evaluates the test polynomial p = T_d(l(x)) through the whitened design, h^2 >= eps c_j^2 T_d(g)^2 / ||Z delta_p||^2, instead of bounding ||Z delta_p||^2 by n / lambda_min. On seed 7 it is 4-20 times tighter and holds at 1000 of 1000 points per law (`src/horizon.py`, `test_poly_denominator`).

G-F3 now applies to both bands (each must hold at >= 99 % of points). The bracket width and the censored-ray count are reported for both. Nothing else changes. 0F-3 code is `src/horizon.py` and `python src/run_0F.py horizon`.

