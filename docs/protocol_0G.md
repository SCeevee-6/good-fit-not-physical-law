# Protocol 0G: a bias-corrected set for noisy regressors (E4)

Written 2026-11-28, frozen at the commit that adds this file, before any run on the gate seeds. Background: G-F1 failed on E4 (0/100 at degrees 2 and 3, `docs/gate_F_result.md`), and protocol 0F item 3 traced it to bias in the estimate, not to the covariance. ROADMAP_v3 Addendum 2 lists this unit as optional.

## What was seen before the freeze

Everything below used development seeds 10-29, outside every gate sample (gate seeds are 100-199, as in 0F). Tables are in `results/0G_dev/`.

1. All of 0F, including the E4 failure and its development notes.
2. First version of the corrected set (`results/0G_dev/dev_{E4,A,B}_b.csv`, seeds 10-29, noise index 1). Median stat/eps and coverage:
   - E4 degree 2: 0F set 2.87, 0/20; corrected 0.74, 20/20.
   - E4 degree 3: 0F set 227, 0/20; corrected 2.52, 3/20, with a non-positive reliability 12 times.
   - A degrees 2/3/4: 18, 17, 18 of 20 (0F set 18, 19, 19). B degrees 2/3: 19, 20 of 20 (0F set 19, 19).
   - Sanity: stored Xs equals H X to 9e-16; the same-row covariance sum_l H_tl D_tl is 3e-16 or smaller.
3. Monte Carlo on seed 10, with the clean trajectories fixed and 400 noise redraws (`tools/dev_0G_mc.py`). Degree 2: bias^2 0.17 against rank 60, Monte Carlo / sandwich variance trace ratio 0.99, coverage 0.96. Degree 3: bias^2 3863 against rank 140, variance ratio median 5.4 with a maximum of 8.6e5, coverage 0.
4. Reliability spectra, seeds 10-12. Degree 2: every direction at 0.64 or above. Degree 3: 2 or 3 of 35 directions below 0.5, the lowest between -0.04 and 0.03.
5. A truncation rule, dropping directions with reliability below kappa (`results/0G_dev/dev_E4_k*.csv`, 200-draw Monte Carlo on seed 10). kappa 0.5: E4 degree 3 0/20, Monte Carlo coverage 0. kappa 0.8: degree 3 0/20, and degree 2 drops to 12/20. Rejected: the eigenbasis of low-reliability directions is too unstable to project on. The frozen rule drops only non-positive reliabilities, which at degree 2 gives the same estimate and set as the first version.
6. Pipeline test of the full runner (seeds 10-11, all noise levels; outputs deleted). Stat/eps at noise index 0 for E4 degree 3: 1.03 and 1.21. Real series seen: hare-lynx (refuted share 0 for R1x, R1y and R2 at degrees 2, 3 and 4; 0F degree 2 mean refuted share was 0.18) and COVID Russia. For Russia at degree 2 the per-law (certified, refuted) shares went from 0F's B1S (0.54, 0.41), B1I (0.11, 0.89), B1R (0.05, 0.93), B2S (0.44, 0.52), B2R (0.32, 0.58), B3p (0.20, 0.50), B3m (0.61, 0.15) to B1S (0.50, 0.35), B1I (0.15, 0.82), B1R (1.00, 0), B2S (0.39, 0.46), B2R (0.91, 0), B3p (0.17, 0), B3m (0.17, 0), with B1R region-certified. So R-G4 below is not blind for hare-lynx or Russia.

## Estimator and set

Code: `src/eiv.py` (`CorrectedSet`). Same noise model as 0F: iid Gaussian state noise e_k with sd sigma_k, Savitzky-Golay smoothing H and derivative D, interior rows only.

- The smoothed state has noise u = H e with per-row variance v_tk = sigma_k^2 sum_l H_tl^2. For probabilists' Hermite polynomials with variance v, E[He_m(x + u; v)] = x^m.
- Corrected Gram and cross term: G~[a, b] = sum_t prod_k He_{a_k + b_k}(xs_tk; v_tk), b~[a, j] = sum_t prod_k He_{a_k}(xs_tk; v_tk) dX_tj. Their expectations equal the noise-free values, because on interior rows of a symmetric window sum_l H_tl D_tl = 0, so (H e)_t and (D e)_t are independent.
- Reliability of a direction: generalized eigenvalue of G~ against the plain Gram Theta' Theta. Directions with reliability <= 1e-6, or outside the range of the plain Gram, are unidentified. A law at a point that depends on them gets infinite width.
- Estimate: xi~ solves G~ xi = b~ on the identified directions (one G~ shared by all equations).
- Covariance: first-order sandwich, V = (I (x) G~)^+ [Z' C Z] (I (x) G~)^+, where C = B B' is 0F's linearised residual model, with the Jacobian at xi~.
- Set: Wald ellipsoid (xi - xi~)' V^-1 (xi - xi~) <= chi2_{rank, 0.95} on the identified directions. Coverage means the truth lies in it. Zones, sample checks and region certificates use the same code as 0F.
- sigma is estimated as in 0F (`lawmap.noise_sd`). Noise index 1 also runs with the true sigma.

## Runs

On Kaggle (4 CPUs, packages pinned as in 0F), from the frozen commit: `python tools/kaggle_launch.py push G <commit> <user>`.

- **Calibration** (`src/run_0G.py calib`): systems E4, A, B; noise indices 1 (2 %, primary), 0 and 2 (sensitivity); seeds 100-199; every library degree of each system (E4 2, 3; A 2, 3, 4; B 2, 3). At noise index 1: zone shares on the 0F grids, sample checks with 200 boundary models on seeds 100-119, and region certificates (30 s budget) for every law certified on its whole grid.
- **Real series** (`src/run_0G.py real`, report only): hare-lynx and the 10 COVID countries at the primary window, same domains and grids as 0F-2.
- Tables: `python src/analyze_0G.py` writes `results/0G_gate.csv`, `0G_calib.csv`, `0G_zones_calib.csv` and `0G_real.csv`.

## Gates

| rule | question | pass rule |
|---|---|---|
| G-G1 | does the corrected set calibrate on E4 at degree 2? | truth covered in >= 90 of 100 seeds (noise index 1, estimated sigma) |
| G-G2 | no harm on the systems that already calibrated | >= 90 of 100 for A at degrees 2, 3, 4 and B at degrees 2, 3, each |
| G-G3 | implementation | 0 sampled models contradicting the zone map over seeds 100-119, all systems, degrees and laws |

Reported, not gating, with predictions where I have them:

- R-G1: E4 degree 3 coverage. Prediction: below 90 of 100 (item 3 above).
- R-G2: mean certified share per system and degree, against the 0F confidence set at the same seeds. No prediction.
- R-G3: coverage at noise indices 0 and 2, and with the true sigma at index 1. No prediction.
- R-G4: COVID median over countries of the mean refuted share at degree 2. Prediction: still >= 0.10, i.e. the conflict between the data and the library at degree 2 is not an artefact of noise in the regressors. Not blind for Russia (item 6).

## What follows from the outcome

- G-G1, G-G2 and G-G3 pass: probability wording extends to E4 at degree 2 at 2 % noise. The corrected set becomes the reference set for synthetic systems in the write-up, with 0F reported next to it.
- G-G1 fails: probability wording stays limited to A and B at 2 % (as after 0F).
- G-G2 fails for a system: the 0F set stays the reference for that system.
- G-G3 fails: bug, stop and fix before any number from 0G is used.
- R-G1 is expected to fail. If so, the write-up says the corrected-score set does not work at degree 3 on four states, and why (reliability near zero in 2-3 directions, variance far above first order).
- Real-series zone maps under the corrected set are reported next to 0F-2. Neither replaces the other, since the iid Gaussian noise model is doubtful for reported cases.

## Limits stated in advance

- Gaussian noise is needed for the Hermite correction and for the independence of H e and D e. The symmetric window is needed for the latter.
- The covariance is first order. Item 3 shows that is enough at degree 2 and not at degree 3.
- The estimated sigma enters the correction. Its error is not propagated.
