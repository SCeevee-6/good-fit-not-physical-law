# Protocol 0R: hare-lynx pelt series (Gate R), frozen 2026-11-24, before any fit on the real series

Question. On a real series, does the model picked by validation break a stated law? Roadmap rule: if not, the real-data claim is dropped.

## Data

Source: Hudson Bay Company trading records of snowshoe hare and Canadian lynx pelts, annual, 1845-1935, all areas of the company. Taken from the CRAN package `tsibbledata` 0.4.1 (`data/pelt.rda`, tarball md5 13d5766ad8bb1cb479e7f5cf478eb59c, checked). Converted with R to `data_raw/hare_lynx_pelts.csv` (91 rows, no missing values). `data_raw/pelt.rda` and `pelt.Rd` are kept next to it. Hare 1,800 to 152,650 pelts, lynx 3,190 to 79,350.

Pelt counts are a proxy for abundance, not a census. This is not the 1900-1920 table that textbooks use for Lotka-Volterra, it is the longer series. A different source would change the numbers here.

Before this freeze I looked at the series length, ranges, the split sizes and the scale numbers below. Nothing was fitted on the real series. The pipeline was tested only on a simulated series (system H sampled annually with 2 % noise), and that output was not kept.

## Setup (all fixed here)

- Split by year: training 1845-1895 (51 points), validation 1896-1915 (20), extrapolation window 1916-1935 (20).
- States: x = hare / mean training hare, y = lynx / mean training lynx. Training maxima x_m, y_m are rounded up to the next 0.05 (2.90 and 2.65 on the real series).
- Derivatives: Savitzky-Golay, order 3, window 7, applied separately inside each of the three segments (no smoothing across a split). Primary window 7, others in the sensitivity table. dt = 1 year.
- Training data is cut into contiguous blocks of 5 smoothed points (the last block takes the remainder). The block bootstrap (8 resamples plus the full set, seed 0) replaces the trajectory bootstrap of the pilot. Pool: pilot STLSQ sweep, degrees 2-4, 12 thresholds, duplicates removed (`gate_b.pool_rows`).
- Accuracy rule: validation derivative RMSE (normalised by the training derivative sd). Accuracy-equivalent (AES) = within 10 % of the best in the pool.
- Primary candidate: lowest validation RMSE in the pool (ties: fewer terms, lower degree, lower index). Second pick: lowest training AIC.

## Laws (stated before any fit)

Domain: triangle x, y >= 0, x/x_m + y/y_m <= c, with c = 3 (primary). The hypotenuse is the trapping shell.
- R1x: x'(0, y) >= -tol for y in [0, c y_m]. R1y: y'(x, 0) >= -tol for x in [0, c x_m]. Structural under arm S.
- R2: dV/dt <= tol on the shell, V = x/x_m + y/y_m = c. Non-structural.
- Tolerances as the pilot: weak 1e-3, strong 2e-2, in normalised units. Verifier: z3 with the 3 s limit and the pilot fallbacks (`verify_spec_robust`), strong tolerance. All three are 1-D queries.

These are modelling constraints, not laws of nature. The shell level c = 3 is my choice. Data check, done before the freeze: the largest V in any year is 1.88 in the training years, 1.22 in validation and 1.01 in the extrapolation window (V with unrounded training maxima; the rounded ones give 1.88 as the maximum over all years), so the data satisfy R2 with room. Positivity is trivial for counts.

A witness counts as extrapolative if its state lies outside the box of all observed (x, y). Positivity witnesses sit on a face (x = 0 or y = 0) and so are always outside it, because the observed minima are 0.03 and 0.11.

## Arms

Same as 0B: pool candidate (unstructured), S, C (z3 witness loop), Rnd (2000 random points per law), budget 15, margin and scoring as in `gate_a.run_arm` and `gate_a.evaluate`. They run for the primary candidate, the AIC pick and the first 60 AES candidates in order of validation RMSE. Rnd stream `default_rng([cand, 5151])`.

Measures: laws broken per candidate with exact depth and volume; validation RMSE ratio of each arm to its own pool candidate; extrapolation-window rollout RMSE ratio (rollout from the first smoothed state of 1916, RK45 as in the pilot, 20 s cap per trajectory, normalised by training state sd, error capped at 10); for each violating primary or AIC candidate, replay of the witness and a 10-year rollout from the witness state (minimum coordinate, largest V).

## Decision rules (frozen, exact fractions)

- **G-R1:** the primary candidate breaks at least one law at strong tolerance.
- **G-R2:** at least 5 AES candidates, and at least 1/5 of them break a law. Candidates share one dataset, so no interval is claimed.
- **G-R3:** the C model of the primary candidate verifies all three laws with validation ratio <= 1.10 and rollout ratio <= 1.10, and over the AES candidates in the arm set f(C) >= 9/10 with median validation ratio <= 1.10 and median rollout ratio <= 1.10.

Verdict. G-R1 and G-R2 pass: the real-data claim is kept. G-R3 also passes: the cheap-fix claim is kept on real data. G-R1 passes and G-R2 fails: inconclusive, a one-model anecdote that is not claimed. G-R1 fails: the real-data claim is dropped.

## Sensitivity (reported, not part of the gate)

SG window in {7, 5, 9} and shell level c in {3, 2, 4}, nine configs, the first being the primary. I report whether the primary violates in each and the AES share that violates. If the verdict only holds for c = 3, I say that.

## Limits

One series, one split, 20 validation points (14 after smoothing). With roughly nine cycles, the validation window holds about two. Annual sampling of a 10-year cycle makes the smoothed derivative crude, and pelt counts are not populations. No truth model exists, so there is no structure-recovery measure. The extrapolation window is used for rollout comparison only, never for selection.

## Code and outputs

`src/real_r.py`, `src/run_0R.py` (one JSON per config, resumable), `src/analyze_0R.py` (CSVs only), `tools/kaggle_0R.py`. Outputs: `results/0R_jobs/`, `results/0R_gate.csv`, `0R_summary.csv`, `0R_candidates.csv`, `0R_witnesses.csv`. `real_r.py` registers system R and patches the system-name map and degree list of `gate_b` at import. `gate_b.py` itself is unchanged and the 0B results stand. Run on Kaggle, same pinned packages as 0B.


## Addendum 1, 2026-11-25: time limits (before any real-series result exists)

The first Kaggle run of this protocol (code at commit 0ad1896c) stalled and was cancelled after about 2.5 hours. Its log stops about 80 seconds in, after numpy overflow warnings, and no config finished. No fit result from the real series was ever written or seen. I think a model with huge coefficients hung a numerical step that had no time cap, but I did not confirm where.

Changes, all about run time. No law, split, tolerance, candidate rule, arm or decision rule changes.
- Per candidate and per law: verification 60 s, volume and depth 60 s, witness replay and rollout 60 s. Per arm run (loop, evaluation and rollout): 180 s. Done with a wall-clock alarm, so it only runs on Linux (Kaggle). A step that hits its limit counts as `timeout`, which means not verified, the same as a z3 timeout in 0B. An arm run that hits its limit counts as not admissible, like an infeasible one. A witness step that fails is recorded as an error and is not used.
- Each config runs in its own process and is killed after 100 minutes. A killed config is written as timed out, listed in `results/0R_lost_configs.csv`, and not retried. If the primary config (window 7, c = 3) is killed or errors, the verdict is undecided and I say so. No config is dropped for any other reason. Total limit for the whole run is 6.5 hours.
- Candidates with a NaN validation error get an infinite one, so a NaN cannot be picked as the best model or break the equivalence rule.
- Progress is written per config to `logs/0R_progress_<config>.log`, and the runner always zips whatever exists.

The new code was tested on the simulated series only (kill path locally, full path on Kaggle). Commit for the real run: see `STATUS.md`.
