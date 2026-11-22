# Protocol 0B: misspecification, needles, 20 seeds (frozen 2026-11-22, before any 0B run)

Status: drafted 2026-11-08 after Gate A was re-decided (`docs/gate_A2_result.md`, outcome 3). Settled and frozen 2026-11-22. No 0B cell, sample or needle on the grid below has run. Code was tested before the freeze only off the grid (seed 99, and a needle at v = 3e-4), and those test outputs were not kept.

Questions.
1. When the true law is outside the library (Holling type II), do validation-equivalent models break the stated laws, and do structure plus a counterexample loop fix that at low cost?
2. When does a random-sampling loop miss a violation that the verifier finds? Answered as a function of violation volume and sampler budget, on needles built on purpose, because natural needles were 1 in 60 in 0A2.
3. Do the Gate A numbers hold over 20 seeds, with confidence intervals?

## Units

### 0B-1. Holling II testbed (system H)

True model, Rosenzweig-MacArthur:
x' = x (1 - x/K) - a x y / (1 + h x), y' = -d y + b x y / (1 + h x),
with a = 1, h = 1, b = 0.8, d = 0.4, K = 4 (`src/systems.py`, `H_PAR`).

Checked by `tools/check_holling.py` (output `results/check_holling.json`, all_ok = True), exact arithmetic plus z3:
- Coexistence point (x*, y*) = (1, 3/2).
- Hopf at K = (b + d h) / (h (b - d h)) = 3. At K = 4 the Jacobian at (x*, y*) has trace 1/8 > 0 and determinant 3/20 > 0, so it is an unstable focus. An RK4 rollout settles on a limit cycle with x in [0.13, 3.17] and y in [0.41, 2.53]. At K = 2.7 it converges to the point.
- Trapping. With V = (b/a) x + y the rational terms cancel: V' = (b/a) x (1 - x/K) - d y. On V = c the maximum of V' is 1.568 - 0.4 c (at x = 2.8), so any c > 98/25 = 3.92 traps. c = 5 gives max V' = -54/125 = -0.432. The cycle stays inside V < 5.
- z3, each law multiplied by (1 + h x) > 0: H1x, H1y, H2 and H2 with the -0.432 margin are unsat for the truth. The tight version (>= the margin) is sat, so the margin is attained. Every query took under 0.01 s.

The polynomial library (degree 2 to 4 in x, y) cannot represent the rational term, so every candidate is misspecified. The positivity structure still holds for the truth (x' is x times something, y' is y times something), so arm S is a correct structure here.

Data. Same generator as the pilot (`data.make_cell`), with T = 40, dt = 0.05, Savitzky-Golay window 21 and order 3, 14 training, 6 held-out and 6 out-of-range trajectories. Training initial states: x ~ U(0.3, 3), y ~ U(0.2, 2). Out-of-range: x ~ U(3.5, 5.5), y ~ U(0.1, 0.6) (high prey, low predator, inside V <= 5). Noise indices 1 and 2 (2 % and 5 %), seeds 0-19. 40 cells.

Laws, domain {x, y >= 0, 0.8 x + y <= 5}, tolerances as in the pilot (weak 1e-3, strong 2e-2):
- H1x: x'(0, y) >= -tol for y in [0, 5]. H1y: y'(x, 0) >= -tol for x in [0, 6.25]. Structural under arm S.
- H2: -V' >= -tol on the segment 0.8 x + y = 5, x in [0, 6.25]. Non-structural.
All three are 1-D queries.

Pool. The pilot pool loop (degrees 2-4, 12 thresholds, full data plus 8 bootstraps, duplicates removed) and the pilot accuracy-equivalence rule (held-out derivative RMSE within 10 % of the cell's best), written without pandas in `src/gate_b.py` (`pool_rows`). The pool code reproduces the pilot pool for cell A_n1_s0 (130 candidates, 114 equivalent, same flags).

Every pool candidate is verified on H1x, H1y, H2 at strong tolerance (z3 first, 3 s limit, the pilot's sub-box fallback). For every violated law:
- volume = share of the law's domain with G < -tol_strong, computed exactly from the real roots of G + tol (all H laws are 1-D);
- depth = exact minimum of G / tol (`verify.exact_min_1d`).

Sample. 4 accuracy-equivalent candidates per cell, without replacement, `numpy.random.default_rng(20261122)`, one stream over the fixed cell order H (noise 1, seeds 0-19), H (noise 2, seeds 0-19), then A, A-mis, B (noise 1, seeds 0-19). Written to `results/0B_sample.csv` before any arm runs, never redrawn. 160 H candidates.

Arms: the pool candidate itself (unstructured STLSQ), S, C, Rnd. Loop, budget (15), witness margin, sampler (2000 points per law) and scoring are those of 0A (`gate_a.run_arm`, `gate_a.evaluate`). The Rnd stream is `default_rng([seed, cand, noise_idx, 5151])`.

Measures, per arm, per noise level and pooled:
- Fraction of pool equivalent candidates that violate at least one law.
- f(pool), f(S), f(C), f(Rnd) on the 160 sampled candidates: a candidate counts if its arm-model verifies all three laws. The H2-only fraction and per-law fractions are reported beside it.
- Held-out derivative RMSE ratio, and out-of-range rollout RMSE ratio, both against the candidate's own pool model. Rollouts use the pilot settings (RK45, rtol 1e-6, blow-up at |x| > 1e3, error capped at 10) plus a 20 s wall-clock cap per trajectory, which counts as a failed rollout.
- Volume and depth of every violation in the pool (`results/0B1_violation_volumes.csv`). This feeds G-B3.

### 0B-2. Needle test (constructed)

Base models. The true A model (degree 4 library), the true B model (degree 3), the first two arm-C models in `results/0A_sample.csv` order from 0A (system A) and from 0A2 B-pos (system B) that verify every law of their system. The A arm-C models are regenerated with `gate_a.run_arm` because 0A did not store coefficients. 6 base models.

Construction (`src/needles.py`). For base model Xi, law L and target (kappa, v): add to one equation i with c_i != 0
  delta f_i(X) = -(t / c_i) ((X_k / M)^p - s),
so that on the law's domain G_new = G_base - t (l^p - s) with l = X_k / M in [0, 1]. X_k is a state coordinate that is a free variable of the query, M its largest value on the domain. The violation set lies in the corner {l > s^(1/p)}. That corner is a scaled copy of the domain (interval, triangle or simplex), so its share is known and points in it can be drawn uniformly.
- t is set so that the minimum of G_new is -kappa * tol_strong. s is found by bisection on log10(1 - s^(1/p)) to hit the target volume.
- The search runs over every corner coordinate k allowed for the law and every power p = 1 .. library degree. It keeps the needle whose measured volume is closest to the target in log10.
- Coefficients are rounded to the 1e-6 grid. Volume and depth are then re-measured on the rounded model with 10^6 corner points. The measured volume v, not the target, is used in all analysis.
- Grid: kappa in {1.05, 1.5, 3}, target v in {1e-2, 1e-3, 1e-4, 1e-5, 1e-6}, every law of the base's system. 6 bases x (3 or 7 laws) x 15 = 450 needles.

For each needle:
- Random detector: at budgets n = 2k, 20k, 200k and 50 sampler seeds each, record whether any uniform domain point (the Rnd arm's sampler, `oracle._sample_domain`) has G < -tol_strong.
- Witness verifier: run once (`verify_spec_robust`, strong, Bernstein first on 3-D laws as in arm C). Record the outcome and time.
- Reference: detection probability 1 - (1 - v)^n.

Dropped from the draft: running the full Rnd loop on the kappa = 1.05 needles. The Rnd loop stops when its 2000-point sampler finds no violation, so for a fixed model the n = 2k miss rate is the rate at which the loop accepts the needle as clean. A separate loop run adds nothing.

### 0B-3. Gate A over 20 seeds

Arms S, C, Rnd on A, A-mis, B-pos, B-fac and H, seeds 0-19, noise index 1, 4 candidates per cell from the 0B sample. The H part is the noise-1 half of the 0B-1 sample (no second draw). B-pos and B-fac run on the same B candidates. Pools for A, A-mis and B are rebuilt with `pool_rows` for all 20 seeds (seeds 0-4 give the pilot pools again). The draw is new, so seeds 0-4 do not reuse the 0A candidates.

Arms, budget and scoring are those of 0A and 0A2 (`gate_a.run_sample`, `gate_a2.run_sample2`, unchanged). f(arm) is the candidate-level fraction verified on every non-structural law (A2; B1S, B2S, B2R; H2). Report f per arm with Wilson 95 % intervals, and f(C) - f(S) and f(C) - f(Rnd) with cluster-bootstrap 95 % intervals (resampling seeds, 2000 replicates, rng 20261122, percentile interval from sorted replicates 50 and 1949). This is the replication of Gate A, not a new gate.

## Decision rules (frozen)

All thresholds are compared as exact fractions of counts.

- **G-B1, misspecification matters:** on H, at least 1/5 of accuracy-equivalent pool candidates violate a law (pooled over the 40 cells), and the lower end of the 95 % cluster-bootstrap interval (resampling cells) is above 1/10. Wilson intervals are reported beside it.
- **G-B2, structure plus loop repairs at low cost:** on H, f(C) >= 9/10 over the 160 sampled candidates (all three laws), and the median out-of-range rollout ratio of C is <= 1.10.
- **G-B3, needles matter in practice:** fit one efficiency eta per budget n by maximum likelihood to the needle detections, with miss probability (1 - v)^(n eta). Set v50(n) = 1 - 0.5^(1/(n eta)). An efficient uniform sampler gives eta = 1, and then v50(20k) is about 3.5e-5. Needles matter if at least 1/20 of the accuracy-equivalent H violators from 0B-1 have volume below v50(20k). A candidate's volume is the largest volume over the laws it violates, since a sampler that hits any of them flags the candidate. If the share is below 1/20, the witness loop is reported as a certificate on top of a random loop that already finds what occurs, not as a better search. Binned observed and reference miss rates are reported beside the fit.
- **Replication:** Gate A outcome 3 is kept if, on the pooled A + A-mis + B-pos data, the 95 % interval for f(C) - f(Rnd) reaches below 1/5 and the interval for f(C) - f(S) lies above 0.

## Answers to the open questions (settled 2026-11-22, before freezing)

1. Holling parameters and trapping constant: confirmed by `tools/check_holling.py` (above). The values a = 1, h = 1, b = 0.8, d = 0.4, K = 4, c = 5 are kept. z3 decides the truth queries in under 0.01 s. Fitted H candidates give polynomial 1-D queries of degree <= 4. In an off-grid test cell (seed 99, 210 candidates) the slowest query took 0.13 s, so the Bernstein fallback is not needed for H. The largest z3 time per cell is reported in `results/0B1_pool_cells.csv`.
2. Needle volume on the 3-D SIR laws: 1e-6 is reachable. The corner region can be made as small as needed, and in tests on the true B model (B2S, B1I) and the true A model (A2, A1x), every target from 1e-2 to 1e-6 at kappa = 3, and every target from 1e-4 down at kappa = 1.05, came within a factor 2. Shallow needles cannot always be large. With kappa = 1.05 the violation set has to sit where the base margin is below 0.05 tol, and for the true B model on B2S that caps the volume near 1e-4. Such needles are kept at their measured volume and flagged `reached = False`. A degree-1 corner power (a planar cut) is included, because higher powers give thinner sets.
3. 0A2 stays the citable Gate A decision. 0B-3 gives the interval next to it.
4. The accuracy-equivalence rule for H is the pilot rule, unchanged. If any H cell has fewer than 4 equivalent candidates, the sample takes all of them, the cell is listed in `results/0B_sample_short_cells.txt`, and a dated addendum is written here before any arm runs.

## Code and outputs

- `src/systems.py` (system H), `src/gate_b.py` (pools, sample, H arms, needle jobs), `src/needles.py`, `src/run_0B.py` (driver, one JSON per job, resumable), `src/analyze_0B.py` (CSVs only).
- `src/savgol_np.py`: numpy copy of `scipy.signal.savgol_filter`, used only when `scipy.signal` will not load, as on this machine. Only interior points reach the fit. Rerunning a 0A2 job (B_s0, B-pos) with it gives the same coefficients as the stored run, with RMSE equal to within 1e-15.
- Phases: `python src/run_0B.py pool`, then `sample`, then `arms` and `needles`, then `python src/analyze_0B.py`.
- Outputs: `results/0B_jobs/{pool,arms,needles}/`, `results/0B_sample.csv`, `results/0B1_*.csv`, `results/0B2_*.csv`, `results/0B3_*.csv`, `results/0B_gate.csv`.

## Compute (estimate)

Pools: about 1 min per H cell (verification and rollouts), seconds for A, A-mis and B. Arms: about 3 s per H candidate, up to about 30 s per B candidate. Needles: about 10 s to build, plus detection (up to 1.1e7 points per needle). Jobs may run on a remote CPU machine (Kaggle) with the same code; the job JSONs are the record either way.
