# Protocol 0C: robustness checks (soft penalty, noise level, independent check), frozen 2026-11-27, before any 0C run

Three checks on the synthetic testbeds, each with a rule fixed here. They come from the roadmap's baseline list (soft penalties) and from gaps the earlier gates left: Gate A and 0B ran at 2 % noise only, and every "verified" verdict rests on one verifier.

Questions.
1. Is a soft penalty on a dense constraint grid, which needs no verifier and no loop, as good as the counterexample loop?
2. Does the Gate A pattern (witness loop far above structure alone, barely above random counterexamples) hold at 0.5 % and 5 % noise?
3. Do models the verifier accepts survive an independent uniform-sampling check?

## 0C-1. Soft-penalty baseline

Sample: the 0B sample at noise index 1 for systems A, A-mis, B and H, 4 candidates per cell, seeds 0-19 (320 candidates, `results/0B_sample.csv`). B uses structure B-pos. Arms S, C and Rnd are taken from the stored 0B arm records, not rerun.

New arms, same structure, same STLSQ threshold and bootstrap variant as the candidate: Soft1, Soft2, Soft3. Each refits the structured library with the pilot's dense constraint grid (about 100 points per law, `ceg.dense_repair` grid), target G(u) >= -0.25 tol_strong (the repair margin of C and Rnd), as a penalty lam * sum_k max(0, target - r_k . xi)^2 added to the least-squares objective, with lam = 1e-2, 1, 1e2. The penalised fit is solved once per lam, with no counterexample loop. Models are judged by the same verifier as every other arm (`gate_a.evaluate`, strong tolerance, Bernstein first, z3 fallback). A fit that comes back infeasible or times out counts as not admissible.

Measures: candidate-level share verified on the non-structural laws (A2; B1S, B2S, B2R; H2), per system and pooled; median validation derivative ratio against the candidate's own pool model.

**R-C1.** The counterexample loop is "needed over soft penalties" if, for every lam, f(C) - f(Soft) >= 1/10 on the pooled 320 candidates, or the median validation ratio of that Soft arm is above 1.10. Otherwise a soft penalty is an adequate cheaper alternative on these testbeds, and I say that. Cluster-bootstrap intervals (seeds, 2000 replicates, rng 20261127, sorted replicates 50 and 1949) are reported for each gap and are not part of the rule.

## 0C-2. Noise sensitivity

Systems A, A-mis, B (B-pos) and H, noise indices 0 (0.5 %) and 2 (5 %), seeds 0-19: 160 cells. Pools are built with the pilot sweep (`gate_b.pool_rows`, no verification needed), 4 accuracy-equivalent candidates per cell drawn without replacement with `default_rng(20261127)`, one stream in the fixed cell order (noise 0, then noise 2; within each A, A-mis, B, H; seeds 0-19), written to `results/0C_sample.csv` before any arm runs. A cell with fewer than 4 equivalent candidates gives all it has and is listed in `results/0C_sample_short_cells.txt`.

Arms S, C, Rnd as in 0A, 0A2 and 0B (`gate_a.run_arm`, `gate_a2.run_arm2` with B-pos), budget 15, same margin and scoring. Rnd stream `default_rng([seed, cand, noise_idx, 5151])`.

**R-C2**, per noise level: on the pooled A + A-mis + B-pos candidates, Gate A outcome 3 holds if the 95 % interval of f(C) - f(S) lies above 0 and the lower end of the interval of f(C) - f(Rnd) is below 1/5. f is the candidate-level share verified on the non-structural laws. Intervals are cluster bootstrap over seeds. For reference the same table is computed at 2 % noise from the 0B arm records. H is reported separately, not in the rule.

## 0C-3. Independent check of accepted models

For every arm model produced in 0C-1 and 0C-2 (status ok): for each law the verifier calls verified, draw 200,000 uniform points from the law's domain (`oracle._sample_domain`, fixed stream) and evaluate G with the independent evaluator `needles.g_vec`.

**R-C3.** The count of verified laws with a sampled point below -tol_strong must be 0. Any other number is reported as a verifier problem, before anything else.

## Not in this protocol

Tolerance sensitivity (the strong tolerance 2e-2 is a modelling choice and stays), library degree beyond 2-4, other noise models, and any new real series.

## Time limits and code

Per arm run 300 s (Linux alarm, a timed-out arm counts as not admissible), as in the addenda of 0R. `src/robust.py`, `src/run_0C.py` (phases pool, sample, arms, soft; one JSON per job, resumable), `src/analyze_0C.py` (CSVs only), `tools/kaggle_0C.py`. Outputs: `results/0C_jobs/`, `results/0C_sample.csv`, `results/0C1_soft.csv`, `results/0C2_noise.csv`, `results/0C3_unsound.csv`, `results/0C_gate.csv`. Run on Kaggle with the pinned packages of 0B. The code was tested before the freeze on the first 0B-sample candidate of A, B and H (soft arms) and on one A candidate at noise index 0 (arms S, C, Rnd). I saw those outputs: the soft arms verified every law at all three lam on the A, B and H candidates, with validation ratios 1.0, 1.0 and 1.23. They are 3 of 320 candidates, nothing was tuned from them, and the files were not kept.
