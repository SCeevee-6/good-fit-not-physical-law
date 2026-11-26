# Protocol 0D: does a dense constraint grid still work at four states? (frozen 2026-11-27, before any 0D run)

Question. 0C found that a soft penalty on a dense constraint grid enforces the laws as well as the counterexample loop on systems with two and three states (`docs/gate_C_result.md`). The open point was dimension: a grid grows exponentially with the number of states and the loop does not. This protocol tests one more dimension, four states, with the same arms and a verifier that works for any number of states.

## Testbed E4 (SEIR)

True model: S' = -b S I, E' = b S I - s E, I' = s E - g I, R' = g I, with b = 0.6, s = 0.4, g = 0.2 (`src/dimtest.py`). Data as the SIR testbed B: 14 training, 6 held-out and 6 out-of-range trajectories, T = 80, dt = 0.25, Savitzky-Golay window 11 order 3, noise index 1 (2 %), seeds 0-19. Initial states: total population N uniform on (0.7, 1.3) (out-of-range (1.3, 1.5)), E0 = U(0.01, 0.1) N, I0 = U(0.01, 0.08) N, R0 = U(0, 0.1) N.

Library: polynomials of degree 2 and 3 in (S, E, I, R) (15 and 35 monomials), the pilot STLSQ sweep with 8 bootstraps, `gate_b.pool_rows`, accuracy rule as before (validation derivative RMSE within 10 % of the cell's best). Sample: 4 equivalent candidates per cell, without replacement, `default_rng(20261128)`, written to `results/0D_sample.csv` before any arm runs. 80 candidates.

Laws, on the 4-simplex S, E, I, R >= 0, S + E + I + R <= 3/2 (faces: the 3-simplex), tolerances weak 1e-3, strong 2e-2 (conservation 5e-3, 2e-2):
- E1S, E1E, E1I, E1R: each rate is >= 0 on the face where its own state is 0 (3 free variables).
- E2S: S' <= 0. E2R: R' >= 0. E3p, E3m: S' + E' + I' + R' within tolerance of 0 (4 free variables).
Structural family: E1E, E1I, E1R, E3p, E3m. Non-structural: E1S, E2S, E2R.

Structure (the four-state version of B-pos): free equations E, I, R, with S' = -(E' + I' + R'), so conservation holds exactly. In the E, I and R equations the coefficients of monomials free of the equation's own state are constrained >= 0, which guarantees positivity at that face. The true model satisfies it. S positivity, S' <= 0 and R' >= 0 are not guaranteed.

## Verifier

Any number of states: Bernstein branch and bound on the d-simplex (`src/bern_nd.py`, longest-edge bisection, 60 s budget), with z3 (3 s) as the fallback on time-out. A counterexample counts only after exact rational re-evaluation of the witness. A law that times out in both counts as not verified. Before the freeze `bern_nd` was checked on 42 random polynomials in 3 and 4 variables (degrees 2 and 3, margins from -0.03 to 0.3 around the tolerance, truth by dense sampling and local optimisation): no wrong verdict, no time-out, at most 0.01 s.

## Arms (all judged by the same verifier, strong tolerance)

- S: structure only. C: structure plus the witness loop (budget 15, margin 0.25 tol, as before). Rnd: structure plus the random loop (2000 points per law, worst point). The loops are `gate_a2.run_arm2` with the d-state verifier (`dimtest.run_arm_nd`).
- SoftA, SoftB: structure plus the soft penalty of 0C (lam = 1, `robust.fit_soft`) on a grid of constraint points with target G >= -0.25 tol. SoftA uses the pilot's budget of about 100 points per law, i.e. points per axis 11, 6 and 4 for laws with 2, 3 and 4 free variables (box grid filtered to the law's domain). SoftB uses 11 points per axis in every law.

Measures per arm: candidate-level share verified on the non-structural laws (all laws as secondary), median validation derivative ratio, median number of constraint rows, time. Rows by number of states are also tabulated from 0B (C rows on A and B-pos) and 0C (pilot-budget grid rows on A and B) next to E4.

## Rules (exact fractions, seeds as clusters, 2000 bootstrap replicates, rng 20261128)

- **R-D1, coverage advantage of the loop at four states:** C verifies at least 1/10 more of the candidates than SoftA (R-D1a), or than SoftB (R-D1b). Either is enough. If neither holds, the dense grid is still adequate at four states in coverage and accuracy, and I say so. Constraint rows are reported, not gated: the grid size is fixed by design, so a row count would pass by construction.
- **R-D2:** on a random quarter of the laws the Bernstein verifier decides on arm models, z3 also runs (3 s). The count of laws where both decide and differ must be 0.
- **R-D3:** for every verified law of every arm model, 200,000 uniform points of the law's domain are evaluated with the independent evaluator. The count of verified laws with a point below -tol_strong must be 0.

Reported, not gated: S, C and Rnd compared as in Gate A (C - S, C - Rnd with cluster intervals).

## Limits

One four-state system, SEIR with the structure that contains the truth. Coefficients, library degree 3 at most, one noise level. A five-state test would need a 5-simplex verifier. A grid coarser or finer than SoftA and SoftB, a different penalty, and hard constraints on the grid were not tried. At four states the simplex is a small share of the box (1/24), so a box grid wastes most points, and a grid built on the simplex would do better. I use the pilot's box construction so the result extends 0C.

## Time limits and code

Per arm run 900 s (Linux alarm; a time-out counts as not admissible). `src/bern_nd.py`, `src/dimtest.py`, `src/run_0D.py` (phases pool, sample, arms; one JSON per job, resumable), `src/analyze_0D.py`, `tools/kaggle_0D.py`. Outputs: `results/0D_jobs/`, `results/0D_sample.csv`, `results/0D_arms.csv`, `results/0D_per_law.csv`, `results/0D_rows_by_dimension.csv`, `results/0D_gate.csv`. Kaggle, pinned packages of 0B.

The pipeline was tested before the freeze on one candidate of an off-grid cell (seed 99). I saw its outputs: all five arms verified all eight laws, with 3 constraint rows for C and 364 and 5148 for SoftA and SoftB. That is one candidate and nothing was tuned from it; the files were not kept.
