# Protocol 0A2: system B rerun with a correct structure (frozen 2026-11-08, before any 0A2 run)

Why. In 0A the B structure kept only R-divisible monomials in the R-equation, so it could not hold the true R' = 0.2 I (`docs/gate_A_result.md`, ROADMAP_v2 Addendum 3). 0A2 reruns B with two structures that do hold the true model. Gate A is re-decided from the 0A results for A and A-mis plus B-pos. Nothing else changes.

## Unchanged from 0A
- Sample: the 20 B rows of `results/0A_sample.csv` (no new draw). A and A-mis are not rerun; their records are read from `results/0A_jobs/`.
- Data: `data.make_cell(S, 0.02, seed, noise_idx=1)`, same bootstrap variants, threshold and library degree from the candidate row.
- Arms S, C, Rnd, budget 15 iterations, 2000 uniform points per spec for Rnd, witness margin 0.25 x tol, strong tolerance, the same RNG stream for Rnd (`[seed, cand, 5151]`).
- Verifier and final judgement: `verify_spec_robust(..., 'strong', prefer='bernstein')` on every spec of the final model.
- Spec families: structural = B1I, B1R, B3p, B3m; non-structural = B1S, B2S, B2R. B1S stays non-structural for both structures so that f means the same thing as in 0A (see B-fac below).
- Scoring note of 0A: f(arm) is the fraction of candidates verified on every non-structural spec of their system; infeasible fits and timeouts count as not admissible and stay in the denominator; median held-out ratio is over arm-models with status ok. The pooled candidate x spec fraction is reported as secondary.

## B structures
Free unknowns are the I- and R-equations; S' = -(I' + R'), so S' + I' + R' = 0 exactly (B3p, B3m).
- **B-pos.** I-equation: only I-divisible monomials (gives B1I exactly). R-equation: every library monomial allowed; coefficients of R-free monomials are constrained >= 0. On the orthant this gives R'(S, I, 0) >= 0 exactly (B1R).
- **B-fac.** B-pos plus S divides S'. As linear equalities: for every S-free monomial m, Xi_I[m] + Xi_R[m] = 0. This gives S'(0, I, R) = 0, so B1S holds exactly. Addendum 3 described the equality as covering the "R-free and S-free" part only; that would give S'(0, I, 0) = 0 and not B1S on the whole face. The full S-free condition is what "S divides S'" means and is what runs here.

Fit details (`src/gate_a2.py`, `fit_free2`). The constrained STLSQ of 0A with the bounds and equalities added as fixed constraints on the active terms. They are not counterexample rows and do not count in n_rows. For B-fac an equality pair (Xi_I[m], Xi_R[m]) is kept or dropped together at each thresholding step, so thresholding one side cannot silently zero the other. After the solve the equalities are restored exactly (average of the pair) and solver residue below zero on bounded terms is set to 0, then coefficients are rounded to the 1e-6 grid as before.

## Check before running (done)
`tools/check_mask.py`, output in `results/0A2_check_mask.txt`. For both structures and both library degrees (2, 3): no true term is masked, no bound is violated by the true coefficients, every equality holds with residual 0, Xi_free M^T reproduces the true Xi exactly, and the true model is verified on all seven B specs. The same script confirms the 0A defect (the old B mask excludes the I term of the R-equation). On noise-free SIR derivatives `fit_free2` recovers the true model exactly for both structures and degrees (max abs error 0.0). No run on the sample data was made before this file was committed.

## Gate A rule (same thresholds as 0A)
Computed on 60 candidates: 40 A and A-mis from 0A, 20 B from 0A2.
- **G-A1:** f(C) - f(S) >= 0.20 AND median held-out ratio of C <= 1.10.
- **G-A2:** f(C) - f(Rnd) >= 0.20.
- Outcome 1 = both pass, 2 = G-A1 fails, 3 = G-A1 passes and G-A2 fails, with the same meaning as in 0A.

The decision is taken with **B-pos**. B-fac is computed the same way and reported next to it. If the two give different outcomes, both are reported and the B-pos outcome is the one used for planning, with the disagreement written into the result.

## Reported alongside (not part of the decision)
- Per system (A, A-mis, B under 0A, B-pos, B-fac) and arm: f, struct_all, median held-out ratio, median rows and iterations, infeasible count, timeouts.
- Per law: verified fraction per spec, arm and structure (`results/0A2_per_law.csv`).
- Rnd sampler misses: Rnd runs that stopped because 2000 points showed no violation while the verifier finds one in the final model.
- Expected checks: arm S verified on all structural B specs for both structures, and on B1S for B-fac. A failure is a finding, reported, not patched.

## Files
Code `src/gate_a2.py`; driver `src/run_gate_a.py` with `main(n, protocol='0A2')` (0A behaviour is the default and unchanged); raw output `results/0A2_jobs/<struct>/<sample_id>.json`; tables from `src/analyze_gate_a2.py`: `results/0A2_runs.csv`, `0A2_summary.csv`, `0A2_per_law.csv`, `0A2_gate.csv`. Result write-up in `docs/gate_A2_result.md`.

## Limits written down before running
20 B candidates, one noise level. A and A-mis are the same 40 records as in 0A, so any change in the gate comes from B only. The B structures still assume the shape of the law: B-pos assumes the R-free part of R' has nonnegative coefficients, B-fac assumes S-loss is proportional to S. Both hold for the true SIR model; neither is claimed beyond it.
