# Protocol 0A: structure first, then counterexamples (frozen 2026-11-08, before any 0A run)

Question. On the pilot testbeds, do counterexample witnesses from the SMT verifier remove law violations that a structure-only fix leaves in place? Does a witness-driven loop beat a loop that uses random counterexamples at the same budget?

This protocol decides Gate A. Its rule and every setting below are fixed now.

## Inputs
- Candidate pool: `results/pilot_candidates.csv` (the pilot grid, unchanged).
- Cells: noise index 1 (rho = 2 %), seeds 0-4, for systems A, A-mis and B. That is 15 cells.
- Candidates: 4 per cell, drawn without replacement from accuracy-equivalent rows (`aes == True`) using `numpy.random.default_rng(20261108)`. The drawn row ids are written to `results/0A_sample.csv` before any arm runs.
- Data: regenerated with `data.make_cell(S, rho, seed, noise_idx=1)`, the same generator as the pilot.

## Structure (the by-construction fix)
- System A: the x-equation keeps only monomials divisible by x; the y-equation keeps only monomials divisible by y. Positivity of both populations then holds exactly, because each equation is x*h or y*h.
- System B: the I- and R-equations keep only monomials divisible by I and R. The S-equation is defined as minus the sum of the I- and R-equations, so S'+I'+R' = 0 holds exactly. Positivity of S is not guaranteed by this and is left to be tested.
- Thresholds, library degree and bootstrap variant are taken from the candidate row. The structural arms re-run constrained STLSQ on the masked library.

## Arms
- **S** (structure only): structured STLSQ, no counterexample rows.
- **C** (witness loop): arm S plus the counterexample loop in `src/ceg.py`. At most 15 iterations. Each violated spec adds one row from the z3 witness (strong tolerance, margin 0.25 x tol).
- **Rnd** (random-counterexample loop): arm S plus a loop of the same shape. Each iteration draws 2000 uniform points per spec domain. If any point has G < -tol_strong, the worst sampled point adds one row. Stop when no sample violates, or after 15 iterations. The final model is judged by z3, not by the sampler.

## Outcome measures (all judged by the verifier, not by the loop)
- Admissible fraction, strong tolerance, per arm, over two spec families:
  - structural family: A1x, A1y, B1S, B1I, B1R, B3p, B3m;
  - non-structural family: A2, B2S, B2R.
- Held-out derivative RMSE, relative to the candidate's own pool RMSE (median ratio).
- Number of constraint rows (median).
- Check: structural-family admissibility of arm S must be 1.0 for A-systems and for B3p/B3m. Any other failure is reported as a finding, not fixed.

## Gate A rule (frozen)
Let f(arm) be the admissible fraction over the non-structural family, pooled across the 60 candidates.
- **G-A1, counterexamples add value:** f(C) - f(S) >= 0.20 AND median held-out ratio of C <= 1.10.
- **G-A2, witnesses beat random counterexamples:** f(C) - f(Rnd) >= 0.20.

Outcomes:
1. G-A1 and G-A2 both pass. Witness counterexamples are needed. Proceed to 0B with the witness loop as the main method.
2. G-A1 fails. Structure is enough on these specs. The project claim is redirected to the rule for when a structural fix is sufficient. The roadmap gets a dated addendum. Proceed to 0B on that question.
3. G-A1 passes and G-A2 fails. Counterexamples help, but the feedback loop explains it, not the witness type. The random-counterexample loop becomes a main arm. Report it as the result.

## Limits written down before running
- 60 candidates, 15 cells, one noise level. This is a go/no-go decision, not a precise effect size.
- The structural fix is an assumption about the true law's shape. The true A and B models satisfy it; that is by design, and the S arm is not claimed to generalise beyond these systems.
- A-mis is included so that the structural fix is tested where the library lacks a true term.

## Clarification (2026-11-08, before any 0A run)
B1S is moved from the structural family to the non-structural family. The S-equation is defined as minus the sum of the I- and R-equations, so S' at S = 0 is -(I h_I + R h_R), whose sign is not guaranteed by the structure. Structural family: A1x, A1y, B1I, B1R, B3p, B3m. Non-structural family: A2, B1S, B2S, B2R.

## Scoring note (2026-11-08, before the rerun; no 0A results exist)
- f(arm) is read at the candidate level: the fraction of the 60 sampled candidates whose arm-model is verified on every non-structural spec of its system. A pooled candidate x spec fraction is reported next to it as a secondary number. The gate is decided on the candidate-level figure.
- A run that ends infeasible, or a spec that ends in timeout, counts as not admissible. Candidates stay in the denominator of 60.
- The median held-out ratio is taken over arm-models that ran (status ok).
- Driver: `src/run_gate_a.py` (one JSON per candidate in `results/0A_jobs/`). Analysis: `src/analyze_gate_a.py`.

## Defect found after the run (2026-11-08)
The structure for system B is wrong. The R-equation keeps only monomials divisible by R, but the true law is R' = 0.2 I, which has no R factor. So the S arm on B cannot represent the true model, and the line above saying the true A and B models satisfy the structure is false for B. A and A-mis are fine (checked with `tools/check_mask.py`: no true term is masked). The frozen verdict is not changed. The correction is a new protocol, `docs/protocol_0A2.md`, to be written before it runs.
