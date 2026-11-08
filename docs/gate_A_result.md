# Gate A result (protocol 0A)

Run 2026-11-08, after the scoring note in `docs/protocol_0A.md`. The first run was cancelled by a kernel reset and left no results; this is the first complete run. Sample `results/0A_sample.csv` (60 candidates, seed 20261108), code `src/gate_a.py`, driver `src/run_gate_a.py`, tables from `src/analyze_gate_a.py`.

## Decision

Outcome 3. G-A1 passes, G-A2 fails. Counterexamples help, but at this budget the feedback loop explains most of it, not the witness type. The random-counterexample loop becomes a main arm.

| criterion | quantity | value | threshold | pass |
|---|---|---|---|---|
| G-A1 | f(C) - f(S) | 0.533 | >= 0.20 | yes |
| G-A1 | median held-out ratio, C | 1.059 | <= 1.10 | yes |
| G-A2 | f(C) - f(Rnd) | 0.100 | >= 0.20 | no |

f is the fraction of the 60 candidates verified on every non-structural spec. The pooled candidate x spec reading gives the same verdict (0.61 and 0.06). All 60 jobs ran, no infeasible fits, no verifier timeouts.

## Per system

| system | f(S) | f(C) | f(Rnd) | median held-out ratio S / C / Rnd | median rows C / Rnd |
|---|---|---|---|---|---|
| A | 0.90 | 1.00 | 1.00 | 1.00 / 1.00 / 1.00 | 0 / 0 |
| A-mis | 0.50 | 1.00 | 1.00 | 1.07 / 1.07 / 1.07 | 2 / 1 |
| B (SIR) | 0.00 | 1.00 | 0.70 | 1.15 / 1.23 / 1.20 | 7.5 / 6 |

## Checks and findings

- Structural family: arm S is verified on every structural spec for all 60 candidates, as required.
- Structure alone leaves violations: B2S fails in 20/20 B candidates, B2R in 19/20, B1S in 10/20, A2 in 12/40 A and A-mis candidates.
- All six Rnd failures are B2S on SIR. In each one the sampler stopped after 3-4 iterations because 2000 random points showed no violation, and z3 then found one. These are narrow violations near the tolerance, the same pattern the pilot found. This is where witnesses beat random points, and it points straight at the 0B needle test.
- On SIR the repaired models cost more accuracy than the pooled median shows: median held-out ratio 1.23 for C and 1.20 for Rnd, against 1.15 for S. The gate uses the pooled median, so this does not change the decision, but it is reported here.

## Limits

60 candidates, 15 cells, one noise level. The 0.10 gap in G-A2 comes from 6 of 20 SIR candidates. This is a go/no-go call, not an effect size.

## Files

`results/0A_jobs/*.json` (raw, one per candidate), `results/0A_runs.csv` (one row per candidate and arm), `results/0A_summary.csv`, `results/0A_gate.csv`, `logs/gate_0A_run.log` (local only, logs are not pushed).

## Defect found after the run (2026-11-08)

The B structure masks out a true term. The R-equation was restricted to monomials divisible by R, but the true law is R' = 0.2 I. The A and A-mis masks keep every true term (`tools/check_mask.py`).

What this does to the result:
- On A and A-mis, where the structure is right, C and Rnd both reach 1.00. The whole G-A2 gap (0.100) comes from B.
- 20 of the 32 candidates that move from failing to passing between S and C are B candidates. Without B, f(C) - f(S) is 0.30 (A 0.10, A-mis 0.50).
- The SIR accuracy figures (S 1.15, C 1.23, Rnd 1.20) are fits under a wrong structure. They say nothing yet about the cost of a correct one.
- B1S failing in 10/20 is partly the same defect: S' = -(I' + R') inherits the bad R-equation.

The verdict stays outcome 3 as frozen, flagged as confounded on B. Protocol 0A2 reruns B with a correct structure under the same rule.

Re-decided with a correct B structure in `docs/gate_A2_result.md` (protocol 0A2): outcome 3 holds; the SIR accuracy cost and most SIR structure-only failures above were the defect.
