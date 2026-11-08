# Gate A re-decided (protocol 0A2)

Run 2026-11-08 after `docs/protocol_0A2.md` was committed (f996151e). 40 jobs (20 B candidates x 2 structures), 6 workers, 31 s, no job errors, no infeasible fits, no verifier timeouts. A and A-mis records are the 40 from 0A, unchanged. Tables from `src/analyze_gate_a2.py`.

## Decision

Outcome 3 with B-pos, the structure the decision is taken on. B-fac gives outcome 3 as well. Counterexamples add value over structure alone, and a random-counterexample loop with 2000 points per spec does as well as the witness loop. The random loop stays a main arm.

| B structure | f(C) - f(S) | median held-out ratio C | G-A1 | f(C) - f(Rnd) | G-A2 | outcome |
|---|---|---|---|---|---|---|
| B-pos (decision) | 0.317 (19/60) | 1.000 | pass | 0.017 (1/60) | fail | 3 |
| B-fac (reported) | 0.200 (12/60) | 1.000 | pass | 0.000 | fail | 3 |
| B under 0A (defective, for reference) | 0.533 | 1.059 | pass | 0.100 | fail | 3 |

Pooled candidate x spec differences (secondary): B-pos 0.21 and 0.01, B-fac 0.12 and 0.00.

B-fac sits exactly on the G-A1 threshold, and all of its 12 moves come from A and A-mis (A 2 of 20, A-mis 10 of 20); B adds nothing because structure alone already passes every B candidate. So with B-fac, G-A1 passes only on the A-systems.

Analysis fix made while reading the result. The first pass of `analyze_gate_a2.py` compared floats and reported B-fac as outcome 2, because 60/60 - 48/60 is 0.19999... in floating point. The gate is now tested on exact fractions of counts. The 0A values (0.533, 0.100) are not near a threshold, so the 0A verdict is unaffected.

## Per system

| system / structure | f(S) | f(C) | f(Rnd) | median (max) held-out ratio S / C / Rnd | candidates needing rows C / Rnd |
|---|---|---|---|---|---|
| A (0A) | 0.90 | 1.00 | 1.00 | 1.000 / 1.000 / 1.000 | median 0 / 0 |
| A-mis (0A) | 0.50 | 1.00 | 1.00 | 1.066 / 1.066 / 1.066 | median 2 / 1 rows |
| B-pos | 0.65 | 1.00 | 0.95 | 1.000 (1.003) / 1.000 (1.000) / 1.000 (1.000) | 7 / 7 of 20, at most 3 rows |
| B-fac | 1.00 | 1.00 | 1.00 | 1.000 / 1.000 / 1.000 | 0 / 0 |

Arm S is verified on every structural spec for both structures (struct_all = 1.00), as required. Under B-fac arm S is also verified on B1S for all 20 candidates, as the construction promises.

## Per law (B, arm S, verified out of 20)

| law | 0A structure | B-pos | B-fac |
|---|---|---|---|
| B1S, S'(0,I,R) >= 0 | 10 | 14 | 20 |
| B2S, S' <= 0 | 0 | 18 | 20 |
| B2R, R' >= 0 | 1 | 19 | 20 |

All B-pos failures are removed by both loops except one (below). Full table: `results/0A2_per_law.csv`.

## Findings

1. The 0A accuracy cost on SIR was the defect. With a correct structure the median held-out ratio on B is 1.000 for every arm and structure, and the largest single value is 1.003. The 1.15 / 1.23 / 1.20 figures in `gate_A_result.md` should not be cited.
2. The 0A structure-only failures on SIR were mostly the defect too. B2S failed in 20/20 and B2R in 19/20 under 0A; with B-pos it is 2/20 and 1/20, and with B-fac 0/20.
3. B-fac is enough on its own. Making S divide S' gives a model that satisfies every B law at strong tolerance on all 20 candidates, with no counterexamples. B-pos needs the loop on 7 of 20.
4. Random sampling missed one violation. Under B-pos the Rnd loop on B_s0_69 (degree 3) stopped after one iteration with 2000 clean points, and the verifier then found a B2S violation. A 121^3 grid gives the worst point at the corner (S, I, R) = (0, 1.5, 0) with G = -0.0204 against tol 0.02, and the violating set is about 7e-6 of the domain. Under 0A there were six such misses; with a correct structure there is one. Natural needles are rare on these testbeds.

## What this changes

- Outcome 3 stands, now on a correct structure. The witness-vs-random question cannot be answered with the needles that occur naturally here (1 in 60). Protocol 0B has to build narrow violations on purpose and measure sampler budget against violation volume, with the witness loop as the reference.
- The useful split for the taxonomy is per law: B1I, B1R, B3p, B3m are enforced by B-pos, B1S by B-fac, and B2S, B2R are not structural in either, yet B-fac candidates satisfied them in this sample. That is an empirical result on 20 candidates, not a guarantee.
- No accuracy cost was measured for structure or repair at this noise level on SIR.

## Limits

20 B candidates, one noise level, the same 40 A-systems records as 0A. The G-A2 gap is now 1 candidate. This is a go/no-go call.

## Files

`results/0A2_jobs/{Bpos,Bfac}/*.json` (raw), `results/0A2_runs.csv`, `results/0A2_summary.csv`, `results/0A2_per_law.csv`, `results/0A2_gate.csv`, `results/0A2_check_mask.txt`. Log `logs/gate_0A2_run.log` (local only).
