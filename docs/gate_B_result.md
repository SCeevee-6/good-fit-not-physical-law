# Protocol 0B result

Run 2026-11-23 on Kaggle (4 CPUs, Python 3.13, pinned numpy 2.4.6, scipy 1.17.1, z3-solver 5.1.0.0, cvxpy 1.9.3, clarabel 0.11.1, sympy 1.14.0), from the frozen commit 9f1cf045. Protocol: `docs/protocol_0B.md`. Runner: `tools/kaggle_0B.py`, `tools/kaggle_launch.py`. Job records: `results/0B_jobs/`. Tables: `src/analyze_0B.py`. 100 pool cells, 480 arm jobs and 450 needles, with no job errors and no cell below 4 equivalent candidates (minimum 48). So no addendum was needed before sampling.

## Gates

| rule | value | threshold | result |
|---|---|---|---|
| G-B1 | 3365 / 3456 = 0.974 of equivalent H candidates break a law; cluster bootstrap [0.966, 0.981] | >= 1/5, lower > 1/10 | pass |
| G-B2 | f(C) = 160/160; median OOD rollout ratio of C 0.588 | >= 9/10, <= 1.10 | pass |
| G-B3 | 0 / 3365 equivalent H violators below v50(20k) = 3.50e-5 | >= 1/20 | fail |
| Replication | pooled A + A-mis + B-pos (240): f(C)-f(Rnd) 0.033 [0.017, 0.050]; f(C)-f(S) 0.354 [0.288, 0.429] | C-Rnd low end < 1/5, C-S low end > 0 | pass |

`results/0B_gate.csv`.

## 0B-1, Holling II (system H)

- Misspecification breaks the laws almost always, and mostly the positivity laws: of the equivalent violators, 2103 break H1x and H1y, 788 only H1x, and 461 involve H2. On the 160 sampled candidates the pool model verifies H1x in 5 %, H1y in 31 %, H2 in 88 %.
- Structure alone fixes H1x and H1y (by construction) but leaves H2 broken more often than the pool did: f(H2) drops from 0.88 (pool) to 0.40 (S). Forcing x' = x h(x, y) and y' = y g(x, y) on a misspecified library pushes the fit to break the trapping law instead.
- C repairs all 160 (at most 8 rows, median 1). Rnd repairs 159. The one Rnd miss, H_n2_s9_c173, still breaks H2.
- Cost: median held-out derivative ratio 1.036 for C (1.056 at 2 %, 1.007 at 5 %). Out-of-range rollouts get better, not worse: median ratio 0.588 for C and 0.591 for S, because many pool models drift or blow up from high-prey starts.
- z3 decided every query, with no time-outs. The slowest query took 0.17 s.

`results/0B1_pool_cells.csv`, `0B1_arms.csv`, `0B1_summary.csv`, `0B1_violation_volumes.csv`.

## 0B-2, needles

- 450 needles were built, with measured volumes from 8.2e-8 to 1.3e-2. 400 came within a factor 2 of their target. The witness verifier found all 450 (270 by z3, 180 by Bernstein, slowest 0.41 s).
- The random detector behaves as a uniform sampler should. The fitted efficiency is 1.04, 0.99 and 0.99 at n = 2k, 20k and 200k, and observed misses match the reference 1 - (1 - v)^n in every volume bin. Example: at n = 20k, needles near 1e-4 are missed 15 % of the time (reference 15 %), and needles near 1e-5 are missed 81 % of the time (reference 81 %).
- v50 is 3.3e-4 at 2k, 3.5e-5 at 20k and 3.5e-6 at 200k.

`results/0B2_needles.csv`, `0B2_detector_fit.csv`.

G-B3 fails by a wide margin. The smallest natural violation volume among the 3365 equivalent H violators is 0.064 (1st percentile 0.44), more than three orders of magnitude above v50(20k). Natural violations in these testbeds are large, so a random loop finds them. As fixed in the protocol, the witness loop is reported as a certificate on top of a random loop that already finds what occurs, not as a better search. The needle curve says where that stops being true: below about 1e-4 of the domain at 2k points.

## 0B-3, Gate A over 20 seeds

| group | n | f(S) | f(C) | f(Rnd) | C-S [95 %] | C-Rnd [95 %] | med iid ratio C |
|---|---|---|---|---|---|---|---|
| A | 80 | 0.925 | 1.000 | 1.000 | 0.075 [0.000, 0.175] | 0.000 [0, 0] | 0.999 |
| A-mis | 80 | 0.525 | 1.000 | 1.000 | 0.475 [0.338, 0.613] | 0.000 [0, 0] | 1.077 |
| B-pos | 80 | 0.488 | 1.000 | 0.900 | 0.513 [0.388, 0.638] | 0.100 [0.050, 0.150] | 0.999 |
| B-fac | 80 | 1.000 | 1.000 | 1.000 | 0.000 [0, 0] | 0.000 [0, 0] | 0.999 |
| H (2 %) | 80 | 0.350 | 1.000 | 1.000 | 0.650 [0.538, 0.775] | 0.000 [0, 0] | 1.056 |
| A + A-mis + B-pos | 240 | 0.646 | 1.000 | 0.967 | 0.354 [0.288, 0.429] | 0.033 [0.017, 0.050] | 1.000 |

Gate A outcome 3 replicates: counterexamples add a lot over structure, and witnesses add little over random counterexamples. The only place Rnd falls short is B-pos (8 of 80 left inadmissible). B-fac again satisfies every SIR law with structure alone. 0A2 stays the citable decision; these are its intervals.

`results/0B3_candidates.csv`, `0B3_summary.csv`.

## Notes

- Kaggle ran Python 3.13; the off-grid tests before the freeze ran locally on 3.11 with the same package versions.
- The pool JSONs hold every candidate's coefficients, so any number here can be recomputed from `results/0B_jobs/` without rerunning the fits.
