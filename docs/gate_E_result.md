# Protocol 0E result (second real series, COVID-19 reported-case SIR panel)

Run 2026-11-26 on Kaggle (same pinned packages as 0B) from commit 942d55e5, the frozen protocol. Protocol: `docs/protocol_0E.md`. Data: JHU CSSE confirmed, recovered and deaths files plus the UID population table, in `data_raw/jhu/`. All 30 configs (10 countries, 3 smoothing windows) finished: 0 errors, 0 killed (`0E_lost_configs.csv` is empty). The run took 28 minutes. Job records: `results/0E_jobs/`.

## Gates (primary window 21 days, 10 countries)

| rule | value | threshold | result |
|---|---|---|---|
| G-E1 | validation-selected model breaks a law in 10 of 10 countries | >= 1/2 | pass |
| G-E2 | 35 of 35 equivalent candidates break a law, bootstrap interval [1.00, 1.00] | >= 1/5, lower > 1/10 | pass |
| G-E3 | B-pos C verifies all seven laws on 34 of 35 (0.971); median over countries of the median validation ratio 1.465, of the rollout ratio 0.641 | f(C) >= 9/10, both medians <= 1.10 | **fail** |

Verdict: the claim holds on a second real series (G-E1, G-E2). The cheap fix is not retained on this series, because of the validation cost. `results/0E_gate.csv`.

## What the numbers say

- Every selected model breaks a law. The laws broken differ by country (B1S in 8 of 10 selected models, B2R in 5, B2S in 4). `0E_candidates.csv`.
- The equivalent sets are tiny: 1 to 10 candidates per country (Turkey 10, Bangladesh 9, Iraq 6, Argentina, Chile, Iran 2 each, the other four 1), out of pools of 94 to 141. So the 100 % pooled share rests on 35 candidates, one third of them from Turkey and Bangladesh, and the bootstrap interval has no width because every country is at 100 %.
- The witness loop (C) repairs 34 of 35 under B-pos and under B-fac. The one failure in each ends at the 15-iteration budget with a violation left. Under B-pos the random loop (Rnd) ends clean in 30 of 35 runs, 4 end infeasible (its random rows contradict the structure) and 1 stops at the budget. The verifier then accepts its models on B1S in 31, B2R in 30 and B2S in 29 of 35. B-pos C - Rnd on the three non-structural laws is 5/35 = 0.14, interval [0.00, 0.29].
- Structure alone: B-pos S verifies the non-structural law B1S in 0 of 35 models and the other two in 19 and 20 of 35. B-fac S verifies B1S in 35, B2S in 34, and B2R in only 11 of 35. B-pos C - S = 34/35 [0.87, 1.00]. B-fac C - S = 23/35 [0.26, 1.00].
- The cost is real here. C changes the validation error by country-median factors of 5.21, 2.32, 3.29, 0.38, 1.03, 1.89, 1.94, 0.51, 1.04, 0.46 (Argentina, Bangladesh, Chile, Indonesia, Iran, Iraq, Russia, South Africa, Turkey, Ukraine). Three countries improve and seven get worse, three of them by more than a factor of 2. Extrapolation-window rollouts get better in most countries (country medians from 0.09 to 4.57, median 0.64).
- Per law, over the 35 candidates (`0E_per_law.csv`): the pool violates B1S in 31, B2R in 16, B2S in 11, B3p and B3m in 13, B1I in 5, B1R in 4.

## Against the earlier results

- The structure-only result from the SIR testbed does not carry over. There, B-fac alone satisfied all seven laws on all 20 and then 80 sampled candidates. Here B-fac alone leaves B2R unverified in 24 of 35. The synthetic truth lies inside the B-fac structure and these series do not, so the guarantee only holds when the structure contains the truth.
- The H-system finding (structure alone can lower the share of models that satisfy a law) shows up again: B2R under B-fac goes from 19 of 35 satisfied in the pool to 11, and B1S under B-pos from 4 to 0. It is not universal (B2S under B-pos goes from 24 to 19, B-fac raises it to 34).
- Gate A outcome 3 (witness loop beats structure, barely beats random counterexamples) is consistent here for B-pos. The difference to the random loop is 5 candidates and its interval includes 0.

## Sensitivity

Windows 15 and 29, same countries: the primary model breaks a law in all 20 configs, and the equivalent candidates break a law in every one of them (1 to 6 candidates per config). `0E_summary.csv`.

## Limits that matter

- Primary validation errors, in training-derivative standard deviations, are above 1 in eight of ten countries (alphabetical, Argentina to Ukraine: 0.69, 2.00, 0.76, 6.97, 1.92, 2.72, 1.04, 1.02, 3.56, 5.83). These fits do not explain the validation derivatives, as in Gate R, so passing G-E1 and G-E2 says that poor fits break laws easily.
- The data live in a corner of the simplex: largest I between 0.0009 and 0.0154, smallest S between 0.89 and 0.99. The laws are posed on the whole simplex, so violations are mostly far from the data.
- Reported cases, ten correlated series, a selection rule that drops many large countries, one pandemic window. Candidates inside a country share the data.
- The validation ratio is a median over only 1 to 10 candidates per country. With one candidate in four countries, the country medians are single numbers.

`0E_summary.csv`, `0E_candidates.csv`, `0E_witnesses.csv`, `0E_per_law.csv`, `0E_gate.csv`.
