# Protocol 0F3b (unit R2, post hoc): the projected lower bound on the 0F-3 horizon rays

Written 2026-11-28, frozen at the commit that adds this file, before the run. **Post hoc**: the 0F-3 rays and their results are already seen (`docs/gate_F_result.md`), so nothing here is a blind test of the horizon. What is new is the lower bound from Theorem 4 of `docs/law_horizon.md`, which did not exist when 0F-3 ran.

## Question

0F-3 found the lower side of the horizon band loose: 53-89 % of rays had no closed bracket by g = 8. Theorem 4 gives a tighter lower bound from the data projected on the ray's direction (the span of T_0(l) .. T_d(l)) instead of the single test polynomial T_d(l). How much of the censoring does it remove on the same rays?

## Run

`src/run_0F3b.py` with `src/horizon2.py`. Same systems (A, B, H, E4), seeds 100-119, sets (0F stacked set at 2 % noise), rays, sup estimates and test points as 0F-3. The random stream is consumed in the same order, so 0F-3's own censoring counts are reproduced as a check. Kaggle, mode F3B. Tables: `src/analyze_0F3b.py` writes `results/0F3b_horizon.csv` and `results/0F3b_gate.csv`.

Seen before the freeze: one local job (A, seed 100), used to confirm the reproduction. 0F-3 counts reproduced in 6 of 6 rows; 0 violations of the projected bound; censored rays for A1x / A1y at degrees 2, 3, 4 went from 167 / 186, 124 / 169, 65 / 115 to 108 / 162, 81 / 148, 0 / 0 (of 200).

## Rules

| rule | kind | pass rule |
|---|---|---|
| projected bound below the exact h^2 at every band point | check | 0 violations (a violation is a bug in the code or in Theorem 4) |
| 0F-3 censoring counts reproduced | check | every row equal |
| censored share, per system and degree | prediction | projected share <= half of the T_d(l) share |
| median log(g_lower / g_true) | report | how far each lower bound's crossing lies beyond the true one |

Corollary 6's closed-form bracket is not compared against these rays. Under the errors-in-variables set Theorem 5 is not exact, and the band itself holds by construction, so a comparison would be circular.
