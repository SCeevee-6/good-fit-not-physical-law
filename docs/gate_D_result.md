# Protocol 0D result (four states)

Run 2026-11-27 on Kaggle (4 CPUs, pinned packages as 0B) from commit 3641866a, the frozen protocol plus a test-only fix noted in it. Protocol: `docs/protocol_0D.md`. 20 pool cells, 80 sampled candidates (every cell had 4 equivalent candidates), 5 arms each, no job errors. The run took 14 minutes. Job records: `results/0D_jobs/`.

## Rules

| rule | value | result |
|---|---|---|
| R-D1a, C verifies >= 1/10 more than the pilot-budget grid (SoftA) | 0.0625 (interval [-0.05, 0.175]) | no |
| R-D1b, C verifies >= 1/10 more than the fine grid (SoftB) | -0.1375 ([-0.238, -0.05]) | no |
| R-D1, coverage advantage of the loop at four states | | **not supported** |
| R-D2, Bernstein vs z3 disagreements | 0 of 703 laws decided by both | holds |
| R-D3, verified laws with a violating point in 200,000 samples | 0 of 400 models | holds |

`results/0D_gate.csv`.

## What happened at four states

Share of the 80 candidates verified on the non-structural laws (E1S, E2S, E2R), median validation ratio, median constraint rows (`0D_arms.csv`):

| arm | verified | validation ratio | rows |
|---|---|---|---|
| S (structure only) | 7 | 1.000 | 0 |
| C (witness loop) | 69 | 1.000 | 15 |
| Rnd (random loop) | 51 | 1.000 | 13 |
| SoftA (grid, pilot budget) | 64 | 1.000 | 364 |
| SoftB (grid, 11 per axis) | 80 | 1.000 | 5148 |

- The fine grid verifies every candidate and beats the loop by 11, at a median of 5148 constraint rows against 15. The pilot-budget grid is at 64, five below the loop with an interval that includes 0. So the finding from 0C extends to four states: a dense grid enforces the laws as well as the loop or better, without accuracy cost.
- The loop did not do well here. 11 of its 80 runs stopped at the 15-iteration budget with a violation left (E2S in 11, E1S in 2). Half of its runs added 15 or more rows (maximum 41). I did not test a larger budget. With 3 states it needed a median of 1 row.
- The random loop is clearly worse at four states. Its sampler of 2000 points found nothing in 77 of 80 runs, yet the verifier rejects 29 of the 80 final models (51 verified). That is the narrow-violation effect from 0B, now at a size that matters: with the simplex a small share of the box and four coordinates, 2000 points are too few. C - Rnd is 0.225 [0.125, 0.338], against 0.03 to 0.10 at two and three states.
- Structure alone is nearly useless on this system: S verifies 7 of 80 (E1S 7, E2S 38, E2R 47). The truth satisfies the structure, but the fitted models do not satisfy the three laws it leaves open. C - S is 0.775 [0.675, 0.863].
- Verifier: Bernstein and z3 agree on all 703 laws both decided, and no accepted model fails the independent sample. The five arms of one candidate took a median of 39 s in total.

## Rows against number of states (`0D_rows_by_dimension.csv`)

| system | states | median rows, loop C | median rows, grid (pilot budget) | rows, fine grid |
|---|---|---|---|---|
| A | 2 | 0 | 300 | |
| B-pos | 3 | 1 | 422 | |
| E4 | 4 | 15 | 364 | 5148 |

The pilot-budget grid stays near 300 to 420 rows because the budget is about 100 points per law. At four states that is 4 points per axis, and the simplex keeps only 35 of the 256 box points per 4-variable law, which is why SoftA drops to 64. The fine grid is 5148 rows at four states and the same grid in five variables would start from 11^5 = 161,051 box points per law. The loop's row count grew from 0 to 15 and began to hit its budget.

## What this does and does not say

- A dense grid has not stopped working at four states. It does need a grid fine enough (SoftB) and then thousands of rows, and the box grid wastes most of its points on the simplex. A grid built on the simplex was not tried.
- The loop's remaining advantage is the row count (24 times fewer than the pilot grid, over 300 times fewer than the fine grid) and the witness it ends with. On coverage it is behind the fine grid and the loop budget is a limit I did not vary.
- Against the random loop the witness loop has a clear coverage advantage at four states. This is the first system in the project where that is large.
- One four-state system whose structure contains the truth, one noise level, one library family. A five-state test needs a verifier for the 5-simplex.

`results/0D_arms.csv`, `0D_per_law.csv`, `0D_rows_by_dimension.csv`.
