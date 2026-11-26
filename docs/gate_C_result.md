# Protocol 0C result (robustness checks)

Run 2026-11-27 on Kaggle (4 CPUs, same pinned packages as 0B) from commit 2109feba, the frozen protocol. Protocol: `docs/protocol_0C.md`. 160 pool cells, 637 sampled candidates (two cells, A-mis at 0.5 % noise with seeds 1 and 6, had fewer than 4 equivalent candidates and gave what they had; `results/0C_sample_short_cells.txt`), 637 x 3 arm runs, 320 x 3 soft-penalty runs. No job errors. The full run took 18 minutes. Job records: `results/0C_jobs/`.

## Rules

| rule | value | result |
|---|---|---|
| R-C1, loop needed over a soft penalty | no: f(C) - f(Soft) = 0 for all three lam, median validation ratio 1.00 for Soft and C | **not needed** on these testbeds |
| R-C2, Gate A outcome 3 at 0.5 % noise | C - S [0.445, 0.529], C - Rnd 0.046 [0.017, 0.080] (n = 237) | holds |
| R-C2, at 5 % noise | C - S [0.329, 0.450], C - Rnd 0.046 [0.021, 0.071] (n = 240) | holds |
| R-C3, accepted models found violating by 200,000-point sampling | 0 of 2871 models | holds |

`results/0C_gate.csv`.

## 0C-1, soft penalty

On the 320 sampled candidates (A, A-mis, B-pos, H at 2 % noise) a penalty on the pilot's dense constraint grid verifies every non-structural law in all 320 models at lam = 1e-2, 1 and 1e2, with the same median validation ratio as C (1.000 pooled). C also verifies 320 of 320. S verifies 183 and Rnd 312 (B 72 of 80). Per system the soft arms match C: A 80, A-mis 80, B 80, H 80 of 80.

The three lam give the same shares and the same median ratios in every system, so the lam range shows no trade-off. A weaker penalty than 1e-2 was not tried, and a larger grid or a higher-dimensional system was not tried either.

What this does and does not say. On three- and two-state testbeds, a dense grid of about 100 points per law is enough to enforce these laws softly, and costs no accuracy. The loop's remaining advantages are the ones the roadmap already named: it adds few constraint rows (median 1 on H in 0B, 3 on the hare-lynx primary config) where the grid adds one row per grid point, and it ends with an exact witness or a certificate. Whether a grid still works as the state dimension grows is not tested here, and a grid grows exponentially with it.

`results/0C1_soft.csv`.

## 0C-2, noise level

Gate A outcome 3 holds at 0.5 %, 2 % and 5 %: structure alone is far behind C, and C is only a few candidates above Rnd.

| noise | pooled A + A-mis + B-pos n | f(S) | f(C) | f(Rnd) |
|---|---|---|---|---|
| 0.5 % | 237 | 0.515 | 1.000 | 0.954 |
| 2 % (0B) | 240 | 0.646 | 1.000 | 0.967 |
| 5 % | 240 | 0.608 | 1.000 | 0.954 |

Three things change with noise.
- At 0.5 % the structure alone fails on A-mis in every candidate (0 of 77) and on H in all but one (1 of 80), against 52 % and 35 % at 2 %. With little noise the missing term is visible in the data.
- At 0.5 % the repair costs accuracy where the library is wrong: median validation ratio of C is 1.49 on A-mis and 1.42 on H, above the 1.10 limit used in the earlier gates. At 2 % it was 1.08 and 1.06, at 5 % 1.02 and 1.01. The earlier "cheap" result was a 2 % result.
- The random loop is behind C on B-pos at every level (0.86 at 0.5 % and 5 %, 0.90 at 2 %), and on H at 0.5 % (0.91, C - Rnd 0.088 [0.038, 0.138]).

`results/0C2_noise.csv`.

## 0C-3, independent check

For 2871 accepted models (1911 arm models, 960 soft models), 200,000 uniform points per verified law found no point below -tol_strong. This is a check against verifier bugs with a sampler that cannot see narrow violations, so it supports the verdicts and does not prove them. `results/0C3_unsound.csv` is empty.

## Limits

- All synthetic, at most three states, one library family (degrees 2 to 4, the pilot's candidate grid).
- The soft-penalty arm uses the same grid and margin as the loop, so it is the hard-grid arm of the pilot with a penalty. A different penalty form, or one weighted by the data fit, was not tried.
- Noise indices 0 and 2 have their own candidate draws, so 2 % (0B) is a separate sample.
