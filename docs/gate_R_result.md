# Protocol 0R result (Gate R, hare-lynx)

Run 2026-11-26 on Kaggle (4 CPUs, Python 3.13, same pinned packages as 0B) from commit fcfb358d, which is the frozen protocol plus addendum 1 (time limits). Protocol: `docs/protocol_0R.md`. Data: `data_raw/hare_lynx_pelts.csv` (Hudson Bay Company pelt counts, CRAN `tsibbledata` 0.4.1). Job records: `results/0R_jobs/`. All 9 configs finished, 0 killed, 0 errors, 0 limit hits in the pool (`0R_lost_configs.csv` is empty). The whole run took 7 minutes. The first run (commit 0ad1896c) stalled and was cancelled with no output, so nothing from it is used.

## Gates (primary config: SG window 7, shell level c = 3)

| rule | value | threshold | result |
|---|---|---|---|
| G-R1 | the model with the lowest validation error (degree 2, 5 terms) breaks R1x and R2 | breaks a law | pass |
| G-R2 | 45 of 47 accuracy-equivalent candidates break a law (0.957) | n >= 5, share >= 1/5 | pass |
| G-R3 | C verifies the selected model (validation ratio 1.051, extrapolation rollout ratio 0.975); over the 47 candidates f(C) = 47/47, medians 0.987 and 0.618 | all <= 1.10, f(C) >= 9/10 | pass |

Verdict: real-data claim kept, cheap fix kept. `results/0R_gate.csv`.

## What the numbers say

- The selected model breaks the positivity law for hare and the trapping law. Its R1x violation covers 0.90 of the face x = 0 (worst point 103 tolerances deep), and its R2 violation covers 0.47 of the shell (60 tolerances deep). Witnesses are confirmed by replay. From the R1x witness (hare 0, lynx 1 in normalised units) a rollout dips to x = -0.06. From the R2 witness (x = 5, y = 3.38) the rollout goes to V = 3.08 and a coordinate reaches -0.80.
- The AIC pick breaks all three laws (R1y on 0.064 of its face).
- Structure alone (S) leaves most models broken: f(S) = 5/47 at c = 3. Structure fixes R1x and R1y but not the trapping law.
- C repairs 47/47 with a median of 3 added rows. Rnd repairs 45/47 (0.957). At c = 4 it is 46/47 for C and 42/47 for Rnd. At c = 2 both are 47/47.
- Cost of C: median validation ratio 0.987 and median extrapolation-window rollout ratio 0.618, so repaired models roll out better over 1916-1935 than the originals.

## Sensitivity (nine configs, not part of the gate)

The primary model and the AIC pick break a law in all nine. The share of equivalent candidates that break a law is 45/47 (window 7), 63/65 (window 5) and 19/20 (window 9), the same at every c. f(C) is 1.0 everywhere except window 7, c = 4 (46/47). f(S) goes from 0.11 to 0.42. Rnd is below C at window 7, c = 3 and 4 and at window 5, c = 4. `results/0R_summary.csv`.

## Limits that matter

- The fits are poor. The best validation error is 0.89 training-derivative standard deviations, so the models explain little of the validation derivative. Annual samples of a roughly 10-year cycle give crude derivatives, and pelt counts are not populations. A model that fits this badly breaks laws easily, so G-R1 and G-R2 passing says less than it did on the simulated systems.
- The positivity witnesses lie on a face (hare 0 with lynx 1 or 2, or lynx 0 at hare 0), and the data never come near a face (smallest observed hare 0.03, lynx 0.11). "Outside the observed box" is true for them by construction, as the protocol said. The trapping witness is at hare 5, far from the data too. So this shows that the fits break the laws where the data say nothing, not that they break them where the data are.
- The 47 candidates share one dataset and are not independent, so no interval is given. The share is a count on one series.
- One series, one split, 20 validation points (14 after smoothing). No true model, so no structure-recovery measure. The shell level c = 3 is my choice.
- The sign of the Rnd versus C difference is small and not tested here.

`results/0R_candidates.csv`, `0R_witnesses.csv`, `0R_summary.csv`, `0R_gate.csv`.
