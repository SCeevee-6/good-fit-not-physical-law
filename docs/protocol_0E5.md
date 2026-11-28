# Protocol 0E5: a five-state testbed for calibration and measurement design

Written 2026-11-28, frozen at the commit that adds this file, before any run on the gate seeds. ROADMAP_v3 Addendum 2 dropped E5 for lack of time. It is back as a side track at the user's request.

## Testbed

E5 (`src/e5.py`): SEIR with two infectious stages, closed population.

    S' = -b S (I1 + I2),  E' = b S (I1 + I2) - s E,  I1' = s E - g I1,  I2' = g I1 - g I2,  R' = g I2

with b = 0.6, s = 0.4, g = 0.4. The time grid, Savitzky-Golay filter, noise model (2 %, noise index 1), initial-condition ranges and law families are E4's, extended to five states:
- face positivity on S, E, I1, I2, R (E1*);
- S' <= 0 (E2S) and R' >= 0 (E2R);
- conservation (E3p, E3m).

Grids: `lawmap.law_grid` gets n5 = 7 points per axis for 5-variable laws. Region certificates use `dimtest.verify_nd` on the 5-simplex.

## What was seen before the freeze

- All of 0F-4 (E4 design: pooled E4 median random/law 1.52; E2S 7 vs 9, E2R 16 vs 31 censored) and 0G (E4 degree 2 calibrated, lowest reliability median 0.78; degree 3 not, median -0.03).
- Development seeds 10-29, corrected set, noise index 1 (`results/0G_dev/dev_E5.csv`). Degree 2: 8/20 covered, median stat/eps 1.11, median lowest reliability 0.091, no non-positive reliability. Degree 3: 1/20 covered, median stat/eps 2.3, a non-positive reliability in 20/20.
- Development seed 10: E2S and E2R start fully undecided (target sets of 116 and 252 points). A degree-3 refit takes 12 s, and 64 candidates about 13 s.
- Pipeline test (calibration seeds 10-11, design seed 0 with cap 2; outputs deleted). Calibration stat/eps 0.62 / 0.71 at degree 2 and 1.51 / 1.85 at degree 3. No region certificate was attempted. Every design arm reached the cap.

## Runs

Kaggle, mode E5, from the frozen commit: `src/run_0E5.py calib`, then `design`.

- calib: the 0G corrected set at degrees 2 and 3, seeds 100-199, estimated and true sigma, plus zones, sample checks (seeds 100-119) and region certificates, as in 0G.
- design: the 0F-4 code unchanged (`src/design.py`), with E5 at degree 3, laws E2S and E2R, seeds 200-219 (the 0F-4 seeds), arms law design / D-optimal / random / Halton, 64 candidates, cap 30 added trajectories, target < 1/20 undecided.
- Tables: `src/analyze_0E5.py` writes `results/0E5_gate.csv`, `0E5_calib.csv`, `0E5_zones.csv` and `0E5_design.csv`.

## Gates and predictions

| rule | question | pass rule | prediction |
|---|---|---|---|
| G-E5a | does the corrected set calibrate E5 at degree 2? | >= 90 of 100 covered | **fails**: development gave 8/20 with lowest reliability about 0.09 |
| G-E5b | implementation | 0 sampled models contradicting the zone map | pass |
| G-E5c | does law design beat random starts at five states? | median random/law >= 1.5 with bootstrap low > 1 | pass (E4 gave 1.52 on the same laws) |
| R-E5a | degree 3 coverage | report | below 90 |
| R-E5c | E5 against E4 | report the difference in median ratio | none |

Post-hoc report, labelled so: coverage against the median lowest reliability for every system and degree in 0G and 0E5. It asks whether reliability predicts when the corrected set can be trusted.

## What follows

- G-E5a fails as predicted: probability wording does not extend to E5. The write-up says the corrected set needs reliability well above 0.1, which at 2 % noise E5 does not have even at degree 2.
- G-E5c passes: the design advantage seen on E4 holds at five states. Whether it grows is answered by R-E5c, descriptively.
- G-E5c fails: the design advantage does not carry to five states, and the design claim stays dropped (as after 0F-4).
- G-E5b fails: bug, stop.

## Limits

One five-state system, one noise level, one library family, 20 design seeds, a coarse 5-D law grid (7 per axis).
