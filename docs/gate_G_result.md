# Protocol 0G result (bias-corrected set for noisy regressors)

Run on Kaggle (4 CPUs, packages pinned as 0F) from the frozen commit 115fdf22, 2026-11-28: 900 calibration jobs in 55 min and 11 real-series jobs in 17 s. No job errors, no time limits. Protocol: `docs/protocol_0G.md`, including what was seen before the freeze. Tables: `src/analyze_0G.py` writes `results/0G_gate.csv`, `0G_calib.csv`, `0G_zones_calib.csv` and `0G_real.csv`.

## Gates

| rule | value | threshold | result |
|---|---|---|---|
| G-G1 E4 degree 2 coverage | 95/100 (median stat/rank 0.974) | >= 90 | pass |
| G-G2 A degree 2 / 3 / 4 coverage | 90 / 95 / 95 of 100 (0F 95 / 96 / 94) | >= 90 each | pass |
| G-G2 B degree 2 / 3 coverage | 97 / 97 of 100 (0F 93 / 92) | >= 90 each | pass |
| G-G3 sampled models contradicting the zone map | 0 (780 system x degree x law checks) | 0 | pass |
| R-G1 E4 degree 3 coverage | 2/100 (0F 0/100), median stat/rank 3.41 | predicted < 90 | held |
| R-G4 COVID degree 2 median refuted share | 0.153 (0F 0.247) | predicted >= 0.10 | held |

So the probability wording now covers E4 at degree 2 at 2 % noise, next to A and B. A at degree 2 passed with no room to spare (exactly 90).

## Calibration

- **E4 degree 2 is fixed.** From 0/100 under the 0F set to 95/100, with a median stat/rank of 0.97 against 1 for an exact chi-square. With the true sigma it is 94/100, so estimating sigma costs nothing here. No seed had a non-positive reliability direction; the median lowest reliability was 0.78.
- **E4 degree 3 is not**, as predicted: 2/100 (3/100 with the true sigma). In 66 of 100 seeds at least one direction had reliability <= 0, and the median lowest reliability was -0.03. The corrected-score estimate exists there, but its error is far from the first-order sandwich (development Monte Carlo: variance about 5 times the sandwich, with heavy tails). A degree-3 library on four states at 2 % noise asks for degree-6 moments the data cannot correct.
- **No harm on A and B.** Coverage 90-97 against 92-96 under 0F. B improves slightly (97 against 92-93).
- **Other noise levels** (R-G3, coverage of 100):

| system, degree | 0.5 % noise, 0G (0F) | 5 % noise, 0G (0F) |
|---|---|---|
| E4, 2 | 93 (74) | 85 (0) |
| E4, 3 | 38 (0) | 7 (0) |
| A, 2 / 3 / 4 | 82 / 86 / 89 (75 / 73 / 79) | 92 / 97 / 94 (94 / 97 / 92) |
| B, 2 / 3 | 98 / 96 (98 / 97) | 96 / 98 (92 / 81) |

  The correction helps at every noise level where 0F under-covered, but it does not fully fix A at 0.5 %. There the problem is not regressor noise (the correction is tiny at that noise level), so the smoothing bias that 0F item 3 called negligible at 2 % is the likely cause at 0.5 %. I have not tested that.

## Zone maps (R-G2)

The corrected set certifies less than the 0F set at every system and degree. Mean certified share of the law domains at 2 %:

| system, degree | 0G | 0F |
|---|---|---|
| E4, 2 / 3 | 0.265 / 0.002 | 0.409 / 0.054 |
| A, 2 / 3 / 4 | 0.609 / 0.371 / 0.011 | 0.699 / 0.429 / 0.026 |
| B, 2 / 3 | 0.894 / 0.489 | 0.953 / 0.641 |

For E4 this is the expected direction: the 0F set was centred on a biased estimate and did not cover the truth, so its certified zones were not backed by anything. For A and B both sets are calibrated and the corrected one is wider, which is the usual price of a corrected-score estimate.

Region certificates: 63 attempted, 63 verified (A2 at degree 2 in 20 seeds and at degree 3 in 15; B1I at degree 2 in 20, B1R and B2R in 4 each). No E4 law was certified on its whole domain, so none was attempted there.

## Real series (report only)

- **Hare-lynx.** Every law is undecided on its whole domain at every degree (0F degree 2: certified 0.20, refuted 0.18). At degree 2 the lowest reliability is 0.32, so the corrected set is wide. At degrees 3 and 4 some directions have reliability <= 0 and every law has infinite width. The 0F refutation of hare positivity and the trapping law therefore depended on treating the noisy pelt counts as exact regressors.
- **COVID, 10 countries.** At degree 2 the median refuted share drops from 0.247 to 0.153 and the median undecided share rises from 0.186 to 0.377. The refutation survives in most countries (Bangladesh 0.34, Iran 0.38, Russia 0.23, Iraq 0.18), and it vanishes in Chile and Turkey (0 from 0.29 and 0.10). Region certificates at degree 2 verified 14 times (B1I in 7 countries, B1R in 5, B1S in 2); at degree 3 all 4 attempts found counterexamples. At degree 3 the median refuted share is 0.086 (0F 0.089).
- The R-G4 prediction held, so the conflict between the data and the degree-2 library is not only an effect of regressor noise. It is smaller than 0F said, though, and the iid Gaussian noise model is still doubtful for reported cases.

## What changes

- Probability wording: A and B at 2 % noise (as before), plus E4 at degree 2 at 2 %. Not E4 at degree 3, and not any system at 0.5 %.
- The corrected set is the reference set for synthetic systems in the write-up, with 0F next to it.
- The hare-lynx law refutation from 0R and 0F is withdrawn as a probability statement. Under the corrected set the data decide nothing about those laws. The real-data claim from 0R (equivalent models break laws, and C repairs them) is about the fitted models and is unaffected.
- For COVID the refuted share at degree 2 is reported under both sets: about a sixth of the law domain under 0G, a quarter under 0F.

## Limits

- Gaussian state noise and a symmetric smoothing window are needed for the correction. Both hold by construction on the synthetic systems and are assumptions on the real series.
- The covariance is first order. It is adequate where the lowest reliability is well above 0 (E4 degree 2 median 0.78) and not where it is near 0.
- The error in the estimated sigma is not propagated. The known-sigma runs suggest this matters little at 2 %.
- 100 seeds per cell. A coverage of 90/100 is consistent with a true rate somewhat below 0.90.
