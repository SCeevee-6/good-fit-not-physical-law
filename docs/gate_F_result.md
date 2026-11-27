# Protocol 0F result (law zones over a confidence set)

Runs on Kaggle (4 CPUs, pinned packages as 0B), 2026-11-27:
- 0F-1 and 0F-2 from the frozen commit 01a24e87 (900 calibration jobs in 3.1 h, 11 real-series jobs in 20 s);
- 0F-3 and 0F-4 from 04350905, which is the frozen protocol plus addenda 1-3 (80 horizon jobs in 2.6 h, 60 design jobs in 9.9 h).

There were no job errors and no time-limit hits. Protocol: `docs/protocol_0F.md`, including what was seen before the freeze. Tables: `src/analyze_0F.py` writes `results/0F_gate.csv`, `0F_calib.csv`, `0F_zones_calib.csv`, `0F_real.csv`, `0F_witnesses.csv`, `0F_horizon.csv` and `0F_design.csv`. Figure: `figures/fig_0F.png`.

## Gates

| rule | value | threshold | result |
|---|---|---|---|
| G-F0 sampled models that contradict the zone map | 0 (780 system x degree x law checks, 200 models each) | 0 | pass |
| G-F0 region certificates (reported) | 93 attempted, 93 verified | - | reported |
| G-F1 coverage, A degree 2 / 3 / 4 | 95 / 96 / 94 of 100 | >= 90 | pass |
| G-F1 coverage, B degree 2 / 3 | 93 / 92 of 100 | >= 90 | pass |
| G-F1 coverage, E4 degree 2 / 3 | 0 / 0 of 100 | >= 90 | **fail** (expected, protocol item 3) |
| G-F2 hare-lynx mean undecided, lowest degree | 0.62 | prediction >= 1/2 | held |
| G-F2 COVID median over countries | 0.19 | prediction >= 1/2 | **did not hold** |
| G-F3 band, lambda_min lower bound | 680,000 / 680,000 points | >= 99 % | pass |
| G-F3 band, exact test-polynomial lower bound | 680,000 / 680,000 points | >= 99 % | pass |
| G-F4 median random / law-design | 1.42 (bootstrap low 1.17); law-design median 2 = D-optimal median 2 | >= 2, low > 1, <= | **fail** |

## 0F-1: the set is right where it should be, except E4

- **Calibration.** Median (r(true) - r(fit)) / q is 0.91-0.99 for A and 1.01-1.02 for B, against 1 for an exact chi-square. Coverage is 92-96 % at a nominal 95 %. Known sigma gives the same (89-96 %), so the noise estimate is not the problem.
- **E4 fails completely.** The median ratio is 3.9 at degree 2 and 286 at degree 3, so no seed covers the truth. This is the errors-in-variables bias in the estimate. It grows like sqrt(n) at a fixed noise level, and E4 has the most rows (protocol item 3).
- **Noise level matters.** At 0.5 % noise, A coverage falls to 73-79 %. At 5 % noise, B degree 3 falls to 81 %. So the 2 % result does not carry to every noise level, and probability wording applies to A and B at 2 % only.
- **The confidence set decides much more than the Rashomon set** (panel b). Mean certified share of the law domains at 2 %:
  - A: 0.38 against 0.11;
  - B: 0.80 against 0.14;
  - E4: 0.23 against 0.005 (E4 not calibrated).

  The refuted share is 0 for every system with the truth in the library, as it must be.
- **Implementation.** No sampled model contradicts the map. All 93 attempted region certificates were verified exactly: A2 at degrees 2 and 3, and B1I, B1R and B2R at degree 2. These are statements that every model in the 95 % set obeys the law on the whole domain.

## 0F-2: real series

Mean over laws at each degree (confidence set, estimated sigma, training data only):

| series | degree | certified | refuted | undecided | region certificates |
|---|---|---|---|---|---|
| hare-lynx | 2 / 3 / 4 | 0.20 / 0 / 0.04 | 0.18 / 0.03 / 0.01 | 0.62 / 0.97 / 0.95 | none attempted |
| COVID, 10 countries, median | 2 | 0.54 | 0.25 | 0.19 | 15 verified, 1 timeout |
| COVID, 10 countries, median | 3 | 0.13 | 0.09 | 0.67 | 2 counterexamples |

- At degree 2, COVID refutes a quarter of the law domain on the plane: median 0.25, from 0.10 (Turkey) to 0.57 (Russia). Refuted means every model in the set breaks the law there. The library and the iid-noise model are both doubtful for reported cases, so this says one of them is wrong; it does not say which. The prediction (undecided >= 1/2) assumed the data would be silent. At degree 2 they are not silent, they are in conflict.
- Region certificates at degree 2 verified B1I in 7 countries, B1R in 6, and B2R and B3m in one each. At degree 3 the two attempted certificates (Iran and Ukraine, B1R) found counterexamples that the grid had missed. The exact check catches what the grid does not.
- Stored witnesses under the confidence set: 16 undecided, 9 refuted and 4 certified, out of 29. The 4 certified ones belong to models outside the set, which were selected on validation error.

## 0F-3: the bound holds, the bracket is wide

- Both bands hold at all 680,000 points. So the upper side, h^2 <= eps N^2 T_d(g)^2, is confirmed on A, B, H and E4 at every degree.
- The lower side is loose. Even the exact test-polynomial bound leaves the bracket open in most rays by g = 8: 53-89 % of rays censored, depending on system and degree.
- Where the bracket closes, its median width log(g_hi / g_lo) is 0.04-1.05, larger at higher degree.
- As a predictor of where a law stops being decided, the bound is one-sided: it gives a guaranteed decided radius, not a guaranteed undecided one.

## 0F-4: design helps only where the problem is hard

Trajectories added until the target set is < 1/20 undecided (median; count of seeds censored at the cap of 31):

| law | start undecided | law design | D-optimal | random | space-filling |
|---|---|---|---|---|---|
| A2 (degree 4) | 1.0 | 2 (1) | 2 (3) | 4 (0) | 4 (0) |
| B2S (degree 3) | 0.12 | 1 | 1 | 1 | 1 |
| B2R (degree 3) | 0.39 | 1 | 1 | 1 | 1.5 |
| E2S (degree 3) | 0.95 | 7 | 7 | 9 | 8 |
| E2R (degree 3) | 1.0 | 16 (1) | 22.5 (0) | 31 (14) | 31 (12) |

- G-F4 fails. The pooled median ratio is 1.42 (bootstrap low 1.17) against a required 2. Without E4 it is 1.00.
- The pairs explain it. The 40 B pairs are settled by one trajectory under every arm, so their ratio is 1 and they are 40 of the 100 pairs.
- Where one trajectory is not enough, law design needs fewer:
  - A2: random / law median 1.58, 2 against 4;
  - E4: 1.52; on E2R, 16 against 31 (censored), with D-optimal at 22.5.

  These are post-hoc subsets and do not change the gate.
- Law design never needed more than D-optimal on median (both 2 pooled). On E2R it needed 6.5 fewer.
- The design claim as stated in the roadmap is dropped. What can be said: on the two laws whose target started fully undecided (A2 at degree 4, E2R), choosing starts for the law roughly halved the trajectories needed against random starts. That is a description of this run, not a passed test.

## What changes

- The headline object works: a calibrated confidence set (A, B at 2 % noise), an exact region certificate, and a zone map that never contradicts sampled models. It certifies laws over whole domains (93 of 93), which is the thing pointwise search cannot do.
- Its limits are now measured. It loses calibration when errors in the regressors dominate (E4, and A at 0.5 % noise). The horizon bound is one-sided. Design is not a general win.
- On real data the confidence set turns undecided regions into refuted ones (COVID degree 2). That is the misspecification signal from F0, now under a set with a probability meaning on synthetic data.

## Limits

- Probability wording holds only for A and B at 2 % noise. Everything on E4 and on real series is a set-membership statement.
- One feasible-GLS step, a linearised noise model, and the regressor-noise bias not corrected. A bias-corrected estimator is the obvious next step and is not done.
- Design used predicted rows with the per-equation filter covariance, not the full errors-in-variables covariance. There were 20 seeds per system, the target set used the truth, and starting data were the 14 training trajectories.
- Real series: one panel, reported cases, and iid noise assumed.
