# F0 result: law zones for the stored fits

Run 2026-11-26, locally (about 13 minutes). It is exploratory and has no gates, as set in `docs/ROADMAP_v3.md`. Code: `src/lawmap.py`, `src/run_F0.py`, `src/run_F0_plane.py`, `src/analyze_F0.py`. Outputs: `results/F0_zones.csv` (419 rows, one per data set x library degree x law), `F0_summary.csv`, `F0_covid_plane.csv`, `F0_witness_zones.csv`, `F0_checks.txt`, `F0_maps.npz`. Figure: `figures/fig_F0_zones.png`.

## What was computed

- Data sets: A, A-mis, B, H and E4 at 2 % noise, seeds 0-4, the same cells as the 0B and 0D pools; hare-lynx at the primary config; the 10 COVID countries at the primary window.
- Set: every model on the full library of a given degree whose training residual sum of squares (each equation scaled by its derivative sd) is within a factor 1.21 of the least-squares minimum. A factor of 1.21 on the sum of squares is the 10 % RMSE rule of the earlier gates. A-mis keeps x^2 out of the library.
- For each law and each point of a grid on its domain: m = law value of the least-squares fit plus tolerance, and h = sqrt(eps K). The zone is certified if m >= h, refuted if m + h < 0, otherwise undecided. The tolerance is the strong one.
- "Near data" means within 0.25 state sd of a training point. Only 171 of the 419 rows have any law-domain point near data.

## Checks

| check | value | required |
|---|---|---|
| stored equivalent models inside the set | 96 % (13913 / 14425 model x law rows) | descriptive |
| largest member deviation / h over all grids | 0.706 | <= 1 |
| models with a stored counterexample on a law the set certifies everywhere | 0 | 0 |
| stored witnesses of in-set models that fall in the certified zone | 0 of 10 | 0 |

The math and the code agree on every stored model. The stored models also use at most 71 % of the half-width, so the set is wider than what STLSQ explores.

## Results

Zone shares at each library degree, averaged over seeds or countries and over laws:

| system | degree | certified | refuted | undecided |
|---|---|---|---|---|
| A | 2 / 3 / 4 | 0.333 / 0 / 0 | 0 / 0 / 0 | 0.667 / 1 / 1 |
| A-mis | 3 / 4 | 0.186 / 0 | 0 / 0 | 0.814 / 1 |
| B (SIR) | 2 / 3 | 0.243 / 0.065 | 0 / 0 | 0.757 / 0.935 |
| H (Holling truth) | 2 / 3 / 4 | 0.398 / 0 / 0 | 0.353 / 0.047 / 0 | 0.249 / 0.953 / 1 |
| E4 (SEIR) | 2 / 3 | 0.010 / 0.001 | 0 / 0 | 0.990 / 0.999 |
| hare-lynx | 2 / 3 / 4 | 0 / 0 / 0 | 0 / 0 / 0 | 1 / 1 / 1 |
| COVID, on S+I+R = 1 | 2 / 3 | 0.076 / 0.032 | 0.013 / 0.013 | 0.911 / 0.955 |

1. **The data decide little.** Away from the data, almost every law is undecided. Hare-lynx is undecided on 100 % of every law domain at every degree, and E4 on 99 %.
2. **They decide near the data.** Certified share near data against the whole domain: SIR 0.58 vs 0.17, SEIR 0.26 vs 0.005. Holling is 0 in both, because near its data the laws are refuted or undecided.
3. **A larger library decides less.** On A the certified share goes 0.333, 0, 0 for degrees 2, 3, 4. The A trapping law A2 is certified everywhere at degree 2 in all 5 seeds and undecided everywhere at degrees 3 and 4. This is the growth in item 7 of ROADMAP_v3: the half-width outside the data grows like T_d.
4. **A refuted zone means the library is wrong.** Refuted means every equally good model in the library breaks the law there. It appears only where the truth is outside the library: H at degree 2 (0.353; on the prey-positivity face H1x it is 0.76) and at degree 3 (0.047), and COVID on the plane (0.013). A, B and E4, whose truth is in the library, show 0 refuted at every degree. A-mis also shows 0, so the absence of a refuted zone does not show the library is right.
5. **Laws that hold with zero margin are almost never certified.** On a face such as x = 0 the true x' is exactly 0, so the margin is only the tolerance, and certification needs h < tol there. A1x and A1y are undecided on 100 % of their faces at every degree. The exceptions are faces the trajectories come close to: B1I is 0.36 certified at degree 2, because SIR runs end near I = 0. This is why these laws are the ones structure enforces in `docs/law_taxonomy.md`. For them the data cannot do the job.
6. **COVID data cannot decide anything off the plane.** The JHU fractions satisfy S + I + R = 1 exactly, so the library matrix has rank 6 of 10 at degree 2. Any multiple of (1 - S - I - R) is invisible to the data, and every law point off the plane has infinite half-width (`F0_zones.csv`, h_inf_share = 1). On the plane (`F0_covid_plane.csv`) the laws are 91 % undecided at degree 2.
7. **Witnesses.** Of the 29 stored witnesses (5 hare-lynx, 24 COVID), 28 sit in the undecided zone and 1 in the refuted zone (Turkey, B1S). None is certified. Only 10 of the 29 come from models inside the training-data set. The others were picked on validation error and fit the training data more than 21 % worse than least squares.

## What this changes for protocol 0F

- The Rashomon set at a factor of 1.21 is wide and leaves little decided. 0F needs the confidence set (G-F1) as its main object, with the Rashomon set reported next to it.
- G-F4 as drafted ("undecided share below 1/20") cannot be met for zero-margin laws far from the data. It should be restricted to laws with positive true margin (A2, H2, B2S, B2R, E2S, E2R), or measured as h relative to the margin.
- When the data lie on an invariant set (the COVID plane), law domains are restricted to that set. Otherwise the result is trivially "undecided".
- For real series, 0F must fix whether the set is defined on training data, validation data or both, since 19 of 29 stored witness models are outside the training-data set.

## Limits

- Grids, not certificates. A zone share is a share of grid points. The exact region certificate is 0F.
- Five seeds per synthetic system, one noise level (2 %), the full library only (no sparse supports), and no probability statement: the set is defined by a fixed ratio.
