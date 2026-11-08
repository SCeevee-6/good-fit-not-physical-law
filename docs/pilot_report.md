# Pilot report

Question: if an equation found from data fits held-out data as well as the best one, can it still break a physical law, and does an exact check catch that better than ordinary validation?

Verdict from the protocol's own rule: REVISE. The gap is real, the checker is sound and the repair works. "Better than ordinary validation" failed its test.

## what was run
- Two polynomial systems with known truth: logistic predator-prey (A) and SIR with varying total population (B). A third arm, A-mis, is system A with the x^2 term removed from the library.
- Candidates: sparse regression (STLSQ), 2-3 library degrees, 12 thresholds, 9 data resamples per cell. 45 cells (3 noise levels, 5 seeds, 3 arms), 5,553 candidates. 3,924 of them are "accuracy-equivalent": held-out error within 10% of the best in the cell.
- Laws are polynomial sign conditions that are linear in the coefficients: populations stay >= 0, trapping region (A), S' <= 0 and R' >= 0 (B), total population conserved (B). Checked at a weak and a strong tolerance. Strong is what I report.
- Checker: z3 with exact rational coefficients, Bernstein branch-and-bound when z3 times out, replay of every counterexample as a trajectory, and an oracle that does not use either (dense grid, differential evolution, long trajectories).

## results

| criterion | result | threshold | |
|---|---|---|---|
| G1 accuracy-equivalent models breaking a law | 42-82% in all 6 groups, 93-100% in A-mis | >= 20% in >= 3 of 6 | pass |
| G2a checker recall / false flags | 1.00 / 0 of 4,477 | >= 0.95 / <= 0.05 | pass |
| G2b best ordinary metric, AUROC | AIC 0.966 | <= 0.70 | fail |
| G3 soundness | 99.51% replayed, 99.97% of queries decided, 0 disagreements in 14,212 cross-checks | >= 99%, >= 95%, 0 | pass |
| G4 repair | 99.25% repaired, median error x0.9993, random constraints 17% | >= 80%, <= +10%, +20 points | pass |

Figures are in `figures/`, per-criterion numbers in `results/criteria_verdicts.csv`.

## what it means
- Violations rise with noise and with the number of terms. Models with no extra terms essentially never violate. Among parsimonious models (k <= k_true+4) the rate is 4-5% at 0.5% noise, so G1 is weaker there (4 of 6 groups reach 20%).
- The model picked by held-out error alone breaks a law in 6 of 15 predator-prey cells and 6 of 15 SIR cells. The AIC pick: 1 of 15 and 0 of 15. In A-mis both picks break a law in all 15 cells.
- AIC works well here because the true law is sparse and inside the library. At fixed complexity (k_true+5 to k_true+12 terms) ordinary metrics are near chance (AUROC 0.50-0.61).
- Random testing is a strong baseline. 200 points catch 97.6% of violations, 2,000 catch 99.4%, the checker 100%. The gap is at the margin: SIR violations within 1.5x the tolerance (n = 73), 200 random points find 27%.
- The 44 counterexamples that failed strict replay (0.49%) are all finite-time blow-ups of the candidate. Under the amended rule replay is 100%. The protocol logs this as a post-hoc change.
- Repaired models are admissible but not the true law (median support F1 0.44). Only 20% of them also pass at the weak tolerance.

## limits
- Two polynomial systems, laws and tolerances chosen by hand. Nothing about trig terms, neural models or real data.
- Truth is inside a sparse library for A and B, which favours AIC. One kind of misspecification only.
- The label is "checker counterexample (replayed) or oracle violation", so checker recall of 1.00 is partly by construction. The oracle agreed with it on all 4,334 candidates it ran on, which is the independent part.
- Not done: a repair loop that takes its counterexamples from random sampling instead of z3. Random testing detects almost as well, so the repair gain may come from the feedback loop and not from the SMT solver.

Protocol and every amendment: `pilot_protocol.md`.
