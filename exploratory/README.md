# Which energy measurements distinguish opposing physical predictions?

Exploratory simulation, independent of the frozen 0F gates. Reproduce from the repository root with Python 3, NumPy and SciPy:

    python exploratory/study.py
    python exploratory/summarize.py

The first command writes 7,200 policy evaluations in trials.csv; the second reproduces the committed summary.csv. Seeds are fixed. These are simulated experiments.

## Physics and observation

Write s=x²+v² and g=a+b(s-1)+cx². Simulate x'=v, v'=-x-0.2v+kvg, and meter'=v²g. The meter reports the integrated auxiliary channel in calibrated units. The channel contributes k times this reading to mechanical energy; k=0 is a virtual probe. The tested sign is that of the auxiliary gain at x=0, s=1.5. It is distinct from net energy change in the damped oscillator. Exclude the indifference region |a+0.5b|<0.08. Integrate by RK4, with one refined-step cross-check per coupling.

Eight training readings at ordinary energy, independently noisy across phases, screen a 150-member grid by a 95% Gaussian residual test with eight degrees of freedom. Each experiment draws a truth with replacement, alternating signs. Each cell has 120 in-grid trials and 120 trials from a separate seeded bank of randomly drawn coefficients. The latter use the same functional form. Seeds are shared across coupling and noise cells. The evaluator retains truth predictions separately; only admitted grid members determine schedules and thresholds.

The design compares three initial energies, eight phases and two schedules. One schedule takes four cumulative readings at 0.25, 0.5, 1 and 2. The other averages four independent half-unit experiments at the same start. Runtime and number of reads match; repeat costs three additional resets. Meter SD is 0.02, with known within-trace equicorrelation rho. Separate resets have independent errors. The detector whitens readings, projects along the difference of the surviving class centroids and places a threshold midway between opposing projected extremes. It optimizes its conditional worst-case accuracy over the permitted schedules. This projection is a heuristic. Overlap of projections triggers a coin flip and cannot prove indistinguishability of the full readings. If training leaves just one sign it is reported directly; an empty set triggers a coin flip. Mean accuracy integrates Gaussian test noise analytically and includes these fallback cases. Conditional floors apply only to admitted grid members.

## Result and interpretation

At k=0.6 and rho=0.6 on off-grid trials, mean accuracy was 0.570741 for a fixed four-reading trace, 0.967432 for a designed four-reading trace, 0.974994 for four fixed repeats, and 0.974998 for unrestricted design. The fixed protocols start at maximum speed at the largest allowed energy. In 113 of 114 still-ambiguous trials the unrestricted policy selected repeats. Across all twelve reported strata its largest absolute gain or loss versus fixed repeats was below 0.001. The absolute comparison depends on the chosen detector, dynamics, allowed starts and noise model.

The positive result concerns phase selection when resets are expensive. The negative result is equally important: adaptive schedule selection has no material average advantage over a well-chosen fixed repeat policy here. The simulation uses specified coefficient ranges, a known coupling, known noise and a finite grid. Coverage of an individual grid truth is a training-set property; off-grid calibration and worst-case continuous-family guarantees require different tests. An appropriate next challenge is a disjoint model family, unknown meter drift shared across resets, and explicit reset cost. Compare against a single best fixed phase selected on development data before interpreting a benefit as adaptation.
