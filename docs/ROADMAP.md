# Roadmap

Where things stand: pilot done, verdict REVISE (`pilot_report.md`). The claim to test next is narrower than "beats validation": exact check plus repair of discovered equations, and what happens when the library cannot express the law.

## prior work this has to sit next to
(references from memory and search snippets, confirm each before citing)
- Counterexample-guided learning with SMT for Lyapunov functions: Chang, Roohi, Gao (NeurIPS 2019); Abate, Ahmed, Giacobbe, Peruffo (2020). Same loop, different object. They synthesize certificates, this repo checks and repairs discovered equations.
- Constrained SINDy by construction: energy-preserving quadratic terms (Loiseau and Brunton 2018), trapping theorem for global stability (Kaptanoglu et al. 2021, `TrappingSR3` in PySINDy).
These are the baselines. If they beat CEG repair everywhere, the contribution becomes the benchmark and the certification result. PySINDy needs scikit-learn, which is blocked in the sandbox I ran in, so the baseline is either my own alternating-SDP version or a run on another machine.

## phase 0, Nov 9 to Dec 6: close the pilot gaps
1. Baselines: energy-preserving constraint (linear equalities) and a trapping-style SDP refit against CEG on the same cells.
2. Repair with random-sample counterexamples (the missing control).
3. Needle test: sweep the volume of the violating region, plot detection of random testing vs the exact check.
4. Misspecification family: drop each true term in turn; truth outside the library (Holling II predator-prey).
5. 20 seeds, bootstrap CIs, margin sweep for weak-tolerance repair.
6. Tests in CI, tag `v0.1-pilot`.
Gate A: CEG repair is at least as good as the baselines on admissibility at equal held-out error, or it wins in the misspecified arm. If not, go benchmark-first.

## phase 1, Dec 7 to Jan 24: fluid-scale polynomial models
- Models: Lorenz-63 (trapping ellipsoid known), the 9-mode shear flow model of Moehlis, Faisst and Eckhardt (2004), 2D Kolmogorov flow Galerkin truncations with 5-9 modes.
- Laws: nonlinear term conserves energy (coefficient identities), dissipation (E' <= -mu E + C outside a ball), shift/reflection symmetry.
- Verifier: z3 does not scale past a few variables. Use SOS (cvxpy with Clarabel or SCS) for degree 4 in up to about 10 variables, linear algebra for the quadratic forms, branch-and-bound for cubic forms on the sphere. Rationalize and re-check certificates exactly.
Gate B: an n = 9 law checked with a certificate in under 60 s. If not, cut the fluid track to n <= 5.

## phase 2, Jan 25 to Apr 4: Navier-Stokes closure track
- DNS: 2D Kolmogorov flow with a pseudo-spectral solver (numpy FFT, 64^2 to 128^2). Validate on Taylor-Green decay and the energy budget (residual < 1e-3) before using any data.
- Reduced model: Fourier truncation to r = 8-16 modes, learned closure for the discarded modes.
- Laws for the closure: net dissipative (a . c(a) <= 0 on the energy shell), bounded trajectories, symmetries. All polynomial sign conditions, so the same spec form applies.
- Compare: Smagorinsky-type closure, unconstrained sparse closure, trapping-constrained, CEG-repaired. Metrics: long-run statistics (energy spectrum, dissipation PDF), short-term error, fraction of initial conditions that blow up, certificate or witness per model, repair cost.
Gate C: DNS validated. Gate D: at least one of (i) a learned closure with baseline-level statistics that is certified where the baseline is not, (ii) a certificate that no model in a library is admissible.

## phase 3, Apr 5 to May 8: freeze and write up
Re-run everything from a tagged commit, archive with a DOI, preprint draft. Stretch (not promised): rigorous bounds on blow-up in learned closures.

## the Navier-Stokes problem itself
The Clay problem (global regularity in 3D) is not a target of this project. The honest connection is certified reduced models and closures of Navier-Stokes. Computer-assisted proofs about fluid singularities are an active area (recent work on Euler and Boussinesq); read it for method, do not promise results.

## discipline
- A protocol file per phase before running, with thresholds and a decision rule. Later changes go in dated addenda, headline criteria stay as written.
- Every table comes from `src/analyze.py`.
- Failed criteria go in the same table as passed ones.

## risks
| risk | plan |
|---|---|
| SOS does not scale or is not exact | rationalize + exact re-check; shrink n; keep Bernstein/z3 for n <= 3 |
| baselines already solve repair | benchmark + certification becomes the result |
| DNS time on a laptop | 64^2 grid, short runs, reuse snapshots |
| pilot tolerance too loose (4-18% of strong-verified models still misbehave) | verify at weak tolerance in phase 0 |
| novelty overlaps prior work | state the delta explicitly in the writeup, cite the baselines |


## Addendum 2026-11-08: revised plan

Why this changed. Two things came up after the pilot. (1) Checking the competition record: the top award goes to one of the category winners, and past winners had a method with a measurable advantage and a clear user. (2) The state of the Navier–Stokes question: in September 2026 an AI-agent claim of finite-time blow-up was reported, and reports disagree on how far the machine check goes. A student project cannot build on that. The pilot also showed that validation metrics and random testing already catch most violations on synthetic data. So the claim has to rest on real data, on a misspecified case, and on a fair comparison with simple by-construction fixes.

Changes to the original plan:

| item | original | revised |
|---|---|---|
| Navier–Stokes | phase 2 closure target | dropped as a target. Motivation only. |
| SOS verifier | phase 1 requirement | optional. The z3 verifier covers up to 3 state variables. Use SOS only for a model with 4 or more. |
| fluid models (Lorenz-63, shear flow, Kolmogorov) | phase 1 | optional stretch |
| real data | none | headline: hare–lynx pelt records (Hudson's Bay), then one epidemic time series. Data source to be confirmed. |
| baselines | constrained SINDy | add by-construction fixes: fitting log-states for positivity, fraction-state parametrisation for conservation, soft-penalty training. A reviewer will ask why these were not used. |
| misspecification | arm A-mis; Holling II in phase 0 | Holling II is the headline synthetic test (true law outside the library), then the real data. |
| Gate A | after all of phase 0 | after 0A (baselines and random-counterexample control) |
| repair target | strong tolerance only | add a weak-tolerance repair arm. The strong tolerance leaves 4–18% of verified models unphysical in trajectory tests. |
| board deliverables | none | live demo (best-fit model goes negative; repaired and by-construction models stay positive); one-page explanation written by the student; tests; reproducible repo |

Decision rules added:
- Gate A. If by-construction fixes match or beat repair on admissibility at equal held-out error, the claim becomes: when verification is needed and when a by-construction fix is enough, with the rule tested on systems not used to derive it. That is a result, not a dead end.
- Gate R. At least one real dataset must show a validation-selected model violating a stated law under extrapolation. If none does, drop the real-data claim and report the synthetic results only.

Timeline. Dates assume an affiliated fair in late February. Set the results freeze at that fair's deadline minus three weeks.
- 0A, Nov 9–Dec 6: by-construction baselines, constrained SINDy, random-counterexample control. Gate A.
- 0B, Dec 7–Jan 10: Holling II misspecification, needle test, 20 seeds with confidence intervals.
- 1, Jan 11–Jan 31: hare–lynx. Gate R.
- 2, Feb 1–Feb 21: epidemic series. Results freeze. Fair write-up, written by the student.
- 3, Feb 22–Apr 4: robustness checks and the board demo. Any new experiment is registered before it runs.
- May 2027: Regeneron ISEF, Los Angeles.

Judging reference. The Society says projects are judged on creativity, innovation and depth of scientific inquiry. The Innovator Award is for innovative research with potential impact in its field and beyond. Recent winners of the top $100,000 award: 2025, a new synthesis route for an antiviral; 2026, a simulation program for protein-like folding. Earlier top awards: 2023 ($75,000), a method for finding exoplanets that orbit close to their stars; 2024 ($75,000), an organic electrochemical transistor.

Odds. About 22 category winners compete for the top prize, and about 1,725 finalists compete overall. Being a finalist does not make a win likely, and no plan makes it near-certain. The levers are a clear test that can fail, a real case where the failure changes a decision, a simple fix that is compared honestly, and a student who can explain every step.
