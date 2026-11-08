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
