# Roadmap v3: where the data decide the law

Date: 2026-11-26. This replaces `docs/ROADMAP_v2.md` for planning. v2, v1 and their addenda stay as history, and every result under them stands. Each new unit gets a frozen protocol before its first run, as before.

## Why change direction

The v2 program answered its questions, and the answers mostly cut against the method it started from:

- Good fits break laws everywhere (pilot G1, G-B1 0.974, G-R2 45/47, G-E2 35/35).
- On 2-4 state systems random sampling and dense grids find and fix what occurs (G-B3 fail, 0C soft penalty 320/320, 0D fine grid 80/80). The exact loop is not needed for search there.
- The violations sit where there is no data (hare-lynx witnesses on the box faces and at hare = 5; COVID data in the corner I <= 0.0154, S >= 0.89 of a simplex the laws are posed on).

The third point is the useful one. If the data say nothing about a region, then whether a fitted model obeys a law there is decided by the fitting procedure, not by the data. Repairing one model does not address that. The better question is which statements about a law the data actually support.

## New question

Given a dataset and a model class, at which states do *all* models that fit the data equally well obey a law, at which states do they all break it, and where can the data not tell? Which new measurement shrinks the third region fastest? How far from the data does the third region start?

## The core objects

Notation: phi(u) is the library evaluated at state u, Theta the library matrix on the data, xi the stacked coefficients, r(xi) the (generalised) residual sum of squares, xi_hat its minimiser on a given support, G = Theta' W Theta.

1. **Data-consistent set.** E = {xi : r(xi) <= tau}, an ellipsoid (xi - xi_hat)' G (xi - xi_hat) <= eps with eps = tau - r(xi_hat). There are two thresholds:
   - Rashomon: tau = (1 + delta) r(xi_hat), delta = 0.1, which matches "accuracy-equivalent" in every earlier gate. Its width does not shrink with n (check 5).
   - Confidence: tau = r(xi_hat) + sigma^2 chi2_{p, 1-alpha}. Its width shrinks like n^-1/2. Savitzky-Golay derivatives have correlated noise, so W comes from the known filter covariance (generalised least squares), not the identity.
2. **Pointwise law bounds.** A law at u is a(u)' xi >= -tol, with a(u) polynomial in u. Over E:
   min, max of a(u)' xi = a(u)' xi_hat -/+ h(u),  h(u) = sqrt(eps K(u)),  K(u) = a(u)' G^-1 a(u).
   (Cauchy-Schwarz. Check 1.)
3. **Three zones.** Let m(u) = a(u)' xi_hat + tol.
   - Certified: m(u) >= h(u). Every equally good model obeys the law at u.
   - Refuted: m(u) + h(u) < 0 (the largest value over E is still below -tol). Every equally good model breaks it at u.
   - Undecided: everything else.

   Refutation is pointwise. "Every model breaks the law somewhere" is weaker and is not what this zone means.
4. **Region certificate.** "Certified on all of D" is the polynomial system m(u) >= 0 and m(u)^2 - eps K(u) >= 0 for all u in D. It has degree 2 deg a, and the existing z3 / Bernstein checker decides it. This is where the exact checker becomes necessary: a grid of points cannot certify a continuum of models over a continuum of states.
5. **Selection does not break it.** With an absolute threshold tau, the full-library set contains every sparse model with r <= tau, so one certificate on the full library covers every STLSQ support (check 2). The price is a wider set when library columns are nearly collinear.
6. **One more measurement.** Adding a data row f changes K by a rank-one update: K_new(u) = K(u) - (a(u)' G^-1 f)^2 / (1 + f' G^-1 f) (check 3). The effect of any candidate measurement on the whole zone map is therefore cheap to compute. That gives measurement design aimed at a law, not at the coefficients.
7. **Law horizon.** Take a single-equation law, a(u) = phi(u), with a full degree-d library and data inside a centrally symmetric convex set K with gauge g(u). Then K(u) is the inverse Christoffel function of the data, and for g(u) >= 1:

   sqrt(eps / n) T_d(g(u)) <= h(u) <= sqrt(eps) N T_d(g(u)),  N^2 = sup over K of K(v).

   The lower side uses p(x) = T_d(l(x)), with l a supporting linear functional, as a test polynomial. The upper side uses the Chebyshev extremal property along the line through u. Check 6 confirms both on a disk (the upper side is tight at the boundary, and the lower side is loose by about 5x in K). Laws that couple equations through a function V, such as trapping, scale this by |grad V(u)|^2. Since T_d(g) grows like (g + sqrt(g^2 - 1))^d / 2, the half-width grows exponentially in d outside the data. The horizon is where h(u) overtakes |m(u)|.

   What I still need to prove or settle:
   - (i) the statement with the noise model written in, as a closed-form horizon radius in sigma, d, n;
   - (ii) whether a library degree above the degree of the law's margin always gives a finite horizon;
   - (iii) a version for libraries that are not full degree (A-mis).

All six numerical checks are in `tools/check_lawmap_math.py`, with output in `results/check_lawmap_math.txt` (ALL_OK True on 2026-11-26).

## Gates

Draft thresholds. The binding values go into `docs/protocol_0F.md` and are frozen before any 0F run.

| gate | question | pass rule (draft) | if it fails |
|---|---|---|---|
| G-F0 | implementation | on synthetic runs, every violating point found by sampling models inside E lies outside the certified zone: 100% | bug, stop |
| G-F1 | calibration of the confidence set (A, B, E4, truth in library) | share of runs with xi_true in E >= 0.90 at alpha = 0.05, over >= 100 runs per system | drop probability wording, keep the Rashomon version |
| G-F2 | real series (0R, 0E): how much the data decide | report certified / refuted / undecided shares of each law's domain. Prediction stated in advance: undecided share >= 1/2 of the law domain and <= 1/5 inside the data hull | descriptive, reported either way |
| G-F3 | horizon (A, B, H, E4 over noise x degree x n) | measured horizon along rays within the proved two-sided band in >= 95% of rays. If a closed-form radius is proved by the freeze: median abs log(measured / predicted) <= log 1.25 | report the band only |
| G-F4 | design: new trajectories chosen for the law vs random starts, space-filling starts and D-optimal starts | trajectories needed until the undecided share of D is < 1/20. Median ratio random / law-design >= 2 with paired bootstrap lower end > 1, and law-design no worse than D-optimal on median. Cap 30 trajectories | drop the design claim |
| G-F5 | COVID retrospective, not gating | map from the early window vs later data: did later states land in zones that were certified, and did the refit agree? | exploratory only |

A measurement in G-F4 is a short trajectory from a chosen initial condition, sampled like the original data. Nobody can measure a derivative at an arbitrary state, but in a lab system or a simulation one can choose where to start.

## Work units

| unit | scope | output | where |
|---|---|---|---|
| F0 | zone maps for stored fits (A, B, H, E4 samples; hare-lynx primary; 10 COVID countries), Rashomon version only, no thresholds | maps plus undecided shares, used to set the 0F design | local, small |
| 0F protocol | write and freeze `docs/protocol_0F.md`: noise model and W, tau, alpha, domains, design arms, gate thresholds | frozen commit | local |
| 0F-1 | G-F0, G-F1 | calibration tables | Kaggle |
| 0F-2 | G-F2 on 0R and 0E stored fits | real-data zone tables and figures | Kaggle |
| 0F-3 | G-F3 sweep | horizon vs band | Kaggle |
| 0F-4 | G-F4 design comparison | trajectories-to-decision table | Kaggle |
| T | proposition 7 written out in full, then (i)-(iii) | `docs/law_horizon.md` with proofs and the numerical checks they rely on | by hand, checked by script |
| E5 (optional) | five-state system with a simplex checker for 5 variables | does law-design's edge grow with dimension? | Kaggle |

## Calendar

Target: project finished on 2027-01-15. Results freeze on 2027-01-03, so the last two weeks are write-up only.

- Nov 26 to Nov 30: F0 (done, see `docs/F0_result.md`), then write and freeze protocol 0F.
- Dec 1 to Dec 10: 0F-1 (G-F0, G-F1) and 0F-2 (real series) on Kaggle. Write out proposition 7 in full.
- Dec 11 to Dec 22: 0F-3 (horizon) and 0F-4 (design) on Kaggle. Work on (i), the closed-form radius.
- Dec 23 to Jan 2: robustness addenda, (ii) and (iii) if time. E5 only if every gate above is done.
- Jan 3: results freeze. Nothing new gets run after this date.
- Jan 4 to Jan 15: figures, write-up, preprint, repo cleanup. Finished on Jan 15.

## Kill rules

- G-F0 fails: stop and fix before anything else runs.
- G-F1 fails even with the filter covariance: probability claims are dropped, and the Rashomon version carries the project.
- G-F4 fails: the result is the zone map and the horizon, without a design claim.
- Proposition 7 does not hold in some case the checks did not cover: say so, keep the counterexample, and narrow the statement.

## What carries over

- The checker (z3 plus Bernstein, `bern_nd.py`), the laws, the testbeds, both real series, the Kaggle runners and the freeze discipline.
- The v2 results become the motivation and the control set. The law taxonomy is one section. 0C and 0D are the evidence that pointwise search is easy here, and the region certificate in item 4 is what search cannot do.
- Not continued: more gates on loop vs grid vs random at 2-4 states. S1 (Navier-Stokes closure) stays parked.

## Derivations I should be able to do without notes

Items 2, 3, 5 and 6, and both sides of the bound in item 7.

## Prior work to read and position against

All of these need confirming against the papers before they are cited:

- Data informativity for analysis and control (van Waarde, Eising, Camlibel, Trentelman), which works with the set of all systems consistent with the data.
- Dissipativity verification for polynomial systems from noisy data with SOS (Martin and Allgower, arXiv 2008.03231).
- Christoffel-Darboux kernels for data analysis (Lasserre, Pauwels, Putinar).
- Optimal design and the Kiefer-Wolfowitz equivalence theorem, the baseline for G-F4.
- Rashomon sets (Semenova, Rudin, Parr).
- Uncertainty in sparse model discovery (ensemble SINDy, Fasel et al.; Bayesian SINDy, Hirsh et al.).
- Counterexample-guided learning (Chang, Roohi, Gao; Abate et al.).
- Constrained SINDy (Loiseau and Brunton; Kaptanoglu et al.).

What I claim as mine has to survive that reading: law zones for discovered equations, the region certificate over the whole consistent set, the horizon bound, and design aimed at deciding a law, tested on synthetic systems and two real series.

## Addendum 1 (2026-11-26): F0 result and changes to the 0F draft

F0 is in `docs/F0_result.md`. It changes the draft gates as follows, before anything in 0F is frozen:

- The confidence set is the main object of 0F. The Rashomon set (factor 1.21 on the residual sum of squares) is reported next to it. F0 showed the Rashomon set is wide: stored models use at most 71 % of its half-width.
- G-F4 is restricted to laws with a positive true margin (A2, H2, B2S, B2R, E2S, E2R). Zero-margin face laws can only be certified where h < tol, which F0 shows the data rarely allow, and structure handles them anyway.
- When the data lie on an invariant set, the law domain is restricted to that set. For COVID this is S + I + R = 1, where the library matrix has rank 6 of 10 at degree 2.
- New item for T: "a refuted zone implies the truth is not in the library, if the law holds" is immediate from the definitions. F0 saw refuted zones only for H and COVID. How often a wrong library shows no refuted zone (A-mis shows none) is a G-F2 measurement.
- 0F fixes the set definition for real series: training data only, as in F0. Validation-selected models outside that set are reported, not dropped.

