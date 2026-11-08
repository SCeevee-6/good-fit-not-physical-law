# Pilot protocol — Counterexample-guided verification of AI-discovered ODEs

**Project question.** Can counterexample-guided mathematical verification identify and eliminate physically inadmissible
differential equations that achieve excellent held-out predictive accuracy, and can it do so more reliably than
conventional validation?

**Status:** pre-registered before any pilot result exists. Written 2026-10-27. Anything changed after data are seen is
logged in section 9 ("Deviations").

---------------------------------------------------------------------------------------------------

## 1. Pilot purpose
The pilot is a feasibility test, not the final study. It asks whether four conditions that the full project depends on
actually hold in a controlled testbed where ground truth is known:

| ID | Hypothesis |
|----|-----------|
| H1 | *The gap exists.* Among candidate equations that are statistically indistinguishable from the best on conventional held-out accuracy, a non-trivial fraction violate a stated physical law. |
| H2 | *The verifier is sound and decisive.* SMT counterexamples replay as real behavioural violations, and nearly all verification queries return a verdict (no time-outs). |
| H3 | *Repair works.* Feeding counterexamples back into a constrained refit removes violations at small cost in accuracy. |
| H4 | *Reliability.* The verifier detects inadmissibility more reliably than conventional validation metrics (and than randomised property testing at a matched evaluation budget). |

## 2. Testbeds (polynomial ODEs with known ground truth)
All parameters are fixed here.

**System A — logistic predator–prey**, state (x, y) = (prey, predator):
x' = 1.0 x − 0.4 x² − 0.5 x y ,  y' = −0.6 y + 0.4 x y. Coexistence equilibrium (1.5, 0.8) is a stable spiral.
Training ICs: x0 ~ U[0.3, 3.0], y0 ~ U[0.2, 2.0]; 14 train, 6 i.i.d. held-out, 6 OOD trajectories.
OOD ICs: x0 ~ U[3.0, 5.0], y0 ~ U[2.0, 3.5]. t in [0, 20], dt = 0.05.

**System B — SIR epidemic**, state (S, I, R):
S' = −0.5 S I , I' = 0.5 S I − 0.2 I , R' = 0.2 I.
Total population N varied per trajectory, N ~ U[0.7, 1.3], to avoid exact collinearity of the library (S+I+R = N would otherwise make
polynomial features linearly dependent). Training ICs: I0 ~ U[0.01, 0.15]·N, R0 ~ U[0, 0.15]·N, S0 = N − I0 − R0.
14 train, 6 i.i.d. held-out, 6 OOD (N ~ U[1.3, 1.5]). t in [0, 60], dt = 0.25.

**Noise.** Gaussian measurement noise on states, sigma = rho · sd(state over the training data), rho ∈ {0.005, 0.02, 0.05}.
**Seeds:** 5 per (system, noise) → 30 cells.

## 3. Physical specifications (the "laws")
Each law is a polynomial sign condition on the vector field f(x; xi) = Theta(x) xi that is **linear in the coefficient vector xi**.
Violation is declared at two tolerances: *weak* and *strong*. **The strong tolerance is the headline** (to avoid the
tautology that any non-zero spurious coefficient "violates" a law). Weak is reported as a sensitivity analysis.

| Spec | System | Statement | Domain D | tol (weak / strong) |
|------|--------|-----------|----------|---------------------|
| A1x | A | Positivity invariance: f_x(0, y) ≥ −tol | y ∈ [0, 6] | 1e-3 / 2e-2 |
| A1y | A | Positivity invariance: f_y(x, 0) ≥ −tol | x ∈ [0, 7.5] | 1e-3 / 2e-2 |
| A2  | A | Trapping region: dV/dt ≤ tol on shell V = 0.4x + 0.5y = 3, x,y ≥ 0 (V = conversion-weighted biomass; true system has dV/dt ≤ −0.5 on the shell) | shell | 1e-3 / 2e-2 |
| B1S, B1I, B1R | B | Positivity invariance on each face of the orthant, f_i ≥ −tol | {x ≥ 0, S+I+R ≤ 1.5} | 1e-3 / 2e-2 |
| B2S | B | Susceptibles never increase: S' ≤ tol | same | 1e-3 / 2e-2 |
| B2R | B | Removed never decrease: R' ≥ −tol | same | 1e-3 / 2e-2 |
| B3  | B | Conservation of total population: abs(S'+I'+R') ≤ tol | same | 5e-3 / 2e-2 |

Domains deliberately extend **beyond** the training data (the physical law holds everywhere, the data do not).
The true systems satisfy every spec at the strong tolerance (checked as a unit test in step 4; any failure is a protocol bug).

## 4. Candidate pool ("AI-discovered" equations)
Own implementation of STLSQ sparse regression on Savitzky–Golay derivative estimates (window/poly order fixed in code, logged).
Per cell: library degree d ∈ {2,3,4} (A) or {2,3} (B) × 12 thresholds (log-spaced 1e-3 … 3e-1) × 8 trajectory-bootstrap resamples,
de-duplicated by (support, rounded coefficients). Roughly 200–300 candidates per cell.

**Conventional validation battery (computed for every candidate)**
1. i.i.d. held-out derivative RMSE (and R², reported but not used as a filter because derivative R² is capped by noise);
2. held-out trajectory rollout RMSE from in-distribution ICs;
3. OOD-IC rollout RMSE (a "stress-test" validation, stronger than usual practice);
4. AIC-style score (n log RSS/n + 2k).

**Accuracy-equivalent set (AES):** candidates whose i.i.d. held-out derivative RMSE is within 10 % of the best RMSE in the same cell.
This is the population in which "excellent held-out accuracy" cannot separate models.

## 5. Verifier and oracle
**Verifier (SMT).** Coefficients → exact rationals; spec polynomial → z3 real arithmetic (nlsat, complete for polynomial real arithmetic);
10 s time-out per query. Output: `verified`, `counterexample` (exact rational point + severity), or `timeout`.

**Replay.** Every counterexample is replayed by integrating the candidate model from (a point shifted into the domain interior by 1e-4 for face specs, or the
counterexample itself on the shell) and checking the *behavioural* consequence (negative population beyond 1e-6, V leaves the trapping set, S increases, total mass drifts).

**Independent oracle (SMT-free).** Dense grid evaluation (≥1e4 points/spec), differential-evolution adversarial search, and long-horizon simulation from boundary-adjacent
and OOD ICs, judged on integrated behaviour (negative states, escape from the sublevel set, S rising, mass drift > tol·T).

**Label.** A candidate is *inadmissible* if (SMT counterexample at strong tolerance AND replay-validated) OR (oracle violation at strong tolerance).
Otherwise *admissible* if SMT verified and oracle finds nothing. Time-outs are excluded from detection metrics and reported separately.
*Circularity note:* spec-level labels depend on the specs I chose; the oracle and replay show that violations are real behaviour of the candidate ODE, not artefacts of the encoding, but cannot show that the spec itself is the right law.

## 6. Counterexample-guided repair (CEGIS)
Learner = constrained STLSQ (cvxpy QP; each counterexample adds one linear inequality row on xi); Verifier = z3. ≤ 15 iterations, support re-selected each iteration.
Arms at **equal constraint budget**: (1) CEG; (2) random-constraint refit (rows at random domain points); (3) post-hoc rejection only; (4) unconstrained, validation-selected best.
Repair target = verified at **strong** tolerance (weak-level verification reported as secondary).

## 7. Detection comparison (H4)
Scores: iid derivative RMSE, rollout RMSE, OOD rollout RMSE, AIC (higher = more suspicious); random-sample falsification (matched number of point evaluations);
verifier (binary). Reported on the AES and on the full pool: AUROC and TPR at 5 % FPR against the label, with bootstrap CIs over cells.

## 8. Pre-registered go / no-go criteria
| Crit. | Statement | Pass threshold |
|-------|-----------|----------------|
| **G1** (H1) | Violating fraction of AES at strong tol., pooled over seeds, per (system × noise) group (6 groups) | ≥ 20 % in ≥ 3 of 6 groups |
| **G2a** (H2) | Verifier vs. label: TPR (recall of label) and FPR (SMT counterexamples that fail replay) | TPR ≥ 0.95 and FPR ≤ 0.05 |
| **G2b** (H4) | Best conventional validation metric (max over the four metrics in §4), AUROC for the label within AES | ≤ 0.70 |
| **G3** (H2) | Replay-validated fraction of SMT counterexamples; decided fraction of queries; disagreements with exact 1-D root isolation | ≥ 99 %; ≥ 95 %; 0 |
| **G4** (H3) | Violators repaired within 15 iterations; median relative increase in iid held-out RMSE; CEG vs. random-constraint arm | ≥ 80 %; ≤ 10 %; CEG ≥ 20 pp higher repaired fraction, paired McNemar p < 0.05 |
Supplementary (not pass/fail): TPR@5 %FPR gap between verifier and best conventional metric; random-falsification TPR; severity-stratified gap; weak-tolerance sensitivity.

**Decision rule.**
* **GO** to the full study if G1, G2a, G2b, G3 pass (G4 failing → GO with the repair method flagged for redesign).
* **REVISE** if G1 and G3 pass but G2a or G2b fail (e.g. OOD validation already detects most violators → sharpen the claim around *which* violations validation misses).
* **NO-GO** (reframe) if G1 fails (no gap to close) or G3 fails (verifier unsound).

## 9. Deviations from the approved plan (logged before any result)
1. The plan's "R² ≥ 0.99" filter for the accuracy-equivalent set is dropped from the definition (derivative R² is capped by measurement noise at rho = 5 %, so it would empty the set); the 10 %-of-best-RMSE rule is the sole criterion and R² is reported.
2. Violation is graded at two tolerances (weak/strong) so that the gap cannot be produced by negligible coefficients; strong is headline.

*(Further deviations are appended below with date and reason.)*

### Addendum A — sensitivity subsets (registered 2026-10-27, after building the candidate pool but **before** any verifier or oracle result existed)
Seeing that the accuracy-equivalent set (AES) contains models with up to 30 (A) / 60 (B) non-zero terms, the headline G1 number could be attributed to
implausible over-parameterised models. The following subsets are therefore reported alongside the primary AES (G1 is still judged on the primary AES):
* **S2 AES-strict:** i.i.d. derivative RMSE within 5 % of best **and** in-distribution rollout RMSE within 1.25× the best in the cell.
* **S3 AES-parsimonious:** AES with k ≤ k_true + 4 non-zero terms (k_true = 5 for A, 4 for B), i.e. models a practitioner might plausibly report.
* **S4 Validation pick:** the single candidate with lowest i.i.d. derivative RMSE per cell (the model conventional validation would select).
* **S5 AIC pick:** the single candidate with lowest AIC score per cell.
For S4 and S5 the quantity reported is the fraction of the 30 cells in which the selected model is inadmissible at the strong tolerance.

### Addendum B — oracle and baseline definitions (registered 2026-10-27 before any oracle output was inspected)
**Local oracle (SMT-free, vector-field level).** For every spec: (i) dense grid of the spec domain (1-D: 4001 pts; 2-D: 101x101; 3-D: 41^3, domain-filtered),
(ii) differential-evolution adversarial minimisation of G (popsize 20, maxiter 80, polished), both evaluating the candidate's right-hand side with an implementation
independent of the verifier (direct monomial products, no sympy/z3/Bernstein). Violation at the strong tolerance if min G < -tol_strong.
**Trajectory oracle (behavioural).** 60 initial states per system (20 interior, 20 on the faces/boundary, 20 on the trapping shell [A] or at large total population [B]),
simulated over the full data horizon (A: 20, B: 60) with blow-up termination. A trajectory is *unphysical* if any of:
blow-up (max|x| > 50); a state variable < -0.01; [A] V = 0.4x+0.5y rises more than 0.05 above 3 when started at V <= 3; [B] S rises > 0.02 above its start, R falls > 0.02 below its start, or S+I+R drifts > 0.05.
The trajectory oracle is used to report the *behavioural consequence rate* of SMT-flagged violations; it is not part of the label.
**Label.** inadmissible = (SMT counterexample at strong tolerance, replay-validated) OR (local-oracle violation at strong tolerance); admissible = SMT verified and local oracle finds nothing.
**Random-falsification baseline.** n random points per spec drawn uniformly from the spec domain (n = 200 and 2000 per spec, i.e. the cheap and the generous budget);
flag if any sample has G < -tol_strong; continuous score = max over specs of (-min sampled G)/tol_strong.

### Addendum C — amendments made after a single-cell preview, before the 30-cell run (registered 2026-10-27)
**What the preview showed (system A and B, 2 % noise, seed 0, well-specified libraries; not part of the pre-registered criteria):**
(i) the SMT verdict and the SMT-free local oracle agreed on all 246 candidates; (ii) random falsification (n = 200) flagged 99.4 % of SMT violators and (n = 2000) 100 %, with no false flags;
(iii) within the AES, AIC and OOD-rollout error had AUROC 0.97-0.99 for the label (iid RMSE 0.92-0.97), i.e. G2b is likely to **fail as originally worded**;
(iv) with a *misspecified* library (system A without the x^2 term) 134/137 candidates, and 80/81 of the AES, violated a law at the strong tolerance although the best held-out error equalled the well-specified best (ratio 1.001).
**Consequences (the original G-criteria are NOT changed; the following are added):**
1. **Arm "A-mis" (misspecified library).** System A, library degrees {3, 4} with the monomial x^2 removed from the x-equation library (zero-padded, same code path), 3 noise levels x 5 seeds = 15 extra cells. Reported separately; not counted in the six well-specified (system x noise) groups of G1.
   *G1-mis:* >= 20 % of AES violating in >= 2 of the 3 noise levels. *V-sel:* fraction of cells in which the validation-selected model (S4) and the AIC-selected model (S5) are inadmissible; verifier-guided selection = the best AES candidate that is verified admissible, or "none exists".
2. **G2b is kept exactly as approved** (max AUROC over iid RMSE, rollout RMSE, OOD-rollout RMSE, AIC <= 0.70). Supplementary G2b': the same maximum over the *practical* metrics only (iid RMSE, in-distribution rollout, AIC), because OOD rollout uses simulated extrapolation data that real studies do not have ("privileged baseline").
3. **Complexity-stratified analysis.** AUROC of each detector within strata of k (k <= k_true+4; k_true+5..k_true+12; > k_true+12) to test whether validation separates violators from admissible models of similar size.
4. **Severity-stratified analysis of random falsification.** Detection rate of random falsification (n = 200, 2000) as a function of worst-case severity (min G / tol), to test whether exact verification adds anything over sampling when violations are narrow.
5. **Tolerance slack.** For SMT-verified models, the rate of trajectory-oracle violations is reported at strong and weak verification level, since a spec tolerance of 0.02 per time unit allows drift/negativity that accumulate over the data horizon.
6. Oracle cost control: the local oracle and trajectory oracle are run on all candidates of the AES and on a 25 % random subsample of non-AES candidates; non-AES labels rely on SMT + replay.

### Addendum D — deviations and clarifications recorded after the 45-cell run (2026-10-28)
1. **Bernstein fallback and 3 s z3 limit.** During verifier development, z3 (nlsat) timed out on 26/812 strong queries in the SIR example cell (tight 3-D bounds on large models). An independent Bernstein-form branch-and-bound on the 3-simplex was added as fallback after a 3 s z3 limit (originally 10 s, no fallback). It was cross-checked against z3 on every query where z3 decided (9,805 checks over the grid, 0 disagreements). Inside repair loops Bernstein is tried first for 3-D specs (speed); the headline verification uses z3 first, and CEG-repaired models were audited with z3 (1,291 audited models, 100 % verified at the strong tolerance).
2. **Replay criterion.** The replay window (0.05 time units; a violation must show c.x(t)-c.x(0) < -1e-7 at some t) was fixed during step 4. Over the grid 44/8,959 strong counterexamples (0.49 %) failed it; a post-hoc diagnosis (`results/replay_failure_diagnosis.csv`) showed all 44 are finite-time blow-ups of the candidate (|x| > 1e6 or integrator failure before t = 0.05) with the violation already visible at t = 1e-5. G3 is judged on the strict figure (99.51 %, passes); labels and the other analyses use the amended rule (blow-up counts as confirmed), under which no SMT-flagged candidate is labelled admissible.
3. **Checkpoint format.** `results/pilot_results.parquet` in the plan was written as CSV (`results/pilot_candidates.csv` and companions) because pyarrow's Parquet module is blocked by the host's application-control policy.
4. **Repair target.** Repairs target the strong tolerance with a 25 % margin. Weak-level verification of repaired models is reported (about 20 %) but was not a repair target.
5. **Grid.** 45 cells = 15 SIR + 15 predator-prey + 15 misspecified-library cells (5 seeds x 3 noise levels each); 5,553 candidates; oracles ran on all accuracy-equivalent candidates and a 25 % subsample of the rest (4,334 candidates).
6. **Not done (flagged for the next phase):** a falsification-guided repair baseline (random-sample counterexamples fed to the same learner), narrow-violation ("needle") stress tests beyond the severity stratification, transcendental (trigonometric) terms, non-sparse or non-library ground truth, real data.
