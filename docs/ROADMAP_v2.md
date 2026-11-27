# Roadmap v2

Superseded for planning by `docs/ROADMAP_v3.md` (2026-11-26). Kept as history.

Date: 2026-11-08. This file supersedes `docs/ROADMAP.md` for planning. The original text and addenda 1 and 2 stay as history. Each unit is governed by its frozen protocol. This file gives the order of work, the gates and the decision rules.

## Question

When a model fits held-out data well, does it break a stated physical law? And what is the right fix: a counterexample-guided repair, a structure built into the model, or nothing?

Claim to earn: a decision rule, backed by a test that can fail. Not "verification beats validation". The pilot already showed that claim fails on synthetic data.

## Current state (handoff)

- Pilot: done. 45 cells, 5,553 candidate ODEs, pre-registered verdict REVISE. Report: `docs/pilot_report.md`.
- Repo: https://github.com/SCeevee-6/good-fit-not-physical-law, commit 3f91a31f. Addenda 1 and 2 and this file are not pushed yet.
- Gate A, protocol 0A (`docs/protocol_0A.md`): frozen 2026-11-08, with one dated clarification (B1S moved to the non-structural family). Sample frozen in `results/0A_sample.csv`: 60 candidates, seed 20261108. Code: `src/gate_a.py`. Reproduction of pilot coefficients matched exactly on three candidates. The full run was cancelled by a kernel reset. No 0A results exist yet.
- Side track S1 (`docs/protocol_S1.md`): frozen, not run.
- Open decisions: competition category (after Gate A), real-data source for hare-lynx, affiliated fair and its deadline.
- Environment: the `ceg-ode` kernel can no longer import pandas, and some SciPy DLLs were blocked after the reset. Use `csv` and numpy. Run long jobs in the background with per-job outputs. ## Gates

| gate | question | pass rule | source |
|---|---|---|---|
| A | Do witness counterexamples add value beyond a structural fix, and beyond random counterexamples at equal budget? | Rule in protocol 0A | `docs/protocol_0A.md` |
| R | Does a real dataset show a validation-selected model breaking a stated law under extrapolation? | If no, drop the real-data claim | this file |
| S1 | Is each learned linear closure energy-stable at truncation K? | Exact quadratic-form test agrees with eigenvalues | `docs/protocol_S1.md` |

Gate A outcomes:
1. Both criteria pass. The witness loop is the main method. Category follows the headline demo.
2. Structure is enough on these specs. The claim becomes the rule for when a structural fix suffices. Recorded as a dated addendum.
3. Counterexamples help, but the feedback loop explains it, not the witness type. The random-counterexample loop becomes a main arm. Report that.

## Work order

Each block opens from files, writes its protocol before its first run, runs long jobs in the background with resumable outputs, and ends with `STATUS.md` and a push via `tools/gh_push.py`.

| block | scope | ends at |
|---|---|---|
| 1 | Gate A: run protocol 0A, analyse, decide | Gate A decision recorded |
| 2 | Phase 0B: Holling II misspecification (true law outside the library), needle test (narrow violations), 20 seeds with confidence intervals | Phase 0 frozen |
| 3 | Side track S1: DNS validation, advection energy-neutrality check, linear closures, exact stability test | Gate S1 results |
| 4 | Phase 1: hare-lynx real data, by-construction baselines on real data | Gate R |
| 5 | Phase 2: epidemic time series, results freeze | Freeze date |
| 6 | Phase 3: robustness checks, board demo, write-up | Board-ready tag |

Block 3 (S1) can start earlier if block 2 is waiting on compute.

## Calendar

Assumes an affiliated fair in late February 2027. Set the freeze at that fair's deadline minus three weeks.

- Nov 9 to Dec 6: block 1
- Dec 7 to Jan 10: block 2
- Dec to Jan: block 3, overlapping block 2 if compute allows
- Jan 11 to 31: block 4
- Feb 1 to 21: block 5, freeze
- Feb 22 to Apr 4: block 6
- May 2027: Regeneron ISEF, Los Angeles (exact dates per the Society)

## Method scope

- Verifier: z3 for up to three state variables, with Bernstein fallback on the simplex. Every counterexample has an exact witness and a trajectory replay.
- SOS only for models with four or more state variables. Optional, not scheduled.
- Fluid-scale polynomial models (Lorenz-63, shear flow): optional stretch, not scheduled.
- Baselines: constrained STLSQ; by-construction fixes (log-states for positivity, fraction parametrisation for conservation, soft penalties); random-counterexample loop.
- Navier-Stokes: side track S1 only. Claims are limited to energy stability of a truncated model with a linear closure. No statements about blow-up, the Millennium problem or the 2026 AI claims.

## Risks and responses

| risk | response |
|---|---|
| Structure is enough (Gate A outcome 2) | Shift the claim. That is a result. |
| Real-data source blocked by the sandbox | Request access, or the user supplies the file. Do not substitute silently. |
| Kernel reset or DLL block mid-run | Per-job outputs, resume from files, `csv` in place of pandas. |
| Overclaiming from two real systems | Report as two cases. No general claim. |
| Fair rules, including AI-use rules | The fair's rulebook governs the submission. Check it before the fair deadline. |
| Context lost between work blocks | Each block starts from the files listed above. This file is the index. |

## Discipline

- Freeze the protocol before the first run. Changes after a run are dated addenda. Headline criteria do not change.
- Every reported number comes from a script in `src/`. Tables are regenerated by script.
- Failed criteria go in the same table as passed ones.
- Public text is plain.

## Change log

- Original plan, before the pilot run: `docs/ROADMAP.md`, original text.
- Addendum 1, 2026-11-08: `docs/ROADMAP.md`, dated section above the original text.
- Addendum 2, 2026-11-08: `docs/ROADMAP_ADDENDUM_2.md`.
- v2, 2026-11-08: this file, a full rewrite of the plan. History kept.


## Addendum, 2026-11-08: Gate A decided

Protocol 0A gave outcome 3 (`docs/gate_A_result.md`). G-A1 passed (f(C) - f(S) = 0.533, median held-out ratio 1.059). G-A2 failed (f(C) - f(Rnd) = 0.100). Changes that follow from the frozen rule:

- The random-counterexample loop is a main arm next to the witness loop from now on, not only a control.
- The claim is that a counterexample loop is needed on top of structure. Witnesses only win where violations are narrow. In 0A all six random-loop misses were B2S on SIR.
- Block 2 (0B) keeps its scope. The needle test is now the main test of witness vs random, and should vary the sample count of the random loop.
- SIR repairs cost more accuracy (median ratio 1.23 for C) than the pooled gate figure shows. 0B reports accuracy cost per system.
- Category decision can now be made. It still follows the headline demo.


## Addendum 3, 2026-11-08: B structure defect, review points, plan for 0A2

**Defect.** In protocol 0A the B structure kept only R-divisible monomials in the R-equation. The true R' = 0.2 I is not R-divisible, so arm S on B was misspecified by construction. A and A-mis are correct. Details in `docs/gate_A_result.md`. The 0A verdict (outcome 3) stays on record, flagged as confounded on B. It is not used for planning until 0A2 is in.

**Wording.** Rule text in 0A says "pooled". The pre-run scoring note fixed the reading per candidate. From 0A2 on, protocols say "fraction of candidates" and report the per law x candidate pooled figure as secondary. On timing: the log of the cancelled first run shows 29 of 60 candidates finished (A and A-mis only) before the reset. Their outputs were held in kernel memory and lost, and none were written to file or recorded.

**Protocol 0A2 (next block, frozen before running).** Same sample, arms, budget, seed and decision rule as 0A. Only B changes, and A and A-mis are not rerun. Two B structures, both checked to contain the true model before running:
- B-pos: I-equation divisible by I. R-equation free, but coefficients of R-free monomials are constrained >= 0, which makes R'(S,I,0) >= 0 on the orthant. S-equation = -(I' + R').
- B-fac (the review's variant): same as B-pos, plus the S-equation must be divisible by S. That is a linear equality: the R-free and S-free part of I' + R' is zero. This enforces S positivity and conservation exactly.
Gate A is re-decided from A, A-mis (0A results reused) and B-pos. B-fac is reported next to it. If B-fac and B-pos give different outcomes, report both. Accuracy cost is reported per system and per law.

**Review points taken.**
- 0B: random-loop sampler budget 2k, 20k and 200k points, plotted against violation volume, with witnesses as the reference.
- Accuracy cost reported per law and per system, not only pooled.
- Working headline, conditional on 0A2: which laws a parametrisation enforces exactly, which need an exact check, what the check costs, and when finite sampling misses a violation.
- S1 (Navier-Stokes) stays small, as one more case of the same question.
- Gate R moves earlier once the hare-lynx source is settled.

**Review points not taken as written.**
- "Accuracy cost is real on SIR": the 0A SIR numbers come from the defective structure. Measure it again in 0A2.
- The demo idea (structured SIR stays positive but S rises): only valid with a correct structure, so it waits for 0A2.
- Category: still decided after 0A2.


## Addendum 4, 2026-11-08: Gate A re-decided (0A2)

Protocol 0A2 ran B-pos and B-fac on the same 20 B candidates. Outcome 3 with B-pos (f(C)-f(S) = 0.317, median held-out ratio C 1.000, f(C)-f(Rnd) = 0.017). B-fac also gives outcome 3, at exactly 0.200 on G-A1, all from A and A-mis. Details in `docs/gate_A2_result.md`.

What changes in the plan:
- The random-counterexample loop is a main arm. This was conditional on 0A2 and now holds.
- No accuracy cost on SIR with a correct structure (median ratio 1.000, max 1.003). The 0A SIR cost is withdrawn.
- Only 1 of 60 candidates has a violation that 2000 random points miss. 0B therefore builds needles on purpose and runs sampler budgets 2k / 20k / 200k against violation volume, witness loop as reference. Natural needles are too rare here to carry the comparison.
- B-fac alone satisfied every SIR law on all 20 candidates. The law taxonomy gets a per-law column for "enforced by which structure" and "satisfied empirically without enforcement".
- Category decision is unblocked; it still follows the headline demo.

## Addendum 5, 2026-11-23: protocol 0B run

Protocol 0B was frozen at commit 9f1cf045 before any run, and ran on Kaggle. Details in `docs/gate_B_result.md`.

- G-B1 passes: 97 % of equivalent Holling II models break a law, mostly positivity.
- G-B2 passes: C repairs 160/160 at a median held-out cost of 3.6 %, and out-of-range rollouts improve.
- G-B3 fails: no natural violation is anywhere near the volume a 20k-point sampler misses (smallest 0.064 vs v50 3.5e-5). The random detector matches 1 - (1 - v)^n on built needles.
- Gate A outcome 3 replicates over 20 seeds: f(C)-f(Rnd) 0.033 [0.017, 0.050], f(C)-f(S) 0.354 [0.288, 0.429].

What changes in the plan:
- Do not claim witness search finds violations that sampling misses on these testbeds. The claim is a certificate (verified admissible, exact witness) plus a measured miss curve that says when sampling stops being enough.
- New observation to follow up: on H, structure alone makes the trapping law fail more often (f(H2) 0.88 pool, 0.40 S). Enforcing one law can push a misspecified fit into breaking another. This is a candidate headline for the law taxonomy.
- Next: hare-lynx real data (network access approved by the user), then the category decision.

## Addendum 6, 2026-11-26: protocol 0R run (Gate R)

Protocol 0R was frozen at 0ad1896c. A first run stalled, and addendum 1 of the protocol added time limits before any real-series result existed. The run that counts is commit fcfb358d. Details in `docs/gate_R_result.md`.

- G-R1 passes: the validation-selected model breaks the hare positivity law and the trapping law.
- G-R2 passes: 45 of 47 accuracy-equivalent candidates break a law.
- G-R3 passes: the witness loop repairs 47 of 47 at a median validation ratio of 0.99.
- Verdict: the real-data claim stays.

What changes in the plan:
- The real-data claim is weaker than the synthetic one. The fits are poor (validation error 0.89 training sd) and the witnesses sit where there is no data. State it as "unconstrained fits break stated laws outside the data", not as a finding about hares and lynx.
- Gate R is passed on one series. A second real series (epidemic, block 5) is what would make it a claim.
- Next: the law taxonomy per law and structure, then the epidemic series.

## Addendum 7, 2026-11-26: law taxonomy

`docs/law_taxonomy.md` and `results/law_taxonomy.csv` tabulate the pilot, 0B and 0R results per system and law. Nothing was rerun.

- Positivity laws and the SIR conservation laws are enforced exactly by structure. Without it, the fit breaks them often (hare positivity in 81 % of equivalent hare-lynx models, H1x in 97 % of H models).
- The non-structural laws (A2, B1S, B2S, B2R, H2, R2) need a check. The witness loop verifies all of them in every sampled model. The random loop is behind on B2S (72/80), B1S (79/80) and R2 (45/47).
- B-fac enforces all seven SIR laws by structure alone on the 80 sampled models. B-pos does not.
- On H, structure alone lowers the share of models that satisfy the trapping law from 0.87 to 0.35 (Wilson 0.25 to 0.46). One system, so it is a candidate finding to test on the epidemic series, not a result.

## Addendum 8, 2026-11-26: protocol 0E run (second real series)

Protocol 0E was frozen at commit 942d55e5 and ran on Kaggle (JHU CSSE panel of 10 countries). Details in `docs/gate_E_result.md`.

- G-E1 passes: the validation-selected model breaks a law in 10 of 10 countries.
- G-E2 passes: all 35 equivalent candidates break a law. The equivalent sets are small (1 to 10 per country), so the interval has no width.
- G-E3 fails: C repairs 34 of 35 but the median validation ratio is 1.46, so the cheap-fix claim is not kept on this series. Extrapolation rollouts improve (median 0.64).
- The B-fac result from the SIR testbed does not carry over: B-fac alone leaves B2R unverified in 24 of 35 models. It held on the testbed because the truth lies inside the structure.
- Structure alone lowers satisfaction of a law again (B2R under B-fac, B1S under B-pos), as on H.

What changes in the plan:
- Two real series now show the same thing for unconstrained fits: they break the stated laws, mostly where there is no data. Both fits are poor, so say it that way.
- The fix has a cost on real data that the synthetic runs did not show. Report the validation ratio per country, not only a pooled median.
- The claim to carry forward is the one in `docs/law_taxonomy.md`: structure enforces what it contains and nothing else, and a structure that does not contain the truth can make other laws worse.
- Next: robustness checks and the board demo (block 6). Competition category still follows the headline demo.

## Addendum 9, 2026-11-27: protocol 0C run (robustness)

Protocol 0C was frozen at commit 2109feba and ran on Kaggle. Details in `docs/gate_C_result.md`.

- R-C1: the counterexample loop is not needed over a soft penalty on a dense grid, on these testbeds. The soft penalty verifies all 320 sampled candidates at all three lam, with the same accuracy as C.
- R-C2: Gate A outcome 3 holds at 0.5 % and 5 % noise.
- R-C3: 0 of 2871 accepted models fail an independent uniform-sampling check.
- At 0.5 % noise the repair is not cheap where the library is wrong: median validation ratio 1.49 (A-mis) and 1.42 (H).

What changes in the plan:
- Do not claim the loop beats soft penalties on accuracy or on share of laws satisfied. The remaining argument is cost in constraints and the certificate, and growth with state dimension. The dimension argument is untested.
- Report the accuracy cost of repair as a function of noise. "Cheap" holds at 2 % and 5 %, not at 0.5 % for a wrong library.
- Next: the choice of headline for the demo. The candidates are the law taxonomy with its structure-containing-the-truth caveat (two real series, three synthetic systems) and the soft-penalty adequacy result. Category follows the headline.

## Addendum 10, 2026-11-27: protocol 0D run (four states)

Protocol 0D was frozen at commit 3cbcb2e9, with a test-only fix before the run (3641866a), and ran on Kaggle. Details in `docs/gate_D_result.md`.

- R-D1: no coverage advantage for the loop at four states. The fine grid (11 points per axis, 5148 rows) verifies 80 of 80 candidates, the loop 69, the pilot-budget grid 64.
- The loop hit its 15-iteration budget in 11 of 80 runs. A larger budget was not tried.
- The random loop is clearly behind: 51 of 80 verified, although its 2000-point sampler stopped clean in 77. C - Rnd is 0.225 [0.125, 0.338].
- Verifier agreement (0 of 703) and the independent sampling check (0 of 400) hold. A d-state Bernstein verifier now exists (`src/bern_nd.py`).

What changes in the plan:
- The soft-penalty result stands up to four states. Do not claim the loop beats a grid on coverage or accuracy. The claims left for the loop are the constraint count (15 against 364 or 5148 rows) and the exact witness.
- The witness loop's advantage over random counterexamples grows with the number of states. That is the dimension result that supports the loop.
- Untested: a loop with a larger budget, a grid built on the simplex, five or more states.
- Headline recommendation stays the taxonomy. The dimension result goes in as one of its controls.
