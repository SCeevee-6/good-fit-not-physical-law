---
name: ceg-constrained-sindy-repair
description: "Repair a sparse-regression (SINDy-style) model that violates a physical law by feeding verifier counterexamples back as linear constraints in a constrained STLSQ refit (cvxpy QP). Use after polynomial-law-verification returns a witness."
---

# ceg-constrained-sindy-repair

## loop
1. Fit: constrained STLSQ. Per equation take the QR of the active, column-scaled library, solve the QP with the constraint rows, threshold small coefficients, repeat until the support is stable (cvxpy, CLARABEL). Round to 1e-6.
2. Verify all specs (strong tolerance).
3. For every violated spec add one row: `r . vec(Xi) >= -0.25 * tol`, r[j,i] = c_i * Theta_j(point(u*)). Linear because every spec is linear in Xi.
4. Refit. Stop at no violation, infeasible QP, or 15 iterations.
Variant "worst": also add the grid argmin of G for each violated spec.

## arms to run at equal budget
CEG witness, CEG + worst point, random rows (same count as CEG used), dense grid (about 300-420 rows), reject only, unconstrained best-by-validation.
Pilot (2,935 violating accuracy-equivalent models): CEG 99.25%, median 3 iterations and 5 rows, median held-out error x0.9993. Random rows 17%. Dense grid about 100% with 60-80x the rows.

## things that went wrong or need care
- Repaired is admissible, not true: median support F1 0.44 (14 terms).
- Repair targets the strong tolerance; only 20% of repaired models also verify at the weak one. Change the margin if weak matters.
- Failures (22 of 2,935) were at 5% noise: iteration cap or solver time-outs.
- Rebuild a candidate from (degree, variant, threshold) rather than storing models; bootstrap indices must come from a seeded generator (`make_variants`), not Python `hash`.
- Missing control: a loop whose counterexamples come from random sampling instead of SMT. Run it before claiming SMT matters for repair.

Code: `src/ceg.py` (learner, arms), `src/run_cell.py` (pipeline per grid cell).
