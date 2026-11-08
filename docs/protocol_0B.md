# Protocol 0B: misspecification, needles, 20 seeds (DRAFT, not frozen)

Status: draft written 2026-11-08 after Gate A was re-decided (`docs/gate_A2_result.md`, outcome 3). Nothing in 0B has run. Values marked *check* are computed by hand and must be confirmed by script before this file is frozen. The open questions at the end are settled before freezing.

Questions.
1. When the true law is outside the library (Holling type II), do validation-equivalent models break the stated laws, and do structure plus a counterexample loop fix that at low cost?
2. When does a random-sampling loop miss a violation that the verifier finds? Answered as a function of violation volume and sampler budget, on needles built on purpose, because natural needles were 1 in 60 in 0A2.
3. Do the Gate A numbers hold over 20 seeds, with confidence intervals?

## Units

### 0B-1. Holling II testbed (system H)

True model, Rosenzweig-MacArthur:
x' = x (1 - x/K) - a x y / (1 + h x), y' = -d y + b x y / (1 + h x).
Proposed values a = 1, h = 1, b = 0.8, d = 0.4, K = 4. Coexistence point x* = d / (b - d h) = 1. Hopf at K = (b + d h) / (h (b - d h)) = 3, so K = 4 gives a limit cycle (*check* numerically).

The polynomial library (degree 2 to 4 in x, y) cannot represent the rational term, so every candidate is misspecified. The positivity structure still holds for the truth (x' is x times something, y' is y times something), so arm S is a correct structure here.

Laws (strong tolerance as in the pilot):
- H1x: x'(0, y) >= 0, H1y: y'(x, 0) >= 0. Structural under arm S.
- H2: trapping on V = (b/a) x + y = c. On V = c, V' = (b/a) x (1 - x/K) - d y. With the values above the maximum over the segment is 1.568 - 0.4 c, so any c > 3.92 traps. Proposed c = 5 (V' <= -0.432 on the segment, *check*). Non-structural.
- Truth check before freezing: multiply each law by (1 + h x) > 0 to get a polynomial sign condition and verify the true model with z3. Script `tools/check_holling.py`.

Data and pool: the pilot generator and pool code with system H added (`src/systems.py`), noise indices 1 and 2 (2 % and 5 %), seeds 0-19, the pilot's candidate pool and accuracy-equivalence rule. 40 cells.

Arms: unstructured STLSQ (the pool), S, C, Rnd (2000 points), all as in 0A.

Measures per cell and pooled:
- Fraction of accuracy-equivalent candidates that violate at least one law (pool).
- f(S), f(C), f(Rnd) on a frozen sample of 4 accuracy-equivalent candidates per cell (160 candidates), drawn as in 0A with a new seed.
- Held-out ratio and out-of-range rollout error ratio (OOD initial conditions), per law and per arm. In 0A2 the in-range cost was nil, so the OOD figure carries the accuracy question here.
- Violation volume of every violating candidate (Monte Carlo, 10^6 points per spec, plus z3 depth min G / tol). This feeds the needle unit.

### 0B-2. Needle test (constructed)

Needles are built from verified models, so depth and volume are known.
- Base models: the true A and B models and the arm-C models from 0A2 that verify on every spec.
- Perturbation: add t * m(u) to one equation, with m a monomial chosen so the violation sits near a corner of the domain (higher degree gives a smaller set), and t set so that min G = -kappa * tol_strong. Depth kappa in {1.05, 1.5, 3}. Degree and corner give the volume.
- Target volume grid: about 1e-2, 1e-3, 1e-4, 1e-5, 1e-6 of the spec domain, measured by Monte Carlo with 10^7 points.
- For each needle: the random detector at budgets n = 2k, 20k, 200k, 50 sampler seeds each, records whether any point has G < -tol. The witness verifier is run once per needle (detection or timeout).
- Reference curve: 1 - (1 - v)^n, with v the measured volume. Report the observed miss rate against it.
- The full Rnd loop is also run on a subset (kappa = 1.05) to see whether a miss ends the loop with an inadmissible model, as in B_s0_69.

### 0B-3. Gate A over 20 seeds

Arms S, C, Rnd on A, A-mis, B-pos, B-fac and H, seeds 0-19, noise index 1, 4 candidates per cell, new sample seed. Report f differences with cluster-bootstrap 95 % intervals (resampling seeds, 2000 replicates) and Wilson intervals for each f. This is the replication of Gate A, not a new gate.

## Decision rules (proposed, to be frozen)

- **G-B1, misspecification matters:** on H, at least 20 % of accuracy-equivalent pool candidates violate a law (pooled over the 40 cells), and the lower 95 % bound is above 10 %.
- **G-B2, structure plus loop repairs at low cost:** on H, f(C) >= 0.90 and median OOD rollout ratio of C <= 1.10.
- **G-B3, needles matter in practice:** let v50(n) be the volume where the observed miss rate of the n-point sampler reaches 50 %. Needles matter if, among violating candidates in the 0B-1 pool, at least 5 % have volume below v50(20k). If not, the witness loop is reported as a certificate on top of a random loop that already finds what occurs, not as a better search.
- **Replication:** Gate A outcome 3 is kept if the 95 % interval for f(C) - f(Rnd) contains values below 0.20 and the interval for f(C) - f(S) lies above 0 on the pooled A + A-mis + B-pos data.

## Compute (estimate)

0B-1: 40 cells at about 10 min each on 6 workers, about 70 min. 0B-2: depends on the volume grid; 200k-point checks are cheap, z3 per needle as in 0A. 0B-3: 100 cells x 4 candidates, about 15 min. Long runs go in the background with one JSON per job.

## Open questions before freezing

1. Holling parameters and trapping constant c: confirm by script, and check that z3 decides the H2 query within the 3 s limit at degree 4; otherwise use the Bernstein fallback (2-D, fine).
2. Needle construction for the 3-D SIR specs: corner monomials are easy on the simplex; check that the measured volume can reach 1e-6 with a degree-4 library.
3. Whether 0B-3 replaces 0A2 as the citable Gate A number, or is reported next to it. Proposal: 0A2 stays the decision, 0B-3 is the interval.
4. Accuracy-equivalence rule for H: reuse the pilot rule unchanged unless H pools come out with very few equivalent candidates (fewer than 4 in a cell). If so, record an addendum before sampling.
