---
name: polynomial-law-verification
description: "Exact check of physical laws (positivity, trapping region, monotonicity, conservation) on a fitted polynomial ODE: z3 nlsat with Bernstein fallback, exact witnesses, trajectory replay, SMT-free oracle. Use when a discovered equation has to be certified or falsified, not just scored."
---

# polynomial-law-verification

Use when a fitted model is f(x) = Theta(x) Xi (polynomial library) and you need yes/no/witness on a law, not an error number.

## spec form
Every law is `G(u) = c . f(point(u)) >= -tol` for all u in a polytope D.
G is linear in Xi, so one spec object feeds the verifier, the oracle and the repair step.
- positivity on a face: point = x with x_i = 0, c = e_i
- trapping region V = const: point on the shell, c = -grad V
- monotone S' <= 0: c = -e_S on the whole domain
- conservation: c = (1,..,1), two specs (+/-)
Two tolerances: weak (1e-3) and strong (2e-2). Report strong. Weak is a sensitivity check.

## verify
1. Round coefficients to 1e-6 and turn them into `Fraction`. The model you verify is the rounded one.
2. z3 `SolverFor('QF_NRA')`, assert D and `G < -tol`. unsat = proof. sat = witness.
3. Re-evaluate the witness in `Fraction` arithmetic. Only an exact violation counts.
4. z3 times out on tight 3-D queries with big models (about 2% of queries). Fall back to Bernstein branch-and-bound on the simplex: homogenise G to degree n in barycentric coordinates, b_alpha = f_alpha * alpha!/n!, bisect the longest edge with precomputed linear maps, witnesses from vertices/centroid. Sound with float slack 1e-9 for "verified"; witnesses are exact-checked.
5. Inside repair loops try Bernstein first for 3-D (ms instead of seconds); audit the final models with z3.
6. Cross-check: 1-D specs against sympy `real_roots` of G'; 3-D against Bernstein. Pilot: 0 disagreements in 14,212 checks.

## replay and oracle
- Replay: integrate the candidate from the witness; c.x has derivative G, so a violation shows c.x(t)-c.x(0) < -1e-7 at some t in (0, 0.05]. Finite-time blow-up counts as confirmed (the pilot's only 44 "failures" were blow-ups; check blow-up before calling a replay a false alarm).
- Oracle that does not use the verifier: dense grid + differential evolution on G, and trajectory tests from boundary/out-of-range starts (negative state, escape, S rising, mass drift). Pilot: 0 disagreements with SMT on 4,334 candidates.

## pitfalls
- A tolerance of 0.02 per time unit still lets 4-18% of "verified" models misbehave over a long horizon. Weak-verified: <= 1.6%.
- Do not label with the verifier alone if you also use it as the detector under test. Keep the oracle and replay in the label.
- nlsat does not scale past a handful of variables. For n >= 6 use SOS (cvxpy + Clarabel/SCS) or quadratic-form linear algebra.

Code: `src/verify.py`, `src/bern.py`, `src/oracle.py`, `src/systems.py` (spec definitions).
