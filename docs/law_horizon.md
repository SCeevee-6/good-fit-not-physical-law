# Law horizon: proofs for ROADMAP_v3 item 7 (unit T)

Date: 2026-11-27. The statements, proofs and the pass rules in section 9 were committed before `tools/check_law_horizon.py` was run for the first time. Results are added at the end, below a line, and nothing above that line is edited after the run except through a dated note.

This writes out item 7 of `docs/ROADMAP_v3.md` and settles its open points (i)-(iii). It also answers the question Addendum 2 left: is there a better lower bound than T_d along one functional? There is, but only one that looks at the data (section 4). Without the data, the T_d(l) test polynomial cannot be beaten (Lemma 2).

## 0. Setting and notation

- State u in R^s. The library is a vector phi(u) of p monomials of degree <= d. "Full library" means every monomial of degree <= d, so p = C(s + d, d). P(phi) is the span: polynomials q(u) = c' phi(u).
- Data: rows phi(u_i), i = 1..n, targets y = Theta xi + e with e ~ N(0, Sigma), Sigma known. G = Theta' Sigma^-1 Theta and xi_hat is the GLS estimate. One equation at a time.
- Confidence set: E = {xi : (xi - xi_hat)' G (xi - xi_hat) <= q}, q = chi2_{p, 1 - alpha}. This is the 0F set in whitened units (eps = q).
- Law at u: phi(u)' xi >= -tol. K(u) = phi(u)' G^-1 phi(u), h(u) = sqrt(q K(u)). Estimated margin m(u) = phi(u)' xi_hat + tol, true margin mu(u) = phi(u)' xi + tol.
- Data norm of a polynomial: ||q||_S^2 = y_q' Sigma^-1 y_q with y_q = (q(u_1), ..., q(u_n)). For iid noise of variance sigma^2, ||q||_S^2 = n ||q||_n^2 / sigma^2 with ||q||_n^2 = (1/n) sum q(u_i)^2, and I write K = sigma^2 Kbar / n, Kbar(u) = phi(u)' (Theta' Theta / n)^-1 phi(u).
- Omega: a centrally symmetric convex body that contains every u_i, centred at 0 (P_d is translation invariant, so this costs nothing for a full library). Its gauge is g(u) = inf{t > 0 : u / t in Omega}. In 0F, Omega is the Mahalanobis ellipsoid of the training states.
- T_d is the Chebyshev polynomial. For t >= 1, T_d(t) = cosh(d arccosh t), and with rho(t) = t + sqrt(t^2 - 1): rho^d / 2 <= T_d(t) <= rho^d. Its inverse on [1, inf) is T_d^-1(y) = cosh(arccosh(y) / d).

## 1. K is an extremal problem

**Lemma 1.** K(u) = max over nonzero q in P(phi) of q(u)^2 / ||q||_S^2, attained at c = G^-1 phi(u).

Proof. For q = c' phi, ||q||_S^2 = c' G c. By Cauchy-Schwarz in the G inner product, (c' phi)^2 = (c' G G^-1 phi)^2 <= (c' G c)(phi' G^-1 phi), with equality at c proportional to G^-1 phi. QED

Three things follow. K depends on the function space and the data norm, not on the basis. A larger library gives a larger K at every point. Any single test polynomial gives a lower bound on K.

## 2. How fast a polynomial can grow outside a symmetric convex body

**Lemma 2.** Let deg q <= d and |q| <= 1 on Omega. If g(u) >= 1 then |q(u)| <= T_d(g(u)). The bound is attained by T_d(l(x)), where l is linear, l(u / g(u)) = 1 and |l| <= 1 on Omega (a supporting functional at the boundary point u / g). So

    M_d(u) := sup{ |q(u)| : deg q <= d, |q| <= 1 on Omega } = T_d(g(u)).

Proof. The segment {lambda u / g : -1 <= lambda <= 1} lies in Omega (u / g is on the boundary, Omega is symmetric and convex). r(lambda) = q(lambda u / g) is a univariate polynomial of degree <= d with |r| <= 1 on [-1, 1]. Chebyshev's inequality gives |r(lambda)| <= T_d(|lambda|) for |lambda| >= 1. Take lambda = g. For attainment, a supporting hyperplane at u / g gives l with l <= 1 on Omega and l(u / g) = 1; symmetry gives |l| <= 1 on Omega. Then |T_d(l)| <= 1 on Omega and T_d(l(u)) = T_d(g). QED

I believe this is a classical fact about convex bodies, but I have not pinned down a reference yet. The consequence is what matters here. The lower side of item 7 used T_d(l) as the test polynomial and bounded its data norm by its sup over Omega. Lemma 2 says no other test polynomial does better *as long as only the sup over Omega is used*. A tighter lower bound has to use where in Omega the data actually are.

## 3. The two-sided band (item 7)

**Theorem 3.** Full library of degree d, iid noise, data in Omega, N2 := sup over Omega of Kbar. For g(u) >= 1:

    T_d(g)^2 <= Kbar(u) <= N2 T_d(g)^2,   and N2 >= p.

So sigma sqrt(q / n) T_d(g) <= h(u) <= sigma sqrt(q / n) sqrt(N2) T_d(g). With a general Sigma the lower side reads K(u) >= T_d(g)^2 / (n lambda_max(Sigma^-1)), which is the lambda_min form used in 0F-3.

Proof. Upper side: for any q and any v in Omega, Lemma 1 at v gives q(v)^2 <= Kbar(v) ||q||_n^2 <= N2 ||q||_n^2. So q / (sqrt(N2) ||q||_n) is bounded by 1 on Omega, and Lemma 2 gives q(u)^2 <= N2 ||q||_n^2 T_d(g)^2. Take the max over q and use Lemma 1 at u. Lower side: q = T_d(l) is in the full library, and ||q||_n <= 1 because every u_i is in Omega. Lemma 1 gives Kbar(u) >= T_d(g)^2. Finally, (1/n) sum_i Kbar(u_i) = tr(M^-1 M) = p with M = Theta' Theta / n, and each u_i is in Omega, so N2 >= p. QED

Notes.
- Only the lower side needs the full library. The upper side needs only degree <= d. This is what section 8 uses.
- The bracket in Kbar is a factor N2 >= p wide. It cannot close in general.
- Laws that weight several equations, a(u) = sum_j w_j(u) e_j (x) phi(u), with the equations fitted on the same rows and independent noise, have K_a(u) = sum_j w_j(u)^2 sigma_j^2 Kbar(u) / n. For trapping, w = grad V, which is the |grad V|^2 factor in the roadmap. This is exact only for block-diagonal G with equal blocks. The stacked errors-in-variables set in 0F is not of that form.

## 4. A better lower bound: the data seen along one direction

**Theorem 4.** Full library. For any linear functional l (not necessarily a supporting one), let t_i = l(u_i), t = l(u), and

    k_d(t) = max over univariate r, deg r <= d, of r(t)^2 / ((1/n) sum r(t_i)^2) = v' M_l^-1 v,

with v_k = T_k(t) and (M_l)_jk = (1/n) sum_i T_j(t_i) T_k(t_i), j, k = 0..d. Then

    Kbar(u) >= sup over l of k_d(l(u)) >= T_d(g)^2 / ((1/n) sum_i T_d(l*(u_i))^2) >= T_d(g)^2,

where l* is a supporting functional at u / g. In one dimension the first inequality is an equality.

Proof. The polynomials r(l(x)) with deg r <= d are a subspace of the full library. Lemma 1 on that subspace gives the first inequality, and the Rayleigh quotient gives v' M^-1 v. The second takes r = T_d and l = l*. The third uses |T_d(l*(u_i))| <= 1. When s = 1, r(l(x)) covers the whole library. QED

The middle term is the "exact test-polynomial" bound already used in 0F-3. The new part is the first term. It is a (d + 1) x (d + 1) solve per direction.

How much it gains: if the projected data have a limit distribution nu on [-1, 1], k_d tends to the inverse Christoffel function of nu. Two closed cases:
- Data at Chebyshev nodes (arcsine nu): the discrete orthogonality of T_k gives M = diag(1, 1/2, ..., 1/2) exactly, so k_d(t) = 1 + 2 sum_{k=1..d} T_k(t)^2. At t = 1 that is 2d + 1, against 1 from the T_d bound.
- Uniform nu: k_d(t) = sum (2k + 1) P_k(t)^2 (Legendre), which is (d + 1)^2 at t = 1.

So near the data the T_d(l) bound undershoots by a factor between 2d + 1 and (d + 1)^2 for these two distributions. Far out, the leading term dominates and the factor tends to a constant. This gap is consistent with 0F-3, where most rays were still censored at g = 8.

## 5. The noise model written in: exact zone probabilities (item (i), part 1)

**Theorem 5.** Assume exact regressors (the u_i carry no error), Gaussian noise with known Sigma, and the law at u as above. Let z(u) = mu(u) / sqrt(K(u)). Then

    P(certified at u) = Phi(z - sqrt q),   P(refuted at u) = Phi(-z - sqrt q),

and P(decided at u) lies in [Phi(|z| - sqrt q), 2 Phi(|z| - sqrt q)].

Proof. m(u) - mu(u) = phi(u)' (xi_hat - xi) ~ N(0, K(u)). The point is certified iff m >= sqrt(q K), that is (m - mu) / sqrt K >= sqrt q - z, with probability Phi(z - sqrt q). It is refuted iff m < -sqrt(q K), with probability Phi(-z - sqrt q). Their sum is at least the larger term Phi(|z| - sqrt q) and at most twice it. QED

So: where |z| >= sqrt q the law is decided with probability at least 1/2. Where |z| <= sqrt q - Phi^-1(1 - beta/2), it is decided with probability at most beta. I call {|z| >= sqrt q} the median decided region and its edge the median horizon.

Where it stops holding. With noise in the regressors, xi_hat is biased and m - mu has a mean b(u) that is not zero. This is the E4 failure in 0F (median calibration ratio 3.9 at degree 2). Then z must be replaced by (mu + b) / sqrt K and the formula is not usable without a bias estimate. With sigma estimated rather than known, Phi becomes a noncentral t; at the n used here the difference is small and I do not pursue it.

## 6. A closed-form horizon radius (item (i), part 2)

With iid noise, z(u) = sqrt(n) mu(u) / (sigma sqrt(Kbar(u))). Put Theorem 3 into it.

**Corollary 6.** Suppose |mu| = mu0 > 0 at the points considered, which have g >= 1 (a margin that does not change outside the data). Let R = mu0 sqrt(n) / (sigma sqrt q). Then

- T_d(g(u)) <= R / sqrt(N2) implies |z| >= sqrt q, so the law is decided at u with probability >= 1/2;
- T_d(g(u)) > R implies |z| < sqrt q;
- T_d(g(u)) >= R sqrt q / (sqrt q - Phi^-1(1 - beta/2)) implies the law is decided at u with probability <= beta.

Along a ray, the median horizon g_med lies in [T_d^-1(R / sqrt(N2)), T_d^-1(R)], where the left end is read as 1 if R / sqrt(N2) < 1.

Proof. Kbar >= T_d^2 gives |z| <= sqrt(q) R / T_d(g). Kbar <= N2 T_d^2 gives |z| >= sqrt(q) R / (sqrt(N2) T_d(g)). Each bullet follows by substitution and Theorem 5. No monotonicity along the ray is needed. QED

Scaling. For large y, T_d^-1(y) is about (2y)^(1/d) / 2. So the horizon radius grows like (n / sigma^2)^(1 / (2d)), and the bracket's ratio tends to N2^(1 / (2d)). At d = 3, four times the data moves the horizon out by a factor 4^(1/6) = 1.26. Doubling the noise pulls it in by 2^(1/3) = 1.26. More data barely moves a high-degree horizon.

Exact example. For data at n Chebyshev nodes on [-1, 1], Kbar(t) = 1 + 2 sum_{k=1..d} T_k(t)^2 exactly (section 4) and N2 = 2d + 1, reached at t = +-1. So the median horizon solves 1 + 2 sum T_k(g)^2 = R^2, and the bracket is [T_d^-1(R / sqrt(2d + 1)), T_d^-1(R)].

## 7. Library degree above the law's degree (item (ii))

**Proposition 7.** Iid noise, full library of degree d.

(a) Suppose that on the law's domain outside Omega, |mu(u)| <= C g(u)^{d_m} with d_m < d. Let R_m = 2 C sqrt(n) / (sigma sqrt q). Then |z(u)| < sqrt q at every such u with g(u) > cosh(log(R_m) / (d - d_m)). So the median horizon is finite at every n and sigma, and its radius is at most that number (when R_m > 1).

(b) Along a ray {g e}, suppose mu(g e) has exact degree d in g with leading coefficient mu_d(e). Let kappa(e) = phi_top(e)' [M^-1]_top phi_top(e), where phi_top holds the degree-d monomials and [M^-1]_top is the matching block of the inverse moment matrix. Then z(g e) -> sqrt(n) mu_d(e) / (sigma sqrt(kappa(e))) as g -> inf. If the absolute value of that limit is above sqrt q, the law is decided with probability >= Phi(|limit| - sqrt q) > 1/2 far out along e. If it is below sqrt q, the median horizon along e is finite.

Proof. (a) Kbar >= T_d(g)^2 >= rho^{2d} / 4 and g <= rho, so |z| <= 2 C sqrt(n) rho^{d_m} / (sigma rho^d) = sqrt(q) R_m / rho^{d - d_m}. This is below sqrt q once rho(g) > R_m^{1 / (d - d_m)}. Since g = (rho + 1/rho) / 2 = cosh(log rho), that is g > cosh(log(R_m) / (d - d_m)). (b) mu(g e) / g^d -> mu_d(e). phi(g e) = D_g phi_e, with D_g scaling each degree-k monomial by g^k, so Kbar(g e) / g^{2d} -> kappa(e). kappa(e) > 0 because [M^-1]_top is the inverse of a Schur complement of a positive definite matrix, and phi_top(e) is not zero for e != 0. QED

Answer to (ii): yes. If the law's margin has degree below the library degree, the horizon is finite in every direction, at every sample size and noise level, with the explicit radius in (a). If the degrees are equal, the horizon can be infinite, and whether it is depends on the top-degree signal-to-noise ratio in (b). The usual case in this project is a truth of degree 2 fitted with a library of degree 3 or 4 (A at d = 3, 4). Every single-equation law there has margin degree <= 2 < d, so its horizon is finite.

## 8. Libraries that are not full degree (item (iii))

Let the library be a set of monomials Mset with top degree d_M, and let d' be the largest k such that every monomial of degree <= k is in Mset. For a direction e, the line degree is d_M(e) = max{|beta| : beta in Mset, e^beta != 0}. If Mset is downward closed (as in A-mis, the degree-2 library without x^2), its span is translation invariant and the centre of Omega does not matter.

**Proposition 8.** Iid noise, data in Omega, Kbar_M and N2_M defined with the library Mset.

(a) Kbar_M <= Kbar_L at every point, for any library L that contains Mset, on the same data.
(b) For u = g e with g >= 1: Kbar_M(u) <= N2_M T_{d_M(e)}(g)^2.
(c) For g(u) >= 1: Kbar_M(u) >= max( T_{d'}(g)^2, max over beta in Mset of (u^beta)^2 / max over Omega of (x^beta)^2 ).

Proof. (a) is Lemma 1, a max over a larger set. (b) is the proof of Theorem 3's upper side. On the line through e, every q in span(Mset) restricts to a polynomial of degree <= d_M(e), because the terms of higher degree all carry a factor e^beta = 0. (c) is Lemma 1 with the test polynomials T_{d'}(l*) (in the span, since every monomial of degree <= d' is there) and x^beta / max |x^beta|. Both have ||.||_n <= 1 on data in Omega. QED

For A-mis, d_M = 2 and d' = 1, and on the prey axis e = (1, 0) the line degree is 1, because y^2 and xy vanish there. So the half-width of the prey equation grows like g along that axis, not like g^2. Along that axis the true prey equation has its x^2 term, which A-mis cannot represent.

## 9. Numerical checks and pass rules (fixed before the run)

Script: `tools/check_law_horizon.py`, output `results/check_law_horizon.txt`, seed 20261127, runs locally in under a few minutes. ALL_OK is True only if every check passes. If one fails, I keep the failing case, narrow the statement, and say so here (ROADMAP_v3 kill rule).

| check | statement | setup | pass rule |
|---|---|---|---|
| H1 | Lemma 1 | 2-D, degree 3, n = 120 in the unit disk, AR(1) noise with correlation 0.6, 200 points with g in [1, 4], 1000 random q each | every ratio q(u)^2 / (||q||^2 K) <= 1 + 1e-9, and the ratio at c = G^-1 phi equals 1 within 1e-9 |
| H2 | Lemma 2 | disk and square, d = 1..5, 3 directions, g in {1.2, 2, 4}; LP max q(u) subject to \|q\| <= 1 on 5000 points of the body (2000 on the boundary, 3000 inside) | LP / T_d(g) in [1 - 1e-7, 1.02] (the grid only relaxes the constraint, so the LP can only be high) |
| H3 | Theorems 3 and 4 | 300 points uniform in the unit disk, d = 2, 3, 4, 400 random points with g in [1.25, 8]; N2 from 24,000 samples plus Nelder-Mead; 1-D, 200 uniform points | the chain T_d^2 <= test bound <= projected bound <= Kbar <= N2 T_d^2 holds at every point (rel. 1e-9), random directions never exceed Kbar, N2 >= p; 1-D projected bound equals Kbar within 1e-8 |
| H4 | Theorem 5 | 2-D, degree 2 truth in the library, n = 200, sigma = 0.05, 20 points with z in [-4, 4], 40,000 replicates | every zone frequency within 4 binomial SE + 0.001 of the formula |
| H5 | Corollary 6 | Chebyshev nodes, n in {50, 200, 800}, d in {2, 3, 4}, mu0 in {0.1, 0.3}, sigma = 0.05 | Kbar = 1 + 2 sum T_k^2 within 1e-8 on g in [1, 20]; g_med in [T_d^-1(R / sqrt(2d + 1)), T_d^-1(R)] in all 18 cases |
| H6 | Proposition 7 | Chebyshev nodes, n = 200, mu = 0.2 + 0.5 t (d_m = 1, C = 0.7), d = 1, 2, 3 | d = 1: z(10^4) within 1e-3 of the limit in (b), and z >= sqrt q on [1, 10^4]; d = 2, 3: z < sqrt q on [bound radius, 10^4] |
| H7 | Proposition 8 | degree-2 library without x^2 vs the full one, 300 points in the unit disk | (a) at 2000 points; (b) on both halves of the prey axis, g in [1.25, 8]; (c) at every sampled point with g >= 1; all within 1e-9 |

Also reported, not gated: N2 against p, median Kbar / projected bound and Kbar / T_d test bound (how much section 4 gains), the actual horizon crossing against the H6 bound radius, and how many prey-axis points the full library puts above the H7 line bound.

What this unit does not do: it does not apply Corollary 6 to the 0F-3 runs. Those used the stacked errors-in-variables set, where Theorem 5 is not exact, and the horizon data are already seen. A post-hoc comparison can be a separate, labelled section of the write-up.

---

Note (2026-11-27, after the first run, before H4b ran). H4 found only 9 usable points, all with z > 0 (0.77 to 3.75), not 20 points over [-4, 4]. So the refuted formula Phi(-z - sqrt q) was only tested where it is close to 0. Addendum H4b (`tools/check_law_horizon_addendum.py`, seed 20261128) repeats H4 with the truth negated so that z < 0. Same n, sigma and replicates, candidates drawn out to g = 6, one point per 0.4-wide bin of z in [-4, 0). Pass rule: at least 8 points, and every zone frequency within 4 binomial SE + 0.001 of the formula.

## Results

(to be filled in after the run)
