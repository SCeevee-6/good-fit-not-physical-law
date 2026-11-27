"""Addendum H4b to docs/law_horizon.md (written after the first run, before this one).

The frozen H4 only found points with z > 0, so the refuted-zone formula Phi(-z - sqrt q) was tested only where it is
close to 0. H4b repeats H4 with the truth negated (z < 0), same n, sigma, replicates and pass rule.

Usage: python tools/check_law_horizon_addendum.py [out_file]
"""
import sys
import numpy as np
from scipy.stats import chi2, norm
from check_law_horizon import expo, phi, disk


def h4b(rng, reps=40000, n=200, sigma=0.05):
    E = expo(2, 2); Xd = disk(rng, n); Th = phi(Xd, E); p = len(E); q = chi2.ppf(0.95, p)
    xi = -np.array([0.3, -0.2, 0.1, 0.05, -0.08, 0.1])
    G = Th.T @ Th / sigma ** 2; A = np.linalg.solve(Th.T @ Th, Th.T)
    cand = disk(rng, 20000) * rng.uniform(0.5, 6, (20000, 1)); P = phi(cand, E)
    K = np.einsum('ij,ij->i', P, np.linalg.solve(G, P.T).T); z = P @ xi / np.sqrt(K)
    pick = [np.nonzero((z >= lo) & (z < lo + 0.4))[0][0] for lo in np.arange(-4.0, 0.0, 0.4) if np.any((z >= lo) & (z < lo + 0.4))]
    P, K, z = P[pick], K[pick], z[pick]
    C = np.zeros(len(pick)); Rf = np.zeros(len(pick))
    for _ in range(reps // 10000):
        M = (Th @ xi + sigma * rng.normal(size=(10000, n))) @ (A.T @ P.T)
        C += (M >= np.sqrt(q * K)).sum(0); Rf += (M < -np.sqrt(q * K)).sum(0)
    C /= reps; Rf /= reps; worst = 0.0
    for f, pth in ((C, norm.cdf(z - np.sqrt(q))), (Rf, norm.cdf(-z - np.sqrt(q)))):
        worst = max(worst, float(np.max(np.abs(f - pth) / (4 * np.sqrt(pth * (1 - pth) / reps) + 1e-3))))
    return dict(points=len(pick), z_range=[round(float(z.min()), 2), round(float(z.max()), 2)],
                refuted_freq_range=[round(float(Rf.min()), 4), round(float(Rf.max()), 4)],
                worst_scaled_dev=round(worst, 3), ok=bool(worst <= 1 and len(pick) >= 8))


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "results/check_law_horizon_addendum.txt"
    r = h4b(np.random.default_rng(20261128))
    open(out, "w").write(f"H4b zone probs, z < 0: {r}\nALL_OK {r['ok']}\n")
    print("ALL_OK", r["ok"], r["points"], r["z_range"])
    sys.exit(0 if r["ok"] else 1)
