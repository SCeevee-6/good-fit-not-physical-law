"""Protocol S1: Galerkin truncation, energy coordinates, closure target and the exact energy-stability test.

Truncation K keeps the Fourier modes with |kx|, |ky| <= K, k != 0. Energy coordinates: with c_k = w_hat_k / N^2 on the
half-plane (kx > 0, or kx = 0 and ky > 0), a = [sqrt2 Re c_k / |k|, sqrt2 Im c_k / |k|], so the resolved energy is
|a|^2 / 2 and the energy rate of any vorticity tendency tau is a . b with b built from tau the same way.
Closure target (per snapshot): tau = -P_K[N(w)] + P_K[N(w_K)], the part of the resolved tendency the truncated system
does not compute itself (N = u . grad w). Linear closure b ~ C a; its energy rate is a' C a = a' Q a, Q = (C + C') / 2.
"""
import numpy as np
from fractions import Fraction


def half_modes(grid, K):
    sel = (np.abs(grid.KX) <= K) & (np.abs(grid.KY) <= K) & (grid.K2 > 0) & ((grid.KX > 0) | ((grid.KX == 0) & (grid.KY > 0)))
    return np.nonzero(sel)


def trunc_mask(grid, K):
    return (np.abs(grid.KX) <= K) & (np.abs(grid.KY) <= K) & (grid.K2 > 0)


def coords(grid, K, zh):
    """Energy coordinates of spectral field zh (vorticity or a vorticity tendency) on the half-plane modes."""
    idx = half_modes(grid, K); c = zh[idx] / grid.N ** 2; kk = np.sqrt(grid.K2[idx])
    return np.concatenate([np.sqrt(2) * c.real / kk, np.sqrt(2) * c.imag / kk])


def from_coords(grid, K, a):
    """Inverse of coords for a real vector a (builds a Hermitian-consistent rfft2 array)."""
    idx = half_modes(grid, K); m = len(idx[0]); kk = np.sqrt(grid.K2[idx])
    c = (a[:m] + 1j * a[m:]) * kk / np.sqrt(2)
    zh = np.zeros(grid.K2.shape, complex); zh[idx] = c * grid.N ** 2
    # kx = 0 column: fill the conjugate (-ky) entries
    ky0 = idx[1] == 0
    if ky0.any():
        rows = idx[0][ky0]; zh[(-rows) % grid.N, 0] = np.conj(zh[rows, 0])
    return zh


def closure_target(grid, K, wh, NLh=None):
    """(a, b): energy coordinates of the resolved state and of tau = -P_K N(w) + P_K N(w_K)."""
    m = trunc_mask(grid, K)
    if NLh is None:
        NLh = grid.nonlinear(wh)
    wK = wh * m
    tau = (-NLh + grid.nonlinear(wK)) * m
    return coords(grid, K, wK), coords(grid, K, tau)


def advection_rate(grid, K, a):
    """Energy rate of the truncated advection -P_K N(w_K) at resolved state a, and the norms for the relative error."""
    wK = from_coords(grid, K, a); n = coords(grid, K, -grid.nonlinear(wK) * trunc_mask(grid, K))
    return float(a @ n), float(np.linalg.norm(a) * np.linalg.norm(n))


def fit_closure(A, B):
    """Least squares C with B ~ A C' (rows are snapshots). Returns C and the fraction of variance explained."""
    Ct, *_ = np.linalg.lstsq(A, B, rcond=None)
    R = B - A @ Ct
    return Ct.T, float(1 - (R ** 2).sum() / ((B - B.mean(0)) ** 2).sum())


def nsd_exact(Q):
    """Exact test that Q (floats, converted exactly to rationals) is negative semidefinite. Gaussian elimination
    without pivoting on M = -Q: M is PSD iff every pivot is >= 0 and every zero pivot has a zero remaining row.
    Returns (is_nsd, step reached, pivot at failure as float or None)."""
    n = len(Q); M = [[-Fraction(float(Q[i][j])) for j in range(n)] for i in range(n)]
    for k in range(n):
        p = M[k][k]
        if p < 0:
            return False, k, float(p)
        if p == 0:
            if any(M[k][j] != 0 for j in range(k + 1, n)):
                return False, k, 0.0
            continue
        for i in range(k + 1, n):
            f = M[i][k] / p
            if f:
                Mi, Mk = M[i], M[k]
                for j in range(k + 1, n):
                    Mi[j] -= f * Mk[j]
    return True, n, None
