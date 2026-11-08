"""Polynomial library utilities shared by every module.

A model is a coefficient matrix Xi with shape (p, n): f_i(x) = sum_j Theta_j(x) * Xi[j, i],
where Theta_j(x) = prod_k x_k ** E[j, k] and E is the exponent matrix from `exponents`.
Coefficients are always rounded to 6 decimals so every candidate is an exact rational model.
"""
import itertools
from fractions import Fraction
import numpy as np

COEF_DECIMALS = 6


def exponents(n, d):
    """All exponent tuples of total degree <= d in n variables, graded-lexicographic order."""
    out = []
    for deg in range(d + 1):
        combos = [e for e in itertools.product(range(deg + 1), repeat=n) if sum(e) == deg]
        combos.sort(reverse=True)
        out.extend(combos)
    return np.array(out, dtype=int)


def theta(X, E):
    """Evaluate library at points X (m, n) -> (m, p)."""
    X = np.atleast_2d(X)
    m = X.shape[0]
    out = np.ones((m, E.shape[0]))
    for k in range(E.shape[1]):
        maxe = int(E[:, k].max())
        if maxe == 0:
            continue
        pw = np.stack([X[:, k] ** e for e in range(maxe + 1)], axis=1)
        out *= pw[:, E[:, k]]
    return out


def term_name(e, names):
    parts = []
    for k, ek in enumerate(e):
        if ek == 1:
            parts.append(names[k])
        elif ek > 1:
            parts.append(f"{names[k]}^{ek}")
    return "*".join(parts) if parts else "1"


def round_xi(Xi):
    return np.round(np.asarray(Xi, dtype=float), COEF_DECIMALS) + 0.0  # +0.0 removes -0.0


def to_fraction(c):
    """Exact rational from a coefficient that is already on the 1e-6 grid."""
    return Fraction(int(round(float(c) * 10 ** COEF_DECIMALS)), 10 ** COEF_DECIMALS)


def make_rhs(Xi, E):
    """Fast RHS for solve_ivp for coefficient matrix Xi (p, n)."""
    Xi = np.asarray(Xi, dtype=float)
    E = np.asarray(E)
    n = E.shape[1]
    maxe = [int(E[:, k].max()) for k in range(n)]

    def rhs(t, x):
        out = np.ones(E.shape[0])
        for k in range(n):
            if maxe[k] == 0:
                continue
            pw = x[k] ** np.arange(maxe[k] + 1)
            out *= pw[E[:, k]]
        return out @ Xi
    return rhs


def embed_true(terms, E, n):
    """terms: list (per state) of {exponent_tuple: coefficient}. Returns Xi (p, n) in library E."""
    index = {tuple(e): j for j, e in enumerate(E.tolist())}
    Xi = np.zeros((E.shape[0], n))
    for i, tdict in enumerate(terms):
        for e, c in tdict.items():
            Xi[index[tuple(e)], i] = c
    return Xi


def support(Xi, E, names):
    """Return list of (state_index, term_name) for nonzero coefficients."""
    return [(i, term_name(E[j], names)) for j in range(Xi.shape[0]) for i in range(Xi.shape[1]) if Xi[j, i] != 0]
