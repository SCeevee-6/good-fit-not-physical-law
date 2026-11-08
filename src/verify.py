"""Exact SMT verifier (z3 nlsat), exact 1-D cross-check (sympy), and trajectory replay of counterexamples.

Spec form:  G(u) = c . f(point(u)) >= -tol  on domain D.  A violation exists iff  exists u in D: G(u) < -tol.
"""
import time
from fractions import Fraction
import numpy as np
import z3
import sympy as sp
from scipy.integrate import solve_ivp

from lib import to_fraction, make_rhs
from systems import spec_row


class ZK:
    def __call__(self, fr):
        fr = Fraction(fr)
        return z3.RealVal(f"{fr.numerator}/{fr.denominator}")


def _is_zero_z3(e):
    try:
        e = z3.simplify(e)
        return z3.is_rational_value(e) and e.numerator_as_long() == 0
    except Exception:
        return False


def _build_G_z3(spec, Xi, E, u):
    K = ZK()
    pt = [z3.simplify(p) if not isinstance(p, (int, float)) else z3.RealVal(p) for p in spec.point(u, K)]
    zero = [_is_zero_z3(p) for p in pt]
    terms = []
    for j in range(E.shape[0]):
        if any(E[j, k] > 0 and zero[k] for k in range(E.shape[1])):
            continue
        mono = None
        for k in range(E.shape[1]):
            for _ in range(int(E[j, k])):
                mono = pt[k] if mono is None else mono * pt[k]
        for i in range(Xi.shape[1]):
            if Xi[j, i] == 0 or spec.c[i] == 0:
                continue
            coef = K(to_fraction(Xi[j, i]) * Fraction(spec.c[i]))
            terms.append(coef if mono is None else coef * mono)
    if not terms:
        return z3.RealVal(0)
    return z3.Sum(terms)


def _exact_G(spec, Xi, E, uF):
    """Exact rational evaluation of G at a rational point."""
    pt = spec.point([Fraction(x) for x in uF], Fraction)
    pt = [Fraction(p) for p in pt]
    G = Fraction(0)
    for j in range(E.shape[0]):
        mono = Fraction(1)
        for k in range(E.shape[1]):
            mono *= pt[k] ** int(E[j, k])
        if mono == 0:
            continue
        for i in range(Xi.shape[1]):
            if Xi[j, i] != 0 and spec.c[i] != 0:
                G += to_fraction(Xi[j, i]) * Fraction(spec.c[i]) * mono
    return G


def _in_domain_exact(spec, uF):
    return all(sum(Fraction(a_k) * u_k for a_k, u_k in zip(a, uF)) >= Fraction(b) for a, b in spec.domain)


def verify_spec(spec, Xi, E, level='strong', timeout_ms=10000, box=None):
    """Return dict(status, u, G, time).  status in {'verified','counterexample','timeout'}."""
    tol = Fraction(spec.tol(level)).limit_denominator(10 ** 6)
    t0 = time.time()
    u = [z3.Real(n) for n in spec.u_names]
    G = _build_G_z3(spec, Xi, E, u)
    s = z3.SolverFor('QF_NRA')
    s.set('timeout', timeout_ms)
    for a, b in spec.domain:
        s.add(z3.Sum([z3.RealVal(f"{Fraction(ak).numerator}/{Fraction(ak).denominator}") * uk for ak, uk in zip(a, u)])
              >= z3.RealVal(f"{Fraction(b).numerator}/{Fraction(b).denominator}"))
    if box is not None:
        for uk, (lo, hi) in zip(u, box):
            s.add(uk >= z3.RealVal(f"{Fraction(lo).numerator}/{Fraction(lo).denominator}"))
            s.add(uk <= z3.RealVal(f"{Fraction(hi).numerator}/{Fraction(hi).denominator}"))
    s.add(G < z3.RealVal(f"{tol.numerator}/{tol.denominator}") * -1)
    res = s.check()
    out = dict(spec=spec.name, level=level, time=time.time() - t0, u=None, G=None, status=None, witness_exact=None)
    if res == z3.unsat:
        out['status'] = 'verified'
    elif res == z3.sat:
        m = s.model()
        uF = []
        for uk in u:
            v = m.eval(uk, model_completion=True)
            if z3.is_rational_value(v):
                uF.append(Fraction(v.numerator_as_long(), v.denominator_as_long()))
            elif z3.is_algebraic_value(v):
                a = v.approx(40)
                uF.append(Fraction(a.numerator_as_long(), a.denominator_as_long()))
            else:
                uF.append(Fraction(0))
        Gx = _exact_G(spec, Xi, E, uF)
        ok = _in_domain_exact(spec, uF) and Gx < -tol
        out.update(status='counterexample', u=[float(x) for x in uF], G=float(Gx), witness_exact=bool(ok))
    else:
        out['status'] = 'timeout'
    return out


def _split_box(box):
    """Bisect every coordinate of a box (list of (lo, hi) Fractions) -> 2^dim sub-boxes."""
    import itertools
    halves = []
    for lo, hi in box:
        mid = (Fraction(lo) + Fraction(hi)) / 2
        halves.append([(Fraction(lo), mid), (mid, Fraction(hi))])
    return [list(c) for c in itertools.product(*halves)]


def verify_spec_robust(spec, Xi, E, level='strong', timeout_ms=3000, sub_timeout_ms=3000, max_depth=2, budget_s=30.0, prefer='z3'):
    """z3 on the whole domain; on time-out fall back to (3-D simplex) Bernstein branch-and-bound, else recursive sub-box splitting.
    Every verdict is exact (witnesses are re-evaluated in rational arithmetic); `method` records which procedure decided."""
    t0 = time.time()
    if prefer == 'bernstein' and len(spec.u_names) == 3:
        from bern import bernstein_verify
        br = bernstein_verify(spec, Xi, E, level, time_budget=10.0)
        if br['status'] != 'timeout':
            br['method'] = 'bernstein'
            return br
    r = verify_spec(spec, Xi, E, level, timeout_ms)
    if r['status'] != 'timeout':
        r['method'] = 'z3'
        return r
    if len(spec.u_names) == 3:
        from bern import bernstein_verify
        br = bernstein_verify(spec, Xi, E, level, time_budget=60.0)
        br.update(method='bernstein_after_z3_timeout', time=time.time() - t0)
        return br
    stack = [(b, 1) for b in _split_box([(Fraction(lo).limit_denominator(1000), Fraction(hi).limit_denominator(1000)) for lo, hi in spec.u_box])]
    unresolved, n_q = False, 0
    while stack:
        box, depth = stack.pop()
        if time.time() - t0 > budget_s:
            unresolved = True
            break
        rr = verify_spec(spec, Xi, E, level, sub_timeout_ms, box=box)
        n_q += 1
        if rr['status'] == 'counterexample':
            rr.update(time=time.time() - t0, method='z3_split', n_sub=n_q)
            return rr
        if rr['status'] == 'timeout':
            if depth < max_depth:
                stack.extend((b, depth + 1) for b in _split_box(box))
            else:
                unresolved = True
    r.update(time=time.time() - t0, method='z3_split', n_sub=n_q, status='timeout' if unresolved else 'verified')
    return r


def verify_model(S, Xi, E, timeout_ms=3000, levels=('strong', 'weak')):
    """Verify all specs of system S. Weak level only queried for specs verified at strong level."""
    res = {}
    for sp_ in S.specs:
        r = verify_spec_robust(sp_, Xi, E, 'strong', timeout_ms)
        res[(sp_.name, 'strong')] = r
        if 'weak' in levels:
            if r['status'] == 'verified':
                res[(sp_.name, 'weak')] = verify_spec_robust(sp_, Xi, E, 'weak', timeout_ms)
            else:
                res[(sp_.name, 'weak')] = dict(r, level='weak')   # strong violation implies weak violation (tol_weak < tol_strong)
    return res


# -------------------------------------------------------------------------------------------
# Exact 1-D cross-check (independent code path: sympy real-root isolation)
# -------------------------------------------------------------------------------------------
def exact_min_1d(spec, Xi, E):
    assert len(spec.u_names) == 1
    u = sp.Symbol('u', real=True)
    Ks = lambda fr: sp.Rational(Fraction(fr).numerator, Fraction(fr).denominator)
    pt = [sp.sympify(p) for p in spec.point([u], Ks)]
    G = sp.Integer(0)
    for j in range(E.shape[0]):
        mono = sp.Integer(1)
        for k in range(E.shape[1]):
            mono *= pt[k] ** int(E[j, k])
        for i in range(Xi.shape[1]):
            if Xi[j, i] != 0 and spec.c[i] != 0:
                fr = to_fraction(Xi[j, i]) * Fraction(spec.c[i])
                G += sp.Rational(fr.numerator, fr.denominator) * mono
    lo, hi = -sp.oo, sp.oo
    for a, b in spec.domain:
        a0, b0 = Fraction(a[0]), Fraction(b)
        bound = sp.Rational((b0 / a0).numerator, (b0 / a0).denominator)
        if a0 > 0:
            lo = max(lo, bound) if lo != -sp.oo else bound
        else:
            hi = min(hi, bound) if hi != sp.oo else bound
    G = sp.expand(G)
    cands = [lo, hi]
    P = sp.Poly(G, u)
    if P.degree() >= 2:
        for r in sp.real_roots(P.diff(u)):
            rv = sp.N(r, 50)
            if lo <= rv <= hi:
                cands.append(r)
    vals = [sp.N(G.subs(u, c), 50) for c in cands]
    return float(min(vals))


# -------------------------------------------------------------------------------------------
# Behavioural replay of a counterexample
# -------------------------------------------------------------------------------------------
def replay(spec, Xi, E, u_star, horizon=0.05, floor=1e-7):
    """Integrate the candidate ODE from the counterexample state. The behaviour functional c . x has derivative exactly G along trajectories,
    so a genuine violation must show c.x(t) - c.x(0) < -floor at some t in (0, horizon].  Returns (confirmed, worst_change)."""
    pt = spec.point([np.asarray(x, dtype=float) for x in u_star], lambda fr: float(fr))
    x0 = np.array([float(p) for p in pt])
    f = make_rhs(Xi, E)
    ts = np.linspace(0, horizon, 201)
    sol = solve_ivp(f, (0, horizon), x0, t_eval=ts, rtol=1e-12, atol=1e-14, method='LSODA')
    c = np.array([float(ci) for ci in spec.c])
    change = c @ sol.y - c @ x0
    worst = float(change.min())
    return worst < -floor, worst
