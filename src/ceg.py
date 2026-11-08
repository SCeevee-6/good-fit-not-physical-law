"""Counterexample-guided repair (CEGIS) with constrained sparse regression, plus comparison arms.

Learner : constrained STLSQ; each constraint row is  r . vec(Xi) >= rhs  (linear because every spec is linear in the coefficients).
Verifier: z3 / Bernstein (verify.verify_spec_robust).
"""
import time
import numpy as np
import cvxpy as cp

from lib import theta, round_xi
from data import stack
from systems import spec_row, in_domain
from verify import verify_spec_robust
from oracle import _sample_domain

MARGIN = 0.25   # repair target: G(u*) >= -MARGIN * tol_strong at every added point


def fit_constrained(Th, dX, th, keep, rows, rhs, max_iter=10):
    """Constrained STLSQ. rows: list of (p, n) arrays; rhs: list of floats. Returns Xi (p, n) or None if infeasible."""
    p, n = Th.shape[1], dX.shape[1]
    m = Th.shape[0]
    act = np.tile(keep[:, None], (1, n))
    scale = np.linalg.norm(Th, axis=0) / np.sqrt(m)
    scale[scale == 0] = 1.0
    Ths = Th / scale
    Xi = np.zeros((p, n))
    for _ in range(max_iter):
        idx = [(j, i) for i in range(n) for j in range(p) if act[j, i]]
        if not idx:
            return np.zeros((p, n))
        z = cp.Variable(len(idx))
        # objective: sum_i ||Ths[:, A_i] z_i - dX_i||^2 / m, via per-equation QR
        obj = 0
        pos = 0
        for i in range(n):
            A = [j for j in range(p) if act[j, i]]
            if not A:
                continue
            Q, R = np.linalg.qr(Ths[:, A])
            qy = Q.T @ dX[:, i]
            obj += cp.sum_squares(R @ z[pos:pos + len(A)] - qy) / m
            pos += len(A)
        obj += 1e-10 * cp.sum_squares(z)
        cons = []
        if rows:
            Rm = np.array([[r[j, i] / scale[j] for (j, i) in idx] for r in rows])
            cons.append(Rm @ z >= np.array(rhs))
        prob = cp.Problem(cp.Minimize(obj), cons)
        try:
            prob.solve(solver=cp.CLARABEL)
        except Exception:
            return None
        if prob.status not in ('optimal', 'optimal_inaccurate') or z.value is None:
            return None
        Xi = np.zeros((p, n))
        for k, (j, i) in enumerate(idx):
            Xi[j, i] = z.value[k] / scale[j]
        new = act & (np.abs(Xi) >= th)
        if (new == act).all():
            break
        act = new
    return round_xi(Xi)


def row_at(spec, E, u):
    """Constraint row (p, n) and rhs for spec at point u."""
    r = spec_row(spec, [np.asarray(x, dtype=float) for x in u], E)
    r = np.asarray(r).reshape(E.shape[0], -1)
    return r, -MARGIN * spec.tol_strong


def worst_point(spec, Xi, E, m=None):
    """Numeric argmin of G over a grid of the spec domain (used by the CEG-worst ablation)."""
    dim = len(spec.u_names)
    m = m or {1: 2001, 2: 81, 3: 25}[dim]
    axes = [np.linspace(lo, hi, m) for lo, hi in spec.u_box]
    mesh = np.meshgrid(*axes, indexing='ij')
    u = [g.ravel() for g in mesh]
    ok = in_domain(spec, u)
    u = [x[ok] for x in u]
    rows = spec_row(spec, u, E)
    G = np.tensordot(rows, Xi, axes=([-2, -1], [0, 1]))
    i = int(G.argmin())
    return [x[i] for x in u], float(G[i])


def _verify_all(S, Xi, E):
    bad, timeout = [], False
    for sp in S.specs:
        r = verify_spec_robust(sp, Xi, E, 'strong', prefer='bernstein')
        if r['status'] == 'counterexample':
            bad.append((sp, r['u']))
        elif r['status'] == 'timeout':
            timeout = True
    return bad, timeout


def ceg_repair(S, Th, dX, E, th, keep, mode='witness', max_iter=15):
    """mode: 'witness' (pure CEG) or 'worst' (witness + numeric worst-case point). Starts from the unconstrained fit."""
    rows, rhs = [], []
    Xi = fit_constrained(Th, dX, th, keep, rows, rhs)
    hist = []
    for it in range(max_iter + 1):
        if Xi is None:
            return dict(status='infeasible', Xi=None, iters=it, n_rows=len(rows))
        bad, timeout = _verify_all(S, Xi, E)
        hist.append(len(bad))
        if not bad:
            return dict(status='repaired' if not timeout else 'timeout', Xi=Xi, iters=it, n_rows=len(rows), hist=hist)
        if it == max_iter:
            break
        for sp, u in bad:
            r, b = row_at(sp, E, u); rows.append(r); rhs.append(b)
            if mode == 'worst':
                uw, _ = worst_point(sp, Xi, E)
                r, b = row_at(sp, E, uw); rows.append(r); rhs.append(b)
        Xi = fit_constrained(Th, dX, th, keep, rows, rhs)
    return dict(status='failed_max_iter', Xi=Xi, iters=max_iter, n_rows=len(rows), hist=hist)


def random_repair(S, Th, dX, E, th, keep, specs_bad, n_rows, rng):
    """One-shot refit with n_rows constraint rows at random domain points (specs drawn from those violated by the initial model)."""
    rows, rhs = [], []
    for k in range(n_rows):
        sp = specs_bad[k % len(specs_bad)]
        u = _sample_domain(sp, 1, rng)
        r, b = row_at(sp, E, [x[0] for x in u]); rows.append(r); rhs.append(b)
    Xi = fit_constrained(Th, dX, th, keep, rows, rhs)
    return finalize(S, Xi, E, len(rows))


def dense_repair(S, Th, dX, E, th, keep, n_per_spec=100):
    """One-shot refit with a dense grid of constraints (about n_per_spec points per spec)."""
    rows, rhs = [], []
    for sp in S.specs:
        dim = len(sp.u_names)
        m = {1: n_per_spec, 2: int(np.sqrt(n_per_spec)) + 1, 3: int(round(n_per_spec ** (1 / 3))) + 1}[dim]
        axes = [np.linspace(lo, hi, m) for lo, hi in sp.u_box]
        mesh = np.meshgrid(*axes, indexing='ij')
        u = [g.ravel() for g in mesh]
        ok = in_domain(sp, u)
        u = [x[ok] for x in u]
        R = spec_row(sp, u, E).reshape(len(u[0]), E.shape[0], -1)
        for r in R:
            rows.append(r); rhs.append(-MARGIN * sp.tol_strong)
    Xi = fit_constrained(Th, dX, th, keep, rows, rhs)
    return finalize(S, Xi, E, len(rows))


def finalize(S, Xi, E, n_rows):
    if Xi is None:
        return dict(status='infeasible', Xi=None, iters=0, n_rows=n_rows)
    bad, timeout = _verify_all(S, Xi, E)
    return dict(status='repaired' if (not bad and not timeout) else ('timeout' if timeout and not bad else 'failed'), Xi=Xi, iters=0, n_rows=n_rows)


def repair_candidate(S, cell, row, seed, drop=None):
    """Run all repair arms for one violating candidate (row = pool row dict with degree, variant, threshold, xi, xi_shape). Returns list of arm records."""
    import pool as _pool
    from sindy import deriv_rmse
    d = int(row['degree']); E = S.E(d)
    drop = [tuple(x) for x in (drop or [])]
    keep = np.array([tuple(e) not in drop for e in E.tolist()])
    variants = dict(_pool.make_variants(len(cell['train']), seed))
    idx = variants[row['variant']]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    Ttrue = S.true_Xi(d)
    th = float(row['threshold'])
    Xi0 = np.array(row['xi']).reshape(row['xi_shape'])
    specs_bad = [sp for sp in S.specs if verify_spec_robust(sp, Xi0, E, 'strong', prefer='bernstein')['status'] == 'counterexample']
    out = []
    vid = 99 if row['variant'] == 'full' else int(str(row['variant'])[4:])
    rng = np.random.default_rng([seed, vid, int(row['degree']), int(round(th * 1e6))])
    def rec(arm, r, t0):
        Xi = r.get('Xi')
        o = dict(arm=arm, status=r['status'], iters=r.get('iters'), n_rows=r.get('n_rows'), seconds=time.time() - t0)
        if Xi is not None:
            rm, _ = deriv_rmse(Xi, E, cell['iid'], cell['deriv_sd'])
            sup, sup_t = Xi != 0, Ttrue != 0
            tp = int((sup & sup_t).sum()); fp = int((sup & ~sup_t).sum()); fn = int((~sup & sup_t).sum())
            o.update(k=int(sup.sum()), iid_rmse=rm, support_f1=2 * tp / max(2 * tp + fp + fn, 1), coef_err=float(np.linalg.norm(Xi - Ttrue) / np.linalg.norm(Ttrue)),
                     xi=Xi.ravel().tolist())
        return o
    t0 = time.time(); r_w = ceg_repair(S, Th, dX, E, th, keep, 'witness'); out.append(rec('ceg_witness', r_w, t0))
    t0 = time.time(); r_x = ceg_repair(S, Th, dX, E, th, keep, 'worst'); out.append(rec('ceg_worst', r_x, t0))
    K = max(r_w.get('n_rows', 1), 1)
    t0 = time.time(); out.append(rec('random_constraints', random_repair(S, Th, dX, E, th, keep, specs_bad, K, rng), t0))
    t0 = time.time(); out.append(rec('dense_constraints', dense_repair(S, Th, dX, E, th, keep), t0))
    return out
