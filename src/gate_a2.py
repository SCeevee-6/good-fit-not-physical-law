"""Protocol 0A2 (docs/protocol_0A2.md): system B rerun with corrected structures B-pos and B-fac.
Arms, budget, sampler, margin and scoring are those of 0A (src/gate_a.py, frozen; reused here unchanged).
Only the B structure and the fit's fixed constraints differ.

Free unknowns are the I- and R-equations, Xi_full = Xi_free @ M_B.T, S-equation = -(I' + R').
  B-pos: I-eq keeps I-divisible monomials. R-eq keeps every monomial; R-free monomials have coefficient >= 0.
  B-fac: B-pos plus S | S', i.e. for every S-free monomial m: Xi_I[m] + Xi_R[m] = 0.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import cvxpy as cp
from lib import theta, round_xi
from systems import SYSTEMS, spec_G
from data import make_cell, stack
from verify import verify_spec_robust
from oracle import _sample_domain
from ceg import row_at
import gate_a
from gate_a import NOISE, M_B, make_variants, library, evaluate

STRUCTS = ('Bpos', 'Bfac')


def structure2(S, struct, E):
    """allowed (p,2), M, free cols, lb (p,2) bool mask for >= 0, eq list of (p,2) rows with row . Xi_free == 0."""
    p = E.shape[0]
    allowed = np.column_stack([E[:, 1] >= 1, np.ones(p, bool)])
    lb = np.zeros((p, 2), bool)
    lb[:, 1] = E[:, 2] == 0                      # R-free monomials in the R-equation
    eq = []
    if struct == 'Bfac':
        for j in np.where(E[:, 0] == 0)[0]:      # S-free monomials: coefficient in I'+R' is zero
            r = np.zeros((p, 2)); r[j, 0] = 1.0; r[j, 1] = 1.0
            eq.append(r)
    elif struct != 'Bpos':
        raise ValueError(struct)
    return allowed, M_B, [1, 2], lb, eq


def project_eq(Xi, eq):
    """Exact projection onto the pairwise equalities (each row is e_{j,I} + e_{j,R})."""
    Xi = Xi.copy()
    for r in eq:
        j = int(np.where(r[:, 0] != 0)[0][0])
        d = 0.5 * (Xi[j, 0] - Xi[j, 1])
        Xi[j, 0], Xi[j, 1] = d, -d
    return Xi


def fit_free2(Th, dXf, allowed, th, rows, rhs, lb, eq, max_iter=10):
    """gate_a.fit_free plus fixed constraints: Xi[lb] >= 0 and eq rows == 0 (on active entries)."""
    p, nf = Th.shape[1], dXf.shape[1]
    m = Th.shape[0]
    act = allowed.copy()
    scale = np.linalg.norm(Th, axis=0) / np.sqrt(m)
    scale[scale == 0] = 1.0
    Ths = Th / scale
    Xi = np.zeros((p, nf))
    for _ in range(max_iter):
        idx = [(j, f) for f in range(nf) for j in range(p) if act[j, f]]
        if not idx:
            return np.zeros((p, nf))
        z = cp.Variable(len(idx))
        obj, pos = 0, 0
        for f in range(nf):
            A = [j for j in range(p) if act[j, f]]
            if not A:
                continue
            Q, R = np.linalg.qr(Ths[:, A])
            qy = Q.T @ dXf[:, f]
            obj = obj + cp.sum_squares(R @ z[pos:pos + len(A)] - qy) / m
            pos += len(A)
        obj = obj + 1e-10 * cp.sum_squares(z)
        cons = []
        if rows:
            Rm = np.array([[r[j, f] / scale[j] for (j, f) in idx] for r in rows])
            cons.append(Rm @ z >= np.array(rhs))
        li = [k for k, (j, f) in enumerate(idx) if lb[j, f]]
        if li:
            cons.append(z[li] >= 0)
        if eq:
            Em = np.array([[r[j, f] / scale[j] for (j, f) in idx] for r in eq])
            Em = Em[np.abs(Em).sum(1) > 0]
            if len(Em):
                cons.append(Em @ z == 0)
        prob = cp.Problem(cp.Minimize(obj), cons)
        try:
            prob.solve(solver=cp.CLARABEL)
        except Exception:
            return None
        if prob.status not in ('optimal', 'optimal_inaccurate') or z.value is None:
            return None
        Xi = np.zeros((p, nf))
        for k, (j, f) in enumerate(idx):
            Xi[j, f] = z.value[k] / scale[j]
        new = act & (np.abs(Xi) >= th)
        # an equality pair is kept or dropped together, so thresholding cannot force a surviving partner to 0
        for r in eq:
            j = int(np.where(r[:, 0] != 0)[0][0])
            keep = new[j, 0] or new[j, 1]
            new[j, 0] = keep and allowed[j, 0]; new[j, 1] = keep and allowed[j, 1]
        if (new == act).all():
            break
        act = new
    Xi = project_eq(Xi, eq)                      # removes solver residual (~1e-9); rounding is sign-symmetric
    Xi[lb & (Xi < 0)] = 0.0                      # same: solver residual only
    return round_xi(Xi)


def run_arm2(arm, S, Th, dXf, allowed, th, M, E, lb, eq, rng, budget=15, samples=2000):
    rows, rhs = [], []
    Xi_free = fit_free2(Th, dXf, allowed, th, rows, rhs, lb, eq)
    iters = 0
    if arm == 'S' or Xi_free is None:
        return Xi_free, iters, 0, None
    stop = None
    for it in range(budget + 1):
        Xi_full = Xi_free @ M.T
        viol = []
        for sp in S.specs:
            if arm == 'C':
                r = verify_spec_robust(sp, Xi_full, E, 'strong', prefer='bernstein')
                if r['status'] == 'counterexample' and r.get('u') is not None:
                    viol.append((sp, r['u']))
            else:
                U = _sample_domain(sp, samples, rng)
                G = spec_G(sp, Xi_full, U, E)
                i = int(np.argmin(G))
                if G[i] < -sp.tol_strong:
                    viol.append((sp, [x[i] for x in U]))
        if not viol:
            stop = 'clean'; break
        if it == budget:
            stop = 'budget'; break
        for sp, u in viol:
            rr, b = row_at(sp, E, u)
            rows.append(rr @ M)
            rhs.append(b)
        iters += 1
        Xi_free = fit_free2(Th, dXf, allowed, th, rows, rhs, lb, eq)
        if Xi_free is None:
            return None, iters, len(rows), 'infeasible'
    return Xi_free, iters, len(rows), stop


def run_sample2(job, struct):
    job = dict(job, seed=int(job['seed']), degree=int(job['degree']), threshold=float(job['threshold']),
               cand=int(job['cand']), pool_iid=float(job['pool_iid']))
    assert job['key'] == 'B'
    rec0 = dict(job, struct=struct)
    seed, degree, variant, th = job['seed'], job['degree'], job['variant'], job['threshold']
    S = SYSTEMS['B']
    cell = make_cell(S, NOISE[1], seed, 1)
    E = library(S, 'B', degree)
    idx = dict(make_variants(len(cell['train']), seed))[variant]
    Xs, dX = stack(cell['train'], idx)
    Th = theta(Xs, E)
    allowed, M, free_cols, lb, eq = structure2(S, struct, E)
    dXf = dX[:, free_cols]
    rng = np.random.default_rng([seed, job['cand'], 5151])   # same stream as 0A
    out = []
    t0 = time.time()
    for arm in ('S', 'C', 'Rnd'):
        Xi_free, iters, nrows, stop = run_arm2(arm, S, Th, dXf, allowed, th, M, E, lb, eq, rng)
        rec = dict(rec0, arm=arm, iters=iters, n_rows=nrows, stop=stop)
        if Xi_free is None:
            rec.update(status='infeasible')
        else:
            Xi_full = Xi_free @ M.T
            rec.update(evaluate(S, 'B', Xi_full, E, cell))
            rec['iid_ratio'] = rec['iid_rmse'] / job['pool_iid']
            rec['xi'] = json.dumps(Xi_full.tolist())
            rec['status'] = 'ok'
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    return out
