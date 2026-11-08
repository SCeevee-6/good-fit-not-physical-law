"""Own STLSQ sparse regression (no hidden dependencies), candidate pool generation and conventional validation battery."""
import numpy as np
from scipy.integrate import solve_ivp
from lib import theta, make_rhs, round_xi, support, exponents
from data import stack


def _ls(Th, y, ridge=1e-10):
    """Column-scaled ridge least squares."""
    s = np.linalg.norm(Th, axis=0)
    s[s == 0] = 1.0
    A = Th / s
    coef = np.linalg.solve(A.T @ A + ridge * np.eye(A.shape[1]), A.T @ y)
    return coef / s


def stlsq(Th, dX, thresh, max_iter=12, ridge=1e-10):
    """Sequentially-thresholded least squares, independently per state equation. Returns Xi (p, n)."""
    p, n = Th.shape[1], dX.shape[1]
    Xi = np.zeros((p, n))
    for i in range(n):
        act = np.ones(p, dtype=bool)
        for _ in range(max_iter):
            c = np.zeros(p)
            if act.sum() == 0:
                break
            c[act] = _ls(Th[:, act], dX[:, i], ridge)
            new = act & (np.abs(c) >= thresh)
            if (new == act).all():
                break
            act = new
        else:
            pass
        c = np.zeros(p)
        if act.sum() > 0:
            c[act] = _ls(Th[:, act], dX[:, i], ridge)
        Xi[:, i] = c
    return round_xi(Xi)


# ------------------------------------------------------------------------------------------
# Conventional validation battery
# ------------------------------------------------------------------------------------------
BLOWUP = 1e3


def rollout(Xi, E, t, x0, tol_blow=BLOWUP):
    """Integrate a candidate model; returns (traj or None if blow-up/failure)."""
    f = make_rhs(Xi, E)
    def ev(t_, x):
        return tol_blow - np.max(np.abs(x))
    ev.terminal = True
    try:
        sol = solve_ivp(f, (t[0], t[-1]), x0, t_eval=t, rtol=1e-6, atol=1e-9, method='RK45', events=ev, max_step=(t[-1] - t[0]) / 50)
    except Exception:
        return None
    if (not sol.success) or sol.y.shape[1] != len(t) or not np.all(np.isfinite(sol.y)):
        return None
    return sol.y.T


def rollout_rmse(Xi, E, trajs, state_sd, penalty=10.0):
    """RMSE of rolled-out states vs. noisy data, normalised by state sd; failures get a fixed penalty."""
    errs = []
    for tr in trajs:
        sl = tr['sl']
        t = tr['t'][sl]
        x0 = tr['Xs'][0]
        Y = rollout(Xi, E, t, x0)
        if Y is None:
            errs.append(penalty)
        else:
            e = np.sqrt((((Y - tr['X'][sl]) / state_sd) ** 2).mean())
            errs.append(min(e, penalty))
    return float(np.mean(errs))


def deriv_rmse(Xi, E, trajs, deriv_sd):
    Xs, dX = stack(trajs)
    pred = theta(Xs, E) @ Xi
    return float(np.sqrt((((pred - dX) / deriv_sd) ** 2).mean())), float(1 - ((pred - dX) ** 2).sum() / ((dX - dX.mean(0)) ** 2).sum())


def aic_score(Xi, E, trajs_train, deriv_sd):
    Xs, dX = stack(trajs_train)
    pred = theta(Xs, E) @ Xi
    m = Xs.shape[0]
    rss = (((pred - dX) / deriv_sd) ** 2).sum()
    k = int((Xi != 0).sum())
    return float(m * np.log(rss / m + 1e-300) + 2 * k)
