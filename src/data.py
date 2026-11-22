"""Data generation: clean ODE trajectories, measurement noise, Savitzky-Golay smoothing/derivatives."""
import numpy as np
from scipy.integrate import solve_ivp
try:
    from scipy.signal import savgol_filter
except ImportError:                      # scipy.signal DLL blocked on some machines; numpy copy of the same filter
    from savgol_np import savgol_filter

SYS_ID = {'A_predator_prey': 0, 'B_sir': 1, 'H_holling': 2}


def simulate_true(sys, ic):
    f = sys.f_true()
    t_eval = np.arange(0.0, sys.T + sys.dt / 2, sys.dt)
    sol = solve_ivp(f, (0, sys.T), ic, t_eval=t_eval, rtol=1e-10, atol=1e-12, method='LSODA')
    if not sol.success:
        raise RuntimeError('true-system integration failed')
    return sol.t, sol.y.T


def smooth(X, sys):
    """SG smoothing + first derivative, edges trimmed. Returns Xs, dX, trimmed index slice."""
    w, p, dt = sys.sg_window, sys.sg_poly, sys.dt
    Xs = savgol_filter(X, w, p, axis=0, mode='interp')
    dX = savgol_filter(X, w, p, deriv=1, delta=dt, axis=0, mode='interp')
    h = w // 2
    sl = slice(h, X.shape[0] - h)
    return Xs[sl], dX[sl], sl


def make_cell(sys, rho, seed, noise_idx=0):
    """Return dict of trajectory lists plus noise scales. Deterministic in (system, noise level index, seed)."""
    rng = np.random.default_rng([seed, SYS_ID[sys.name], noise_idx, 12345])
    def gen(ic_fn, k):
        out = []
        for _ in range(k):
            ic = ic_fn(rng)
            t, X = simulate_true(sys, ic)
            out.append({'ic': ic, 't': t, 'X_clean': X})
        return out
    train = gen(sys.ic_train, sys.n_train)
    iid = gen(sys.ic_train, sys.n_iid)
    ood = gen(sys.ic_ood, sys.n_ood)
    pooled = np.vstack([tr['X_clean'] for tr in train])
    sd = pooled.std(axis=0)
    sigma = rho * sd
    for grp in (train, iid, ood):
        for tr in grp:
            tr['X'] = tr['X_clean'] + rng.normal(size=tr['X_clean'].shape) * sigma
            Xs, dX, sl = smooth(tr['X'], sys)
            tr['Xs'], tr['dX'], tr['sl'] = Xs, dX, sl
            # clean derivative on the same trimmed window (used only for diagnostics)
            f = sys.f_true()
            tr['dX_clean'] = np.array([f(0, x) for x in tr['X_clean'][sl]])
    dstd = np.vstack([tr['dX'] for tr in train]).std(axis=0)
    return {'train': train, 'iid': iid, 'ood': ood, 'sigma': sigma, 'state_sd': sd, 'deriv_sd': dstd, 'rho': rho, 'seed': seed}


def stack(trajs, idx=None):
    """Stack smoothed states and derivatives over trajectories (optionally a bootstrap index list)."""
    if idx is None:
        idx = range(len(trajs))
    Xs = np.vstack([trajs[i]['Xs'] for i in idx])
    dX = np.vstack([trajs[i]['dX'] for i in idx])
    return Xs, dX
