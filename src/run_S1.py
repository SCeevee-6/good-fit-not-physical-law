"""Protocol S1 run (docs/protocol_S1.md with Addendum 1): one DNS, budget check, advection neutrality, closure data.
python src/run_S1.py  -> results/S1_run.json, results/S1_samples.npz
Env S1_QUICK=1: spin-up 1, sampling 1 time unit (pipeline test only; outputs are discarded).
Analysis (fits, eigenvalues, exact test, figure): python src/analyze_S1.py
"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from ns_dns import Grid, initial, run          # noqa: E402
import ns_galerkin as G                         # noqa: E402

QUICK = os.environ.get('S1_QUICK') == '1'
N, NU, DT, SEED = 64, 0.02, 1e-3, 0
T_SPIN, T_SAMPLE, EVERY = (1.0, 1.0, 100) if QUICK else (50.0, 100.0, 100)
KS = (1, 2, 3)
N_NEUTRAL = 1000


def main():
    t0 = time.time(); g = Grid(N, NU); wh = initial(g, SEED)
    n_spin = int(round(T_SPIN / DT)); n_samp = int(round(T_SAMPLE / DT))
    E1, P1, D1, wh = run(g, wh, DT, n_spin, EVERY)
    samples = {K: ([], []) for K in KS}; times = []
    def on_sample(i, w):
        NLh = g.nonlinear(w); times.append(T_SPIN + i * DT)
        for K in KS:
            a, b = G.closure_target(g, K, w, NLh); samples[K][0].append(a); samples[K][1].append(b)
    E2, P2, D2, wh = run(g, wh, DT, n_samp, EVERY, on_sample)
    E = np.concatenate([E1, E2[1:]]); P = np.concatenate([P1, P2[1:]]); D = np.concatenate([D1, D2[1:]])
    # S-G1: change in E against the trapezoid of P - D, whole run and every window of one time unit
    src = 0.5 * DT * ((P - D)[1:] + (P - D)[:-1]); scale = 0.5 * DT * ((np.abs(P) + np.abs(D))[1:] + (np.abs(P) + np.abs(D))[:-1])
    res_all = abs((E[-1] - E[0]) - src.sum()) / scale.sum()
    w = int(round(1.0 / DT)); cs, cr, cE = np.cumsum(np.r_[0, src]), np.cumsum(np.r_[0, scale]), E
    idx = np.arange(0, len(E) - w, w)
    res_win = float(np.max(np.abs((cE[idx + w] - cE[idx]) - (cs[idx + w] - cs[idx])) / (cr[idx + w] - cr[idx])))
    # S-G2: truncated advection is energy-neutral at random resolved states
    rng = np.random.default_rng([SEED, 7]); neutral = {}
    for K in KS:
        m = 2 * len(G.half_modes(g, K)[0]); worst = 0.0
        for _ in range(N_NEUTRAL):
            r, nn = G.advection_rate(g, K, rng.normal(size=m)); worst = max(worst, abs(r) / nn)
        neutral[str(K)] = worst
    # energy-coordinate consistency: |a|^2 / 2 equals the energy of the truncated field
    ec = {}
    for K in KS:
        wK = wh * G.trunc_mask(g, K); Ek, _, _ = g.diagnostics(wK); a = G.coords(g, K, wK)
        ec[str(K)] = abs(0.5 * a @ a - Ek) / max(Ek, 1e-300)
    rec = dict(N=N, nu=NU, dt=DT, seed=SEED, t_spin=T_SPIN, t_sample=T_SAMPLE, n_samples=len(times),
               budget_residual_whole=float(res_all), budget_residual_max_window=res_win,
               E_mean_sampling=float(E2.mean()), P_mean_sampling=float(P2.mean()), D_mean_sampling=float(D2.mean()),
               advection_neutrality_max_rel=neutral, energy_coord_rel_err=ec, seconds=round(time.time() - t0, 1))
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    json.dump(rec, open(os.path.join(ROOT, 'results', 'S1_run.json'), 'w'), indent=1)
    np.savez_compressed(os.path.join(ROOT, 'results', 'S1_samples.npz'), t=np.array(times),
                        **{f'A{K}': np.array(samples[K][0]) for K in KS}, **{f'B{K}': np.array(samples[K][1]) for K in KS},
                        E=E[::EVERY], P=P[::EVERY], D=D[::EVERY])
    print(json.dumps(rec))


if __name__ == '__main__':
    main()
