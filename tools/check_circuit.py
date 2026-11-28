"""Worked example for the paper: a one-port whose terminal data cannot decide an energy law (docs/example_circuit.md).
Series R-L branch in parallel with series R-C branch, L = tau R, C = tau / R, starting at rest.

Checks (exact where possible):
 C1 (sympy) with x the capacitor voltage and y = R i_L, both obey tau z' = v - z, so y = x from rest and the
    port current is i = y/R + (v - x)/R = v/R for every drive v(t) and every tau.
 C2 (sympy) stored energy E = L i_L^2 / 2 + C x^2 / 2 = tau x^2 / R and dE/dt = 2 x (v - x) / R.
 C3 the 2 V -> 1 V step at t = 1: x(1) = 2 (1 - exp(-1/tau)); storage rate just after the step for tau = 0.2 and 3,
    and the threshold tau* = 1 / ln 2 where it changes sign. Cross-checked by integrating the two branches separately.
Writes results/check_circuit.txt and figures/fig_circuit.png.
Usage: python tools/check_circuit.py
"""
import os, sys
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def symbolic():
    t, R, tau = sp.symbols('t R tau', positive=True)
    v = sp.Function('v')(t); x = sp.Function('x')(t); y = sp.Function('y')(t)
    L, C = tau * R, tau / R
    # branch equations: L i_L' + R i_L = v (i_L = y / R);  C x' = (v - x) / R
    eq_y = sp.simplify(L * sp.diff(y / R, t) + R * (y / R) - v)       # = tau y' + y - v
    eq_x = sp.simplify(C * sp.diff(x, t) - (v - x) / R)                # = (tau x' + x - v) / R
    same_ode = sp.simplify(eq_y - R * eq_x.subs(x, y)) == 0
    i_port = sp.simplify((y / R + (v - x) / R).subs(y, x))
    E = sp.simplify(L * (x / R) ** 2 / 2 + C * x ** 2 / 2)
    dE = sp.simplify(sp.diff(E, t).subs(sp.diff(x, t), (v - x) / tau))
    return dict(same_ode=bool(same_ode), port_current_is_v_over_R=bool(sp.simplify(i_port - v / R) == 0),
                energy_is_tau_x2_over_R=bool(sp.simplify(E - tau * x ** 2 / R) == 0),
                rate_is_2x_v_minus_x_over_R=bool(sp.simplify(dE - 2 * x * (v - x) / R) == 0))


def drive(t):
    return 2.0 if t < 1.0 else 1.0


def simulate(tau, R=1.0, T=3.0):
    """Integrate capacitor voltage and inductor current separately (no use of y = x)."""
    L, C = tau * R, tau / R
    f = lambda t, s: [(drive(t) - s[0]) / (R * C), (drive(t) - R * s[1]) / L]
    s1 = solve_ivp(f, (0, 1), [0.0, 0.0], rtol=1e-11, atol=1e-13, dense_output=True)
    s2 = solve_ivp(f, (1, T), s1.y[:, -1], rtol=1e-11, atol=1e-13, dense_output=True)
    return s1, s2, L, C


def main():
    out = []; ok = True
    sym = symbolic(); out.append(f'C1/C2 symbolic: {sym}'); ok &= all(sym.values())
    rates = {}
    for tau in (0.2, 3.0):
        x1 = 2 * (1 - np.exp(-1 / tau)); rate = 2 * x1 * (1 - x1)
        s1, s2, L, C = simulate(tau)
        xc, il = s1.y[:, -1]
        rate_num = (1.0 - xc) * xc / 1.0 + (1.0 - il) * il * 1.0            # C x x' + L i i' at t = 1+
        port = il + (1.0 - xc)                                               # i_L + i_C with R = 1
        rates[tau] = rate
        out.append(f'C3 tau={tau}: x(1)={x1:.6f}, storage rate at 1+ closed form {rate:.4f} W, separate-branch integration {rate_num:.4f} W, port current {port:.10f} A, port power {port * 1.0:.10f} W')
        ok &= abs(rate - rate_num) < 1e-7 and abs(port - 1.0) < 1e-8
    tstar = 1 / np.log(2)
    out.append(f'C3 sign change at tau* = 1/ln 2 = {tstar:.6f}: tau < tau* releases energy, tau > tau* stores it')
    ok &= rates[0.2] < 0 < rates[3.0]
    out.append(f'ALL_OK {ok}')
    os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
    open(os.path.join(ROOT, 'results', 'check_circuit.txt'), 'w').write('\n'.join(out) + '\n')
    figure()
    print('ALL_OK', ok)
    return ok


def figure():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'legend.fontsize': 7, 'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'DejaVu Sans'})
    cols = {0.2: '#1f5fa8', 3.0: '#d08a1c'}
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3), constrained_layout=True)
    tt = np.linspace(0, 3, 1201)
    for tau in (0.2, 3.0):
        s1, s2, L, C = simulate(tau)
        S = np.where(tt < 1, 0, 1)
        st = np.array([s1.sol(t) if t < 1 else s2.sol(t) for t in tt])
        xc, il = st[:, 0], st[:, 1]; v = np.where(tt < 1, 2.0, 1.0)
        ax[0].plot(tt, il + (v - xc), color=cols[tau], lw=1.6 if tau == 0.2 else 1.0, ls='-' if tau == 0.2 else '--', label=f'τ = {tau:g}')
        ax[1].plot(tt, 2 * xc * (v - xc), color=cols[tau], lw=1.2, label=f'τ = {tau:g}')
    ax[0].set_title('Port current is identical', loc='left'); ax[0].set_xlabel('time'); ax[0].set_ylabel('port current (A)')
    ax[0].legend(frameon=False, loc='upper right')
    ax[1].axhline(0, color='0.6', lw=0.6); ax[1].axvline(1, color='0.8', lw=0.6, ls=':')
    ax[1].set_title('Internal storage rate is not', loc='left'); ax[1].set_xlabel('time'); ax[1].set_ylabel('dE/dt (W)')
    for tau, yy in ((0.2, -3.9195), (3.0, 0.4910)):
        ax[1].annotate(f'{yy:+.2f} W', xy=(1.0, yy), xytext=(1.35, yy * 0.55 + (0.6 if tau == 3.0 else -0.4)), color=cols[tau], fontsize=7,
                       arrowprops=dict(arrowstyle='-', color=cols[tau], lw=0.6))
    taus = np.geomspace(0.05, 20, 400); x1 = 2 * (1 - np.exp(-1 / taus)); r = 2 * x1 * (1 - x1)
    ax[2].semilogx(taus, r, color='0.25', lw=1.2); ax[2].axhline(0, color='0.6', lw=0.6)
    ax[2].axvline(1 / np.log(2), color='0.5', lw=0.8, ls='--')
    ax[2].text(1 / np.log(2) * 1.15, -2.6, 'τ* = 1/ln 2', fontsize=7, color='0.3')
    ax[2].set_title('Same port data, either sign', loc='left'); ax[2].set_xlabel('time constant τ'); ax[2].set_ylabel('dE/dt at t = 1⁺ (W)')
    for a, l in zip(ax, 'abc'):
        a.text(-0.2, 1.08, l, transform=a.transAxes, fontweight='bold', fontsize=9, va='bottom'); a.margins(0.04)
    os.makedirs(os.path.join(ROOT, 'figures'), exist_ok=True)
    fig.savefig(os.path.join(ROOT, 'figures', 'fig_circuit.png'), dpi=300)


if __name__ == '__main__':
    sys.exit(0 if main() else 1)
