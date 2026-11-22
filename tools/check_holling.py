"""Truth checks for system H (Rosenzweig-MacArthur, Holling II) before protocol 0B is frozen.

(a) z3: each law, multiplied by (1 + h x) > 0, holds for the true model. Solve times vs the 3 s limit.
(b) coexistence point, Hopf value of K, Jacobian at K, RK4 rollouts at K and at 0.9 K_hopf (numpy only).
(c) trapping constant: exact max of V' on {V = c, x, y >= 0}, smallest c that traps.

Writes results/check_holling.json, prints one line.
usage: python tools/check_holling.py [--a 1 --h 1 --b 0.8 --d 0.4 --K 4 --c 5 --zt 3000]
"""
import argparse, json, os, time
from fractions import Fraction as F
import numpy as np
import z3

NAMES = ('a', 'h', 'b', 'd', 'K', 'c')


def zq(f):
    return z3.RealVal(f"{f.numerator}/{f.denominator}")


def z3_check(cons, ms):
    s = z3.Solver(); s.set('timeout', ms); s.add(*cons)
    t0 = time.perf_counter(); r = s.check()
    return str(r), time.perf_counter() - t0


def rhs(u, a, h, b, d, K):
    x, y = u
    fr = x * y / (1 + h * x)
    return np.array([x * (1 - x / K) - a * fr, -d * y + b * fr])


def rk4(u0, p, dt, n):
    u = np.array(u0, float); out = np.empty((n, 2))
    for i in range(n):
        k1 = rhs(u, *p); k2 = rhs(u + dt / 2 * k1, *p)
        k3 = rhs(u + dt / 2 * k2, *p); k4 = rhs(u + dt * k3, *p)
        u = u + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); out[i] = u
    return out


def main():
    ap = argparse.ArgumentParser()
    for k, v in zip(NAMES, ('1', '1', '0.8', '0.4', '4', '5')):
        ap.add_argument('--' + k, default=v)
    ap.add_argument('--zt', type=int, default=3000)
    ap.add_argument('--out', default='results/check_holling.json')
    A = ap.parse_args()
    a, h, b, d, K, c = (F(getattr(A, k)) for k in NAMES)
    res = {'params': {k: str(v) for k, v in zip(NAMES, (a, h, b, d, K, c))}}

    # (b) exact equilibrium, Hopf, Jacobian
    assert b - d * h > 0
    xs = d / (b - d * h)
    ys = (1 + h * xs) * (1 - xs / K) / a
    KH = (b + d * h) / (h * (b - d * h))
    J11 = xs * (-1 / K + a * h * ys / (1 + h * xs) ** 2)
    J12 = -a * xs / (1 + h * xs)
    J21 = b * ys / (1 + h * xs) ** 2
    tr, det = J11, -J12 * J21
    p = tuple(float(v) for v in (a, h, b, d, K))
    tail = rk4((float(xs) * 1.05, float(ys)), p, 0.01, 200000)[100000:]
    pb = p[:4] + (0.9 * float(KH),)
    xb = float(xs); yb = float((1 + h * xs) * (1 - xs / F(pb[4]).limit_denominator(10**9)) / a)
    tb = rk4((xb * 1.05, yb), pb, 0.01, 200000)[100000:]
    res['b'] = dict(x_star=str(xs), y_star=str(ys), K_hopf=str(KH), K_hopf_eq_2xs_plus_1_over_h=bool(KH == 2 * xs + 1 / h),
                    trace=str(tr), det=str(det), unstable_focus=bool(tr > 0 and det > 0),
                    rk4_tail_x=[float(tail[:, 0].min()), float(tail[:, 0].max())],
                    rk4_tail_y=[float(tail[:, 1].min()), float(tail[:, 1].max())],
                    limit_cycle=bool(np.ptp(tail[:, 0]) > 0.05),
                    K_below=pb[4], below_converges=bool(np.ptp(tb[:, 0]) < 1e-3),
                    traj_inside_trap=bool(float((tail[:, 0] * float(b / a) + tail[:, 1]).max()) < float(c)))

    # (c) exact trapping: V = (b/a) x + y, V' = (b/a) x (1 - x/K) - d y  (rational terms cancel)
    r = b / a
    q2, q1, q0 = -r / K, r + d * r, -d * c           # V' on the segment as a function of x
    xr = c / r; xm = -q1 / (2 * q2)
    pts = [F(0), xr] + ([xm] if 0 <= xm <= xr else [])
    gmax = max(q2 * t * t + q1 * t + q0 for t in pts)
    xm_c = min(max(xm, F(0)), xr)
    c_min = (q2 * xm_c ** 2 + q1 * xm_c) / d          # valid when xm <= c/r
    res['c'] = dict(Vdot_max=str(gmax), Vdot_max_float=float(gmax), argmax_x=str(xm_c), c_min_trap=str(c_min),
                    argmax_inside=bool(0 <= xm <= xr), traps=bool(gmax < 0))

    # (a) z3 on polynomial forms: (1 + h x) * law
    x, y = z3.Reals('x y')
    zA, zH, zB, zD, zK, zC = (zq(v) for v in (a, h, b, d, K, c))
    P = x * (1 - x / zK) * (1 + zH * x) - zA * x * y          # (1+hx) x'
    Q = -zD * y * (1 + zH * x) + zB * x * y                  # (1+hx) y'
    Vd = (zB / zA) * P + Q                                    # (1+hx) V'
    tests = {
        'H1x': [x == 0, y >= 0, y <= zC, P < 0],
        'H1y': [y == 0, x >= 0, x <= zC * zA / zB, Q < 0],
        'H2': [x >= 0, y >= 0, (zB / zA) * x + y == zC, Vd > 0],
        'H2_margin': [x >= 0, y >= 0, (zB / zA) * x + y == zC, Vd > zq(gmax) * (1 + zH * x)],
        'H2_tight': [x >= 0, y >= 0, (zB / zA) * x + y == zC, Vd >= zq(gmax) * (1 + zH * x)],
    }
    res['a'] = {}
    for k, cons in tests.items():
        st, t = z3_check(cons, A.zt)
        res['a'][k] = dict(z3=st, sec=round(t, 4))
    res['a']['within_limit'] = all(v['sec'] * 1000 < A.zt for v in res['a'].values())
    ok = (all(res['a'][k]['z3'] == 'unsat' for k in ('H1x', 'H1y', 'H2', 'H2_margin'))
          and res['a']['H2_tight']['z3'] == 'sat' and res['b']['unstable_focus'] and res['b']['limit_cycle']
          and res['b']['below_converges'] and res['c']['traps'])
    res['all_ok'] = bool(ok)
    os.makedirs(os.path.dirname(A.out) or '.', exist_ok=True)
    with open(A.out, 'w') as fh:
        json.dump(res, fh, indent=1)
    print(f"check_holling all_ok={ok} x*={xs} y*={ys} K_hopf={KH} Vdot_max={gmax} c_min={c_min} -> {A.out}")


if __name__ == '__main__':
    main()
