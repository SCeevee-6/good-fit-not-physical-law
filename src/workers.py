"""Module-level worker functions for multiprocessing (spawn) on Windows."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np


def verify_model_record(args):
    """args = (tag, sysname, xi, shape, degree). Returns list of per-(spec, level) records incl. replay of strong counterexamples."""
    tag, sysname, xi, shape, d = args
    from systems import SYSTEMS
    from verify import verify_model, replay
    S = SYSTEMS[sysname]
    E = S.E(d)
    Xi = np.array(xi).reshape(shape)
    t0 = time.time()
    res = verify_model(S, Xi, E)
    recs = []
    for (spn, lv), r in res.items():
        sp_ = [s for s in S.specs if s.name == spn][0]
        rec = dict(tag=tag, spec=spn, level=lv, status=r['status'], method=r.get('method', 'z3'), time=r.get('time'), G=r.get('G'),
                   u=r.get('u'), witness_exact=r.get('witness_exact'), copied_from_strong=(lv == 'weak' and res[(spn, 'strong')]['status'] == 'counterexample'))
        if lv == 'strong' and r['status'] == 'counterexample' and r.get('u') is not None:
            rec['replay_ok'] = bool(replay(sp_, Xi, E, r['u'])[0])
        recs.append(rec)
    return recs


def oracle_record(args):
    """args = (tag, sysname, xi, shape, degree, seed). Local oracle + trajectory oracle + random falsification (n=200, 2000)."""
    tag, sysname, xi, shape, d, seed = args
    from systems import SYSTEMS
    import oracle
    S = SYSTEMS[sysname]
    E = S.E(d)
    Xi = np.array(xi).reshape(shape)
    t0 = time.time()
    lo = oracle.local_oracle(S, Xi, E, seed=seed)
    tj = oracle.trajectory_oracle(S, Xi, E, seed=seed)
    rng = np.random.default_rng([seed, 99])
    f200, s200 = oracle.random_falsify(S, Xi, E, 200, rng)
    f2000, s2000 = oracle.random_falsify(S, Xi, E, 2000, rng)
    return dict(tag=tag, local_min={k: v['min'] for k, v in lo.items()},
                local_viol_strong=any(v['viol_strong'] for v in lo.values()), local_viol_weak=any(v['viol_weak'] for v in lo.values()),
                local_viol_specs=[k for k, v in lo.items() if v['viol_strong']],
                traj=tj, traj_unphysical=tj['unphysical'], rf200=f200, rf200_score=s200, rf2000=f2000, rf2000_score=s2000, seconds=time.time() - t0)
