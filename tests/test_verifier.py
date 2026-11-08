"""Verifier checks on hand-built models. Run: python tests/test_verifier.py   (or pytest)."""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src'))
import numpy as np
from systems import SYSTEMS
from verify import verify_model, replay


def _model(key, d, mods):
    S = SYSTEMS[key]
    E = S.E(d)
    Xi = S.true_Xi(d).copy()
    idx = [tuple(e) for e in E.tolist()]
    for (i, e), delta in mods.items():
        Xi[idx.index(tuple(e)), i] += delta
    return S, E, np.round(Xi, 6)


def _status(S, E, Xi):
    res = verify_model(S, Xi, E)
    return {k: v['status'] for k, v in res.items()}, res


def test_true_models_pass():
    for key, d in (('A', 3), ('B', 3)):
        S, E, Xi = _model(key, d, {})
        st, _ = _status(S, E, Xi)
        assert all(v == 'verified' for v in st.values()), (key, st)


def test_predator_leak_strong():
    S, E, Xi = _model('A', 3, {(0, (0, 1)): -0.05})
    st, res = _status(S, E, Xi)
    assert st[('A1x', 'strong')] == 'counterexample'
    assert st[('A1y', 'strong')] == 'verified' and st[('A2', 'strong')] == 'verified'


def test_weak_only_violation():
    S, E, Xi = _model('A', 3, {(0, (0, 1)): -0.002})
    st, _ = _status(S, E, Xi)
    assert st[('A1x', 'strong')] == 'verified' and st[('A1x', 'weak')] == 'counterexample'


def test_cubic_escape():
    S, E, Xi = _model('A', 3, {(0, (3, 0)): 0.05})
    st, _ = _status(S, E, Xi)
    assert st[('A2', 'strong')] == 'counterexample'


def test_sir_mass_leak_and_replay():
    S, E, Xi = _model('B', 3, {(2, (0, 0, 1)): 0.05})
    st, res = _status(S, E, Xi)
    assert st[('B3p', 'strong')] == 'counterexample' and st[('B2S', 'strong')] == 'verified'
    sp = [s for s in S.specs if s.name == 'B3p'][0]
    r = res[('B3p', 'strong')]
    assert r['witness_exact'] and replay(sp, Xi, E, r['u'])[0]


def test_allowed_constant_gain():
    S, E, Xi = _model('A', 3, {(1, (0, 0)): 0.05})
    st, _ = _status(S, E, Xi)
    assert all(v == 'verified' for k, v in st.items() if k[1] == 'strong'), st


if __name__ == '__main__':
    fns = [v for k, v in sorted(globals().items()) if k.startswith('test_')]
    for f in fns:
        f()
        print('ok', f.__name__)
    print(len(fns), 'passed')
