"""Merge results/cells/*.pkl into flat CSV files (cand, spec, repair, audit). Run in a fresh interpreter."""
import os, sys, glob, pickle, json
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, '..', 'results')
C, SP, R, A = [], [], [], []
files = sorted(glob.glob(os.path.join(RES, 'cells', '*.pkl')))
for f in files:
    d = pickle.load(open(f, 'rb'))
    c = d['cand'].copy(); c['xi'] = c['xi'].apply(json.dumps); c['xi_shape'] = c['xi_shape'].apply(lambda s: json.dumps(list(s)))
    c['cell_seconds'] = d['seconds']; C.append(c)
    SP.append(d['spec'])
    if len(d['repair']):
        r = d['repair'].copy(); r['xi'] = r['xi'].apply(lambda v: json.dumps(v) if isinstance(v, list) else ''); R.append(r)
    if len(d['audit']): A.append(d['audit'])
pd.concat(C, ignore_index=True).to_csv(os.path.join(RES, 'pilot_candidates.csv'), index=False)
pd.concat(SP, ignore_index=True).to_csv(os.path.join(RES, 'pilot_spec_queries.csv'), index=False)
pd.concat(R, ignore_index=True).to_csv(os.path.join(RES, 'pilot_repairs.csv'), index=False)
pd.concat(A, ignore_index=True).to_csv(os.path.join(RES, 'pilot_repair_audit.csv'), index=False)
print('combined', len(files), 'cells')
