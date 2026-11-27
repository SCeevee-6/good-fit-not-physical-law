"""Summarise F0: results/F0_summary.csv (mean zone shares per system x degree x law) and results/F0_checks.txt."""
import os, csv
from collections import defaultdict
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); RES = os.path.join(ROOT, 'results')


def f(x):
    return float(x) if x not in ('', None) else np.nan


def main():
    Z = list(csv.DictReader(open(os.path.join(RES, 'F0_zones.csv'))))
    W = list(csv.DictReader(open(os.path.join(RES, 'F0_witness_zones.csv'))))
    grp = defaultdict(list)
    for r in Z:
        sysname = r['key'] if r['group'] == 'synthetic' else ('R' if r['key'] == 'R' else 'COVID')
        grp[(sysname, int(r['degree']), r['law'])].append(r)
    out = []
    for (s, d, law), rs in sorted(grp.items()):
        g = lambda k: float(np.nanmean([f(r[k]) for r in rs]))
        out.append(dict(system=s, degree=d, law=law, n_sets=len(rs), certified=round(g('certified'), 4),
                        refuted=round(g('refuted'), 4), undecided=round(g('undecided'), 4),
                        any_region_certified=sum(r['region_certified'] == 'True' for r in rs),
                        h_inf_share=round(g('h_inf_share'), 4), median_cond=float(np.median([f(r['cond']) for r in rs])),
                        n_aes=sum(int(r['n_aes']) for r in rs), n_aes_in_E=sum(int(r['n_aes_in_E']) for r in rs)))
    with open(os.path.join(RES, 'F0_summary.csv'), 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out[0].keys())); wr.writeheader(); wr.writerows(out)
    ratios = [f(r['max_member_ratio']) for r in Z]
    lines = [
        f"zone rows {len(Z)}, summary rows {len(out)}",
        f"max member ratio over all rows {np.nanmax(ratios):.4f} (must be <= 1)",
        f"contradictions {sum(int(r['contradictions']) for r in Z)} (must be 0)",
        f"stored equivalent models inside E {sum(int(r['n_aes_in_E']) for r in Z)} / {sum(int(r['n_aes']) for r in Z)} (law rows, counted once per law)",
        f"rows with near-data grid points {sum(int(r['n_near']) > 0 for r in Z)} / {len(Z)}",
        f"witnesses {len(W)}; of candidates in E {sum(w['cand_in_E'] == 'True' for w in W)}; zones of in-E witnesses: "
        + str({z: sum(1 for w in W if w['cand_in_E'] == 'True' and int(w['zone']) == z) for z in (1, 0, -1)}) + " (1 must be 0)",
        f"witness zones, all: " + str({z: sum(1 for w in W if int(w['zone']) == z) for z in (1, 0, -1)}),
    ]
    open(os.path.join(RES, 'F0_checks.txt'), 'w').write('\n'.join(lines) + '\n')
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
