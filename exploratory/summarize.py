"""Summarize independent trial blocks, preserving the paired design."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

root = Path(__file__).resolve().parent
groups = defaultdict(list)
with (root/'trials.csv').open() as f:
    for row in csv.DictReader(f):
        group = 'grid' if int(row['trial']) < 120 else 'off_grid'
        groups[row['coupling'],row['rho'],group,row['policy']].append(row)

with (root/'summary.csv').open('w',newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['coupling','rho','truth_set','policy','n','ambiguous','mean_accuracy',
                     'gain_vs_fixed_repeat','paired_standard_error','repeat_selected'])
    for (k,r,g,p),rows in groups.items():
        ref = groups[k,r,g,'fixed_repeat']
        acc = [float(x['true_accuracy']) for x in rows]
        diff = [a-float(b['true_accuracy']) for a,b in zip(acc,ref)]
        ambiguous = [x for x in rows if min(int(x['npos']),int(x['nneg']))>0]
        writer.writerow([k,r,g,p,len(rows),len(ambiguous),round(mean(acc),6),
                         round(mean(diff),6),round(stdev(diff)/len(rows)**.5,6),
                         sum(x['arm']=='repeat' for x in ambiguous)])
