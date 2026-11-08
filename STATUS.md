# STATUS (2026-11-08, chat 2)

## Decided
- Gate A re-decided under protocol 0A2: outcome 3 with B-pos (f(C)-f(S) 0.317, med ratio C 1.000, f(C)-f(Rnd) 0.017); B-fac also outcome 3 (0.200 exactly, 0.000). `docs/gate_A2_result.md`, ROADMAP_v2 Addendum 4.
- 0A result kept on record, superseded for B by 0A2.

## Frozen
- 0A: `docs/protocol_0A.md`, `results/0A_sample.csv`, `src/gate_a.py`, `results/0A_jobs/`.
- 0A2: `docs/protocol_0A2.md` (commit f996151e, before the run), `src/gate_a2.py`, `results/0A2_jobs/`, `results/0A2_*.csv`.

## Next
- Protocol 0B drafted: `docs/protocol_0B.md` (DRAFT, not frozen). Before freezing: write `tools/check_holling.py` (truth check of H laws, trapping c, limit cycle at K=4), settle the four open questions, then commit the frozen version. No 0B run before that.

## Resume
- ceg-ode kernel, cwd repo root: `sys.path.insert(0, 'src'); import run_gate_a; run_gate_a.main(6, protocol='0A2')` (skips finished JSONs; 31 s for 40 jobs). Tables: `import analyze_gate_a2; analyze_gate_a2.main()`.
- Mask check: `runpy.run_path('tools/check_mask.py')`.
- No pandas in ceg-ode; use csv and numpy. Gate thresholds are compared as exact fractions.
