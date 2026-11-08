# STATUS (2026-11-08, end of chat 1)

## Decided
- Gate A (protocol 0A): outcome 3 on record (G-A1 pass, G-A2 fail). Confounded on B: the B structure masked the true term R' = 0.2 I. See `docs/gate_A_result.md` and Addendum 3 in `docs/ROADMAP_v2.md`.

## Frozen
- `docs/protocol_0A.md`, `results/0A_sample.csv`, `src/gate_a.py`, raw outputs `results/0A_jobs/`.

## Next (chat 2, before 0B)
1. Write `docs/protocol_0A2.md` from Addendum 3 (B-pos and B-fac structures, same rule) and commit it.
2. Before running, check both B structures contain the true model (`tools/check_mask.py` pattern).
3. Run B only (20 candidates x 2 structures). The full 0A took about 1 min on 6 workers.
4. Re-decide Gate A, then write protocol 0B (Holling II, needle test with 2k/20k/200k sampler budgets, 20 seeds).

## Resume
- In the ceg-ode kernel: `sys.path.insert(0, 'src'); import run_gate_a; run_gate_a.main(6)` (skips finished JSONs). Tables: `import analyze_gate_a; analyze_gate_a.main()`.
- No pandas in ceg-ode; use csv and numpy.
