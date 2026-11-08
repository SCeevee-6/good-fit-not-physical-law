# STATUS (2026-11-08, end of chat 1)

## Passed / decided
- Gate A (protocol 0A): outcome 3. G-A1 pass, G-A2 fail. See `docs/gate_A_result.md` and `results/0A_gate.csv`.
- Roadmap: `docs/ROADMAP_v2.md`, with a dated addendum recording the Gate A decision.

## Frozen
- `docs/protocol_0A.md` (with the B1S clarification and the pre-rerun scoring note), `results/0A_sample.csv`, `src/gate_a.py`.
- Raw 0A outputs: `results/0A_jobs/`.

## Next
- Chat 2: write protocol 0B (Holling II misspecification, needle test with random-loop sample counts, 20 seeds with CIs) and commit it before the first run.
- Side track S1: `docs/protocol_S1.md`, frozen, not run.
- Open: competition category, hare-lynx data source, affiliated fair and deadline.

## Resume
- Rerun 0A (about 1 min on 6 workers): in the ceg-ode kernel, `sys.path.insert(0, 'src'); import run_gate_a; run_gate_a.main(6)`. It skips candidates whose JSON exists; delete `results/0A_jobs/` for a clean rerun.
- Tables: `import analyze_gate_a; analyze_gate_a.main()`.
- Environment: no pandas in ceg-ode; use csv and numpy.
