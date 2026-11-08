---
name: prereg-pilot-workflow
description: "Run a pilot with hypotheses, thresholds and a decision rule written down before results, then log every later change as a dated addendum. Use when a pilot must tell you whether a project is worth scaling, and when peeking at a preview happens."
---

# prereg-pilot-workflow

1. Write `docs/pilot_protocol.md` first: question, hypotheses, systems, specs, metrics, pass thresholds (G1..Gn), decision rule (GO / REVISE / NO-GO). Save it as an artifact and commit it.
2. Smoke-test one cell end to end before the grid. A single-cell preview is allowed, but anything it changes is an addendum written before the full run, with what the preview showed.
3. Never edit a headline criterion after seeing data. Add a supplementary version next to it (G2b stayed as approved; G2b' added).
4. Name privileged baselines (anything using data a real study would not have) and report with and without.
5. Add sensitivity subsets before the run (parsimonious, strict accuracy set, top pick per cell).
6. After the run append one more addendum for every deviation (new fallback, amended rule, format change) with the reason.
7. Report failed criteria in the same table as passed ones. In the pilot G2b failed, and that changed the claim.
8. Label construction: if the method under test is part of the label, say so and give the independent check.
9. Keep per-cell results resumable and combine to CSV; keep `analyze.py` the only place numbers come from.

Template files: `docs/pilot_protocol.md`, `src/analyze.py`.
