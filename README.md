# good-fit-not-physical-law

Two questions. If an equation found from data predicts held-out data as well as the best one, can it still break a physical law? And if it does, can we catch it exactly and fix it from the counterexamples?

Pilot: 45 cells, 5,553 candidate ODEs (predator-prey, SIR, and predator-prey with a term missing from the library).

- 42-82% of the equally accurate candidates broke a stated law (positivity, trapping region, S' <= 0, conservation). 93-100% when the library lacks a term.
- z3 plus a Bernstein fallback decided every candidate. 0 disagreements with two independent checks (14,212 comparisons).
- Feeding counterexamples back into a constrained refit repaired 99.25% with no median accuracy cost. Random constraints, same number of rows: 17%.
- It does not beat AIC at ranking violators (AUROC 0.966) and 200 random test points already catch 97.6% of violations. Pilot verdict was REVISE, not GO. Details in `docs/pilot_report.md`.

Next steps are in `docs/ROADMAP.md`.

## layout

| path | what |
|---|---|
| `src/` | systems and laws (`systems.py`), STLSQ pool (`pool.py`), verifier (`verify.py`, `bern.py`), oracle (`oracle.py`), repair (`ceg.py`), one grid cell (`run_cell.py`), tables (`analyze.py`), figures (`figures.py`) |
| `docs/` | protocol with the addenda, pilot report, roadmap |
| `results/` | per-candidate CSVs and summary tables |
| `figures/` | fig0-fig5 |
| `skills/` | the four skill write-ups I made from this (verification, repair, Windows sandbox notes, prereg workflow) |
| `tests/` | verifier checks on hand-built good and bad models |
| `tools/gh_push.py` | pushes a folder to GitHub through the API |

## run

```
conda env create -f environment.yml
conda activate ceg-ode
python tests/test_verifier.py
python src/run_grid.py 7 B:1:0 A:1:0 Amis:1:0     # key:noise_idx:seed, workers first
python src/combine.py
```

`run_grid.py` writes one pickle per cell in `results/cells/` and skips cells that exist. Noise index 0, 1, 2 is 0.5%, 2%, 5%.

Only tested on Windows. joblib, scikit-learn and Parquet were blocked in the sandbox this ran in, so nothing here uses them.

## AI use

See `AI_USE.md`.
