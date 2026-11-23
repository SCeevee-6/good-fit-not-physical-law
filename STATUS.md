# STATUS (2026-11-23, chat 3)

## Decided
- Gate A: protocol 0A2, outcome 3 (`docs/gate_A2_result.md`). Replicated over 20 seeds in 0B-3.
- Protocol 0B ran (`docs/gate_B_result.md`, ROADMAP_v2 Addendum 5): G-B1 pass, G-B2 pass, G-B3 fail, replication pass.

## Frozen
- 0A, 0A2 as before.
- 0B: `docs/protocol_0B.md` (commit 9f1cf045, before the run), `src/gate_b.py`, `src/needles.py`, `src/run_0B.py`, `results/0B_jobs/`, `results/0B_sample.csv`.

## Next
- Hare-lynx real data. The user approved network access; pick the source and record it before fitting.
- Law taxonomy per law and structure, including the H result that structure alone breaks the trapping law more often.
- Competition category follows the headline demo.

## Resume
- Kaggle (heavy runs): `python tools/kaggle_launch.py push <main|needles> <commit> <user>`, then `status` / `fetch`. The token's account is shreyanchethan. The runner pulls the given commit from GitHub, so the repo must be public or the code uploaded another way.
- Local: `python src/run_0B.py <pool|sample|arms|needles> [workers]` (resumable), `python src/analyze_0B.py` (CSVs only).
- Holling truth check: `python tools/check_holling.py`.
- This machine blocks `scipy.signal`; `src/savgol_np.py` stands in, and gives the same fits.
- No pandas in analysis code. Gate thresholds are compared as exact fractions.
