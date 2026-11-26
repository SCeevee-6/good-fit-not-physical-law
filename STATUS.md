# STATUS (2026-11-26, chat 3)

## Decided
- Gate A: protocol 0A2, outcome 3 (`docs/gate_A2_result.md`). Replicated over 20 seeds in 0B-3.
- Protocol 0B ran (`docs/gate_B_result.md`, ROADMAP_v2 Addendum 5): G-B1 pass, G-B2 pass, G-B3 fail, replication pass.

- Protocol 0R ran (`docs/gate_R_result.md`, ROADMAP_v2 Addendum 6): G-R1, G-R2, G-R3 pass, real-data claim kept, with the limits listed there.

## Frozen
- 0A, 0A2 as before.
- 0B: `docs/protocol_0B.md` (commit 9f1cf045, before the run), `src/gate_b.py`, `src/needles.py`, `src/run_0B.py`, `results/0B_jobs/`, `results/0B_sample.csv`.

- 0R: `docs/protocol_0R.md` (frozen at 0ad1896c, addendum 1 for time limits), run at commit fcfb358d, `src/real_r.py`, `src/run_0R.py`, `src/analyze_0R.py`, `data_raw/`, `results/0R_jobs/`.

## Next
- Epidemic series (chat 5): choose the source, write the protocol and freeze it before the first fit.
- Law taxonomy per law and structure, including the H result that structure alone breaks the trapping law more often.
- Competition category follows the headline demo.

## Resume
- 0R on Kaggle: `python tools/kaggle_launch.py push R <commit> <user> [test]`, then `status R <user>` and `fetch R <user> <dir>`. Unzip into the repo, then `python src/analyze_0R.py`. Each config runs in its own process with a time limit.
- Kaggle (heavy runs): `python tools/kaggle_launch.py push <main|needles> <commit> <user>`, then `status` / `fetch`. The token's account is shreyanchethan. The runner pulls the given commit from GitHub, so the repo must be public or the code uploaded another way.
- Local: `python src/run_0B.py <pool|sample|arms|needles> [workers]` (resumable), `python src/analyze_0B.py` (CSVs only).
- Holling truth check: `python tools/check_holling.py`.
- This machine blocks `scipy.signal`; `src/savgol_np.py` stands in, and gives the same fits.
- No pandas in analysis code. Gate thresholds are compared as exact fractions.
