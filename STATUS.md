# STATUS (2026-11-26, chat 3)

## Decided
- Gate A: protocol 0A2, outcome 3 (`docs/gate_A2_result.md`). Replicated over 20 seeds in 0B-3.
- Protocol 0B ran (`docs/gate_B_result.md`, ROADMAP_v2 Addendum 5): G-B1 pass, G-B2 pass, G-B3 fail, replication pass.

- Protocol 0R ran (`docs/gate_R_result.md`, ROADMAP_v2 Addendum 6): G-R1, G-R2, G-R3 pass, real-data claim kept, with the limits listed there.

- Law taxonomy done (`docs/law_taxonomy.md`, `results/law_taxonomy.csv`, Addendum 7). The H result (structure lowers trapping-law satisfaction) is untested elsewhere.

- Protocol 0E ran (`docs/gate_E_result.md`, ROADMAP_v2 Addendum 8): G-E1 and G-E2 pass, G-E3 fails on the validation cost (median ratio 1.46). The cheap-fix claim is not kept on this series.

## Frozen
- 0A, 0A2 as before.
- 0B: `docs/protocol_0B.md` (commit 9f1cf045, before the run), `src/gate_b.py`, `src/needles.py`, `src/run_0B.py`, `results/0B_jobs/`, `results/0B_sample.csv`.

- 0R: `docs/protocol_0R.md` (frozen at 0ad1896c, addendum 1 for time limits), run at commit fcfb358d, `src/real_r.py`, `src/run_0R.py`, `src/analyze_0R.py`, `data_raw/`, `results/0R_jobs/`.

- 0E: `docs/protocol_0E.md` (frozen at 942d55e5), `src/real_e.py`, `src/run_0E.py`, `src/analyze_0E.py`, `data_raw/jhu/`, `results/0E_jobs/`.

## Next
- Epidemic series (chat 5): choose the source, write the protocol and freeze it before the first fit.
- Robustness checks and the board demo (chat 6). The category choice follows the headline demo.
- Competition category follows the headline demo.

## Resume
- 0E on Kaggle: `python tools/kaggle_launch.py push E <commit> <user> [test]`, then `status E <user>` and `fetch E <user> <dir>`. Unzip into the repo, then `python src/analyze_0E.py`. 30 configs take about 30 minutes on 4 CPUs.
- 0R on Kaggle: `python tools/kaggle_launch.py push R <commit> <user> [test]`, then `status R <user>` and `fetch R <user> <dir>`. Unzip into the repo, then `python src/analyze_0R.py`. Each config runs in its own process with a time limit.
- Kaggle (heavy runs): `python tools/kaggle_launch.py push <main|needles> <commit> <user>`, then `status` / `fetch`. The token's account is shreyanchethan. The runner pulls the given commit from GitHub, so the repo must be public or the code uploaded another way.
- Local: `python src/run_0B.py <pool|sample|arms|needles> [workers]` (resumable), `python src/analyze_0B.py` (CSVs only).
- Holling truth check: `python tools/check_holling.py`.
- This machine blocks `scipy.signal`; `src/savgol_np.py` stands in, and gives the same fits.
- No pandas in analysis code. Gate thresholds are compared as exact fractions.
