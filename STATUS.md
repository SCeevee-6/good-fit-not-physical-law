# STATUS (2026-11-26)

## Decided
- Gate A: protocol 0A2, outcome 3 (`docs/gate_A2_result.md`). Replicated over 20 seeds in 0B-3.
- Protocol 0B ran (`docs/gate_B_result.md`, ROADMAP_v2 Addendum 5): G-B1 pass, G-B2 pass, G-B3 fail, replication pass.

- Protocol 0R ran (`docs/gate_R_result.md`, ROADMAP_v2 Addendum 6): G-R1, G-R2, G-R3 pass, real-data claim kept, with the limits listed there.

- Law taxonomy done (`docs/law_taxonomy.md`, `results/law_taxonomy.csv`, Addendum 7). The H result (structure lowers trapping-law satisfaction) is untested elsewhere.

- Protocol 0E ran (`docs/gate_E_result.md`, ROADMAP_v2 Addendum 8): G-E1 and G-E2 pass, G-E3 fails on the validation cost (median ratio 1.46). The cheap-fix claim is not kept on this series.

- Protocol 0C ran (`docs/gate_C_result.md`, ROADMAP_v2 Addendum 9): soft penalty matches the loop on the synthetic testbeds (R-C1: loop not needed), Gate A outcome 3 holds at 0.5 % and 5 % noise, independent check finds 0 of 2871 accepted models violating. Repair costs accuracy at 0.5 % noise for a wrong library (median ratio 1.4-1.5).

- Protocol 0D ran (`docs/gate_D_result.md`, ROADMAP_v2 Addendum 10): at four states a fine dense grid verifies 80/80, the loop 69/80 (budget stops), pilot-budget grid 64/80, random loop 51/80. R-D1 (loop coverage advantage) not supported; C - Rnd 0.225.

## Frozen
- 0A, 0A2 as before.
- 0B: `docs/protocol_0B.md` (commit 9f1cf045, before the run), `src/gate_b.py`, `src/needles.py`, `src/run_0B.py`, `results/0B_jobs/`, `results/0B_sample.csv`.

- 0R: `docs/protocol_0R.md` (frozen at 0ad1896c, addendum 1 for time limits), run at commit fcfb358d, `src/real_r.py`, `src/run_0R.py`, `src/analyze_0R.py`, `data_raw/`, `results/0R_jobs/`.

- 0E: `docs/protocol_0E.md` (frozen at 942d55e5), `src/real_e.py`, `src/run_0E.py`, `src/analyze_0E.py`, `data_raw/jhu/`, `results/0E_jobs/`.

- 0C: `docs/protocol_0C.md` (frozen at 2109feba), `src/robust.py`, `src/run_0C.py`, `src/analyze_0C.py`, `results/0C_jobs/`, `results/0C_sample.csv`.

- 0D: `docs/protocol_0D.md` (frozen at 3cbcb2e9, run at 3641866a), `src/bern_nd.py`, `src/dimtest.py`, `src/run_0D.py`, `src/analyze_0D.py`, `results/0D_jobs/`.

## Next
- Direction changed to `docs/ROADMAP_v3.md`: law zones over the data-consistent set, region certificate, law horizon, measurement design. v2 results are the motivation and controls.
- Identities used by v3 are checked in `tools/check_lawmap_math.py` (`results/check_lawmap_math.txt`, ALL_OK True).
- F0 done (`docs/F0_result.md`, ROADMAP_v3 Addendum 1): zone maps for all stored fits; checks pass (member ratio <= 0.706, 0 contradictions).
- 0F done (`docs/gate_F_result.md`, ROADMAP_v3 Addendum 2): G-F0 pass, G-F1 pass on A and B, fail on E4, G-F2 hare-lynx held, COVID did not, G-F3 pass, G-F4 fail.
- Next: T (proofs for item 7, `docs/law_horizon.md`), optional bias-corrected set for E4 as a new frozen unit, then the write-up. Finish target 2027-01-15, freeze 2027-01-03.
- README rewritten as a summary of every unit and its outcome.

## Resume
- Law-map checks: `python tools/check_lawmap_math.py` (about 20 s).
- 0F on Kaggle: `python tools/kaggle_launch.py push F|F3|F4 <commit> <user>`, then `status` and `fetch` with the same mode. Unzip 0F.zip into the repo, then `python src/analyze_0F.py`. Times: calib + real 3.1 h, horizon 2.6 h, design 9.9 h.
- F0: `python src/run_F0.py 5` (about 13 min locally), `python src/run_F0_plane.py`, `python src/analyze_F0.py`.
- 0D on Kaggle: `python tools/kaggle_launch.py push D <commit> <user> [test]`, then `status D <user>` and `fetch D <user> <dir>`. Unzip into the repo, then `python src/analyze_0D.py`. The full run takes about 15 minutes on 4 CPUs.
- 0C on Kaggle: `python tools/kaggle_launch.py push C <commit> <user> [test]`, then `status C <user>` and `fetch C <user> <dir>`. Unzip into the repo, then `python src/analyze_0C.py`. The full run takes about 20 minutes on 4 CPUs.
- 0E on Kaggle: `python tools/kaggle_launch.py push E <commit> <user> [test]`, then `status E <user>` and `fetch E <user> <dir>`. Unzip into the repo, then `python src/analyze_0E.py`. 30 configs take about 30 minutes on 4 CPUs.
- 0R on Kaggle: `python tools/kaggle_launch.py push R <commit> <user> [test]`, then `status R <user>` and `fetch R <user> <dir>`. Unzip into the repo, then `python src/analyze_0R.py`. Each config runs in its own process with a time limit.
- Kaggle (heavy runs): `python tools/kaggle_launch.py push <main|needles> <commit> <user>`, then `status` / `fetch`. The token's account is shreyanchethan. The runner pulls the given commit from GitHub, so the repo must be public or the code uploaded another way.
- Local: `python src/run_0B.py <pool|sample|arms|needles> [workers]` (resumable), `python src/analyze_0B.py` (CSVs only).
- Holling truth check: `python tools/check_holling.py`.
- This machine blocks `scipy.signal`; `src/savgol_np.py` stands in, and gives the same fits.
- No pandas in analysis code. Gate thresholds are compared as exact fractions.
