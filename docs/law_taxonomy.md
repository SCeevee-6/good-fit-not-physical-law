# Law taxonomy

One row per system and law. It tabulates results that already exist (pilot pools, 0B, 0R). Nothing was rerun, and no gate is involved. Script: `tools/law_table.py`, table: `results/law_taxonomy.csv`.

Columns: how many accuracy-equivalent pool models violate the law at strong tolerance; how many of the sampled models verify it after structure only (S), the witness loop (C), and the random-counterexample loop (Rnd). Class rule, in exact fractions: pool violation share <= 1/20 is "rarely violated"; else S verifies all is "enforced by structure"; else C verifies >= 9/10 is "needs the counterexample loop".

| system | law | family | pool violating | S | C | Rnd | class |
|---|---|---|---|---|---|---|---|
| A | A1x | structural | 288/571 (0.50) | 80/80 | 80/80 | 80/80 | enforced by structure |
| A | A1y | structural | 245/571 (0.43) | 80/80 | 80/80 | 80/80 | enforced by structure |
| A | A2 | non-structural | 90/571 (0.16) | 74/80 | 80/80 | 80/80 | needs the counterexample loop |
| Amis | A1x | structural | 353/444 (0.80) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Amis | A1y | structural | 116/444 (0.26) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Amis | A2 | non-structural | 191/444 (0.43) | 42/80 | 80/80 | 80/80 | needs the counterexample loop |
| Bpos | B1S | non-structural | 147/472 (0.31) | 53/80 | 80/80 | 79/80 | needs the counterexample loop |
| Bpos | B1I | structural | 5/472 (0.01) | 80/80 | 80/80 | 80/80 | rarely violated, no fix needed |
| Bpos | B1R | structural | 141/472 (0.30) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bpos | B2S | non-structural | 136/472 (0.29) | 61/80 | 80/80 | 72/80 | needs the counterexample loop |
| Bpos | B2R | non-structural | 155/472 (0.33) | 72/80 | 80/80 | 80/80 | needs the counterexample loop |
| Bpos | B3p | structural | 282/472 (0.60) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bpos | B3m | structural | 236/472 (0.50) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B1S | non-structural | 147/472 (0.31) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B1I | structural | 5/472 (0.01) | 80/80 | 80/80 | 80/80 | rarely violated, no fix needed |
| Bfac | B1R | structural | 141/472 (0.30) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B2S | non-structural | 136/472 (0.29) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B2R | non-structural | 155/472 (0.33) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B3p | structural | 282/472 (0.60) | 80/80 | 80/80 | 80/80 | enforced by structure |
| Bfac | B3m | structural | 236/472 (0.50) | 80/80 | 80/80 | 80/80 | enforced by structure |
| H | H1x | structural | 3338/3456 (0.97) | 80/80 | 80/80 | 80/80 | enforced by structure |
| H | H1y | structural | 2361/3456 (0.68) | 80/80 | 80/80 | 80/80 | enforced by structure |
| H | H2 | non-structural | 461/3456 (0.13) | 28/80 | 80/80 | 80/80 | needs the counterexample loop |
| R | R1x | structural | 38/47 (0.81) | 47/47 | 47/47 | 47/47 | enforced by structure |
| R | R1y | structural | 37/47 (0.79) | 47/47 | 47/47 | 47/47 | enforced by structure |
| R | R2 | non-structural | 42/47 (0.89) | 5/47 | 47/47 | 45/47 | needs the counterexample loop |

Sources and sizes. A, A-mis, B pools are the pilot pools at noise 2 % (all seeds), H is the 0B pool (40 cells), R is the primary 0R config. Arm columns are the 0B sample (80 candidates per system, 2 % noise, 20 seeds), for R all 47 equivalent candidates. B-pos and B-fac use the same 80 B candidates. The pool shares count every equivalent model, the arm columns only the sample, so they differ by sampling.

## What the table shows

1. Positivity laws are broken by the unconstrained fit very often: the hare positivity law in 81 % of equivalent hare-lynx models, A1x in 50 % of A models, H1x in 97 % of H models. Writing x' = x h(x, y) and y' = y g(x, y) removes all of them, at no cost in the sample (S verifies 80/80 and 47/47).
2. The conservation laws of the SIR model (B3p, B3m) and the removed-never-decreases law B1R are also enforced by the structure, because the S-equation is defined as minus the sum of the other two.
3. The non-structural laws are the ones that need a check: A2, B1S, B2S, B2R, H2, R2. The witness loop verifies every one of them in every sampled model. The random loop is behind on B2S (72/80 under B-pos, so 8 models are still inadmissible), on B1S (79/80) and on R2 (45/47).
4. B-fac, the structure where the S-equation is divisible by S, enforces all seven SIR laws on all 80 sampled models, including the three that B-pos leaves to the loop. For this system one structural choice replaces the loop, but only if it is the right one (B-pos is not enough).
5. Structure alone can make a non-structural law fail more often. On H the share of pool models that satisfy the trapping law H2 is 0.87, and after structure only 28/80 = 0.35 (95 % Wilson 0.25 to 0.46). No other law moves that far down: A-mis A2 -0.04, B1S -0.03, R2 0.00, with 80 models each these are within noise. In the other direction B2R gains 0.23 under B-pos and the B-fac laws gain around 0.3.
6. Violation volume, where measured: median share of the domain violated is 0.85 for H1x, 0.36 for H1y, 0.17 for H2, and 0.89, 0.16, 0.67 for R1x, R1y, R2. These are large, which matches the 0B needle result: a sampler finds them.
7. The class rule gives no row in the "not fixed by C" class. The one row marked "rarely violated" is B1I (5/472 in the pool); it is also enforced by structure.

## Limits

- Six systems or structures, one noise level in the arm columns, and a pool that depends on the pilot's candidate grid. The classes are a description of these runs, not a theorem about the laws.
- R is a poor fit on a single series, so its pool shares say little about hares and lynx (see `docs/gate_R_result.md`).
- Item 5 rests on one system (H). It is a candidate finding to test on the epidemic series, not a result.
