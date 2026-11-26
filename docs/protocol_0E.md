# Protocol 0E: COVID-19 reported-case SIR panel (second real series), frozen 2026-11-26, before any fit on these series

Question. On a second real dataset, does the model picked by validation break stated SIR laws, and does the witness loop repair it at low cost? Gate R passed on one series with weak fits (`docs/gate_R_result.md`). A panel of countries gives clusters, so this one can carry intervals.

## Data

Source: Johns Hopkins CSSE COVID-19 repository, time series files as given by the user: confirmed and recovered (global). I also use the deaths file and the UID lookup table (population) from the same repository and folder, because active cases need deaths and the SIR fractions need N. Downloaded 2026-11-26 from raw.githubusercontent.com, saved unchanged in `data_raw/jhu/`:

| file | sha256 |
|---|---|
| time_series_covid19_confirmed_global.csv | e6234a59eec4359d2577358b5220e1a7e3da74c162913cdb7d882db1413f98c2 |
| time_series_covid19_recovered_global.csv | 381bb7527a52114d3f07b40c25dc8aba0b6283aa535e97cba97152ce3cbdb526 |
| time_series_covid19_deaths_global.csv | 4e87757a3e059c45650a1e1856614f8e339ee3c70653c378a7ba6f7b0ee8c72e |
| UID_ISO_FIPS_LookUp_Table.csv | 80bbd734ae7547dac8915d0960d01f43975c817bcddebebd19c026718a1368db |

The recovered series stops being updated on 2021-08-04 (134 of the 201 countries that have both files have their last increase on that date), so every window ends on 2021-08-04. The files run to 2023-03-09 and the later part is not used. Cite JHU CSSE and see the repository for the data terms before reusing the files.

These are reported cases, not infections. Ascertainment, testing policy and the definition of "recovered" differ by country and over time, and nothing here corrects for that.

## States

For a country with population N (UID lookup, country row): C, R, D = cumulative confirmed, recovered, deaths.
S = 1 - C/N, I = (C - R - D)/N, R_state = (R + D)/N. S + I + R_state = 1 by construction. "S" is the population not yet reported, not a measured susceptible count. Units are population fractions, so the data sit in a small corner of the simplex (I is far below 1). I keep that on purpose, as in the earlier testbeds where data cover little of the domain, and I record the largest I and smallest S per country.

## Panel (selection rule fixed here, applied to the files only, no fits)

A country is eligible if all of these hold:
1. It has a single row (blank Province/State) in the confirmed, recovered and deaths files, and a population in the lookup table.
2. Population >= 10 million.
3. Over the window (first day with confirmed >= 100, through 2021-08-04): no missing values, confirmed never decreases, recovered never decreases, recovered last increases on or after 2021-07-29, and I >= 0 on every day.

No hand edits to the data. 37 countries pass; 103 fail on population, 26 on a decreasing recovered series, 15 on recovered stopping early, 6 on decreasing confirmed, 5 have no population, 1 never reaches 100. The panel is the 10 eligible countries with the most confirmed cases on 2021-08-04 (`data_raw/jhu/panel.json`):

| country | population | window start | days | confirmed at end |
|---|---|---|---|---|
| Russia | 146 M | 3/17/20 | 506 | 6.27 M |
| Turkey | 84 M | 3/19/20 | 504 | 5.82 M |
| Argentina | 45 M | 3/20/20 | 503 | 4.98 M |
| Iran | 84 M | 2/26/20 | 526 | 4.02 M |
| Indonesia | 274 M | 3/15/20 | 508 | 3.53 M |
| South Africa | 59 M | 3/18/20 | 505 | 2.48 M |
| Ukraine | 44 M | 3/25/20 | 498 | 2.34 M |
| Iraq | 40 M | 3/13/20 | 510 | 1.67 M |
| Chile | 19 M | 3/15/20 | 508 | 1.62 M |
| Bangladesh | 165 M | 4/6/20 | 486 | 1.31 M |

Many large countries fail rule 3 because their recovered series has downward corrections. The panel is therefore not a random sample of countries, and these ten series share global waves, so they are not independent.

## Setup

- Split by time, per country: first 60 % of days train, next 20 % validation, last 20 % extrapolation window.
- Derivatives: Savitzky-Golay, order 3, window 21 days (primary), applied separately inside each segment. Sensitivity windows 15 and 29. dt = 1 day. Training data are cut into blocks of one window length for the block bootstrap (8 resamples plus the full set, seed 0).
- Library: polynomial degrees 2 and 3 in (S, I, R_state), no structure, pilot STLSQ sweep (12 thresholds, duplicates removed), `gate_b.pool_rows`. Accuracy rule: validation derivative RMSE, normalised by the training derivative sd, equivalent = within 10 % of the best in the pool.
- Primary candidate per country: lowest validation RMSE (ties: fewer terms, lower degree, lower index). Second pick: lowest training AIC. NaN errors count as infinite.

## Laws

The seven SIR laws of system B, unchanged (`src/systems.py`): S, I, R positivity at the faces (B1S, B1I, B1R), S' <= 0 (B2S), R' >= 0 (B2R), and S' + I' + R' within tolerance of 0 (B3p, B3m), on the domain S, I, R >= 0, S + I + R <= 3/2. Tolerances weak 1e-3, strong 2e-2. Verifier: Bernstein branch and bound first, z3 fallback (3 s limit), strong tolerance, as in arm C. Structural family B1I, B1R, B3p, B3m. Non-structural family B1S, B2S, B2R.

These are modelling constraints. The conservation laws hold for the data by construction here, because S + I + R = 1.

## Arms

B-pos and B-fac as in 0A2 (`gate_a2.structure2`, `gate_a2.run_arm2`), arms S, C (witness loop) and Rnd (2000 random points per law), budget 15, margin and scoring as in 0A and 0A2, Rnd stream `default_rng([cand, 5151])`. For the primary, the AIC pick and the first 20 equivalent candidates in order of validation RMSE. Measures: laws verified per arm model, validation RMSE ratio to the candidate's own pool model, extrapolation-window rollout RMSE ratio (rollout from the first smoothed state of the window, pilot settings, 20 s cap), violation volume (200,000 uniform points of the law's domain) and depth for every violated law, replay and a 30-day rollout from each witness of the primary and AIC picks.

## Decision rules (frozen, exact fractions, country is the cluster)

A decision needs at least 8 of the 10 panel countries to finish at the primary window. Otherwise the verdict is undecided.
- **G-E1:** the validation-selected model breaks at least one law in at least half of the finished countries.
- **G-E2:** pooled over countries, at least 1/5 of the accuracy-equivalent candidates break a law, and the lower end of the 95 % cluster bootstrap (resampling countries, 2000 replicates, rng 20261126, sorted replicates 50 and 1949) is above 1/10.
- **G-E3:** B-pos arm C verifies all seven laws on at least 9/10 of the arm-set candidates (pooled), and the median over countries of the median validation ratio and of the median rollout ratio are both <= 1.10.

Verdict. G-E1 and G-E2 pass: the claim holds on a second real series. G-E3 also passes: the cheap fix holds. One of G-E1, G-E2 passes: inconclusive. Neither passes: the claim is dropped for epidemic data.

Reported, not gated: candidate-level C - S and C - Rnd on the non-structural family for B-pos and B-fac, with country-cluster bootstrap intervals; per-country tables; windows 15 and 29.

## Limits

- Reported cases, one pandemic, ten correlated series, a filter that drops many large countries. The window covers 2020-2021 only, and the extrapolation window contains the waves of 2021 that the training part never saw.
- Differencing cumulative counts gives noisy daily derivatives even with a 21-day smoother, and weekday reporting effects remain. Model fits may be poor, as in Gate R. I report validation errors so that can be read off.
- The laws are posed on the whole simplex while the data cover a corner. A violation far from the data says the fit is unconstrained there, not that it is wrong where the data are. The witness table records whether each witness lies outside the observed box.
- The H result (structure alone lowers satisfaction of one law, `docs/law_taxonomy.md`) is a hypothesis to look at here, not a gate. I will report B-pos and B-fac satisfaction per law next to the pool.

## Time limits and code

The limits of protocol 0R addendum 1 apply (60 s per verification step, 180 s per arm run, each config in its own process killed after 100 minutes, 6.5 hours in total). A killed or failed config is listed in `results/0E_lost_configs.csv` and not retried. Jobs are ordered with the primary window first, so any config lost to the total limit is a sensitivity config.

Code: `src/real_e.py`, `src/run_0E.py`, `src/analyze_0E.py`, `tools/kaggle_0E.py`. Outputs: `results/0E_jobs/`, `results/0E_gate.csv`, `0E_summary.csv`, `0E_candidates.csv`, `0E_witnesses.csv`, `0E_lost_configs.csv`. 30 configs (10 countries, 3 windows). Run on Kaggle with the pinned packages of 0B. The pipeline is tested on a simulated SIR series first, and that output is not kept.
