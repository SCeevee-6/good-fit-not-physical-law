# Protocol 0H: is A's under-coverage at 0.5 % noise caused by smoothing bias?

Written 2026-11-28, frozen at the commit that adds this file, before any run on the gate seeds. Background: under the corrected set (0G), A at 0.5 % noise covered the truth in 82, 86 and 89 of 100 seeds at degrees 2, 3 and 4 (`docs/gate_G_result.md`). In that result I guessed the cause was the deterministic bias of the Savitzky-Golay filter, which 0F item 3 had called negligible at 2 %. This unit tests that guess.

## What was seen before the freeze

- Everything in 0F and 0G.
- Development seeds 30-49, noise index 0 (`results/0G_dev/dev_0H_A_n0.csv`). Covered out of 20 at degrees 2 / 3 / 4: base 14 / 16 / 18, oracle 18 / 19 / 20, known sigma 14 / 15 / 17. Median bias noncentrality 3.1 / 3.4 / 3.7, maximum 10.0 / 10.5 / 11.2, against ranks 12 / 20 / 30.
- A pipeline test on seeds 30-31 at both noise indices, including the two window arms (outputs deleted). Stat/eps for the window arms ranged from 0.21 to 1.71.

## Arms

All arms use the 0G corrected set (`src/eiv.py`), system A, degrees 2, 3, 4, and fresh seeds 200-299 at noise indices 0 (0.5 %, primary) and 1 (2 %, control). Code: `src/run_0H.py`.

- base: as in 0G (window 21, estimated sigma).
- oracle: the filter bias delta = D x_clean - Theta(H x_clean) xi_true, computed from the clean trajectories, is subtracted from the derivatives before fitting. It removes every deterministic error between the smoothed data and the true law, filter bias and integrator error alike. It is an oracle, not a method.
- known: base with the true sigma.
- w15, w31: window 15 or 31 instead of 21, sigma re-estimated with that window.
- For the base arm, the bias noncentrality nu = c' V^-1 c, where c is the shift of the identified coordinates caused by delta alone.

## Gates

| rule | question | pass rule |
|---|---|---|
| G-H1 | does removing the filter bias restore coverage at 0.5 %? | oracle covers >= 90 of 100 at each of degrees 2, 3, 4 (noise index 0) |
| G-H2 | is the base arm under-covering on these fresh seeds? | base covers < 90 of 100 at one or more of degrees 2, 3, 4 (noise index 0) |
| G-H3 | is the bias larger relative to the set at low noise? | median nu at noise index 0 is >= 4 times the median at noise index 1, at each degree |

Reported, no prediction: coverage of every arm at both noise levels; w15 and w31 as a practical lever.

## What follows

- G-H1 and G-H2 pass: smoothing bias is the cause of A's low-noise under-coverage, and the write-up says so. Then G-H3 explains why it shows at 0.5 % and not at 2 %: the bias is fixed while the set shrinks with the noise.
- G-H1 fails: something other than deterministic bias is involved at low noise, and the write-up says the cause is open.
- G-H2 fails (base already >= 90 at every degree on fresh seeds): the 0G number was partly sampling. It is reported, and G-H1 is still read.
- The window arms are descriptive. No window is adopted after seeing them.

## Limits

One system and one filter family. The oracle uses the true coefficients, so it shows the size of the problem, not a fix for it.
